/*
============================================================================
Haplotype-resolved data analysis methods
============================================================================
Copyright (C) 2018

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
============================================================================
Contact: Gear genomics (gear_genomics@embl.de)
============================================================================
*/

#ifndef UTIL_H
#define UTIL_H

#include <boost/algorithm/string.hpp>
#include <boost/tokenizer.hpp>
#include <boost/math/distributions/chi_squared.hpp>
#include <boost/math/distributions/hypergeometric.hpp>


namespace halo
{

  inline bool
  getSMTag(std::string const& header, std::string const& fileName, std::string& sampleName) {
    std::set<std::string> smIdentifiers;
    std::string delimiters("\n");
    typedef std::vector<std::string> TStrParts;
    TStrParts lines;
    boost::split(lines, header, boost::is_any_of(delimiters));
    TStrParts::const_iterator itH = lines.begin();
    TStrParts::const_iterator itHEnd = lines.end();
    bool rgPresent = false;
    for(;itH!=itHEnd; ++itH) {
      if (itH->find("@RG")==0) {
	std::string delim("\t");
	TStrParts keyval;
	boost::split(keyval, *itH, boost::is_any_of(delim));
	TStrParts::const_iterator itKV = keyval.begin();
	TStrParts::const_iterator itKVEnd = keyval.end();
	for(;itKV != itKVEnd; ++itKV) {
	  size_t sp = itKV->find(":");
	  if (sp != std::string::npos) {
	    std::string field = itKV->substr(0, sp);
	    if (field == "ID") {
	      rgPresent = true;
	      std::string rgSM = itKV->substr(sp+1);
	      smIdentifiers.insert(rgSM);
	    }
	  }
	}
      }
    }
    if (!rgPresent) {
      sampleName = fileName;
      return true;
    } else if (smIdentifiers.size() == 1) {
      sampleName = *(smIdentifiers.begin());
      return true;
    } else {
      sampleName = "";
      return false;
    }
  }

  inline unsigned hash_string(const char *s) {
    unsigned h = 37;
    while (*s) {
      h = (h * 54059) ^ (s[0] * 76963);
      s++;
    }
    return h;
  }
  
  inline std::size_t hash_pair(bam1_t* rec) {
    std::size_t seed = hash_string(bam_get_qname(rec));
    boost::hash_combine(seed, rec->core.tid);
    boost::hash_combine(seed, rec->core.pos);
    boost::hash_combine(seed, rec->core.mtid);
    boost::hash_combine(seed, rec->core.mpos);
    return seed;
  }

  inline std::size_t hash_pair_mate(bam1_t* rec) {
    std::size_t seed = hash_string(bam_get_qname(rec));
    boost::hash_combine(seed, rec->core.mtid);
    boost::hash_combine(seed, rec->core.mpos);
    boost::hash_combine(seed, rec->core.tid);
    boost::hash_combine(seed, rec->core.pos);
    return seed;
  }


  template<typename TAlignedReads>
  inline bool
  _firstPairObs(bam1_t* rec, TAlignedReads const& lastAlignedPosReads) {
    if (rec->core.tid == rec->core.mtid) return ((rec->core.pos < rec->core.mpos) || ((rec->core.pos == rec->core.mpos) && (lastAlignedPosReads.find(hash_string(bam_get_qname(rec))) == lastAlignedPosReads.end())));
    else return (rec->core.tid < rec->core.mtid);
  }
  
  
  inline uint32_t alignmentLength(bam1_t const* rec) {
    uint32_t* cigar = bam_get_cigar(rec);
    uint32_t alen = 0;
    for (std::size_t i = 0; i < rec->core.n_cigar; ++i)
      if ((bam_cigar_op(cigar[i]) == BAM_CMATCH) || (bam_cigar_op(cigar[i]) == BAM_CEQUAL) || (bam_cigar_op(cigar[i]) == BAM_CDIFF) || (bam_cigar_op(cigar[i]) == BAM_CDEL) || (bam_cigar_op(cigar[i]) == BAM_CREF_SKIP)) alen += bam_cigar_oplen(cigar[i]);
    return alen;
  }

  inline uint32_t
  lastAlignedPosition(bam1_t const* rec) {
    return rec->core.pos + alignmentLength(rec);
  }

  inline uint32_t halfAlignmentLength(bam1_t const* rec) {
    return (alignmentLength(rec) / 2);
  }

  template<typename TIterator, typename TValue>
  inline void
  _getMedian(TIterator begin, TIterator end, TValue& median) {
    std::nth_element(begin, begin + (end - begin) / 2, end);
    median = *(begin + (end - begin) / 2);
  }

  template<typename TIterator, typename TValue>
  inline void
  _getMAD(TIterator begin, TIterator end, TValue median, TValue& mad) {
    std::vector<TValue> absDev;
    for(;begin<end;++begin)
      absDev.push_back(std::abs((TValue)*begin - median));
    _getMedian(absDev.begin(), absDev.end(), mad);
  }


  template<typename TCounter>
  inline uint32_t
  _getReadSupportPercentile(TCounter const& watsonCount, TCounter const& crickCount, uint32_t perc) {
    TCounter support(watsonCount.size(), 0);
    typename TCounter::iterator itSupport = support.begin();
    typename TCounter::const_iterator itCrick =crickCount.begin();
    for(typename TCounter::const_iterator itWatson = watsonCount.begin(); itWatson != watsonCount.end(); ++itWatson, ++itCrick, ++itSupport)
      *itSupport = *itWatson + *itCrick;
    std::sort(support.begin(), support.end());
    return support[(int32_t) ((perc * watsonCount.size()) / 100)];
  }

  template<typename TWRatioVector>
  inline std::pair<float, float>
  _biModalSynchronizedMinima(TWRatioVector& wRatio) {
    std::sort(wRatio.begin(), wRatio.end());
    uint32_t binCount = (std::pow(wRatio.size(), 0.1/0.3) * (wRatio[wRatio.size()-1] - wRatio[0])) / (wRatio[(int) (3*wRatio.size()/4)] - wRatio[(int) (wRatio.size()/4)]);
    std::vector<uint32_t> histogram(binCount + 1, 0);
    for(typename TWRatioVector::iterator itW = wRatio.begin(); itW != wRatio.end(); ++itW) ++histogram[(uint32_t)((*itW)*binCount)];
    uint32_t smallestBinVal = histogram[0] + histogram[binCount];
    float crickCut = 0;
    float watsonCut = 1;
    for(std::size_t i = 0; i<binCount/2; ++i) {
      if (histogram[i] + histogram[binCount - i] < smallestBinVal) {
	smallestBinVal = histogram[i] + histogram[binCount - i];
	crickCut = (float) i / (float) binCount;
	watsonCut = (float) (binCount - i) / (float) binCount;
      }
    }
    return std::make_pair(crickCut, watsonCut);
  }

  template<typename TWRatioVector>
  inline std::pair<float, float>
  _biModalMinima(TWRatioVector& wRatio) {
    std::sort(wRatio.begin(), wRatio.end());
    uint32_t binCount = (std::pow(wRatio.size(), 0.1/0.3) * (wRatio[wRatio.size()-1] - wRatio[0])) / (wRatio[(int) (3*wRatio.size()/4)] - wRatio[(int) (wRatio.size()/4)]);
    std::vector<uint32_t> histogram(binCount + 1, 0);
    for(typename TWRatioVector::iterator itW = wRatio.begin(); itW != wRatio.end(); ++itW) ++histogram[(uint32_t)((*itW)*binCount)];
    uint32_t smallestBinVal = histogram[0];
    float crickCut = 0;
    for(std::size_t i = 0; i<binCount/2; ++i) {
      if (histogram[i] < smallestBinVal) {
	smallestBinVal = histogram[i];
	crickCut = (float) i / (float) binCount;
      }
    }
    smallestBinVal = histogram[binCount];
    float watsonCut = 1;
    for(std::size_t i = binCount; i>binCount/2; --i) {
      if (histogram[i] < smallestBinVal) {
        smallestBinVal = histogram[i];
        watsonCut = (float) i / (float) binCount;
      }
    }
    return std::make_pair(crickCut, watsonCut);
  }


  template<typename TCount, typename TPrecision>
  inline void
  _fisher(TCount const a, TCount const b, TCount const c, TCount const d, TPrecision& pval) {
    TCount N = a + b + c + d;
    TCount r = a + c;
    TCount s = c + d;
    TCount max_for_k = std::min(r, s);
    TCount min_for_k = (TCount) std::max(0, int(r + s - N));
    boost::math::hypergeometric_distribution<> hgd(r, s, N);
    TPrecision cutoff = pdf(hgd, c);
    pval = 0.0;
    for(int k = min_for_k; k <= max_for_k; ++k) {
      TPrecision p = pdf(hgd, k);
      if (p <= cutoff) pval += p;
    }
  }

  template<typename TVector>
  inline std::pair<float, float>
  _getMedianMAD(TVector& vec) {
    float median = 0;
    _getMedian(vec.begin(), vec.end(), median);
    float mad = 0;
    _getMAD(vec.begin(), vec.end(), median, mad);
    return std::make_pair(median, mad);
  }

}

#endif
