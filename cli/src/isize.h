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
Contact: Gear Genomics (gear_genomics@embl.de)
============================================================================
*/

#ifndef ISIZE_H
#define ISIZE_H

#define _SECURE_SCL 0
#define _SCL_SECURE_NO_WARNINGS
#include <iostream>
#include <vector>
#include <fstream>

#define BOOST_DISABLE_ASSERTS
#include <boost/unordered_map.hpp>
#include <boost/program_options/cmdline.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/date_time/gregorian/gregorian.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/filesystem.hpp>
#include <boost/progress.hpp>
#include <boost/dynamic_bitset.hpp>

#include <htslib/faidx.h>
#include <htslib/sam.h>
#include <htslib/vcf.h>

namespace halo
{

  struct ISize {
    uint8_t layout;
    int32_t minisize;
    int32_t median;
    int32_t maxisize;
  };    
  
  template<typename TIterator, typename TValue>
  inline void
  getMedian(TIterator begin, TIterator end, TValue& median) {
    std::nth_element(begin, begin + (end - begin) / 2, end);
    median = *(begin + (end - begin) / 2);
  }

  template<typename TIterator, typename TValue>
  inline void
  getMAD(TIterator begin, TIterator end, TValue median, TValue& mad) {
    std::vector<TValue> absDev;
    for(;begin<end;++begin) absDev.push_back(std::abs((TValue)*begin - median));
    getMedian(absDev.begin(), absDev.end(), mad);
  }

  // F+ 0
  // F- 1
  // R+ 2
  // R- 3
  inline uint8_t layout(bam1_t const* rec) {
    if (rec->core.flag & BAM_FREAD1) {
      if (!(rec->core.flag & BAM_FREVERSE)) {
	if (!(rec->core.flag & BAM_FMREVERSE)) return (rec->core.pos < rec->core.mpos) ? 0 : 1;
	else return (rec->core.pos < rec->core.mpos) ? 2 : 3;
      } else {
	if (!(rec->core.flag & BAM_FMREVERSE)) return (rec->core.pos > rec->core.mpos) ? 2 : 3;
	else return (rec->core.pos > rec->core.mpos) ? 0 : 1;
      }
    } else {
      if (!(rec->core.flag & BAM_FREVERSE)) {
	if (!(rec->core.flag & BAM_FMREVERSE)) return (rec->core.pos < rec->core.mpos) ? 1 : 0;
	else return (rec->core.pos < rec->core.mpos) ? 2 : 3;
      } else {
	if (!(rec->core.flag & BAM_FMREVERSE)) return (rec->core.pos > rec->core.mpos) ? 2 : 3;
	else return (rec->core.pos > rec->core.mpos) ? 1 : 0;
      }
    }
  }

  
  template<typename TConfig>
  inline void
    estimateISize(TConfig const& c, std::vector<samFile*> const& samfile, std::vector<hts_idx_t*> const& idx, std::vector<bam_hdr_t*> const& hdr, std::vector<ISize>& isvec) {

    typedef std::vector<int32_t> TISize;
    typedef std::vector<TISize> TISizeOrientation;
    typedef std::vector<TISizeOrientation> TSampleISize;
    TSampleISize si(c.files.size(), TISizeOrientation());
    for(uint32_t i = 0; i < si.size(); ++i) si[i].resize(4, TISize());
    
    boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] " << "Estimate Insert Size" << std::endl;
    boost::progress_display show_progress( hdr[0]->n_targets );
    for (int refIndex = 0; refIndex<hdr[0]->n_targets; ++refIndex) {
      ++show_progress;
      for(unsigned int file_c = 0; file_c < c.files.size(); ++file_c) {
	// Mate map
	typedef boost::unordered_map<std::size_t, bool> TMate;
	TMate mateMap;
	
	hts_itr_t* iter = sam_itr_queryi(idx[file_c], refIndex, 0, hdr[0]->target_len[refIndex]);
	bam1_t* rec = bam_init1();
	int32_t lastAlignedPos = 0;
	std::set<std::size_t> lastAlignedPosReads;
	while (sam_itr_next(samfile[file_c], iter, rec) >= 0) {
	  if (rec->core.flag & (BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP | BAM_FSUPPLEMENTARY | BAM_FUNMAP)) continue;
	  if ((rec->core.qual < c.minMapQual) || (rec->core.tid<0)) continue;
	  
	  // Clean-up the read store for identical alignment positions
	  if (rec->core.pos > lastAlignedPos) {
	    lastAlignedPosReads.clear();
	    lastAlignedPos = rec->core.pos;
	  }
	  
	  if (rec->core.flag & BAM_FPAIRED) {
	    if ((rec->core.mtid<0) || (rec->core.flag & BAM_FMUNMAP)) continue;
	    if (rec->core.mtid != rec->core.tid) continue;
	    if (_firstPairObs(rec, lastAlignedPosReads)) {
	      // First read
	      lastAlignedPosReads.insert(hash_string(bam_get_qname(rec)));
	      std::size_t hv = hash_pair(rec);
	      mateMap[hv] = true;
	    } else {
	      // Second read
	      std::size_t hv = hash_pair_mate(rec);
	      if ((mateMap.find(hv) == mateMap.end()) || (!mateMap[hv])) continue; // Mate discarded
	      mateMap[hv] = false;
	      
	      // Insert size filter
	      int32_t isize = (rec->core.pos + alignmentLength(rec)) - rec->core.mpos;
	      si[file_c][layout(rec)].push_back(isize);
	    }
	  }
	}
	bam_destroy1(rec);
	hts_itr_destroy(iter);
      }
    }

    // Summary stats
    for(unsigned int file_c = 0; file_c < c.files.size(); ++file_c) {    
      // Layout
      uint8_t lyt = 0;
      uint32_t maxcount = si[file_c][lyt].size();
      for(uint32_t i = 1; i < 4; ++i) {
	if (si[file_c][i].size() > maxcount) {
	  maxcount = si[file_c][i].size();
	  lyt = i;
	}
      }

      // Insert size cutoffs
      int32_t med = 0;
      getMedian(si[file_c][lyt].begin(), si[file_c][lyt].end(), med);
      int32_t mad = 0;
      getMAD(si[file_c][lyt].begin(), si[file_c][lyt].end(), med, mad);
      if (mad < 10) mad = 10;
      int32_t minisize = med - 5 * mad;
      if (minisize < 25) minisize = 25;
      int32_t maxisize = med + 5 * mad;

      // Median should be even
      med = (int32_t) (med / 2) * 2;

      std::cout << c.sampleName[file_c] << ",F+=" << si[file_c][0].size() << ",F-=" << si[file_c][1].size() << ",R+=" << si[file_c][2].size() << ",R-=" << si[file_c][3].size() << ",Default=" << (int32_t) lyt << ",Median=" << med << ",Mad=" << mad << ",minISize=" << minisize << ",maxISize=" << maxisize << std::endl;

      // Assign library parameters
      isvec[file_c].layout = lyt;
      isvec[file_c].minisize = minisize;
      isvec[file_c].median = med;
      isvec[file_c].maxisize = maxisize;
    }
  }

	    
	
}

#endif
