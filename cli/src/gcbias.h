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

#ifndef GCBIAS_H
#define GCBIAS_H

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

#include "util.h"
#include "isize.h"

namespace halo
{

  struct GcBiasConfig {
    uint16_t minMapQual;
    uint32_t minchrsize;
    int32_t percentid;
    std::vector<std::string> sampleName;
    boost::filesystem::path outfile;
    boost::filesystem::path genome;
    std::vector<boost::filesystem::path> files;
  };


  template<typename TConfig, typename TSampleGc>
  inline void
  loadGcBias(TConfig const& c, TSampleGc& sgc) {
    if (c.gcbiasprof) {
      // Parse GC bias profile
      std::ifstream file(c.gcbias.string().c_str(), std::ios_base::in | std::ios_base::binary);
      boost::iostreams::filtering_streambuf<boost::iostreams::input> dataIn;
      dataIn.push(boost::iostreams::gzip_decompressor());
      dataIn.push(file);
      std::istream instream(&dataIn);
      std::string gline;
      while(std::getline(instream, gline)) {
	if ((gline.size()) && (gline[0] == '#')) continue;
	typedef boost::tokenizer< boost::char_separator<char> > Tokenizer;
	boost::char_separator<char> sep("\t");
	Tokenizer tokens(gline, sep);
	Tokenizer::iterator tokIter = tokens.begin();
	if (tokIter!=tokens.end()) {
	  std::string sample = *tokIter++;
	  int32_t gccount = boost::lexical_cast<int32_t>(*tokIter++);
	  ++tokIter; // Skip GC fraction
	  int32_t refcount = boost::lexical_cast<int32_t>(*tokIter++);
	  int32_t samplecount = boost::lexical_cast<int32_t>(*tokIter++);
	  double obsexp = boost::lexical_cast<double>(*tokIter++);
	  for(uint32_t i = 0; i < c.sampleName.size(); ++i) {
	    if (c.sampleName[i] == sample) {
	      if (obsexp > 1.75) sgc[i][gccount] = 1.75;
	      else if (obsexp < 0.25) sgc[i][gccount] = 0.25;
	      else sgc[i][gccount] = obsexp;
	    }
	  }
	}
      }
      dataIn.pop();
    }
  }
  
  int gcbias(int argc, char **argv) {

#ifdef PROFILE
    ProfilerStart("gcbias.prof");
#endif
    
    GcBiasConfig c;

    // Parameter
    boost::program_options::options_description generic("Generic options");
    generic.add_options()
      ("help,?", "show help message")
      ("genome,g", boost::program_options::value<boost::filesystem::path>(&c.genome), "genome fasta file")
      ("map-qual,q", boost::program_options::value<uint16_t>(&c.minMapQual)->default_value(1), "min. mapping quality")
      ("percentid,p", boost::program_options::value<int32_t>(&c.percentid)->default_value(98), "min. required percent identity")
      ("minchrsize,m", boost::program_options::value<uint32_t>(&c.minchrsize)->default_value(10000000), "min. chr size")
      ("outfile,o", boost::program_options::value<boost::filesystem::path>(&c.outfile)->default_value("gc.prof"), "output file")
      ;
    
    boost::program_options::options_description hidden("Hidden options");
    hidden.add_options()
      ("input-file", boost::program_options::value< std::vector<boost::filesystem::path> >(&c.files), "input bam files")
      ;
    
    boost::program_options::positional_options_description pos_args;
    pos_args.add("input-file", -1);
    
    boost::program_options::options_description cmdline_options;
    cmdline_options.add(generic).add(hidden);
    boost::program_options::options_description visible_options;
    visible_options.add(generic);
    boost::program_options::variables_map vm;
    boost::program_options::store(boost::program_options::command_line_parser(argc, argv).options(cmdline_options).positional(pos_args).run(), vm);
    boost::program_options::notify(vm);
    
    // Check command line arguments
    if ((vm.count("help")) || (!vm.count("input-file")) || (!vm.count("genome"))) {
      std::cout << "Usage: halo " << argv[0] << " [OPTIONS] -g <ref.fa> <align1.bam> <align2.bam> ... <alignN.bam>" << std::endl;
      std::cout << visible_options << "\n";
      return 1;
    } 
    
    // Show cmd
    boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] ";
    std::cout << "halo ";
    for(int i=0; i<argc; ++i) { std::cout << argv[i] << ' '; }
    std::cout << std::endl;

    // Check BAM files
    for(unsigned int file_c = 0; file_c < c.files.size(); ++file_c) {
      if (!(boost::filesystem::exists(c.files[file_c]) && boost::filesystem::is_regular_file(c.files[file_c]) && boost::filesystem::file_size(c.files[file_c]))) {
	std::cerr << "Input BAM file is missing: " << c.files[file_c].string() << std::endl;
	return 1;
      }
    }

    // Check reference
    if (!(boost::filesystem::exists(c.genome) && boost::filesystem::is_regular_file(c.genome) && boost::filesystem::file_size(c.genome))) {
      std::cerr << "Reference file is missing: " << c.genome.string() << std::endl;
      return 1;
    } else {
      faidx_t* fai = fai_load(c.genome.string().c_str());
      if (fai == NULL) {
	if (fai_build(c.genome.string().c_str()) == -1) {
	  std::cerr << "Fail to open genome fai index for " << c.genome.string() << std::endl;
	  return 1;
	} else fai = fai_load(c.genome.string().c_str());
      }
      fai_destroy(fai);
    }
    
    // Load bam files
    typedef std::vector<samFile*> TSamFile;
    typedef std::vector<hts_idx_t*> TIndex;
    typedef std::vector<bam_hdr_t*> THeader;
    TSamFile samfile(c.files.size());
    TIndex idx(c.files.size());
    THeader hdr(c.files.size());
    c.sampleName.resize(c.files.size());
    for(unsigned int file_c = 0; file_c < c.files.size(); ++file_c) {
      samfile[file_c] = sam_open(c.files[file_c].string().c_str(), "r");
      if (samfile[file_c] == NULL) {
	std::cerr << "Fail to open file " << c.files[file_c].string() << std::endl;
	return 1;
      }
      idx[file_c] = sam_index_load(samfile[file_c], c.files[file_c].string().c_str());
      if (idx[file_c] == NULL) {
	std::cerr << "Fail to open index for " << c.files[file_c].string() << std::endl;
	return 1;
      }
      hdr[file_c] = sam_hdr_read(samfile[file_c]);
      if (hdr[file_c] == NULL) {
	std::cerr << "Fail to open header for " << c.files[file_c].string() << std::endl;
	return 1;
      }
      faidx_t* fai = fai_load(c.genome.string().c_str());
      for(int32_t refIndex=0; refIndex < hdr[file_c]->n_targets; ++refIndex) {
	std::string tname(hdr[file_c]->target_name[refIndex]);
	if (!faidx_has_seq(fai, tname.c_str())) {
	  std::cerr << "BAM file chromosome " << hdr[file_c]->target_name[refIndex] << " is NOT present in your reference file " << c.genome.string() << std::endl;
	  return 1;
	}
      }
      fai_destroy(fai);
      std::string sampleName;
      if (!getSMTag(std::string(hdr[file_c]->text), c.files[file_c].stem().string(), sampleName)) {
	std::cerr << "Only one sample (@RG:SM) is allowed per input BAM file " << c.files[file_c].string() << std::endl;
	return 1;
      } else c.sampleName[file_c] = sampleName;
    }

    // Estimate insert size
    typedef std::vector<ISize> TISizeInfo;
    TISizeInfo isize(c.files.size(), ISize());
    estimateISize(c, samfile, idx, hdr, isize);

    // GC-Bias
    typedef std::vector<int32_t> TGcCount;
    typedef std::vector<TGcCount> TGenomicGc;
    TGenomicGc gcbias(c.files.size(), TGcCount()); // Sample GC at given ISize
    TGenomicGc refgc(c.files.size(), TGcCount()); // Reference GC at given ISize
    for(unsigned int file_c = 0; file_c < c.files.size(); ++file_c) {
      refgc[file_c].resize(isize[file_c].median + 2, 0);
      gcbias[file_c].resize(isize[file_c].median + 2, 0);
    }
    
    // Parse bam (contig by contig)
    now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] " << "Estimate GC bias" << std::endl;
    boost::progress_display show_progress( hdr[0]->n_targets );
    for (int refIndex = 0; refIndex<hdr[0]->n_targets; ++refIndex) {
      //for (int refIndex = 20; refIndex<21; ++refIndex) {
      ++show_progress;
      if (hdr[0]->target_len[refIndex] < c.minchrsize) continue;

      // Ignore sex chromosomes
      std::string cn(hdr[0]->target_name[refIndex]);
      if ((cn == "chrX")  || (cn == "X") || (cn == "chrY") || (cn == "Y")) continue;
      
      // Load chromosome
      char* seq = NULL;
      faidx_t* fai = fai_load(c.genome.string().c_str());
      int32_t seqlen = -1;
      std::string tname(hdr[0]->target_name[refIndex]);
      seq = faidx_fetch_seq(fai, tname.c_str(), 0, hdr[0]->target_len[refIndex], &seqlen);
      
      // GC- and N-content
      typedef boost::dynamic_bitset<> TBitSet;
      TBitSet nrun(hdr[0]->target_len[refIndex], false);
      TBitSet gcref(hdr[0]->target_len[refIndex], false);
      for(uint32_t i = 0; i < hdr[0]->target_len[refIndex]; ++i) {
	if ((seq[i] == 'c') || (seq[i] == 'C') || (seq[i] == 'g') || (seq[i] == 'G')) gcref[i] = 1;
	if ((seq[i] == 'n') || (seq[i] == 'N')) nrun[i] = 1;
      }

      // Reference GC
      for(unsigned int file_c = 0; file_c < c.files.size(); ++file_c) {
	if (( (int32_t) hdr[0]->target_len[refIndex] > isize[file_c].median + 2) && (hdr[0]->target_len[refIndex] >= c.minchrsize)) {
	  int32_t halfwin = (int32_t) (isize[file_c].median / 2);
	  int32_t nsum = 0;
	  int32_t gcsum = 0;
	  for(int32_t pos = halfwin; pos < (int32_t) hdr[0]->target_len[refIndex] - halfwin; ++pos) {
	    if (pos == halfwin) {
	      for(int32_t i = pos - halfwin; i<pos+halfwin+1; ++i) {
		nsum += nrun[i];
		gcsum += gcref[i];
	      }
	    } else {
	      nsum -= nrun[pos - halfwin - 1];
	      gcsum -= gcref[pos - halfwin - 1];
	      nsum += nrun[pos + halfwin];
	      gcsum += gcref[pos + halfwin];
	    }
	    //std::string refslice = boost::to_upper_copy(std::string(seq + pos - halfwin, seq + pos + halfwin + 1));
	    //std::cerr << refslice << ',' << refslice.size() << ',' << gcsum << ',' << nsum << std::endl;
	    if (!nsum) ++refgc[file_c][gcsum];
	  }
	}
      }

      // Sample GC
      for(unsigned int file_c = 0; file_c < c.files.size(); ++file_c) {
	// Qualities and alignment length
	typedef boost::unordered_map<std::size_t, bool> TMateMap;
	TMateMap mateMap;

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

	  // Sequence
	  std::string sequence;
	  sequence.resize(rec->core.l_qseq);
	  uint8_t* seqptr = bam_get_seq(rec);
	  for (int i = 0; i < rec->core.l_qseq; ++i) sequence[i] = "=ACMGRSVTWYHKDBN"[bam_seqi(seqptr, i)];

	  // Reference slice
	  std::string refslice = boost::to_upper_copy(std::string(seq + rec->core.pos, seq + lastAlignedPosition(rec)));
	      
	  // Percent identity
	  uint32_t rp = 0; // reference pointer
	  uint32_t sp = 0; // sequence pointer
	  uint32_t* cigar = bam_get_cigar(rec);
	  int32_t matchCount = 0;
	  int32_t mismatchCount = 0;
	  for (std::size_t i = 0; i < rec->core.n_cigar; ++i) {
	    if ((bam_cigar_op(cigar[i]) == BAM_CMATCH) || (bam_cigar_op(cigar[i]) == BAM_CEQUAL) || (bam_cigar_op(cigar[i]) == BAM_CDIFF)) {
	      // match or mismatch
	      for(std::size_t k = 0; k<bam_cigar_oplen(cigar[i]);++k) {
		if (sequence[sp] == refslice[rp]) ++matchCount;
		else ++mismatchCount;
		++sp;
		++rp;
	      }
	    } else if (bam_cigar_op(cigar[i]) == BAM_CDEL) {
	      ++mismatchCount;
	      rp += bam_cigar_oplen(cigar[i]);
	    } else if (bam_cigar_op(cigar[i]) == BAM_CINS) {
	      ++mismatchCount;
	      sp += bam_cigar_oplen(cigar[i]);
	    } else if (bam_cigar_op(cigar[i]) == BAM_CSOFT_CLIP) {
	      sp += bam_cigar_oplen(cigar[i]);
	    } else if(bam_cigar_op(cigar[i]) == BAM_CHARD_CLIP) {
	    } else if (bam_cigar_op(cigar[i]) == BAM_CREF_SKIP) {
	      rp += bam_cigar_oplen(cigar[i]);
	    } else {
	      std::cerr << "Unknown Cigar options" << std::endl;
	      return 1;
	    }
	  }
	  double percid = 0;
	  if (matchCount + mismatchCount > 0) percid = (double) matchCount / (double) (matchCount + mismatchCount);
	  if (percid < ((double) c.percentid / (double) 100)) continue;
	  
	  
	  // Paired-end data
	  if (rec->core.flag & BAM_FPAIRED) {
	    if ((rec->core.mtid<0) || (rec->core.flag & BAM_FMUNMAP)) continue;
	    // Ignore translocation paired-ends
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

	      // Read-pair orientation
	      if (layout(rec) != isize[file_c].layout) continue;
	      
	      // Insert size filter
	      int32_t is = (rec->core.pos + alignmentLength(rec)) - rec->core.mpos;
	      if ((is < isize[file_c].minisize) || (is > isize[file_c].maxisize)) continue;

	      // Count fragment mid-points
	      int32_t halfwin = (int32_t) (isize[file_c].median / 2);
	      int32_t pos = rec->core.mpos + (int32_t) (is/2);
	      int32_t fragstart = pos - halfwin;
	      int32_t fragend = pos + halfwin + 1;
	      if ((fragstart >= 0) && (fragend < (int32_t) hdr[0]->target_len[refIndex])) {
		int32_t ncount = 0;
		for(int32_t i = fragstart; i < fragend; ++i) {
		  if (nrun[i]) ++ncount;
		}
		if (!ncount) {
		  int32_t gccount = 0;
		  for(int32_t i = fragstart; i < fragend; ++i) {
		    if (gcref[i]) ++gccount;
		  }
		  ++gcbias[file_c][gccount];
		}
	      }
	    }
	  }
	}
	bam_destroy1(rec);
	hts_itr_destroy(iter);
      }
      if (seq != NULL) free(seq);
      fai_destroy(fai);
    }

    
    // Output
    now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] " << "Output GC bias" << std::endl;
    boost::progress_display sp( c.files.size() );

    boost::iostreams::filtering_ostream rfile;
    rfile.push(boost::iostreams::gzip_compressor());
    rfile.push(boost::iostreams::file_sink(c.outfile.string().c_str(), std::ios_base::out | std::ios_base::binary));

    // Library statistics
    for(unsigned int file_c = 0; file_c < c.files.size(); ++file_c) {
      rfile << "#isize," << c.sampleName[file_c] << ',' << (int32_t) isize[file_c].layout << ',' << isize[file_c].minisize << ',' << isize[file_c].median << ',' << isize[file_c].maxisize << std::endl;
    }

    // GC bias
    rfile << "#sample\tgcsum\tgcfrac\treference\tobservation\tobsexp" << std::endl;
    for(unsigned int file_c = 0; file_c < c.files.size(); ++file_c) {
      ++sp;
      double totalSample = 0;
      double totalRef = 0;
      for(uint32_t i = 0; i < refgc[file_c].size(); ++i) {
	totalSample += gcbias[file_c][i];
	totalRef += refgc[file_c][i];
      }
      double norm = totalRef / totalSample;
      double fragsize = refgc[file_c].size();
      for(uint32_t i = 0; i < refgc[file_c].size(); ++i) {
	double frac = (double) (i) / fragsize;
	// Debug: Raw counts, Header: std::cerr << "sample\tgcsum\tgcfrac\tcount\ttype" << std::endl;
	//std::cerr << c.sampleName[file_c] << "\t" << i << "\t" << frac << "\t" << refgc[file_c][i] << "\tReference" << std::endl;
	//std::cerr << c.sampleName[file_c] << "\t" << i << "\t" << frac << "\t" << norm * gcbias[file_c][i] << "\tSample" << std::endl;
	// Debug: Observed/Expected
	rfile << c.sampleName[file_c] << "\t" << i << "\t" << frac << "\t" << refgc[file_c][i] << "\t" << gcbias[file_c][i] << "\t" << (norm * (double) gcbias[file_c][i]) / (double) refgc[file_c][i] << std::endl;
      }
    }
    rfile.pop();
    
    // End
    now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] Done." << std::endl;

    // Close bam
    for(unsigned int file_c = 0; file_c < c.files.size(); ++file_c) {
      hts_idx_destroy(idx[file_c]);
      sam_close(samfile[file_c]);
      bam_hdr_destroy(hdr[file_c]);
    }

#ifdef PROFILE
    ProfilerStop();
#endif
    return 0;
  }

}

#endif
