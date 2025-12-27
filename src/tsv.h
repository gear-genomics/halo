#ifndef TSV_H
#define TSV_H

#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/device/file.hpp>

namespace halo
{
  

  template<typename TConfig, typename THeader, typename TSampleWC>
  inline void
  haloTsvOut(TConfig const& c, THeader const& hdr, TSampleWC const& sWC) {
    boost::iostreams::filtering_ostream rfile;
    rfile.push(boost::iostreams::gzip_compressor());
    rfile.push(boost::iostreams::file_sink(c.outfile.string().c_str(), std::ios_base::out | std::ios_base::binary));
    rfile << "sample\tchr\tpos\twatson\tcrick" << std::endl;
    for(uint32_t i = 0; i < c.sampleName.size(); ++i) {
      for (int32_t refIndex = 0; refIndex<hdr[i]->n_targets; ++refIndex) {
	if (!sWC[0][refIndex].size()) continue;
	int32_t pos = 0;
	for (uint32_t k = 0; k < sWC[i][refIndex].size(); ++k) {
	  rfile << c.sampleName[i] << "\t" << hdr[i]->target_name[refIndex] << "\t" << pos << "\t" << sWC[i][refIndex][k].first << "\t" << sWC[i][refIndex][k].second << std::endl;
	  pos += c.window;
	}
      }
    }
    rfile.pop();
  }
 
}

#endif
