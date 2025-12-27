[![C/C++ CI](https://github.com/gear-genomics/halo/workflows/C/C++%20CI/badge.svg)](https://github.com/gear-genomics/halo/actions)
[![Docker CI](https://github.com/gear-genomics/halo/workflows/Docker%20CI/badge.svg)](https://hub.docker.com/r/geargenomics/halo/)
[![GitHub Releases](https://img.shields.io/github/release/gear-genomics/halo.svg)](https://github.com/gear-genomics/halo/releases)

# Halo: Haplotype processor

[GEAR](https://www.gear-genomics.com/) hosts the `halo` web application at [https://www.gear-genomics.com/halo/](https://www.gear-genomics.com/halo/). The `halo` command-line tool is needed to compute the input file for the `halo` web application.

## Installing halo

The `halo` command-line tool is available as a pre-compiled statically linked binary from [Halo's github release page](https://github.com/gear-genomics/halo/releases/), as a singularity container [SIF file](https://github.com/gear-genomics/halo/releases/) or as a minimal [Docker container](https://hub.docker.com/r/geargenomics/halo/).

```bash
git clone --recurse-submodules https://github.com/gear-genomics/halo.git
cd halo/cli
make all
make install
```

## Usage

Halo is used to create the input JSON file for the [web application](https://www.gear-genomics.com/halo/).

```bash
./bin/halo -h
```

This command-line application was designed for [Strand-Seq]() single-cell DNA data but the viewer can also be used for any haplotype-resolved data. To convert the Strand-Seq alignment data into JSON format please provide one input BAM file for each input cell, e.g.:

```bash
halo count -g hg38.fa -o sc.json.gz sc1.bam sc2.bam sc3.bam
```

The output file `sc.json.gz` can then be uploaded and visualized at [https://www.gear-genomics.com/halo/](https://www.gear-genomics.com/halo/).

License
-------
Halo is distributed under the GPLv3 license. Consult the accompanying [LICENSE](https://github.com/gear-genomics/halo/blob/main/LICENSE) file for more details.

Credits
-------
[HTSlib](https://github.com/samtools/htslib) is heavily used for all genomic alignment variant processing. [Boost](https://www.boost.org/) for various data structures and algorithms.
