[![C/C++ CI](https://github.com/gear-genomics/halo/workflows/C/C++%20CI/badge.svg)](https://github.com/gear-genomics/halo/actions)
[![Docker CI](https://github.com/gear-genomics/halo/workflows/Docker%20CI/badge.svg)](https://hub.docker.com/r/geargenomics/halo/)
[![GitHub license](https://img.shields.io/badge/License-BSD%203--Clause-blue.svg)](https://github.com/gear-genomics/halo/blob/master/LICENSE)

# Halo: Haplotype processor

[GEAR](https://www.gear-genomics.com/) hosts the `halo` web application at [https://www.gear-genomics.com/halo/](https://www.gear-genomics.com/halo/). The `halo` command-line tool is needed to compute the input file for the `halo` web application.

### Installation

```bash
git clone --recurse-submodules https://github.com/gear-genomics/halo.git
cd halo/cli
make all
make install
```

### Usage

```bash
./bin/halo -h
```
