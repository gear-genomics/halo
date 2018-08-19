[![Build Status](https://travis-ci.org/gear-genomics/halo.svg?branch=master)](https://travis-ci.org/gear-genomics/halo)

# halo: haplotype processor

[GEAR](https://gear.embl.de) hosts the `halo` web application at https://gear.embl.de/halo.

## Command line tool

The `halo` command line tool is needed to compute the input file for the `halo` web application.

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
