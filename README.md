# Newmap

## Introduction

Newmap is a software package that efficiently identifies uniquely mappable
regions of any genome. It accomplishes this task by outputting read lengths at
every position that are unique to that genome. From the range of unique read
lengths produced, the single-read mappability and the multi-read mappability
for a specific read length can be generated.

Newmap can search for unique k-mer/read lengths on specific values, or entire
continuous ranges using a binary search method allowing for finding the
minimum possible unique k-mer/read length.

## Prerequisites
Newmap requires a CPU that supports the AVX2 instruction set.

Newmap requires Python 3.9 or later. It also requires numpy and the AvxWindowFMIndex
library both of which are installed automatically when using the methods below.

Currently only Linux is supported, but it may be possible to build and run on
other operating systems. Notably a compiler with OpenMP support is required
for parallel processing. See the
[documentation](https://newmap.readthedocs.io/en/latest/usage.html#installation)
for more details on building from source.

## Documentation
The latest documentation for Newmap is available on [Read the
Docs](https://newmap.readthedocs.io).

All commands have a `--help` option to provide additional usage information.

## Quick start

### Installation

#### Python Package Index (PyPI)
```python
pip install newmap
```

#### Bioconda
```bash
conda install bioconda::newmap
```

### Usage

#### 0. (Optional) Follow along with an example test genome
You can download a test genome from the Newmap repository to follow along with
the usage example below.
```bash
curl -sL https://raw.githubusercontent.com/hoffmangroup/newmap/refs/heads/master/tests/data/genome.fa > genome.fa
```
To speed up creating the index in the following step, it is recommended to use
options `--seed-length=1` and `--compression-ratio=1` specifically for the very
small test genome. Otherwise it would be recommended to use the default values.

#### 1. Create an index for a genome
```bash
newmap index genome.fa
```
By default this will create a `genome.awfmi` file in the current directory.

#### 2. Find the minimum unique k-mer lengths for the genome using the index
Searching the entire genome, using 20 threads, printing status information, and
searching lengths ranging from 20 to 200 bp:
```bash
newmap search --verbose --num-threads=20 --search-range=20:200 --output-directory=unique_lengths genome.fa
```
This will create `*.unique.uint8` files (one for each sequence ID) in the `unique_lengths` directory.

#### 3. Convert the unique lengths to mappability tracks
To output single-read and multi-read mappability for a 24 bp read length:
```bash
newmap track --single-read=24.bed --multi-read=24.wig 24 unique_lengths/*.unique.uint8
```
For both single-read and multi-read mappability, this will generate a single
file that contains the mappability for all sequences listed in the
`unique_lengths` directory.
The resulting BED file will be the single read mappability, and the WIG file
will be the multi-read mappability.


Credits
-------
Newmap is a reimplementation of the output of Umap. Umap was developed by Mehran Karimzadeh.
The repository for that implmemention is found at [https://www.github.com/hoffmangroup/umap](https://www.github.com/hoffmangroup/umap).
Umap in turn was originally developed by Anshul Kundaje and was written in MATLAB.
The original repository is available [https://sites.google.com/site/anshulkundaje/projects/mappability](https://sites.google.com/site/anshulkundaje/projects/mappability).

This project uses the excellent [AvxWindowFMIndex
library](https://github.com/TravisWheelerLab/AvxWindowFmIndex). Read their
published article here
([https://doi.org/10.1186/s13015-021-00204-6](https://doi.org/10.1186/s13015-021-00204--6))
