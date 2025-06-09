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

Newmap requires a CPU that supports the AVX2 instruction set.

OpenMP is required for parallel processing.

## Documentation / Read The Docs
The latest for Newmap is available on [Read the Docs](https://newmap.readthedocs.io).

## Quick start

### Installation

#### Python Package Index (PyPI)
```python
pip install newmap
```

#### Bioconda
```bash
conda install newmap
```

### Usage

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
