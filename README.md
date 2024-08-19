# Newmap

## Introduction

Newmap is a software package that efficiently identifies uniquely mappable
regions of any genome. It accomplishes this task by outputting the minimum
length of a k-mer that uniquely maps to a given position in the genome. From
these ranges of values, the single-read mappability and the multi-read
mappability for a k-length can be generated.

Newmap can search for minimum k-mer lengths for specific values, or entire
continuous ranges using a binary search method allowing for computing all
possible minimum unique k-mer lengths for any sequence.

Newmap requires a CPU that supports the AVX2 instruction set.

## Quick start

### Installation
To install Newmap:
```python
pip install newmap
```
OpenMP is required for parallel processing.

### Create an index for a genome
```bash
newmap generate-index genome.fa
```
By default this creates an `index.awfmi` file in the current directory.

### Find the minimum unique k-mer lengths for a chromosome
For chromosome 1, and lengths ranging from 20 to 200:
```bash
newmap unique-lengths --kmer-batch-size 20000000 --thread-count 4 20:200 index.awfmi chr1.fna.gz
```
This will create a `chr1.unique.uint8` file in the current directory.

### Convert the unique lengths to mappability output files
For k-mer lengths of 24:
```bash
newmap generate-mappability -k 24 -m k24_multiread_mappability.wig -s k24_singleread_mappability.bed chr1.unique.uint8
```
The resulting BED file will be the single read mappbility, and the WIG file will be the multi-read mappability.


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
