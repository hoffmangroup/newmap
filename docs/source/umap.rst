Umap
====

Newmap is intended to be a simple, faster, and easy to use replacement for the
orignal `Umap <https://bismap.hoffmanlab.org>`_ project.

------------
Reproduction
------------

To reproduce the ``unique.uint8`` datasets that would have been generated from
Umap, the following criteria must be met:

1. K-mers overlapping specifically with dinucleotide sequence ``N`` are not
   considered unique. This is the default.
2. Only the following k-mer lengths should be included: 24,36,50,100,150,200
