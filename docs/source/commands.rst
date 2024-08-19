.. _commands:

Commands
========

All commands start with the ``newmap`` prefix. Help with is available for each
command by running ``newmap <command> --help``.

.. _generate-index:

--------------
generate-index
--------------
Generates an index file for a given sequence file.

Options
-------
- `index-file`: The name of the index file to generate. Defaults to index.awfmi.

FM-index parameters
-------------------
- `suffix-array-compression-ratio`: The compression ratio of the suffix array. Defaults to 8.
- `kmer-length-in-seed-table`: The length of the k-mer in the seed table. Defaults to 12.

Positional Arguments
--------------------
- `sequence-file`: The name of the input sequence file to generate an index for. Required.

Example:

.. code-block:: console

    $ newmap generate-index --index-file=hg38.awfmi hg38.fa

This will generate an index file named `hg38.awfmi` for the sequence `hg38.fa`.

FM-index technical details
^^^^^^^^^^^^^^^^^^^^^^^^^^
The default parameters are the recommended set to be used for matching
dinucleotide sequences and likely do not need to be changed. The parameters may
be changed for technical reasons trading off disk space and/or memory available
to adjust performance. Each increase in the `suffix-array-compression-ratio`
reduces the index file size at the cost of number of operations to get a count
on the occurances of a given k-mer. Each increase in the
`kmer-length-in-seed-table` increases the memory required to speed up k-mer
searches in the index. Each increase by 1, multiplies the memory usage of the
index by 4.


.. _unique-lengths:

--------------
unique-lengths
--------------
Generates a binary file containing the minimum unique lengths of sequence found
to be unique at each position from a given range of k-mer lengths. See the
:ref:`unique-file-format` for details on the file output.

Options
-------
- `initial-search-length`: The initial k-mer length to search for unique sequences.
  Only valid when the set of lengths of k-mer lengths is a continuous range
  with the ``kmer-lengths`` postional argument (which is a pair of values
  separated by a colon). Useful to use when the majority of largest minimum
  unique lengths are likely to be much smaller the maximum search length from
  your specified range.
- `kmer-batch-size`: The maximum number of sequence positions to search for at
  a time per sequence ID. Useful for controlling memory requirements. Default
  is 1000000.
- `thread-count`: The number of threads to use for counting on the index.
  Default is 1.
- `verbose`: Print verbose output. Includes summary statistics at end of each
  sequence. Default is False.

Positional Arguments
--------------------
- `kmer-lengths`: The range of k-mer lengths to search for unique sequences. A
  colon seperated pair of values specifies a continuous range. A comma
  seperated list specifies specific lengths to search.
- `index-file`: The name of the index file to use for searching for unique sequences.
- `fasta-file`: The name of the fasta file containing sequence(s) where each
  sequence ID will have a ``unique`` file generated. Must be equal to or a
  subset of the sequence used to generate the index used for ``index-file``.

Example:

.. code-block:: console

    $ newmap unique-lengths --kmer-lengths 20:200 hg38.awfmi chr1.fna

This will generate a "unique" binary file from the sequence with it's id (e.g.
``chr1``) with the suffix of the underyling data type (``chr1.unique.uint8``)
containing the minimum unique length found from the given range of k-mer
lengths of 20 to 200 bp long.

K-mer search ranges
^^^^^^^^^^^^^^^^^^^

The `kmer-lengths` parameter can be a comma seperated list of k-mer lengths or
a colon seperated range. A comma seperated list will be linearly searched and
is assumed to be ordered from smallest to largest. It is recommended to use
this method when only a few k-mer lengths are needed. A colon seperated range
will have `all` lengths inclusively searched for using a binary search method.
As a result the range of k-mer lengths can increase significantly with only a
logarithmic increase in compute time.

The verbose output will print statistics such as the minimum and maximum k-mer
lengths that were found to be unique from the specified range. This can be
useful as a guideline for future search ranges on other sequences.
Notably if your the largest k-mer length found is much smaller than the maximum
length and your minimum is larger than your (colon seperated) range, it
signifies that the sequence has likely been exhaustively searched.

Ambiguous bases
^^^^^^^^^^^^^^^

Due to the implementation of the AWFM-index, `all non-ACGT bases are treated as
an equivalent base
<https://almob.biomedcentral.com/articles/10.1186/s13015-021-00204-6/tables/1>`_.
Newmap takes the approach of only permitting ACGT bases and their lowercase
soft-masked equivalent conventionally introduced by software such as
`RepeatMasker <https://www.repeatmasker.org>`_. All other character codes are
treated as ambiguous bases and are excluded from the search for unique minimum
length k-mers.

Threading
^^^^^^^^^

The threading option only applies to the counting of the k-mers in the index.
It has `close to linear performance on counting up to 20
<https://almob.biomedcentral.com/articles/10.1186/s13015-021-00204-6#Sec23>`_
with some diminishing returns afterwards.


.. _generate-mappability:

--------------------
generate-mappability
--------------------
Generates mappability files from a given ``unique`` file (see
:ref:`unique-file-format`). There are two types of mappability files that can
be generated:

1. Single-read mappability (see :ref:`single-read-mappability`)
2. Multi-read mappability (see :ref:`multi-read-mappability`)

Options
-------

- `kmer-length`: The length of the k-mer to use for mappability. Defaults to 24.
- `single-read-bed-file`: The name of the BED file to write the single-read mappability to. Specify ``-`` for ``stdout``.
- `multi-read-wig-file`: The name of the WIG file to write the multi-read mappability to. Specify ``-`` for ``stdout``.
- `verbose`: Print verbose output. Default is False.

.. note::

    Only ``single-read-bed-file`` or ``multi-read-wig-file`` can output to ``stdout`` when both are specified on the command line.


Mappability datasets
^^^^^^^^^^^^^^^^^^^^
The mappability datasets are generated from the minimum unique length dataset
and defined for a given k-mer length.

.. _single-read-mappability:

Single-read mappability
^^^^^^^^^^^^^^^^^^^^^^^
Single-read mappability is a binary value (0 or 1) for each position in the
sequence where a 1 signifies that there exists for a length k, at least 1
unique k-mer that overlaps that position and 0 otherwise.

The resulting BED file from this command will place the resulting binary value
in the "score" column of the BED file.

.. _multi-read-mappability:

Multi-read mappability
^^^^^^^^^^^^^^^^^^^^^^
Multi-read mappability is a floating point value between 0 and 1 for each
position in the sequence. Each value represents the fraction of sequence
positions that have a unique k-mer length which overlap that sequence position.
For example, for a given sequence position for a k-mer length of 24, if all
24-mers that overlap that position are also unique at their respective
positions, the resulting value will be 1. If only 12 24-mers (half the amount)
are unique at their respective positions, the resulting value will be 0.5.
All values are put into a WIG file. The WIG file will have a "fixedStep" format
and may be very large.

Example:

.. code-block:: console

    $ newmap generate-mappability -k 24 -m k24_multiread_mappability.wig -s k24_singleread_mappability.bed chr1.unique.uint8
