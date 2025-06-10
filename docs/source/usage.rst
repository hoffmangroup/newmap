Usage
=====

Newmap is set of tools for quantifying mappability for any sequence. It
accomplishes this by finding unique read/k-mer lengths at every position in
your sequence. From these unique lengths, both single-read and multi-read
mappability can be generated for a specific read length.

To produce mappability tracks, the following steps are performed:

1. Create an FM-index for the sequence of interest.
2. Search and record unique read/k-mer lengths from a specified range for each
   position in your sequence. This search is done by counting the number of
   "aligned reads" to the indexed sequence.
3. Generate mappability tracks from the unique read/k-mer lengths for a
   specified read length.

These unique counts are stored in a ``.unique`` file which is a binary file
containing integers of the minimum unique length found from the search range.
The single-read mappability is stored in a BED file, while the multi-read
mappability is stored in a WIG file.

Newmap requires a CPU with AVX2 support, and has only been tested on Linux.

For instructions on how to install see the `Installation`_ section.
For a quick example on how to get started, see the `Quickstart`_ section.
For details on the command line options, see the :ref:`Commands` section.
For details on the ``unique`` file format, see the :ref:`unique-file-format` section.

.. _installation:

Installation
============

To install Newmap, you can install using pip:

.. code-block:: console

   $ pip install newmap

You can also install from Bioconda:

.. code-block:: console

   $ conda install -c bioconda newmap

To install from source (on linux), clone the repository and run:

.. code-block:: console

   $ ./build_libraries.sh
   $ pip install .

The first command will download and compile the required libraries (from a git
submodule), and the second will install the package.


.. _quickstart:

Quickstart
==========

Create an index for a reference sequence (genome.fa)
----------------------------------------------------
.. code-block:: console

   $ newmap index genome.fa

By default this creates an ``genome.awfmi`` file in the current directory
which is the index to count occurances of k-mers on.


Find unique k-mer lengths for a chromosome (chr1.fna.gz)
--------------------------------------------------------------------
.. code-block:: console

    $ newmap search --thread-count 4 --search-range=20:200 genome.awfmi chr1.fna.gz

This will create a ``.unique`` (e.g. ``chr1.unique.uint8``) binary file
containing integers of the minimum unique k-mer length found in the lengths
ranged from 20 to 200 for each position in chr1.

.. warning:: The sequence in the FASTA file **must** be the same or a subset of the sequence used to generate the index.

Convert the minimum unique k-mer lengths to mappability output files
--------------------------------------------------------------------
To create mappability data for 24-mers:

.. code-block:: console

    $ newmap search --multi-read=k24_multiread_mappability.wig --single-read=k24_singleread_mappability.bed 24 chr1.unique.uint8

The resulting BED file will be the single-read mappability for chr1, and the
WIG file will be the multi-read mappability.
