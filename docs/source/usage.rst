Usage
=====

Newmap allows you to create a FM-index for your sequence of interest, then
generate k-mers from the same sequence to count their occurances on the index.
Specifically it attempts to find a minimum unique count 1 for ranges of k.
These unique counts are stored in a ``.unique`` file which is a binary file
containing integers of the minimum unique length found. These files can be used
to generate mappability data for a give k-mer length.

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

   $ newmap generate-index genome.fa

By default this creates an ``index.awfmi`` file in the current directory which
is the index to count occurances of k-mers on.


Find the minimum unique k-mer lengths for a chromosome (chr1.fna.gz)
--------------------------------------------------------------------
.. code-block:: console

    $ newmap unique-lengths --kmer-batch-size 20000000 --thread-count 4 20:200 index.awfmi chr1.fna.gz

This will create a ``.unique`` (e.g. ``chr1.unique.uint8``) binary file
containing integers of the minimum unique k-mer length found in the lengths
ranged from 20 to 200 for each position in chr1.

.. warning:: The sequence in the FASTA file **must** be the same or a subset of the sequence used to generate the index.

Convert the minimum unique k-mer lengths to mappability output files
--------------------------------------------------------------------
To create mappability data for 24-mers:

.. code-block:: console

    $ newmap generate-mappability -k 24 -m k24_multiread_mappability.wig -s k24_singleread_mappability.bed chr1.unique.uint8

The resulting BED file will be the single-read mappability for chr1, and the
WIG file will be the multi-read mappability.
