.. _unique-file-format:

Unique File Format
==================
The ``unique`` file format is a naive binary file format that stores at each
position the minimum length k, a k-mer sequence is unique from the sequence ID
(e.g. chr1) matched against a generated index (see :ref:`generate-index`). The
range of values for k depends on the parameters of :ref:`unique-lengths`.
Values of 0 are positions in the sequence where no minimum length was found. It
can be used to generate mappability datasets for a given k-mer size since a
k-mer that is unique to that sequence will also be unique for a "k+1"-mer. The
conversion to mappability datasets for a given k-mer length is done with the
:ref:`generate-mappability` tool.

----------
Data types
----------
The suffix of the filename specifies the type of the underyling binary data.
The suffix ``uint8`` specifies that each minimum length is represented with a
single unsigned 8 bit integer (1 byte each), and ``uint16`` likewise has each
length represented by unsigned 16 bit integer (2 bytes each). No other data is
stored in the file. The data type is chosen based on the maximum length
specified in range specified to :ref:`unique-lengths`. For example, in a search
range from 20 to 255, the maximum unique minimum length is less than 256 (which
is the maximum value that can be represented with an unsigned byte), therefore
the ``uint8`` format will be used.

-----
Usage
-----
Since the file format is simple, it can be easily processed with any language.
Below is an example of how to read a ``unique.uint8`` file in Python using
numpy.

.. code-block:: python

    import numpy as np

    minimum_unique_lengths = np.fromfile('chr1.unique.uint8', dtype=np.uint8)

    # Print fraction of unique minimum lengths found
    print(np.count_nonzero(minimum_unique_lengths) / len(minimum_unique_lengths)


