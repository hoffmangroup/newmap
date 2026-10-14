.. _bisulfite:

Bisulfite Mappability Example
=============================

Newmap can be used to generate bisulfite mappability tracks for a given reference genome.
Here we provide a protocol for generating non-directional, and directional including Post-Bisulfite adaptor tagging (PBAT) bisulfite mappability.
All examples shown below will be for 50 bp reads and using 16 threads as an example.
For this example, verbosity will be enabled by default for any troubleshooting purposes and should not slow down output by any significant margin.
Adjust the read lengths and threads as necessary.

-------------------------------
Steps for bisulfite mappability
-------------------------------

For all bisulfite mappability, the first step is to convert the reference two different ways:

1. Create a converted genome reference sequence where all C's are converted to T's.
2. Create a converted genome reference sequence where all G's are converted to A's.

These sequences will be referred to as ``genome.CtoT.fa`` and ``genome.GtoA.fa`` respectively.

Using :ref:`newmap index <index>` and the converted sequences:

3. Generate an index where the reference genome has all C's converted to T's (``newmap index genome.CtoT.fa``).
4. Generate an index where the reference genome has all G's converted to A's (``newmap index genome.GtoA.fa``).

These indices will be referred to as ``genome.CtoT.awfmi`` and ``genome.GtoA.awfmi`` respectively.

At this point you should have a C to T converted reference sequence and its index and the same likewise for the G to A conversion.

Directional bisulfite mappability (C to T converted reads only)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Run :ref:`newmap search <search>` using the C to T converted reference sequence on both the C to T index and G to A index.
Match only the C to T converted reads against the C to T index, excluding its reverse complement.
Match only the C to T converted reads against the G to A index, excluding its forward strand (or its reverse complement only).

5. Run :ref:`newmap search <search>`

.. code-block:: bash

    newmap search \
        --verbose \
        --search-range=50 \
        --num-threads=16 \
        --exclude-strands=norc,nofw \
        genome.CtoT.fa \
        genome.CtoT.awfmi,genome.GtoA.awfmi

Run :ref:`newmap track <track>` to generate the 50 bp mappability tracks.


PBAT bisulfite mappability (G to A converted reads only)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Identical to the directional bisulfite mappability above except only matching with the G to A converted sequence:

5. Run :ref:`newmap search <search>`

.. code-block:: bash

    newmap search \
        --verbose \
        --search-range=50 \
        --num-threads=16 \
        --exclude-strands=norc,nofw \
        genome.GtoA.fa \
        genome.GtoA.awfmi,genome.CtoT.awfmi

Run :ref:`newmap track <track>` to generate the 50 bp mappability tracks.


Non-directional Bisulfite Mappability (C to T and G to A converted reads)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
A combination of both search examples above.

- C to T sequence to the C to T index (excluding reverse complement)
- C to T sequence to the G to A index (excluding forward strand)
- G to A sequence to the C to T index (excluding forward strand)
- G to A sequence to the G to A index (excluding reverse complement)

5. Run :ref:`newmap search <search>`

.. code-block:: bash

    newmap search \
        --verbose \
        --search-range=50 \
        --num-threads=16 \
        --exclude-strands=norc,nofw,nofw,norc \
        genome.CtoT.fa,genome.GtoA.fa \
        genome.CtoT.awfmi,genome.GtoA.awfmi

Run :ref:`newmap track <track>` to generate the 50 bp mappability tracks.
