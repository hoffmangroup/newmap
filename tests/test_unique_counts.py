import unittest

from pathlib import Path
from tempfile import NamedTemporaryFile
from util import TEST_DATA_PATH

import numpy as np

from newmap.main import (DEFAULT_SUFFIX_ARRAY_COMPRESSION_RATIO,
                         DEFAULT_KMER_LENGTH_IN_SEED_TABLE)
from newmap.generate_index import generate_fm_index
from newmap.unique_counts import write_unique_counts

# Expected minimum unique lengths at each position
# NB: In order to manually count correctly, it is important to remember to
# count the reverse complement of that k-mer as well.
EXPECTED_CHR1_COUNTS = [0, 10, 9, 8, 7, 6, 5, 4, 4, 4, 6, 5, 4, 4,
                        4, 0, 0, 0, 0, 0]
EXPECTED_CHR2_COUNTS = [10, 10, 9, 8, 7, 6, 5, 4, 4, 4, 0, 0, 0, 0,
                        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                        0, 0]


class TestCountKmers(unittest.TestCase):
    genome_index_file = NamedTemporaryFile(mode="w")
    genome_index_filename = genome_index_file.name
    fasta_filename = str(TEST_DATA_PATH / 'genome.fa')
    num_threads = 1

    @classmethod
    def setUpClass(cls):
        generate_fm_index(cls.fasta_filename,
                          cls.genome_index_filename,
                          DEFAULT_SUFFIX_ARRAY_COMPRESSION_RATIO,
                          DEFAULT_KMER_LENGTH_IN_SEED_TABLE)

    @classmethod
    def tearDownClass(cls):
        cls.genome_index_file.close()

    def test_binary_search(self):
        self.search(use_binary_search=True)

    def test_linear_search(self):
        self.search(use_binary_search=False)

    def search(self, use_binary_search):
        write_unique_counts(Path(self.fasta_filename),
                            Path(self.genome_index_filename),
                            15,  # Batch size
                            list(range(4, 11)),  # Kmer lengths 4 to 10
                            0,  # Initial search length
                            [],  # Include chr ids
                            [],  # Exclude chr ids
                            self.num_threads,
                            use_binary_search)

        # Check the results in chr1.unique.uint8 and chr2.unique.uint8
        chr1_results = np.fromfile('chr1.unique.uint8', dtype=np.uint8)
        chr2_results = np.fromfile('chr2.unique.uint8', dtype=np.uint8)

        self.assertTrue(chr1_results.size == 20)
        self.assertTrue(chr2_results.size == 30)

        self.assertTrue(np.array_equal(chr1_results, EXPECTED_CHR1_COUNTS))
        self.assertTrue(np.array_equal(chr2_results, EXPECTED_CHR2_COUNTS))


if __name__ == '__main__':
    unittest.main()
