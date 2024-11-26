import unittest

from tempfile import NamedTemporaryFile
from util import TEST_DATA_PATH

from newmap._c_newmap_count_kmers import count_kmers
from newmap.main import (DEFAULT_SUFFIX_ARRAY_COMPRESSION_RATIO,
                         DEFAULT_KMER_LENGTH_IN_SEED_TABLE)
from newmap.generate_index import generate_fm_index


class TestCountKmers(unittest.TestCase):
    def setUp(self):
        self.genome_index_filename = NamedTemporaryFile(mode="w").name
        generate_fm_index(str(TEST_DATA_PATH / 'genome.fa'),
                          self.genome_index_filename,
                          DEFAULT_SUFFIX_ARRAY_COMPRESSION_RATIO,
                          DEFAULT_KMER_LENGTH_IN_SEED_TABLE)
        self.num_threads = 1

    # @unittest.skip("Test relies on large data file not in respository")
    def test_count_kmers(self):
        counts = count_kmers(self.genome_index_filename,
                             [b'AAAA', b'AT', b'TAT', b'CCC', b'NNN', b'TCGT'],
                             self.num_threads)
        self.assertEqual(counts, [9, 3, 1, 8, 0, 0])

    # @unittest.skip("Test relies on large data file not in respository")
    def test_count_wrong_type(self):
        with self.assertRaises(TypeError):
            count_kmers(self.genome_index_filename, ["AAAA"],
                        self.num_threads)

    # @unittest.skip("Test relies on large data file not in respository")
    def test_empty_byte_string(self):
        with self.assertRaises(ValueError):
            count_kmers(self.genome_index_filename, [b'AAAA', b'', b'TAT'],
                        self.num_threads)
