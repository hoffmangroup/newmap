import unittest

from tempfile import NamedTemporaryFile
from util import TEST_DATA_PATH

from newmap._c_newmap_count_kmers import count_kmers, count_kmers_from_sequence
from newmap.main import (DEFAULT_COMPRESSION_RATIO,
                         DEFAULT_SEED_LENGTH)
from newmap.index import generate_fm_index


class TestCountKmers(unittest.TestCase):
    def setUp(self):
        self.genome_index_filename = NamedTemporaryFile(mode="w").name
        generate_fm_index(str(TEST_DATA_PATH / 'genome.fa'),
                          self.genome_index_filename,
                          DEFAULT_COMPRESSION_RATIO,
                          DEFAULT_SEED_LENGTH)
        self.num_threads = 1

    def test_count_kmers(self):
        counts = count_kmers(self.genome_index_filename,
                             [b'AAAA', b'AT', b'TAT', b'CCC', b'NNN', b'TCGT'],
                             self.num_threads)
        self.assertEqual(counts, [9, 3, 1, 8, 0, 0])

    def test_count_wrong_type(self):
        with self.assertRaises(TypeError):
            count_kmers(self.genome_index_filename, ["AAAA"],
                        self.num_threads)

    def test_empty_byte_string(self):
        with self.assertRaises(ValueError):
            count_kmers(self.genome_index_filename, [b'AAAA', b'', b'TAT'],
                        self.num_threads)

    def test_sequence_counting(self):
        counts = count_kmers_from_sequence(
                    self.genome_index_filename,
                    b'AAAAATTTTTATCGAATCGA',
                    [0, 4, 9],
                    [4, 2, 3],
                    self.num_threads)
        self.assertEqual(counts, [9, 3, 1])
