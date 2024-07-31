import unittest

from util import TEST_DATA_PATH

from newmap._c_newmap_count_kmers import count_kmers


class TestCountKmers(unittest.TestCase):
    def setUp(self):
        self.genome_index_filename = str(TEST_DATA_PATH / 'genome.awfmi')
        self.num_threads = 1

    @unittest.skip("Test relies on large data file not in respository")
    def test_count_kmers(self):
        counts = count_kmers(self.genome_index_filename,
                             [b'AAAA', b'AT', b'TAT', b'CCC', b'NNN', b'TCGT'],
                             self.num_threads)
        self.assertEqual(counts, [9, 3, 1, 8, 0, 0])

    @unittest.skip("Test relies on large data file not in respository")
    def test_count_wrong_type(self):
        with self.assertRaises(TypeError):
            count_kmers(self.genome_index_filename, ["AAAA"],
                        self.num_threads)

    @unittest.skip("Test relies on large data file not in respository")
    def test_empty_byte_string(self):
        with self.assertRaises(ValueError):
            count_kmers(self.genome_index_filename, [b'AAAA', b'', b'TAT'],
                        self.num_threads)
