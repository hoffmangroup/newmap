import filecmp
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

from tests.util import TEST_DATA_PATH

from newmap._c_newmap_generate_index import generate_fm_index
from newmap.main import (DEFAULT_SUFFIX_ARRAY_COMPRESSION_RATIO,
                         DEFAULT_KMER_LENGTH_IN_SEED_TABLE)


class TestGenerateIndex(unittest.TestCase):
    def setUp(self):
        self.test_genome_filename = "test.awfmi"

    @unittest.skip("Test relies on large data file not in respository")
    def test_generate_fm_index(self):
        self.reference_sequence = str(TEST_DATA_PATH / 'genome.fa')

        with TemporaryDirectory() as temp_dir_name:
            test_genome_index_filename = \
                str(Path(temp_dir_name) / self.test_genome_filename)

            generate_fm_index(self.reference_sequence,
                              test_genome_index_filename,
                              DEFAULT_SUFFIX_ARRAY_COMPRESSION_RATIO,
                              DEFAULT_KMER_LENGTH_IN_SEED_TABLE)

            # Compare our created index with the expected index already made
            self.assertTrue(filecmp.cmp(test_genome_index_filename,
                                        str(TEST_DATA_PATH / 'genome.awfmi')))

    def test_no_fasta(self):
        self.fake_reference_sequence = str(TEST_DATA_PATH / 'genome_foo.fasta')

        with TemporaryDirectory() as temp_dir_name:

            test_genome_index_filename = \
                str(Path(temp_dir_name) / self.test_genome_filename)

            with self.assertRaises(FileNotFoundError):
                generate_fm_index(self.fake_reference_sequence,
                                  test_genome_index_filename,
                                  DEFAULT_SUFFIX_ARRAY_COMPRESSION_RATIO,
                                  DEFAULT_KMER_LENGTH_IN_SEED_TABLE)
