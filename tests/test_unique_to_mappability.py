import unittest

from pathlib import Path
from tempfile import NamedTemporaryFile
from util import TEST_DATA_PATH

from newmap.main import (DEFAULT_COMPRESSION_RATIO,
                         DEFAULT_SEED_LENGTH)
from newmap.index import generate_fm_index
from newmap.search import SearchConfig, write_unique_counts
from newmap.track import write_mappability_files


class TestCountKmers(unittest.TestCase):
    genome_index_file = NamedTemporaryFile(mode="w")
    genome_index_filename = genome_index_file.name
    fasta_filename = str(TEST_DATA_PATH / 'genome.fa')
    num_threads = 1

    @classmethod
    def setUpClass(cls):
        generate_fm_index(cls.fasta_filename,
                          cls.genome_index_filename,
                          DEFAULT_COMPRESSION_RATIO,
                          DEFAULT_SEED_LENGTH)

        write_unique_counts(SearchConfig(
            fasta_filepaths=[Path(cls.fasta_filename)],
            fmindex_filepaths=[Path(cls.genome_index_filename)],
            kmer_lengths=list(range(4, 11)),
            kmer_batch_size=15,
            is_binary_search=True,
            num_threads=cls.num_threads,
        ))

    @classmethod
    def tearDownClass(cls):
        cls.genome_index_file.close()

    def test_mappability_output(self):
        unique_count_filenames = [Path('chr1.unique.uint8'),
                                  Path('chr2.unique.uint8')]
        kmer_mappability_length = 10

        write_mappability_files(unique_count_filenames,
                                kmer_mappability_length,
                                "genome.10.bed",
                                "genome.10.wig",
                                verbose=False)


if __name__ == '__main__':
    unittest.main()
