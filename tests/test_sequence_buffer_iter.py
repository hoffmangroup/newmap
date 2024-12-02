from io import BytesIO
import unittest

from util import TEST_DATA_PATH
from newmap.fasta import sequence_segments


class TestSequenceBufferIterator(unittest.TestCase):

    def setUp(self) -> None:
        self.chr2_fasta_file = open(TEST_DATA_PATH / 'chr2.fa', 'rb')
        self.genome_fasta_file = open(TEST_DATA_PATH / 'genome.fa', 'rb')

    def tearDown(self) -> None:
        self.chr2_fasta_file.close()
        self.genome_fasta_file.close()

    def test_entire_sequence(self):
        buffer_iter = sequence_segments(self.chr2_fasta_file, 1000)
        sequence_buffer = next(buffer_iter)
        self.assertEqual(b'chr2', sequence_buffer.id)
        self.assertEqual(b'CGCANCAGAGCANCGNCG', sequence_buffer.data)
        self.assertTrue(sequence_buffer.epilogue)
        self.assertRaises(StopIteration, next, buffer_iter)

    def test_sequence_overlap(self):
        buffer_iter = sequence_segments(self.chr2_fasta_file, 5, 2)

        sequence_buffer = next(buffer_iter)
        self.assertEqual(b'chr2', sequence_buffer.id)
        self.assertEqual(b'CGCAN', sequence_buffer.data)
        self.assertFalse(sequence_buffer.epilogue)

        sequence_buffer = next(buffer_iter)
        self.assertEqual(b'chr2', sequence_buffer.id)
        self.assertEqual(b'ANCAG', sequence_buffer.data)
        self.assertFalse(sequence_buffer.epilogue)

        sequence_buffer = next(buffer_iter)
        self.assertEqual(b'chr2', sequence_buffer.id)
        self.assertEqual(b'AGAGC', sequence_buffer.data)
        self.assertFalse(sequence_buffer.epilogue)

        sequence_buffer = next(buffer_iter)
        self.assertEqual(b'chr2', sequence_buffer.id)
        self.assertEqual(b'GCANC', sequence_buffer.data)
        self.assertFalse(sequence_buffer.epilogue)

        sequence_buffer = next(buffer_iter)
        self.assertEqual(b'chr2', sequence_buffer.id)
        self.assertEqual(b'NCGNC', sequence_buffer.data)
        self.assertFalse(sequence_buffer.epilogue)

        sequence_buffer = next(buffer_iter)
        self.assertEqual(b'chr2', sequence_buffer.id)
        self.assertEqual(b'NCG', sequence_buffer.data)
        self.assertTrue(sequence_buffer.epilogue)

        self.assertRaises(StopIteration, next, buffer_iter)

    def test_long_sequence_length(self):
        buffer_iter = sequence_segments(self.chr2_fasta_file, 100, 4)
        sequence_buffer = next(buffer_iter)

        self.assertTrue(sequence_buffer.epilogue)
        self.assertEqual(len(sequence_buffer.data), 18)

    def test_long_lookahead_length(self):
        buffer_iter = sequence_segments(self.chr2_fasta_file, 16, 10)

        sequence_buffer = next(buffer_iter)
        self.assertFalse(sequence_buffer.epilogue)
        self.assertEqual(len(sequence_buffer.data), 16)
        self.assertEqual(b'CGCANCAGAGCANCGN', sequence_buffer.data)

        sequence_buffer = next(buffer_iter)
        self.assertTrue(sequence_buffer.epilogue)
        self.assertEqual(len(sequence_buffer.data), 12)
        self.assertEqual(b'AGAGCANCGNCG', sequence_buffer.data)

    def test_multiline_sequence(self):
        buffer_iter = sequence_segments(self.chr2_fasta_file, 12, 4)
        sequence_buffer = next(buffer_iter)
        self.assertEqual(b'chr2', sequence_buffer.id)
        self.assertEqual(b'CGCANCAGAGCA', sequence_buffer.data)

        sequence_buffer = next(buffer_iter)
        self.assertEqual(b'chr2', sequence_buffer.id)
        self.assertEqual(b'AGCANCGNCG', sequence_buffer.data)

    def test_multiple_sequence(self):
        buffer_iter = sequence_segments(self.genome_fasta_file, 25, 2)
        sequence_buffer = next(buffer_iter)

        self.assertEqual(b'chr1', sequence_buffer.id)
        self.assertEqual(b'AAAAATTTTTATCGAATCGA', sequence_buffer.data)

        sequence_buffer = next(buffer_iter)
        self.assertEqual(b'chr2', sequence_buffer.id)
        self.assertEqual(b'AAAAAAAAAAGGGGGGGGGGCCCCC', sequence_buffer.data)

        sequence_buffer = next(buffer_iter)
        self.assertEqual(b'chr2', sequence_buffer.id)
        self.assertEqual(b'CCCCCCC', sequence_buffer.data)

        self.assertRaises(StopIteration, next, buffer_iter)

    def test_multiple_sequences_multiline(self):
        buffer_iter = sequence_segments(self.genome_fasta_file, 11, 2)

        sequence_buffer = next(buffer_iter)
        self.assertEqual(b'chr1', sequence_buffer.id)
        self.assertEqual(b'AAAAATTTTTA', sequence_buffer.data)

        sequence_buffer = next(buffer_iter)
        self.assertEqual(b'chr1', sequence_buffer.id)
        self.assertEqual(b'TATCGAATCGA', sequence_buffer.data)

        sequence_buffer = next(buffer_iter)
        self.assertEqual(b'chr2', sequence_buffer.id)
        self.assertEqual(b'AAAAAAAAAAG', sequence_buffer.data)

        sequence_buffer = next(buffer_iter)
        self.assertEqual(b'chr2', sequence_buffer.id)
        self.assertEqual(b'AGGGGGGGGGG', sequence_buffer.data)

        sequence_buffer = next(buffer_iter)
        self.assertEqual(b'chr2', sequence_buffer.id)
        self.assertEqual(b'GGCCCCCCCCC', sequence_buffer.data)

        sequence_buffer = next(buffer_iter)
        self.assertEqual(b'chr2', sequence_buffer.id)
        self.assertEqual(b'CCC', sequence_buffer.data)

        self.assertRaises(StopIteration, next, buffer_iter)

    def test_exact_alignment_epilogue(self):
        fasta_file = BytesIO(b'>chr1\nAAAAATTTTTATCGAATCGA\n')
        buffer_iter = sequence_segments(fasta_file, 11, 2)

        sequence_buffer = next(buffer_iter)
        self.assertEqual(b'chr1', sequence_buffer.id)
        self.assertEqual(b'AAAAATTTTTA', sequence_buffer.data)

        sequence_buffer = next(buffer_iter)
        self.assertEqual(b'chr1', sequence_buffer.id)
        self.assertEqual(b'TATCGAATCGA', sequence_buffer.data)

        self.assertRaises(StopIteration, next, buffer_iter)
        fasta_file.close()

    def test_single_nucleotide_epilogue(self):
        fasta_file = BytesIO(b'>chr1\nATCGATCGA\n')
        buffer_iter = sequence_segments(fasta_file, 5, 2)

        sequence_buffer = next(buffer_iter)
        self.assertEqual(b'chr1', sequence_buffer.id)
        self.assertEqual(b'ATCGA', sequence_buffer.data)

        sequence_buffer = next(buffer_iter)
        self.assertEqual(b'chr1', sequence_buffer.id)
        self.assertEqual(b'GATCG', sequence_buffer.data)

        sequence_buffer = next(buffer_iter)
        self.assertEqual(b'chr1', sequence_buffer.id)
        self.assertEqual(b'CGA', sequence_buffer.data)

        self.assertRaises(StopIteration, next, buffer_iter)
        fasta_file.close()

    def test_sizes_and_epilogues(self):
        buffer_iter = sequence_segments(self.genome_fasta_file, 20, 0)

        sequence_buffer = next(buffer_iter)
        self.assertEqual(b'chr1', sequence_buffer.id)
        self.assertEqual(len(sequence_buffer.data), 20)
        self.assertTrue(sequence_buffer.epilogue)

        sequence_buffer = next(buffer_iter)
        self.assertEqual(b'chr2', sequence_buffer.id)
        self.assertEqual(len(sequence_buffer.data), 20)
        self.assertFalse(sequence_buffer.epilogue)

        sequence_buffer = next(buffer_iter)
        self.assertEqual(b'chr2', sequence_buffer.id)
        self.assertEqual(len(sequence_buffer.data), 10)
        self.assertTrue(sequence_buffer.epilogue)

        self.assertRaises(StopIteration, next, buffer_iter)


if __name__ == "__main__":
    unittest.main()
