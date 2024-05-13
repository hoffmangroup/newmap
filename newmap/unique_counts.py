from argparse import ArgumentParser
from pathlib import Path
from sys import stderr

import numpy as np

from newmap._c_newmap_count_kmers import count_kmers
from newmap.util import optional_gzip_open
from newmap.fasta import sequence_segments


UNIQUE_COUNT_FILENAME_FORMAT = "{}.unique.uint8"

DEFAULT_KMER_BATCH_SIZE = 100000
DEFAULT_THREAD_COUNT = 1
DEFAULT_MINIMUM_KMER_LENGTH = 20
DEFAULT_MAXIMUM_KMER_LENGTH = 200

UMAP_KMER_LENGTHS = (24, 36, 50, 100, 150, 200)


def get_args():
    parser = ArgumentParser(
        description="Creates a binary file with minimum kmer length"
                    " for uniqueness at each sequence position")

    parser.add_argument(
        "--fasta-file", "-f",
        help="Filename of fasta file for kmer generation")

    parser.add_argument(
        "--index-file", "-i",
        help="Filename of reference index file for kmer counting.")

    parser.add_argument(
        "--reverse-index-file", "-r",
        help="Filename of reversed complement reference index file for kmer counting.")

    parser.add_argument(
        "--kmer-batch-size", "-s",
        default=DEFAULT_KMER_BATCH_SIZE,
        type=int,
        help="Maximum number of kmers to batch per reference sequence from "
             "given fasta file."
             "Use to control memory usage. "
             "Default is {}" .format(DEFAULT_KMER_BATCH_SIZE))

    parser.add_argument(
        "--umap-kmer-lengths", "-k",
        action="store_true",
        help="Use Umap kmer lengths to generate/reproduce unique counts."
             "Overrides minimum and maximum kmer lengths."
    )

    parser.add_argument(
        "--minimum-kmer-length", "-l",
        type=int,
        default=DEFAULT_MINIMUM_KMER_LENGTH,
        help="Minimum kmer length to consider for uniqueness. "
             "Default is {}" .format(DEFAULT_MINIMUM_KMER_LENGTH))

    parser.add_argument(
        "--maximum-kmer-length", "-u",
        type=int,
        default=DEFAULT_MAXIMUM_KMER_LENGTH,
        help="Minimum kmer length to consider for uniqueness. "
             "Default is {}" .format(DEFAULT_MAXIMUM_KMER_LENGTH))

    parser.add_argument(
        "--thread-count", "-t",
        default=DEFAULT_THREAD_COUNT,
        type=int,
        help="Number of threads to parallelize kmer counting. "
             "Default is {}" .format(DEFAULT_THREAD_COUNT))

    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Print additional information to standard error",)

    args = parser.parse_args()

    fasta_filename = args.fasta_file
    index_filename = args.index_file
    reverse_index_filename = args.reverse_index_file
    kmer_batch_size = args.kmer_batch_size
    umap_kmer_lengths = args.umap_kmer_lengths
    min_kmer_length = args.minimum_kmer_length
    max_kmer_length = args.maximum_kmer_length
    num_threads = args.thread_count
    verbose = args.verbose

    out_list = [fasta_filename,
                index_filename,
                reverse_index_filename,
                kmer_batch_size,
                umap_kmer_lengths,
                min_kmer_length,
                max_kmer_length,
                num_threads,
                verbose]

    return out_list


def write_unique_counts(fasta_filename: Path,
                        index_filename: Path,
                        reverse_index_filename: Path,
                        kmer_batch_size: int,
                        use_umap_kmer_lengths: bool,
                        min_kmer_length: int,
                        max_kmer_length: int,
                        num_threads: int,
                        verbose: bool = False):

    if use_umap_kmer_lengths:
        kmer_lengths = (24, 36, 50, 100, 150, 200)
        max_kmer_length = max(kmer_lengths)
    else:
        kmer_lengths = range(min_kmer_length, max_kmer_length + 1)

    # NB: We open the file in binary mode to get read-only bytes
    with optional_gzip_open(fasta_filename, "rb") as fasta_file:
        # Allow a lookahead on the last kmer to include the final dinucleotide
        # in the sequence buffer
        lookahead_length = max_kmer_length - 1

        # Get sequence lengths long enough to process the longest kmer of
        # interest
        # Requires a lookahead so that last dinucleotide has enough to
        # form a kmer of at least the maximum length
        requested_sequence_length = kmer_batch_size + lookahead_length


        current_sequence_id = b''
        # For each sequence buffer from the fasta file
        for sequence_segment in sequence_segments(fasta_file,
                                                  requested_sequence_length,
                                                  lookahead_length):

            # If we are on a new sequence
            if current_sequence_id != sequence_segment.id:
                # Truncate any existing file for the new sequence
                open(UNIQUE_COUNT_FILENAME_FORMAT.format(
                     sequence_segment.id.decode()), "wb").close()
                # Update the current working sequence id
                current_sequence_id = sequence_segment.id

                if verbose:
                    print("Writing unique counts for sequence ID: {}".format(
                          sequence_segment.id.decode()), file=stderr)

            # If this is the last sequence segment
            if sequence_segment.epilogue:
                # The number of kmers is equal to the entire length of the
                # sequence segment
                num_kmers = len(sequence_segment.data)
            # Otherwise
            else:
                # The number of kmers is equal to kmer batch size
                # (or the sequence segment length minus the lookahead)
                num_kmers = kmer_batch_size

            # Initialize the unique counts for this sequence segment to 0
            segment_unique_counts = np.zeros(num_kmers,
                                             dtype=np.uint8)

            if verbose:
                print("Processing {} kmers".format(num_kmers), file=stderr)

            # For each kmer length
            for kmer_length in kmer_lengths:
                # Create all kmers of that length from the sequence buffer
                # NB: At the epilogue of the sequence, out of bounds indexing
                # after the end of the data will automatically be truncated by
                # numpy
                if verbose:
                    print("Counting {}-mers".format(kmer_length), file=stderr)

                # Create the list of kmers for this kmer length and sequence
                kmers = [sequence_segment.data[i:i+kmer_length]
                         for i in range(num_kmers)]

                # Count the occurances of kmers on the forward strand
                count_list = np.array(count_kmers(
                                        str(index_filename),
                                        kmers,
                                        num_threads),
                                      dtype=np.uint32)

                # Count the occurances of kmers on the reverse strand
                reverse_count_list = np.array(
                                     count_kmers(
                                       str(reverse_index_filename),
                                       kmers,
                                       num_threads),
                                     dtype=np.uint32)

                if verbose:
                    print("{} unique counts before merging and ambiguity "
                          "filtering".format(
                              np.count_nonzero(
                                np.where(
                                    (count_list + reverse_count_list) == 1,
                                    1, 0))),
                          file=stderr)

                # Create ambiguity filter
                # TODO: Find a better method of handling this?
                ambigious_kmer_locations = np.array([b'N' in kmer for kmer in kmers])

                # Remove counts of kmers where there is ambiguity
                count_list = np.where(ambigious_kmer_locations,
                                      0, count_list)
                reverse_count_list = np.where(ambigious_kmer_locations,
                                              0, reverse_count_list)

                count_list = count_list + reverse_count_list

                if verbose:
                    print("{} unique counts after ambiguity filtering".format(
                      np.count_nonzero(np.where(count_list == 1, 1, 0))),
                          file=stderr)

                # If we have a unique count (count of 1), record the kmer
                # length if no previous kmer length has been recorded
                # Otherwise keep existing kmer length
                segment_unique_counts = \
                    np.where((count_list == 1) & (segment_unique_counts == 0),
                             kmer_length,  # NB: value on truth
                             segment_unique_counts).astype(
                             dtype=np.uint8, casting="unsafe")

                if verbose:
                    print("{} unique counts found after merging with previous "
                          "counts".format(
                          np.count_nonzero(
                            np.where(segment_unique_counts == kmer_length,
                                     kmer_length, 0))),
                          file=stderr)

            # Append the unique counts to a unique count file per sequence
            with open(UNIQUE_COUNT_FILENAME_FORMAT.format(
              sequence_segment.id.decode()), "ab") as unique_count_file:
                segment_unique_counts.tofile(unique_count_file)


def main():
    fasta_filename, \
     index_filename, \
     reverse_index_filename, \
     kmer_batch_size, \
     use_umap_kmer_lengths, \
     min_kmer_length, \
     max_kmer_length, \
     num_threads, \
     verbose = get_args()

    write_unique_counts(Path(fasta_filename),
                        Path(index_filename),
                        Path(reverse_index_filename),
                        kmer_batch_size,
                        use_umap_kmer_lengths,
                        min_kmer_length,
                        max_kmer_length,
                        num_threads,
                        verbose)


if __name__ == "__main__":
    main()
