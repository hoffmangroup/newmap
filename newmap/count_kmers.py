from argparse import ArgumentParser
from pathlib import Path
from typing import BinaryIO, Generator
from sys import stdout, stderr

from newmap._c_newmap_count_kmers import count_kmers
from newmap.util import optional_gzip_open


DEFAULT_KMER_SIZE = 100
DEFAULT_KMER_BATCH_SIZE = 100000
DEFAULT_THREAD_COUNT = 1
FASTA_FILE_IGNORE_DELIMITERS = (b'>', b';')


class KmerList:
    def __init__(self, sequence_name: bytes):
        self.kmers = []
        self.sequence_name = sequence_name


# TODO: Add required arguments for fasta and index
def get_args():
    parser = ArgumentParser(
        description="Prints all unique kmers per line from all sequences in a"
        " given FASTA file to standard output")

    parser.add_argument(
        "--kmer-length", "-k",
        default=DEFAULT_KMER_SIZE,
        type=int,
        help="Kmer length. Default is {}".format(
            DEFAULT_KMER_SIZE))

    parser.add_argument(
        "--fasta-file", "-f",
        help="Filename of fasta file for kmer generation")

    parser.add_argument(
        "--index-file", "-i",
        help="Filename of reference index file for kmer counting.")

    parser.add_argument(
        "--kmer-batch-size", "-s",
        default=DEFAULT_KMER_BATCH_SIZE,
        type=int,
        help="Maximum number of kmers to batch per reference sequence from "
             "given fasta file."
             "Use to control memory usage. "
             "Default is {}" .format(DEFAULT_KMER_SIZE))

    parser.add_argument(
        "--thread-count", "-t",
        default=DEFAULT_THREAD_COUNT,
        type=int,
        help="Number of threads to parallelize kmer counting. "
             "Default is {}" .format(DEFAULT_KMER_SIZE))

    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Print additional information to standard error",)

    args = parser.parse_args()

    kmer_length = args.kmer_length
    fasta_filename = args.fasta_file
    index_filename = args.index_file
    kmer_batch_size = args.kmer_batch_size
    num_threads = args.thread_count
    verbose = args.verbose

    out_list = [kmer_length,
                fasta_filename,
                index_filename,
                kmer_batch_size,
                num_threads,
                verbose]

    return out_list


# TODO: Add option for sequence selection
def generate_sequence_kmers(
        fasta_file: BinaryIO,
        kmer_length: int,
        kmer_batch_size: int) -> Generator[KmerList, None, None]:

    """
    Yields a KmerList for each (optionally chosen) sequence in the fasta file
    """

    # Create an empty set
    current_reference_sequence_name = b''  # NB: dummy unused sequence name
    kmer_list = KmerList(current_reference_sequence_name)

    # Create a line-based buffer for the current set of nucleotides
    nucleotide_buffer = b''
    kmer_current_index = 0  # NB: index into the nucleotide_buffer

    # For each line in the fasta
    for line in fasta_file:
        # If we are on a new sequence
        # NB: Assume that either of the delimiters are indicators of a
        # new sequence, notably including comments
        if line.startswith(FASTA_FILE_IGNORE_DELIMITERS):
            # Yield the current kmer set if there are elements
            if len(kmer_list.kmers) > 0:
                yield kmer_list

            # Create a new kmer list
            # Get the new reference sequence name
            current_reference_sequence_name = \
                line.split()[0][1:]  # NB: remove leading '>'
            kmer_list = KmerList(current_reference_sequence_name)

            # Empty the nucleotide buffer
            nucleotide_buffer = b''
            # Reset the current window index
            kmer_current_index = 0
        # Otherwise if it does not start with '>' or ';'
        else:
            # Slide a kmer window across the current nucleotide list
            nucleotide_buffer += line.rstrip()

            # Record the current kmer window indices
            kmer_start_index = kmer_current_index
            kmer_end_index = kmer_current_index + kmer_length

            # While there are enough nucleotides in the buffer
            while kmer_end_index <= len(nucleotide_buffer):
                # Add the current sliding window kmer to the set
                kmer_list.kmers.append(
                    nucleotide_buffer[kmer_current_index:kmer_end_index])

                # If we hit our batch size, yield the kmer list
                if len(kmer_list.kmers) >= kmer_batch_size:
                    yield kmer_list
                    # Create a new kmer list
                    # Preserve the working reference sequence name
                    kmer_list = KmerList(current_reference_sequence_name)

                # Advance the current window slice indexes by one
                kmer_current_index += 1
                kmer_end_index += 1

            # After removing all possible kmers from the buffer
            # Calculate amount the sliding window moved by
            kmer_window_shift = kmer_current_index - kmer_start_index

            # Truncate the beginning of nucleotide buffer by the shift amount
            nucleotide_buffer = nucleotide_buffer[kmer_window_shift:]
            # And move back the current window index by the shift amount
            kmer_current_index -= kmer_window_shift

    yield kmer_list

def write_counts(kmer_length: int,
                 fasta_filename: Path,
                 index_filename: Path,
                 kmer_batch_size: int,
                 num_threads: int,
                 verbose: bool = False):

    # NB: We open the file in binary mode to get read-only bytes
    with optional_gzip_open(fasta_filename, "rb") as fasta_file:
        # NB: The fasta_file type is BinaryIO but mypy can't infer
        # from the "rb" parameter to optional_gzip_open
        # For every kmer list per sequence in the fasta file

        # TODO: Fix epilogue at the end of sequence to count kmers < kmer_length
        # e.g. prepend 'N's to the remaining kmers per sequence

        # TODO: Time difference between sequence generation and counting
        # I suspect sequence generation is the bottleneck
        for kmer_list in generate_sequence_kmers(fasta_file,  # type: ignore
                                                 kmer_length,
                                                 kmer_batch_size):

            if verbose:
                stderr.buffer.write(
                    b"Counting occurances of " +
                    bytes(str(len(kmer_list.kmers)), encoding="ascii") +
                    b" kmers for sequence " + kmer_list.sequence_name + b"\n")
                stderr.flush()

            # Count the list of kmers
            count_list = count_kmers(str(index_filename),
                                     kmer_list.kmers, num_threads)

            # Print the count and kmer sequence to standard output
            for i, (count, kmer) in enumerate(zip(
              count_list, kmer_list.kmers)):
                if verbose:
                    byte_string = bytes("{}\t{}\t{}\n".format(
                        i, kmer.decode(), count), encoding='utf-8')
                    stderr.buffer.write(byte_string)
                else:
                    stdout.buffer.write(
                        "{}\n".format(count).encode())


def main():
    kmer_length, \
     fasta_filename, \
     index_filename, \
     kmer_batch_size, \
     num_threads, \
     verbose = get_args()

    write_counts(kmer_length,
                 Path(fasta_filename),
                 Path(index_filename),
                 kmer_batch_size,
                 num_threads,
                 verbose)


if __name__ == "__main__":
    main()
