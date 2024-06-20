from pathlib import Path

import numpy as np

from newmap._c_newmap_count_kmers import count_kmers
from newmap.util import optional_gzip_open, verbose_print
from newmap.fasta import sequence_segments


COMPLEMENT_TRANSLATE_TABLE = bytes.maketrans(b'ACGT', b'TGCA')
UMAP_KMER_LENGTHS = (24, 36, 50, 100, 150, 200)
UNIQUE_COUNT_FILENAME_FORMAT = "{}.unique.uint8"


# KMERQueryState = namedtuple("KMERQueryState", ???)
# TODO: Need to track state for each kmer for each position in
# the sequence segment across all kmer lengths
# Should track:
# - If the current position has already found a minimum kmer length
# - If the current position should be ignored due to an ambiguous base
#   - All kmers longer than this at this position should also be ignored by definition
# - Otherwise if no unique minimum kmer length has been found
#   - These are the only positions that should be submitted for counting
class KMERQueryState:
    # Current minimum kmer length found at this position
    # Where 0 indicates no minimum kmer length found
    min_kmer_length: int
    # If this position has been processed due to minimum length or ambiguity
    processed: bool

    def __init__(self):
        self.min_kmer_length = 0
        self.processed = False


def write_unique_counts(fasta_filename: Path,
                        index_filename: Path,
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

                verbose_print(verbose,
                              "Writing unique counts for sequence ID: "
                              "{}".format(sequence_segment.id.decode()))

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
            # segment_unique_counts = np.zeros(num_kmers,
            #                                  dtype=np.uint8)

            verbose_print(verbose, "Processing {} kmers".format(num_kmers))
            # $ newmap unique-lengths --fasta-file chr19.fna --index-file mm39-index.awfmi --umap-kmer-lengths --thread-count 8
            kmer_query_states = [KMERQueryState() for _ in range(num_kmers)]

            # For each kmer length
            for kmer_length in kmer_lengths:
                # Create all kmers of that length from the sequence buffer
                verbose_print(verbose, "Counting {}-mers".format(kmer_length))

                # Create our working list of kmers to count
                working_kmers = []
                # For each kmer in the sequence segment
                for i in range(num_kmers):
                    # If a previous kmer length was not found to be ambiguous
                    # Or a previous kmer length was not found to be unique
                    if not kmer_query_states[i].processed:
                        # Create the kmer from the sequence segment
                        # NB: At the epilogue of the sequence, out of bounds
                        # indexing after the end of the data will automatically
                        # be truncated by numpy
                        kmer = sequence_segment.data[i:i+kmer_length]
                        # If it contains an ambiguous base
                        if b'N' in kmer:  # TODO: Add option for which bases
                            # Ignore it for all longer kmer lengths (i.e. all
                            # future iterations)
                            kmer_query_states[i].processed = True
                        # Otherwise
                        else:
                            # Add it to our working kmer lengths
                            working_kmers.append(kmer)

                # If there are no kmers to count due to ambiguity
                if not working_kmers:
                    verbose_print(verbose,
                                  "No {}-mers remaining to be found, skipping "
                                  "to next kmer batch".format(kmer_length))
                    # Skip to the next kmer batch
                    break
                else:
                    verbose_print(verbose,
                                  "{} {}-mers remaining to be counted"
                                  .format(len(working_kmers), kmer_length))

                # Count the occurances of kmers on the forward strand
                count_list = np.array(count_kmers(
                                        str(index_filename),
                                        working_kmers,
                                        num_threads),
                                      dtype=np.uint32)

                # TODO: Add no reverse complement option to skip this to
                # support bisulfite treated kmer counting
                count_list += np.array(
                  count_kmers(str(index_filename),
                              [kmer.translate(COMPLEMENT_TRANSLATE_TABLE)[::-1]
                               for kmer in working_kmers],
                              num_threads), dtype=np.uint32)

                # Update the minimum kmer length found for each kmer
                count_index = 0
                min_kmer_lengths_found = 0
                # For each kmer in our current batch
                for kmer_query_state in kmer_query_states:
                    # Where our working query list is unprocessed
                    if not kmer_query_state.processed:
                        # And the count is 1
                        if count_list[count_index] == 1:
                            # Save the minimum kmer length found
                            kmer_query_state.min_kmer_length = kmer_length
                            # Mark the query state as processed
                            kmer_query_state.processed = True
                            min_kmer_lengths_found += 1
                        # Increment the count indexing
                        count_index += 1

                verbose_print(verbose, "{} unique {}-mers found".format(
                    min_kmer_lengths_found, kmer_length))

            # Append the unique counts to a unique count file per sequence
            segment_unique_counts = np.array(
                [kmer_query_state.min_kmer_length
                 for kmer_query_state in kmer_query_states],
                dtype=np.uint8)

            with open(UNIQUE_COUNT_FILENAME_FORMAT.format(
              sequence_segment.id.decode()), "ab") as unique_count_file:
                segment_unique_counts.tofile(unique_count_file)


def main(args):
    fasta_filename = args.fasta_file
    index_filename = args.index_file
    kmer_batch_size = args.kmer_batch_size
    use_umap_kmer_lengths = args.umap_kmer_lengths
    min_kmer_length = args.minimum_kmer_length
    max_kmer_length = args.maximum_kmer_length
    num_threads = args.thread_count
    verbose = args.verbose

    write_unique_counts(Path(fasta_filename),
                        Path(index_filename),
                        kmer_batch_size,
                        use_umap_kmer_lengths,
                        min_kmer_length,
                        max_kmer_length,
                        num_threads,
                        verbose)
