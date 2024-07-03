from math import ceil, log2
from pathlib import Path

import numpy as np
import numpy.typing as npt

from newmap._c_newmap_count_kmers import count_kmers
from newmap.util import ceil_div, optional_gzip_open, verbose_print
from newmap.fasta import sequence_segments

KMER_RANGE_SEPARATOR = ":"

COMPLEMENT_TRANSLATE_TABLE = bytes.maketrans(b'ACGT', b'TGCA')
UMAP_KMER_LENGTHS = (24, 36, 50, 100, 150, 200)
UNIQUE_COUNT_FILENAME_FORMAT = "{}.unique.uint8"


class KMERQueryState:
    # Should track:
    # - If the current position has already found a minimum kmer length
    # - If the current position should be ignored due to an ambiguous base
    #   - All kmers longer than this at this position should also be ignored by definition
    # - Otherwise if no unique minimum kmer length has been found
    #   - These are the only positions that should be submitted for counting

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
                        kmer_lengths: list[int],
                        num_threads: int,
                        use_binary_search=False,
                        verbose: bool = False):

    if (verbose and
            use_binary_search):
        max_kmer_length = max(kmer_lengths)
        min_kmer_length = min(kmer_lengths)
        verbose_print(verbose, "Max {} iterations over range {}-{}".format(
                      ceil(log2(max_kmer_length - min_kmer_length) + 1),
                      min_kmer_length, max_kmer_length))

    # NB: We open the file in binary mode to get read-only bytes
    with optional_gzip_open(fasta_filename, "rb") as fasta_file:
        # Allow a lookahead on the last kmer to include the final dinucleotide
        # in the sequence buffer
        lookahead_length = max(kmer_lengths) - 1

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
                              "Writing minimum unique lengths for sequence "
                              "ID: {}".format(sequence_segment.id.decode()))

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

            verbose_print(verbose, "Processing {} kmers".format(num_kmers))

            if use_binary_search:
                segment_unique_counts = binary_search(index_filename,
                                                      sequence_segment,
                                                      kmer_lengths,
                                                      num_kmers,
                                                      num_threads,
                                                      verbose)
            else:
                segment_unique_counts = linear_search(index_filename,
                                                      sequence_segment,
                                                      kmer_lengths,
                                                      num_kmers,
                                                      num_threads,
                                                      verbose)

            # Append the unique counts to a unique count file per sequence
            with open(UNIQUE_COUNT_FILENAME_FORMAT.format(
              sequence_segment.id.decode()), "ab") as unique_count_file:
                segment_unique_counts.tofile(unique_count_file)

        # TODO: Output summary statistics per sequence:
        # Amount of unique lengths found
        # Amount of positions skipped due to ambiguity
        # Amount of positions that had no unique length found


def get_kmer_counts(index_filename: Path,
                    kmers: list[bytes],
                    num_threads: int) -> npt.ArrayLike:  # uint32

    # Count the occurances of kmers on the forward strand
    count_list = np.array(count_kmers(
                            str(index_filename),
                            kmers,
                            num_threads),
                          dtype=np.uint32)

    # TODO: Add no reverse complement option to skip this to
    # support bisulfite treated kmer counting
    # TODO: Add option for a complement table
    count_list += np.array(
      count_kmers(str(index_filename),
                  [kmer.translate(COMPLEMENT_TRANSLATE_TABLE)[::-1]
                   for kmer in kmers],
                  num_threads), dtype=np.uint32)

    return count_list


def binary_search(index_filename,
                  sequence_segment,
                  kmer_lengths,
                  num_kmers,
                  num_threads,
                  verbose) -> npt.ArrayLike:  # uint8

    max_kmer_length = max(kmer_lengths)
    min_kmer_length = min(kmer_lengths)
    # NB: Floor division for midpoint
    # NB: Avoid an overflow error by dividing first before sum
    starting_kmer_length = (max_kmer_length // 2) + (min_kmer_length // 2)

    # Track which kmer positions have finished searching,
    # skipping any kmers starting with an ambiguous base
    # NB: Iterating over bytes returns ints
    finished_search = np.array([c == ord(b'N') for c in
                                sequence_segment.data[:num_kmers]])
    # Track current kmer query (for minimum) length
    current_length_query = np.full(num_kmers, starting_kmer_length,
                                   dtype=np.uint32)
    upper_length_bound = np.full(num_kmers, max_kmer_length, dtype=np.uint32)
    lower_length_bound = np.full(num_kmers, min_kmer_length, dtype=np.uint32)
    # previous_length_query = np.zeros(num_kmers, dtype=np.uint8)

    # List of minimum lengths (where 0 is nothing was found)
    unique_lengths = np.zeros(num_kmers, dtype=np.uint8)

    verbose_print(verbose, "Skipping {} ambigious positions".format(
                  finished_search.sum()))

    iteration_count = 1

    # While there are still kmers to search
    while not np.all(finished_search):
        # Create our working list of kmers to count
        working_kmers = []
        # Track which kmer positions need to be counted on the index
        should_count = np.copy(~finished_search)

        verbose_print(verbose, "Iteration {}".format(iteration_count))
        verbose_print(verbose, "{} k-mer positions remaining".format(
                      np.sum(should_count)))

        ambigious_overlap_count = 0
        # For each unfinished searched position
        for i in np.nonzero(~finished_search)[0]:
            # Create the kmer from the sequence segment
            # NB: At the epilogue of the sequence, out of bounds
            # indexing after the end of the data will automatically
            # be truncated by numpy
            current_kmer_length = current_length_query[i]
            kmer = sequence_segment.data[i:i+current_kmer_length]

            # If it contains an ambiguous base
            if b'N' in kmer:  # NB: TODO: Add option for which bases
                # We can only attempt to get shorter on this kmer

                # If we are already at the minimum length
                if current_kmer_length == min_kmer_length:
                    # Mark this position as finished
                    finished_search[i] = True
                # Else if we have reduced our bounds next to each other
                # we have finished searching on this position as well
                elif upper_length_bound[i] - lower_length_bound[i] <= 1:
                    # Mark this position as finished
                    finished_search[i] = True
                # Otherwise attempt to reduce the length further
                else:
                    update_query_lengths(i,
                                         current_length_query,
                                         upper_length_bound,
                                         lower_length_bound,
                                         decrease=True)

                    ambigious_overlap_count += 1
                # Don't count this kmer on the index
                should_count[i] = False
            # Otherwise
            else:
                # Add it to our kmers to count
                working_kmers.append(kmer)

        verbose_print(verbose,
                      "{} k-mer positions skipped due to ambiguity".format(
                          ambigious_overlap_count))

        # Get the occurances of the kmers on both strands
        count_list = get_kmer_counts(index_filename,
                                     working_kmers,
                                     num_threads)

        # Get the indices of the kmers that should be counted
        # NB: nonzero returns a tuple for each dimension and we only have 1
        counted_positions = should_count.nonzero()[0]

        # XXX: Assert equal lengths of counted_positions and count_list

        # Binary search for the minimum kmer length based on the counts
        # For every k-mer submitted for counting (i.e. non-ambiguous k-mer)
        for i, count in zip(counted_positions, count_list):
            current_kmer_length = current_length_query[i]

            # If we have searched at the minimum or maximum lengths
            # This kmer position is finished searching
            if (current_kmer_length == min_kmer_length or
                current_kmer_length == max_kmer_length):
                # Mark as this position as finished
                finished_search[i] = True
            # Otherwise if we have reduced our bounds next to each other
            # we have finished searching on this position as well
            elif upper_length_bound[i] - lower_length_bound[i] <= 1:
                # Mark this position as finished
                finished_search[i] = True

            # If the result is 1, we have found a unique kmer length
            if count == 1:
                # If there is no current count for this positions
                # Or if this is a smaller length than what we have recorded
                # currently
                if (unique_lengths[i] == 0 or
                   current_kmer_length < unique_lengths[i]):
                    # Record the current length as the current minimum
                    # TODO: Parameterize or fix this
                    unique_lengths[i] = \
                        current_kmer_length.astype(np.uint8, casting='safe')

            # Update the query lengths for the next iteration if this position
            # is not finished searching
            if not finished_search[i]:
                if count == 1:
                    # This position could still be shorter, so continue
                    # searching at a smaller length
                    update_query_lengths(
                        i,
                        current_length_query,
                        upper_length_bound,
                        lower_length_bound,
                        decrease=True)
                # Otherwise (if the number of occurances is > 0)
                else:
                    # Try to find a longer length that is unique
                    update_query_lengths(
                        i,
                        current_length_query,
                        upper_length_bound,
                        lower_length_bound,
                        decrease=False)

        iteration_count += 1

    return unique_lengths


def update_query_lengths(index: int,
                         current_length_query_list: npt.NDArray[np.uint32],
                         upper_length_bound_list: npt.NDArray[np.uint32],
                         lower_length_bound_list: npt.NDArray[np.uint32],
                         decrease=True):
    # NB: There is a theoretical chance when calculating the midpoint, we might
    # overflow whatever integer type we are using (32 bit) if we sum then
    # divide, however ceiling/floor division divide then sum is not equivalent
    # We are assuming that lengths will be less than half the max of a u32

    current_length = current_length_query_list[index]
    new_length = current_length

    # If we are decreasing the length
    if decrease:
        # Set the new upper length bound to the current length
        upper_length_bound_list[index] = current_length
        # Calculate the new current length as the floor mid point between the
        # current and lower length bound
        lower_length_bound = lower_length_bound_list[index]
        new_length = (current_length + lower_length_bound) // 2
    # Otherwise we are increasing the length
    else:
        # Set the new lower length bound to the current length
        lower_length_bound_list[index] = current_length
        # Calculate the new current length as the ceil mid point between the
        # current and upper length bound
        upper_length_bound = upper_length_bound_list[index]
        new_length = ceil_div(current_length + upper_length_bound, 2)

    # Save the newly calculated length
    current_length_query_list[index] = new_length


def linear_search(index_filename,
                  sequence_segment,
                  kmer_lengths,
                  num_kmers,
                  num_threads,
                  verbose) -> npt.NDArray[np.uint8]:
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

        # Count the occurances of the kmers on both strands
        count_list = get_kmer_counts(index_filename,
                                     working_kmers,
                                     num_threads)

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

    segment_unique_counts = np.array(
        [kmer_query_state.min_kmer_length
         for kmer_query_state in kmer_query_states],
        dtype=np.uint8)

    return segment_unique_counts



def main(args):
    fasta_filename = args.fasta_file
    index_filename = args.index_file
    kmer_batch_size = args.kmer_batch_size
    kmer_lengths_arg = args.kmer_lengths
    num_threads = args.thread_count
    verbose = args.verbose

    # Parse the kmer lengths
    # NB: Either comma seperated or a range seperated by a colon
    kmer_lengths = []
    use_binary_search = False
    # If a colon is present in the kmer lengths argument
    if KMER_RANGE_SEPARATOR in kmer_lengths_arg:
        try:
            min_kmer_length, max_kmer_length = map(
                int, kmer_lengths_arg.split(KMER_RANGE_SEPARATOR))
        except ValueError:
            raise ValueError("Only one colon allowed for length range")

        if min_kmer_length > max_kmer_length:
            raise ValueError("K-mer range start length is larger than the end "
                             "length")
        kmer_lengths = list(range(min_kmer_length, max_kmer_length + 1))
        use_binary_search = True
    # Otherwise assume a single digit or a list of digits seperated by a comma
    # NB: Assume in sorted order?
    else:
        kmer_lengths = list(map(int, kmer_lengths_arg.split(",")))

    write_unique_counts(Path(fasta_filename),
                        Path(index_filename),
                        kmer_batch_size,
                        kmer_lengths,
                        num_threads,
                        use_binary_search,
                        verbose)
