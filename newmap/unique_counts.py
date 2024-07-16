from math import ceil, log2
from pathlib import Path

import numpy as np
import numpy.typing as npt

from newmap._c_newmap_count_kmers import count_kmers
from newmap.util import optional_gzip_open, verbose_print
from newmap.fasta import sequence_segments

KMER_RANGE_SEPARATOR = ":"

COMPLEMENT_TRANSLATE_TABLE = bytes.maketrans(b'ACGT', b'TGCA')
UNIQUE_COUNT_FILENAME_FORMAT = "{}.unique.uint8"


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

        # Keep track of ambigious positions in the sequence for statistics
        total_ambiguous_positions = 0
        total_unique_lengths_count = 0
        total_no_unique_lengths_count = 0

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

                print_summary_statisitcs(verbose,
                                         total_unique_lengths_count,
                                         total_ambiguous_positions,
                                         total_no_unique_lengths_count)

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

            verbose_print(verbose, "Processing {} k-mers".format(num_kmers))

            if use_binary_search:
                segment_unique_counts, ambiguous_count = \
                    binary_search(index_filename,
                                  sequence_segment,
                                  kmer_lengths,
                                  num_kmers,
                                  num_threads,
                                  verbose)
            else:
                segment_unique_counts, ambiguous_count = \
                    linear_search(index_filename,
                                  sequence_segment,
                                  kmer_lengths,
                                  num_kmers,
                                  num_threads,
                                  verbose)

            # Update summary statistics
            unique_lengths_count = np.count_nonzero(segment_unique_counts)
            total_ambiguous_positions += ambiguous_count
            total_unique_lengths_count += unique_lengths_count
            total_no_unique_lengths_count += (num_kmers -
                unique_lengths_count - ambiguous_count)

            # Append the unique counts to a unique count file per sequence
            with open(UNIQUE_COUNT_FILENAME_FORMAT.format(
              sequence_segment.id.decode()), "ab") as unique_count_file:
                segment_unique_counts.tofile(unique_count_file)

        print_summary_statisitcs(verbose,
                                 total_unique_lengths_count,
                                 total_ambiguous_positions,
                                 total_no_unique_lengths_count)


def print_summary_statisitcs(verbose,
                             total_unique_lengths_count,
                             total_ambiguous_positions,
                             total_no_unique_lengths_count):
    if (verbose and
            total_unique_lengths_count):
        verbose_print(verbose,
                      f"{total_unique_lengths_count} unique lengths found")
        verbose_print(verbose,
                      f"{total_ambiguous_positions} positions skipped due to "
                      "ambiguity")
        verbose_print(verbose,
                      f"{total_no_unique_lengths_count} positions with no "
                      "unique length found")


def get_kmer_counts(index_filename: Path,
                    kmers: list[bytes],
                    num_threads: int) -> npt.NDArray[np.uint32]:

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


# TODO: Add types
def binary_search(index_filename,
                  sequence_segment,
                  kmer_lengths,
                  num_kmers,
                  num_threads,
                  verbose) -> tuple[npt.NDArray[np.uint8], int]:

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

    ambiguous_positions_skipped = finished_search.sum()
    # TODO: Switch to f-strings
    verbose_print(verbose, f"Skipping {ambiguous_positions_skipped} ambiguous "
                  "positions")

    # Track current kmer query (for minimum) length
    current_length_query = np.full(num_kmers, starting_kmer_length,
                                   dtype=np.uint32)

    # Inclusive bounds for remaining possible lengths
    lower_length_bound = np.full(num_kmers, min_kmer_length, dtype=np.uint32)

    # The upper search length is bounded by the minimum of:
    # The maximum length in our search query range
    # Or the maximum length that does not overlap with an ambiguous base
    upper_length_bound = np.full(num_kmers, max_kmer_length, dtype=np.uint32)
    upper_bound_change_count = 0
    short_kmers_discarded_count = 0

    # NB: The following is effecitvely O(n^2) where the max k-mer length and
    # sequence length are similar in size/magnitude
    # There might be a better way based on finding the ambiguous bases in the
    # sequence buffer, and then setting the max lengths of the previous
    # positions based on their location up to the max k away

    # For every non-ambiguous starting position
    for i in np.nonzero(~finished_search)[0]:
        # Get what would the maximum length kmer for this position
        max_length_kmer = sequence_segment.data[i:i+max_kmer_length]
        # Search for the first occurance of a base that is ambiguous
        # NB: The index is 0-based, so the index is equal to the length
        # that excludes its own position
        # NB: This is very slow for large sequences
        maximum_non_ambiguous_length = max_length_kmer.find(b'N')
        # If we found an index where an ambiguous base is
        if maximum_non_ambiguous_length != -1:
            # If the found length is longer than the minimum kmer length
            if maximum_non_ambiguous_length >= min_kmer_length:
                # Set the maximum length (to the index of the ambiguous base)
                upper_length_bound[i] = maximum_non_ambiguous_length
                # Recalculate the current query length (as a midpoint)
                current_length_query[i] = np.floor_divide(
                    upper_length_bound[i] + lower_length_bound[i], 2)
                upper_bound_change_count += 1
            # Otherwise
            else:
                # Cannot find a length that would be smaller than the minimum
                # Mark this position as finished
                short_kmers_discarded_count += 1
                finished_search[i] = True

    upper_bound_change_count = np.count_nonzero(
        upper_length_bound[(~finished_search).nonzero()] < max_kmer_length)

    if (verbose and
       upper_bound_change_count):
        verbose_print(verbose, f"{upper_bound_change_count} k-mer search "
                               "ranges truncated due to ambiguity")
    if (verbose and
       short_kmers_discarded_count):
        verbose_print(verbose, f"{short_kmers_discarded_count} k-mers shorter "
                      "than the minimum length discarded due to ambiguity")

    # List of minimum lengths (where 0 is nothing was found)
    unique_lengths = np.zeros(num_kmers, dtype=np.uint32)

    iteration_count = 1

    # While there are still kmers to search
    while not np.all(finished_search):
        # Create our working list of kmers to count
        working_kmers = []

        verbose_print(verbose, "Iteration {}".format(iteration_count))
        verbose_print(verbose, "{} k-mer positions remaining".format(
                      np.sum(~finished_search)))

        # Track which kmer positions need to be counted on the index
        # Create a list of indices where each index refers to the corresponding
        # position in the given sequence segment
        counted_positions = np.nonzero(np.copy(~finished_search))[0]

        # Create a list of kmers to count on the index
        for i in counted_positions:
            current_kmer_length = current_length_query[i]
            kmer = sequence_segment.data[i:i+current_kmer_length]
            working_kmers.append(kmer)

        # Get the occurances of the kmers on both strands
        count_list = get_kmer_counts(index_filename,
                                     working_kmers,
                                     num_threads)

        # TODO: Assert for 0 counts since there must be a mismatch between
        # sequence and index

        # Assert that the number of indices to count and the number of counts
        # are equal
        assert counted_positions.size == count_list.size, \
            "Number of counted positions ({}) and number of counts ({}) " \
            "do not match".format(len(counted_positions), len(count_list))

        # Where we have counts of 1
        unique_lengths[counted_positions] = np.where(
            (count_list == 1) &
            # And if there is no current unique length recorded
            ((unique_lengths[counted_positions] == 0) |
            # Or there is a smaller length found than the current unique length
             (current_length_query[counted_positions] <
              unique_lengths[counted_positions])),
            # Record the minimum kmer length found if it less than the current
            current_length_query[counted_positions],
            # Otherwise keep the current unique length
            unique_lengths[counted_positions])

        # Update the query length and bounds for the next iteration

        # Lower the upper bounds of our search range on positions where
        # We need to decrease our k-mer length (i.e. counts == 1)
        # Set the new upper (inclusive) bound to the current query length - 1
        upper_length_bound[counted_positions] = np.where(
            count_list == 1,
            current_length_query[counted_positions] - 1,
            upper_length_bound[counted_positions])

        # Raise the lower bounds of our search range on positions where
        # We need to increase our k-mer length (i.e. counts > 1)
        # Set the new lower (inclusive) bound to the current query length + 1
        lower_length_bound[counted_positions] = np.where(
            count_list > 1,
            current_length_query[counted_positions] + 1,
            lower_length_bound[counted_positions])

        # Calculate the new query length as the midpoint between the updated
        # upper and lower bounds
        current_length_query[counted_positions] = np.floor_divide(
            upper_length_bound[counted_positions] +
            lower_length_bound[counted_positions], 2)

        # If we have reduced our bounds to overlapping we have finished
        # searching on this position
        finished_search[counted_positions] = \
            finished_search[counted_positions] | \
            (upper_length_bound[counted_positions] <
                lower_length_bound[counted_positions])

        iteration_count += 1

    verbose_print(verbose,
                  f"Finished searching in {iteration_count-1} iterations")

    # TODO: Parameterize the type or change depending on lengths found
    return unique_lengths.astype(np.uint8), ambiguous_positions_skipped


def linear_search(index_filename,
                  sequence_segment,
                  kmer_lengths,
                  num_kmers,
                  num_threads,
                  verbose) -> npt.NDArray[np.uint8]:
    # Track which kmer positions have finished searching,
    # skipping any kmers starting with an ambiguous base
    # NB: Iterating over bytes returns ints
    finished_search = np.array([c == ord(b'N') for c in
                                sequence_segment.data[:num_kmers]])

    ambiguous_positions_skipped = finished_search.sum()
    verbose_print(verbose, f"Skipping {ambiguous_positions_skipped} ambiguous "
                  "positions")

    # List of minimum lengths (where 0 is nothing was found)
    unique_lengths = np.zeros(num_kmers, dtype=np.uint32)

    # For each kmer length
    for kmer_length in kmer_lengths:
        verbose_print(verbose, "{} k-mers remaining".format(
                      (~finished_search).sum()))
        verbose_print(verbose, "Counting {}-mers".format(kmer_length))

        # Create our working list of kmers to count
        working_kmers = []
        for i in np.nonzero(~finished_search)[0]:
            # Create the kmer from the sequence segment
            # NB: At the epilogue of the sequence, out of bounds
            # indexing after the end of the data will automatically
            # be truncated by numpy
            kmer = sequence_segment.data[i:i+kmer_length]
            # If it contains an ambiguous base
            if b'N' in kmer:
                # Ignore it for all longer kmer lengths (i.e. all
                # future iterations)
                finished_search[i] = True
            else:
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

        counted_positions = np.nonzero(np.copy(~finished_search))[0]

        # Assert that the number of indices to count and the number of counts
        # are equal
        assert counted_positions.size == count_list.size, \
            "Number of counted positions ({}) and number of counts ({}) " \
            "do not match".format(len(counted_positions), len(count_list))

        unique_lengths[counted_positions] = np.where(
            # Where the count is 1
            (count_list == 1) &
            # And if there is no current unique length recorded
            (unique_lengths[counted_positions] == 0),
            # Record the minimum kmer length found if it less than the current
            kmer_length,
            # Otherwise keep the current unique length
            unique_lengths[counted_positions])

        # We are finished searching at this position if the count was 1
        finished_search[counted_positions] = (count_list == 1)

        verbose_print(verbose, "{} unique {}-mers found".format(
            np.count_nonzero(unique_lengths == kmer_length), kmer_length))

    return unique_lengths.astype(np.uint8), ambiguous_positions_skipped


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
