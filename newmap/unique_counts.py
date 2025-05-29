from math import ceil, log2
from pathlib import Path
from typing import Union

import numpy as np
import numpy.typing as npt

from newmap._c_newmap_count_kmers import count_kmers_from_sequence
from newmap.util import INDEX_EXTENSION
from newmap.util import optional_gzip_open, verbose_print
from newmap.fasta import SequenceSegment, sequence_segments

KMER_RANGE_SEPARATOR = ":"
SEQUENCE_ID_SEPARATOR = ","

COMPLEMENT_TRANSLATE_TABLE = bytes.maketrans(b'ACGTacgt', b'TGCAtgca')
ALLOWED_BASES = b'ACGTactg'
UNIQUE_COUNT_FILENAME_FORMAT = "{}.unique.{}"


def write_unique_counts(fasta_filename: Path,
                        index_filename: Path,
                        kmer_batch_size: int,
                        kmer_lengths: list[int],
                        initial_search_length: int,
                        include_sequence_ids: list[bytes],
                        exclude_sequence_ids: list[bytes],
                        num_threads: int,
                        use_binary_search=False,
                        output_directory: Path = Path("."),
                        verbose: bool = False):

    max_kmer_length = max(kmer_lengths)
    min_kmer_length = min(kmer_lengths)

    # Use the maximum k-mer length to determine the data type used to store
    # lengths
    if max_kmer_length <= np.iinfo(np.uint8).max:
        data_type = np.uint8
        unique_count_suffix = "uint8"
    elif max_kmer_length <= np.iinfo(np.uint16).max:
        data_type = np.uint16
        unique_count_suffix = "uint16"
    else:
        data_type = np.uint32
        unique_count_suffix = "uint32"

    if (verbose and
            use_binary_search):
        verbose_print(verbose, "Max {} iterations over range {}-{}".format(
                      ceil(log2(max_kmer_length - min_kmer_length) + 1),
                      min_kmer_length, max_kmer_length))

    # Set a flag to ensure that some sequence was processed in case the include
    # or exclude sequence IDs are too strict
    sequence_processed = False

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

        # Keep of summary statistics
        total_ambiguous_positions = 0
        total_unique_lengths_count = 0
        total_no_unique_lengths_count = 0
        # Useful to the real upper and lower bounds for a sequence
        max_length_found = min_kmer_length
        min_length_found = max_kmer_length

        current_sequence_id = b''
        current_unique_filepath = Path("")

        # For each sequence buffer from the fasta file
        for sequence_segment in sequence_segments(fasta_file,
                                                  requested_sequence_length,
                                                  lookahead_length):

            # If we are on a new sequence
            if current_sequence_id != sequence_segment.id:
                # If we are only including specific sequences IDs
                if (include_sequence_ids or
                   exclude_sequence_ids):
                    # And this sequence IDs is not in the list
                    # Or if this sequence is in the exclude list
                    if ((sequence_segment.id not in include_sequence_ids) or
                       (sequence_segment.id in exclude_sequence_ids)):
                        # Skip the current buffer
                        continue

                # NB: At this point we are guaranteed to be on a new valid
                # sequence ID to be processed

                # Mark at least one sequence will be processed (for warnings)
                sequence_processed = True

                # Print out the summary statistics from the previously
                # processed chromosome
                print_summary_statisitcs(verbose,
                                         current_sequence_id,
                                         total_unique_lengths_count,
                                         total_ambiguous_positions,
                                         total_no_unique_lengths_count,
                                         max_length_found,
                                         min_length_found)

                # Update the current working sequence id
                current_sequence_id = sequence_segment.id
                # Update the current unique count filename
                current_unique_filepath = \
                    output_directory / UNIQUE_COUNT_FILENAME_FORMAT.format(
                        sequence_segment.id.decode(), unique_count_suffix)

                # Truncate any existing file for the new sequence
                open(current_unique_filepath, "wb").close()

                verbose_print(verbose,
                              "Writing minimum unique lengths for sequence "
                              "ID: {}".format(sequence_segment.id.decode()))

            num_kmers = get_num_kmers(sequence_segment, max_kmer_length)

            verbose_print(verbose, "Processing {} k-mers".format(num_kmers))

            if use_binary_search:
                segment_unique_counts, ambiguous_count = \
                    binary_search(index_filename,
                                  sequence_segment,
                                  min_kmer_length,
                                  max_kmer_length,
                                  initial_search_length,
                                  num_threads,
                                  data_type,  # type: ignore
                                  verbose)
            else:
                segment_unique_counts, ambiguous_count = \
                    linear_search(index_filename,
                                  sequence_segment,
                                  kmer_lengths,
                                  num_kmers,
                                  num_threads,
                                  data_type,  # type: ignore
                                  verbose)

            # Update summary statistics
            unique_lengths_count = np.count_nonzero(segment_unique_counts)
            total_ambiguous_positions += ambiguous_count
            total_unique_lengths_count += unique_lengths_count
            total_no_unique_lengths_count += (num_kmers -
                                              unique_lengths_count -
                                              ambiguous_count)

            # If any amount of unique lengths were found
            if unique_lengths_count > 0:
                max_length_found = max(max_length_found,
                                       segment_unique_counts.max())
                min_length_found = min(
                    min_length_found,
                    segment_unique_counts[
                        segment_unique_counts.nonzero()
                     ].min()
                )
            # Otherwise
            # NB: If no unique lengths were found
            # neither the max or the minimum will not change
            else:
                verbose_print(verbose, "No unique lengths found for this "
                                       "sequence segment")

            # Append the unique counts to a unique count file per sequence
            with open(current_unique_filepath, "ab") as unique_count_file:
                segment_unique_counts.tofile(unique_count_file)

        print_summary_statisitcs(verbose,
                                 current_sequence_id,
                                 total_unique_lengths_count,
                                 total_ambiguous_positions,
                                 total_no_unique_lengths_count,
                                 max_length_found,
                                 min_length_found)

        # If somehow we skipped all the sequences
        if not sequence_processed:
            # Print out a warning
            # If include sequences were used
            if include_sequence_ids:
                # None of the included sequences were found
                raise ValueError("None of the included sequences were found: "
                                 f"{include_sequence_ids}")
            # If excluded sequences were used
            elif exclude_sequence_ids:
                # Excluded sequences removed all possilibities
                raise ValueError("The excluded sequences were too strict and "
                                 "nothing was processed: "
                                 f"{exclude_sequence_ids}")


def binary_search(index_filename: Path,
                  sequence_segment: SequenceSegment,
                  min_kmer_length: int,
                  max_kmer_length: int,
                  initial_search_length: int,
                  num_threads: int,
                  data_type: Union[np.uint8, np.uint16, np.uint32],
                  verbose: bool) -> tuple[npt.NDArray[np.uint], int]:

    num_kmers = get_num_kmers(sequence_segment, max_kmer_length)

    # Track which kmer positions have finished searching,
    # skipping any kmers starting with an ambiguous base
    finished_search = get_ambiguous_sequence_mask(sequence_segment,
                                                  num_kmers)

    # Print out the number of ambiguous positions skipped
    ambiguous_positions_skipped = finished_search.sum()
    verbose_print(verbose, f"Skipping {ambiguous_positions_skipped} ambiguous "
                  "positions")

    # Track minimum and maximum kmer lengths for each position
    # NB: Inclusive bounds for search lengths
    lower_length_bound = np.full(num_kmers, min_kmer_length, dtype=data_type)

    # The upper search length is bounded by the minimum of:
    # The maximum length in our search query range
    upper_length_bound = np.full(num_kmers, max_kmer_length, dtype=data_type)

    # Update the upper search bounds based on ambiguous bases in the sequence
    # and sequence length
    update_upper_search_bound(upper_length_bound,
                              finished_search,
                              max_kmer_length,
                              len(sequence_segment.data))

    # Create the current length query as the midpoint between the lower and
    # upper search bounds
    current_length_query = np.floor(
        (upper_length_bound / 2) + (lower_length_bound / 2)
    ).astype(data_type)

    # If there was an initial search length specified
    if initial_search_length:
        # Set all current search lengths to at most the initial search length
        current_length_query[current_length_query > initial_search_length] = \
            initial_search_length

    # Mark any positions whose upper search bound is less than the minimum
    # search length as complete
    finished_search[upper_length_bound < min_kmer_length] = True

    # If verbosity is on
    if verbose:
        # Calculate out the number of kmer upper search ranges truncated
        upper_bound_change_count = np.count_nonzero(
            upper_length_bound[(~finished_search).nonzero()] < max_kmer_length)
        # And print it out
        verbose_print(verbose, f"{upper_bound_change_count} k-mer maximum "
                               "search ranges truncated due to ambiguity")

        # Calculate the number of of kmers too short to be counted
        short_kmers_discarded_count = (finished_search.sum() -
                                       ambiguous_positions_skipped)
        # And print it out
        verbose_print(verbose, f"{short_kmers_discarded_count} k-mers shorter "
                      "than the minimum length discarded due to ambiguity")

    # List of unique minimum lengths (where 0 is nothing was found)
    unique_lengths = np.zeros(num_kmers, dtype=data_type)

    iteration_count = 1

    # While there are still kmers to search
    while not np.all(finished_search):
        verbose_print(verbose, "Iteration {}".format(iteration_count))
        verbose_print(verbose, "{} k-mer positions remaining".format(
                      np.sum(~finished_search)))

        # Track which kmer positions need to be counted on the index
        # Create a list of indices where each index refers to the corresponding
        # position in the given sequence segment
        kmer_indices = np.nonzero(~finished_search)[0]

        # Get the occurences of the kmers on both strands
        count_list = get_kmer_counts(index_filename,
                                     sequence_segment.data,
                                     kmer_indices.tolist(),
                                     current_length_query[
                                        kmer_indices].tolist(),
                                     num_threads)

        # Assert that the number of indices to count and the number of counts
        # are equal
        assert kmer_indices.size == count_list.size, \
            "Number of counted positions ({}) and number of counts ({}) " \
            "do not match".format(len(kmer_indices), len(count_list))

        # Where we have counts of 1
        unique_lengths[kmer_indices] = np.where(
            (count_list == 1) &
            # And if there is no current unique length recorded
            ((unique_lengths[kmer_indices] == 0) |
             # Or there is a smaller length found than the current min length
             (current_length_query[kmer_indices] <
              unique_lengths[kmer_indices])),
            # Record the minimum kmer length found if it less than the current
            current_length_query[kmer_indices],
            # Otherwise keep the current unique length
            unique_lengths[kmer_indices])

        # If we have a k-mer count of 1 and the current length queried is the
        # same as the lower bound (i.e. can't get smaller for a unique count)
        # This position has finished searching
        finished_search[kmer_indices] = np.where(
            count_list == 1,
            current_length_query[kmer_indices] ==
            lower_length_bound[kmer_indices],
            finished_search[kmer_indices])

        # If we have a k-mer count > 1 and the current length queried is the
        # same as the uppper bound (i.e. can't find a unique length or larger)
        # This position has finished searching
        finished_search[kmer_indices] = np.where(
            count_list > 1,
            current_length_query[kmer_indices] ==
            upper_length_bound[kmer_indices],
            finished_search[kmer_indices])

        # Update the query length and bounds for the next iteration

        # Lower the upper bounds of our search range on positions where
        # we need to decrease our k-mer length (i.e. counts == 1)
        # Set the new upper (inclusive) bound to the current query length - 1
        upper_length_bound[kmer_indices] = np.where(
            count_list == 1,
            current_length_query[kmer_indices] - 1,
            upper_length_bound[kmer_indices])

        # Raise the lower bounds of our search range on positions where
        # we need to increase our k-mer length (i.e. counts > 1)
        # Set the new lower (inclusive) bound to the current query length + 1
        lower_length_bound[kmer_indices] = np.where(
            count_list > 1,
            current_length_query[kmer_indices] + 1,
            lower_length_bound[kmer_indices])

        # Calculate the new query length as the midpoint between the updated
        # upper and lower bounds
        # NB: Avoid overflow by dividing first before sum
        current_length_query[kmer_indices] = np.floor(
            (upper_length_bound[kmer_indices] / 2) +
            (lower_length_bound[kmer_indices] / 2)).astype(data_type)

        iteration_count += 1

    verbose_print(verbose,
                  f"Finished searching in {iteration_count-1} iterations")

    return unique_lengths, ambiguous_positions_skipped


def linear_search(index_filename: Path,
                  sequence_segment: SequenceSegment,
                  kmer_lengths: list[int],
                  num_kmers: int,
                  num_threads: int,
                  data_type: Union[np.uint8, np.uint16, np.uint32],
                  verbose: bool) -> tuple[npt.NDArray[np.uint], int]:
    # Track which kmer positions have finished searching,
    # skipping any kmers starting with an ambiguous base
    # NB: Iterating over bytes returns ints
    finished_search = get_ambiguous_sequence_mask(sequence_segment,
                                                  num_kmers)

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

        max_kmer_query_lengths = []

        # For every position that does not have an ambiguous base
        for i in np.nonzero(~finished_search)[0]:
            # Create the kmer from the sequence segment
            # NB: At the epilogue of the sequence, out of bounds
            # indexing after the end of the data will automatically
            # be truncated by numpy
            kmer = sequence_segment.data[i:i+kmer_length]
            # TODO: Fix this by only permitted allowed bases in ALLOWED_BASES
            # If it contains an ambiguous base
            if b'N' in kmer:
                # Ignore it for all longer kmer lengths (i.e. all
                # future iterations)
                finished_search[i] = True
            # Otherwise:
            else:
                # Record the length of this k-mer
                max_kmer_query_lengths.append(len(kmer))

        kmer_indices = np.nonzero(~finished_search)[0]

        # If there are no kmers to count due to ambiguity
        if kmer_indices.size == 0:
            verbose_print(verbose,
                          f"No {kmer_length}-mers remaining to be found, "
                          f"skipping to next kmer batch")
            # Skip to the next kmer batch
            break
        else:
            verbose_print(verbose,
                          f"{kmer_indices.size} {kmer_length}-mers remaining "
                          f"to be counted")

        # Count the kmer occurences on the index
        count_list = get_kmer_counts(index_filename,
                                     sequence_segment.data,
                                     kmer_indices.tolist(),
                                     max_kmer_query_lengths,
                                     num_threads)

        # Assert that the number of indices to count and the number of counts
        # are equal
        assert kmer_indices.size == count_list.size, \
            "Number of counted positions ({}) and number of counts ({}) " \
            "do not match".format(len(kmer_indices), len(count_list))

        unique_lengths[kmer_indices] = np.where(
            # Where the count is 1
            (count_list == 1) &
            # And if there is no current unique length recorded
            (unique_lengths[kmer_indices] == 0),
            # Record the minimum kmer length found if it less than the current
            kmer_length,
            # Otherwise keep the current unique length
            unique_lengths[kmer_indices])

        # We are finished searching at this position if the count was 1
        finished_search[kmer_indices] = (count_list == 1)

        verbose_print(verbose, "{} unique {}-mers found".format(
            np.count_nonzero(unique_lengths == kmer_length), kmer_length))

    return unique_lengths.astype(data_type), ambiguous_positions_skipped


def get_kmer_counts(index_filename: Path,
                    sequence_data: bytes,
                    index_list: list[int],
                    kmer_lengths: list[int],
                    num_threads: int) -> npt.NDArray[np.uint32]:

    # Count the occurences of kmers on the forward strand
    count_list = np.array(count_kmers_from_sequence(
                            str(index_filename),
                            sequence_data,
                            index_list,
                            kmer_lengths,
                            num_threads),
                          dtype=np.uint32)

    # Count the occurrences of kmers on the reverse complement strand
    # TODO: Add no reverse complement option to skip this to
    # support bisulfite treated kmer counting
    # TODO: Add option for a complement table
    reverse_complement_sequence = \
        sequence_data.translate(COMPLEMENT_TRANSLATE_TABLE)[::-1]

    sequence_data_length = len(sequence_data)
    reverse_index_list = [sequence_data_length - i - kmer_length
                          for i, kmer_length in zip(index_list, kmer_lengths)]

    count_list += np.array(count_kmers_from_sequence(
                             str(index_filename),
                             reverse_complement_sequence,
                             reverse_index_list,
                             kmer_lengths,
                             num_threads),
                           dtype=np.uint32)

    # If any element in the count list is 0
    if np.any(count_list == 0):
        # There is very likely a mismatch between sequence and the index
        # There is also a chance there is a bug in the index

        # NB: Get first element of tuple, then first element in numpy array
        first_problem_kmer_index = (count_list == 0).nonzero()[0][0]

        problem_sequence_index = index_list[first_problem_kmer_index]
        problem_sequence_length = kmer_lengths[first_problem_kmer_index]

        first_problem_kmer = sequence_data[
            problem_sequence_index:
            problem_sequence_index+problem_sequence_length].decode("utf-8")

        raise RuntimeError(f"The following generated k-mer was not found in "
                           f"the index:\n{first_problem_kmer}\n"
                           f"Possibly a mismatch between the sequence and the "
                           f"index.")

    return count_list


def get_num_kmers(sequence_segment: SequenceSegment,
                  max_kmer_length: int):
    """Returns the maximum number of kmers from a SequenceSegment"""
    sequence_length = len(sequence_segment.data)
    # If this is the last sequence segment
    if sequence_segment.epilogue:
        # The number of kmers is equal to the entire length of the
        # sequence segment
        return sequence_length
    # Otherwise
    else:
        # The number of kmers is equal to kmer batch size
        # (or the sequence segment length minus the lookahead)
        lookahead_length = max_kmer_length - 1
        return sequence_length - lookahead_length


def get_ambiguous_sequence_mask(sequence_segment: SequenceSegment,
                                num_positions: int):
    """Returns a boolean array of ambiguous positions in a sequence segment
       Where True is an ambiguous position and False is a non-ambiguous"""

    # Track which kmer positions have finished searching,
    # skipping any kmers starting with an ambiguous base
    sequence_buffer = np.frombuffer(sequence_segment.data,
                                    dtype=np.uint8,
                                    count=num_positions)

    # Initially set all positions to non-ambiguous (False)
    # ambiguous_array_positions = np.full(num_positions, False, dtype=bool)

    # Initially mark every position as ambiguous
    allowed_positions = np.full(num_positions, False, dtype=bool)
    # For every allowed base
    for base in ALLOWED_BASES:
        # Set this position to non-ambiguous
        allowed_positions |= (sequence_buffer == base)

    # Return mask
    return ~allowed_positions


def update_upper_search_bound(upper_search_bound: npt.NDArray,
                              ambiguity_mask: npt.NDArray,
                              max_search_length: int,
                              sequence_buffer_length: int):
    """Modifies upper_search_bound to the maximum k-mer search length possible
    at each position based on the ambiguity mask and start and end of the
    working sequence positions"""

    # We assume that there is at most (max_search_length - 1) positions after
    # the length of the ambiguity mask to allow searching at the final
    # position at the maximum length
    excess_buffer_length = sequence_buffer_length - ambiguity_mask.size
    assert excess_buffer_length < max_search_length, \
        "Excess sequence buffer length is greater than the maximum search " \
        "length"

    # If all positions non-ambiguous
    # NB: Likely the most common case
    if (~ambiguity_mask).all():
        # If there is not enough sequence buffer to search at the maximum
        # NB: There is a possibility that the maximum search range is larger
        # than our entire buffer
        if excess_buffer_length < (max_search_length - 1):
            # Set the last position's maximum search length to the excess
            # buffer length and all positions before it increasing to the
            # maximum search length
            start_value = max_search_length - 1
            end_value = excess_buffer_length
            new_values = np.arange(start_value, end_value, -1)
            assign_remaining_array_values(upper_search_bound, new_values)
    # Else if all positions are ambiguous
    elif ambiguity_mask.all():
        # Nothing to do
        return
    # Otherwise there is a mix of ambiguous and non-ambiguous positions
    else:
        # Find the 0-based start-inclusive and end-exclusive intervals (for
        # indexing) of non-ambiguous areas including the start and end in the
        # interval list

        # Compare the ambiguity list with a shifted version of itself
        left_shift = ambiguity_mask[1:]
        reference_mask = ambiguity_mask[:-1]  # NB: Keep same dimensions

        # Where our reference mask is True but left shift is False
        # We have the start of an non-ambiguous region
        # NB: Add 1 due to the shift
        start_indices = np.where(reference_mask & ~left_shift)[0] + 1

        # Where our reference mask is False but left shift is True
        # We have the end of an non-ambiguous region
        end_indices = np.where(~reference_mask & left_shift)[0] + 1

        # If we have end indices
        if end_indices.size:
            # If the first non-ambiguous region has an ending before the first
            # start position
            # Or if there is no starting index
            if not start_indices.size or \
              end_indices[0] < start_indices[0]:
                # The first start position must be at the beginning (0)
                start_indices = np.insert(start_indices, 0, 0)

        # Figure out if we have reached end of the current buffer
        reached_end_of_buffer = False

        # If we have start indices
        if start_indices.size:
            # If the last non-ambiguous start is after the last end position
            # Or if there is no ending index
            if not end_indices.size or \
               start_indices[-1] > end_indices[-1]:
                # The last end position must be at the end
                end_indices = np.append(end_indices, ambiguity_mask.size)
                # We have reached the end of the buffer
                reached_end_of_buffer = True

        # Iterate over each start and end non-ambiguous interval
        for start, end in zip(start_indices, end_indices):
            # If we are the on the last interval and are at the end of the
            # buffer
            end_value = 0  # NB: Exclusive in np.arange, will go down to 1
            if end == end_indices[-1] and \
               reached_end_of_buffer:
                # Update the max search length at the end of this interval
                end_value = excess_buffer_length

            non_ambiguous_interval = upper_search_bound[start:end]

            # A descending array from max search length to our largest possible
            # max search length
            start_value = max_search_length - 1
            new_values = np.arange(max_search_length - 1, end_value, -1)

            # Assign the new values to the end of the non-ambiguous interval
            assign_remaining_array_values(non_ambiguous_interval, new_values)


def assign_remaining_array_values(input_array: npt.NDArray,
                                  new_values: npt.NDArray):
    """Modifies the input array so the last n values are assigned the contents
    of the new values array regardless of dimension of either array."""
    num_values = min(new_values.size, input_array.size)
    # Only update if there is anything to assign
    if num_values > 0:
        input_array[-num_values:] = new_values[-num_values:]


def print_summary_statisitcs(verbose: bool,
                             sequence_id: bytes,
                             total_unique_lengths_count: int,
                             total_ambiguous_positions: int,
                             total_no_unique_lengths_count: int,
                             max_length_found: int,
                             min_length_found: int):
    if (verbose and
            total_unique_lengths_count):
        verbose_print(verbose,
                      f"Finished writing minimum unique lengths for sequence "
                      f"ID: {sequence_id.decode()}")
        verbose_print(verbose,
                      f"{total_unique_lengths_count} unique lengths found")
        verbose_print(verbose,
                      f"{total_ambiguous_positions} positions skipped due to "
                      "ambiguity")
        verbose_print(verbose,
                      f"{total_no_unique_lengths_count} positions with no "
                      "unique length found")
        verbose_print(verbose,
                      f"{max_length_found}-mer maximum unique length found")
        verbose_print(verbose,
                      f"{min_length_found}-mer minimum unique length found")


def main(args):
    fasta_filename = args.fasta_file
    index_filename = args.index_file
    kmer_batch_size = args.kmer_batch_size
    kmer_lengths_arg = args.search_range
    output_directory = args.output_directory
    initial_search_length = args.initial_search_length
    include_sequences_arg = args.include_sequences
    exclude_sequences_arg = args.exclude_sequences
    num_threads = args.num_threads
    verbose = args.verbose

    # If an output directory was specified
    if output_directory:
        # Create the output directory (and parents) if it does not exist
        output_directory = Path(output_directory)
        # Do not raise an error if it already exists
        output_directory.mkdir(parents=True, exist_ok=True)
    # Otherwise use the current working directory
    else:
        output_directory = Path(".")

    # If no index filename was specified
    if not index_filename:
        # Use the basename of the fasta file
        index_filename = Path(fasta_filename).with_suffix("."+INDEX_EXTENSION)

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
            raise ValueError("Could not parse k-mer search range format")

        if min_kmer_length > max_kmer_length:
            raise ValueError("K-mer range start length is larger than the end "
                             "length")
        kmer_lengths = list(range(min_kmer_length, max_kmer_length + 1))
        use_binary_search = True
    # Otherwise assume a single digit or a list of digits seperated by a comma
    # NB: Assume in sorted order?
    else:
        kmer_lengths = list(map(int, kmer_lengths_arg.split(",")))

    if (initial_search_length and
       not use_binary_search):
        raise ValueError("Initial search length only valid when a range of "
                         "k-mer lengths is given")

    if (include_sequences_arg and
       exclude_sequences_arg):
        raise ValueError("Cannot include and exclude sequences at the same "
                         "time")

    include_sequence_ids = []
    if include_sequences_arg:
        include_sequence_ids = \
            [s.encode() for s in
             include_sequences_arg.split(SEQUENCE_ID_SEPARATOR)]

    exclude_sequence_ids = []
    if exclude_sequences_arg:
        exclude_sequence_ids = \
            [s.encode() for s in
             exclude_sequences_arg.split(SEQUENCE_ID_SEPARATOR)]

    write_unique_counts(Path(fasta_filename),
                        Path(index_filename),
                        kmer_batch_size,
                        kmer_lengths,
                        initial_search_length,
                        include_sequence_ids,
                        exclude_sequence_ids,
                        num_threads,
                        use_binary_search,
                        output_directory,
                        verbose)
