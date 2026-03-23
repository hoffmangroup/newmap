from contextlib import ExitStack
from collections.abc import Callable
from dataclasses import dataclass
from functools import partial
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
INDEX_FILE_SEPARATOR = ","
FASTA_FILE_SEPARATOR = ","

COMPLEMENT_TRANSLATE_TABLE = bytes.maketrans(b'ACGTacgt', b'TGCAtgca')
ALLOWED_BASES = b'ACGTactg'
UNIQUE_COUNT_FILENAME_FORMAT = "{}.unique.{}"


# NB: This is the default function for logging when no verbosity is specified
def nil_search_log(*_: str):
    pass


@dataclass(frozen=True)
class SearchConfig:
    # Position args
    fasta_filepaths: list[Path]
    fmindex_filepaths: list[Path]

    # Output arguments
    kmer_lengths: list[int]
    is_binary_search: bool

    use_reverse_complement: bool
    output_directory: Path

    include_sequence_ids: list[bytes]
    exclude_sequence_ids: list[bytes]

    # Performance arguments
    num_threads: int
    kmer_batch_size: int
    initial_search_length: int

    # Verbosity
    # NB: Only used when additional calculations are needed for logging
    verbose: bool
    # Logging function
    log: Callable[[str], None]

    @classmethod
    def from_args(cls, args):
        fasta_file_specification = args.fasta_file
        index_file_specification = args.index_file
        kmer_lengths_arg = args.search_range
        output_directory = args.output_directory
        initial_search_length = args.initial_search_length
        include_sequences_arg = args.include_sequences
        exclude_sequences_arg = args.exclude_sequences

        # If an output directory was specified
        if output_directory:
            # Create the output directory (and parents) if it does not exist
            output_directory = Path(output_directory)
            # Do not raise an error if it already exists
            output_directory.mkdir(parents=True, exist_ok=True)
        # Otherwise use the current working directory
        else:
            output_directory = Path(".")

        # Process the fasta filename specificaiton into a list of Paths
        # even if only a single fasta file is specified
        if FASTA_FILE_SEPARATOR in fasta_file_specification:
            fasta_filenames = [Path(filename) for filename in
                               fasta_file_specification
                               .split(FASTA_FILE_SEPARATOR)]
        else:
            fasta_filenames = [Path(fasta_file_specification)]

        # Process the index filename specificaiton into a list of Paths
        # even if only a single index file is specified

        # If no index filename was specified
        if not index_file_specification:
            # Use the basename of the (first) fasta file and cwd
            fasta_filename = fasta_filenames[0]
            index_filename = Path(fasta_filename).stem + "." + INDEX_EXTENSION
            # NB: Explicit list of length 1
            index_filenames = [Path(index_filename)]
        elif INDEX_FILE_SEPARATOR in index_file_specification:
            index_filenames = [Path(filename) for filename in
                               index_file_specification
                               .split(INDEX_FILE_SEPARATOR)]
        else:
            # NB: A single index filename in a list
            index_filenames = [Path(index_file_specification)]

        # Check that the index files exists
        for index_filename in index_filenames:
            if not index_filename.is_file():  # type: ignore
                raise FileNotFoundError(f"Index file not found: "
                                        f"{index_filename}")

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
                raise ValueError("K-mer range start length is larger than the "
                                 "end length")
            kmer_lengths = list(range(min_kmer_length, max_kmer_length + 1))
            use_binary_search = True
        # Otherwise assume a single digit or a list of digits seperated by a
        # comma
        # NB: Assume in sorted order?
        else:
            kmer_lengths = list(map(int, kmer_lengths_arg.split(",")))

        if (initial_search_length and
           not use_binary_search):
            raise ValueError("Initial search length only valid when a range "
                             "of k-mer lengths is given")

        if (include_sequences_arg and
           exclude_sequences_arg):
            raise ValueError("Cannot specify both include and exclude "
                             "sequences")

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

        if args.verbose:
            logging_function = partial(verbose_print, True)
        else:
            logging_function = nil_search_log

        return cls(
                   # Position args
                   fasta_filepaths=fasta_filenames,
                   fmindex_filepaths=index_filenames,

                   # Output arguments
                   kmer_lengths=kmer_lengths,
                   is_binary_search=use_binary_search,

                   use_reverse_complement=not args.norc,
                   output_directory=output_directory,

                   include_sequence_ids=include_sequence_ids,
                   exclude_sequence_ids=exclude_sequence_ids,

                   log=logging_function,
                   verbose=args.verbose,

                   # Performance arguments
                   num_threads=args.num_threads,
                   kmer_batch_size=args.kmer_batch_size,
                   initial_search_length=initial_search_length
        )


def write_unique_counts(config: SearchConfig):

    max_kmer_length = max(config.kmer_lengths)
    min_kmer_length = min(config.kmer_lengths)

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

    if (config.is_binary_search):
        config.log("Max {} iterations over range {}-{}".format(
                   ceil(log2(max_kmer_length - min_kmer_length) + 1),
                   min_kmer_length, max_kmer_length))

    # Set a flag to ensure that some sequence was processed in case the include
    # or exclude sequence IDs are too strict
    sequence_processed = False

    with ExitStack() as stack:
        fasta_files = [stack.enter_context(optional_gzip_open(filename, "rb"))
                       for filename in config.fasta_filepaths]

        # Allow a lookahead on the last kmer to include the final dinucleotide
        # in the sequence buffer
        lookahead_length = max(config.kmer_lengths) - 1

        # Get sequence lengths long enough to process the longest kmer of
        # interest
        # Requires a lookahead so that last dinucleotide has enough to
        # form a kmer of at least the maximum length
        requested_sequence_length = config.kmer_batch_size + lookahead_length

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

        # Create a list iterable sequences for each fasta file
        sequences = [sequence_segments(
                     fasta_file,  # type: ignore
                     requested_sequence_length,
                     lookahead_length)
                     for fasta_file in fasta_files]

        # For each list of sequence segments from each corresponding fasta file
        # NB: Where each list will contain the nth sequence segment from each
        # fasta file
        for sequence_segment_list in zip(*sequences):
            # NB: We assume that all sequences match in length and are in order
            # in each fasta file
            # NB: We only take the meta information from the first sequence
            # file in order to process sequence IDs correctly
            sequence_segment = sequence_segment_list[0]

            # If we are on a new sequence
            if current_sequence_id != sequence_segment.id:
                # If we are only including specific sequences IDs
                if (config.include_sequence_ids or
                   config.exclude_sequence_ids):
                    # And this sequence IDs is not in the list
                    # Or if this sequence is in the exclude list
                    if ((sequence_segment.id not in
                         config.include_sequence_ids) or
                       (sequence_segment.id in config.exclude_sequence_ids)):
                        # Skip the current buffer
                        continue

                # NB: At this point we are guaranteed to be on a new valid
                # sequence ID to be processed

                # Mark at least one sequence will be processed (for warnings)
                sequence_processed = True

                # Print out the summary statistics from the previously
                # processed chromosome
                print_summary_statisitcs(config,
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
                    config.output_directory / \
                    UNIQUE_COUNT_FILENAME_FORMAT.format(
                        sequence_segment.id.decode(), unique_count_suffix)

                # Truncate any existing file for the new sequence
                open(current_unique_filepath, "wb").close()

                config.log("Writing unique lengths for sequence "
                           "ID: {}".format(sequence_segment.id.decode()))

            num_kmers = get_num_kmers(sequence_segment, max_kmer_length)

            config.log("Processing {} k-mers".format(num_kmers))

            if config.is_binary_search:
                segment_unique_counts, ambiguous_count = \
                    binary_search(config,
                                  sequence_segment_list,
                                  min_kmer_length,
                                  max_kmer_length,
                                  data_type  # type: ignore
                                  )
            else:
                segment_unique_counts, ambiguous_count = \
                    linear_search(config,
                                  sequence_segment_list,
                                  num_kmers,
                                  data_type  # type: ignore
                                  )

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
                config.log("No unique lengths found for this "
                           "sequence segment")

            # Append the unique counts to a unique count file per sequence
            with open(current_unique_filepath, "ab") as unique_count_file:
                segment_unique_counts.tofile(unique_count_file)

        print_summary_statisitcs(config,
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
            if config.include_sequence_ids:
                # None of the included sequences were found
                raise ValueError("None of the included sequences were found: "
                                 f"{config.include_sequence_ids}")
            # If excluded sequences were used
            elif config.exclude_sequence_ids:
                # Excluded sequences removed all possilibities
                raise ValueError("The excluded sequences were too strict and "
                                 "nothing was processed: "
                                 f"{config.exclude_sequence_ids}")


def binary_search(config: SearchConfig,
                  sequence_segments: tuple[SequenceSegment],
                  min_kmer_length: int,
                  max_kmer_length: int,
                  data_type: Union[np.uint8, np.uint16, np.uint32]
                  ) -> tuple[npt.NDArray[np.uint], int]:

    # NB: Use the first (and maybe only) sequence to determine the number of
    # kmers and search bounds. We assume all sequences specified are equal in
    # length and ambiguous positions
    first_sequence_segment = sequence_segments[0]

    num_kmers = get_num_kmers(first_sequence_segment, max_kmer_length)

    # Track which kmer positions have finished searching,
    # skipping any kmers starting with an ambiguous base
    finished_search = get_ambiguous_sequence_mask(first_sequence_segment,
                                                  num_kmers)

    # Print out the number of ambiguous positions skipped
    ambiguous_positions_skipped = finished_search.sum()
    config.log(f"Skipping {ambiguous_positions_skipped} ambiguous "
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
                              len(first_sequence_segment.data))

    # Create the current length query as the midpoint between the lower and
    # upper search bounds
    current_length_query = np.floor(
        (upper_length_bound / 2) + (lower_length_bound / 2)
    ).astype(data_type)

    # If there was an initial search length specified
    if config.initial_search_length:
        # Set all current search lengths to at most the initial search length
        current_length_query[
            current_length_query > config.initial_search_length] = \
            config.initial_search_length

    # Mark any positions whose upper search bound is less than the minimum
    # search length as complete
    finished_search[upper_length_bound < min_kmer_length] = True

    # If verbosity is on
    if config.verbose:
        # Calculate out the number of kmer upper search ranges truncated
        upper_bound_change_count = np.count_nonzero(
            upper_length_bound[(~finished_search).nonzero()] < max_kmer_length)
        # And print it out
        config.log(f"{upper_bound_change_count} k-mer maximum "
                   "search ranges truncated due to ambiguity")

        # Calculate the number of of kmers too short to be counted
        short_kmers_discarded_count = (finished_search.sum() -
                                       ambiguous_positions_skipped)
        # And print it out
        config.log(f"{short_kmers_discarded_count} k-mers shorter "
                   "than the minimum length discarded due to ambiguity")

    # List of unique minimum lengths (where 0 is nothing was found)
    unique_lengths = np.zeros(num_kmers, dtype=data_type)

    iteration_count = 1

    # While there are still kmers to search
    while not np.all(finished_search):
        config.log("Iteration {}".format(iteration_count))
        config.log("{} k-mer positions remaining".format(
                   np.sum(~finished_search)))

        # Track which kmer positions need to be counted on the index
        # Create a list of indices where each index refers to the corresponding
        # position in the given sequence segment
        kmer_indices = np.nonzero(~finished_search)[0]

        # Get the occurences of the kmers on both strands
        count_list = get_kmer_counts(config,
                                     [segment.data for segment in
                                      sequence_segments],
                                     kmer_indices.tolist(),
                                     current_length_query[
                                        kmer_indices].tolist())

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

    config.log(f"Finished searching in {iteration_count-1} iterations")

    return unique_lengths, ambiguous_positions_skipped


def linear_search(config: SearchConfig,
                  sequence_segments: tuple[SequenceSegment],
                  num_kmers: int,
                  data_type: Union[np.uint8, np.uint16, np.uint32]
                  ) -> tuple[npt.NDArray[np.uint], int]:

    # NB: Use the first (and maybe only) sequence to determine the number of
    # kmers and search bounds. We assume all sequences specified are equal in
    # length and ambiguous positions
    first_sequence_segment = sequence_segments[0]

    # Track which kmer positions have finished searching,
    # skipping any kmers starting with an ambiguous base
    # NB: Iterating over bytes returns ints
    finished_search = get_ambiguous_sequence_mask(first_sequence_segment,
                                                  num_kmers)

    ambiguous_positions_skipped = finished_search.sum()
    config.log(f"Skipping {ambiguous_positions_skipped} ambiguous "
               "positions")

    # List of minimum lengths (where 0 is nothing was found)
    unique_lengths = np.zeros(num_kmers, dtype=np.uint32)

    # For each kmer length
    for kmer_length in config.kmer_lengths:
        config.log("{} k-mers remaining".format((~finished_search).sum()))
        config.log("Counting {}-mers".format(kmer_length))

        kmer_query_lengths = []

        # For every position that does not have an ambiguous base
        for i in np.nonzero(~finished_search)[0]:
            # Create the kmer from the sequence segment
            # NB: At the epilogue of the sequence, out of bounds
            # indexing after the end of the data will automatically
            # be truncated by numpy
            kmer = first_sequence_segment.data[i:i+kmer_length]
            # TODO: Fix this by only permitted allowed bases in ALLOWED_BASES
            # If it contains an ambiguous base
            if b'N' in kmer:
                # Ignore it for all longer kmer lengths (i.e. all
                # future iterations)
                finished_search[i] = True
            # Otherwise:
            else:
                # Record the length of this k-mer
                kmer_query_lengths.append(len(kmer))

        kmer_indices = np.nonzero(~finished_search)[0]

        # If there are no kmers to count due to ambiguity
        if kmer_indices.size == 0:
            config.log(f"No {kmer_length}-mers remaining to be found, "
                       f"skipping to next kmer batch")
            # Skip to the next kmer batch
            break
        else:
            config.log(f"{kmer_indices.size} {kmer_length}-mers remaining "
                       f"to be counted")

        # Count the kmer occurences on the index
        count_list = get_kmer_counts(config,
                                     [segment.data for segment in
                                      sequence_segments],
                                     kmer_indices.tolist(),
                                     kmer_query_lengths)

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

        config.log("{} unique {}-mers found".format(
            np.count_nonzero(unique_lengths == kmer_length), kmer_length))

    return unique_lengths.astype(data_type), ambiguous_positions_skipped


def get_kmer_counts(config: SearchConfig,
                    sequences: list[bytes],
                    sequence_indices: list[int],
                    kmer_lengths: list[int]) -> npt.NDArray[np.uint32]:

    count_list = np.zeros(len(sequence_indices), dtype=np.uint32)
    # For each index (usually just one)
    # NB: We assume a guarantee that we have > 0 indexes at this point
    for fmindex_filepath in config.fmindex_filepaths:

        if len(config.fmindex_filepaths) > 1:
            config.log(f"Counting kmers from index: {fmindex_filepath}")

        # For every sequence segment in the list of sequence segments for this
        # index (usually just one)
        for i, sequence in enumerate(sequences):

            if len(sequences) > 1:
                # TODO: Print out sequence ID?
                config.log(f"Counting sequence #{i+1}")

            # Count the occurences of kmers on the forward strand
            count_list += np.array(count_kmers_from_sequence(
                                   str(fmindex_filepath),
                                   sequence,
                                   sequence_indices,
                                   kmer_lengths,
                                   config.num_threads),
                                   dtype=np.uint32)

            # If we are not skipping the reverse complement strand (default)
            if config.use_reverse_complement:
                # Count the occurrences of kmers on the reverse complement
                # strand
                # TODO: Add option for a complement table
                # TODO: The reverse complement should only be calculated once
                reverse_complement_sequence = \
                    sequence.translate(COMPLEMENT_TRANSLATE_TABLE)[::-1]

                sequence_length = len(sequence)
                # TODO: Should ideally only do this once
                reverse_index_list = [sequence_length - i - kmer_length
                                      for i, kmer_length in
                                      zip(sequence_indices, kmer_lengths)]

                count_list += np.array(count_kmers_from_sequence(
                                         str(fmindex_filepath),
                                         reverse_complement_sequence,
                                         reverse_index_list,
                                         kmer_lengths,
                                         config.num_threads),
                                       dtype=np.uint32)

    # If any element in the count list is 0
    # NB: np.all will evalute to False if any element is 0 and it is more
    # likely to short-circuit evaluate than count_nonzero
    contains_a_zero_count = not np.all(count_list)
    if contains_a_zero_count:
        # There is very likely a mismatch between sequence and the index
        # There is also a chance there is a bug in the index

        # NB: Get first element of tuple, then first element in numpy array
        first_problem_kmer_index = (count_list == 0).nonzero()[0][0]

        problem_sequence_index = sequence_indices[first_problem_kmer_index]
        problem_sequence_length = kmer_lengths[first_problem_kmer_index]

        # Assume first sequence in the list of sequence is the problem sequence
        first_problem_kmer = sequences[0][
            problem_sequence_index:
            problem_sequence_index+problem_sequence_length].decode("utf-8")
        # TODO: Check if I can just use ascii above

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
    # TODO: Remove assertion
    assert excess_buffer_length < max_search_length, \
        "Excess sequence buffer length is greater than the maximum search " \
        "length"

    # If all positions non-ambiguous
    # NB: Likely the most common case
    # XXX?: Change condition ordering?
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
        # TODO: Consider putting this in into its own function, e.g.
        # `find_non_ambiguous_intervals`
        # Find the 0-based start-inclusive and end-exclusive intervals (for
        # indexing) of non-ambiguous areas including the start and end in the
        # interval list

        # Compare the ambiguity list with a shifted version of itself
        left_shift = ambiguity_mask[1:]
        reference_mask = ambiguity_mask[:-1]  # NB: Keep same dimensions

        # Where our reference mask is True but left shift is False
        # We have the start of an non-ambiguous region
        # NB: Add 1 due to the shift
        # TODO: Add ascii art
        # TODO: Look at shifting so addition of 1 to indicies is not necessary
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
        # XXX: numpy dstack?
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


# TODO: Change to from `assign` to `set`
def assign_remaining_array_values(input_array: npt.NDArray,
                                  new_values: npt.NDArray):
    """Modifies the input array so the last n values are assigned the contents
    of the new values array regardless of dimension of either array."""
    num_values = min(new_values.size, input_array.size)
    # Only update if there is anything to assign
    if num_values > 0:
        input_array[-num_values:] = new_values[-num_values:]


def print_summary_statisitcs(config: SearchConfig,
                             sequence_id: bytes,
                             total_unique_lengths_count: int,
                             total_ambiguous_positions: int,
                             total_no_unique_lengths_count: int,
                             max_length_found: int,
                             min_length_found: int):

    if (total_unique_lengths_count):
        config.log(f"Finished writing unique lengths for sequence "
                   f"ID: {sequence_id.decode()}")
        config.log(f"{total_unique_lengths_count} unique lengths found")
        config.log(f"{total_ambiguous_positions} positions skipped due to "
                   "ambiguity")
        config.log(f"{total_no_unique_lengths_count} positions with no "
                   "unique length found")
        config.log(f"{max_length_found}-mer maximum unique length found")
        config.log(f"{min_length_found}-mer minimum unique length found")


def main(args):
    config = SearchConfig.from_args(args)
    write_unique_counts(config)
