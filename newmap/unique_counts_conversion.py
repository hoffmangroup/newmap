from pathlib import Path

import numpy as np

DEFAULT_KMER_SIZE = 24
CHROMOSOME_FILENAME_DELIMITER = ".unique"
# chr_name, start, end, k-mer length, value
BED_FILE_LINE_FORMAT = "{}\t{}\t{}\tk{}\t{}\t.\n"

def create_multiread_mappability_from_unique_file(
     unique_lengths_filename: Path,
     kmer_length: int):

    # Read the unique k-mer lengths from the unique length file
    unique_kmer_lengths = np.fromfile(str(unique_lengths_filename),
                                      dtype=np.uint8)
    # Create a boolean array of all values that are not 0 and less than or
    # equal to the k-mer length
    unique_mappability = np.logical_and(unique_kmer_lengths <= kmer_length,
                                        unique_kmer_lengths != 0)

    # Create a zero-array of length of the unique length chromosome file
    multiread_mappability = np.zeros(len(unique_mappability), dtype=np.float64)
    # For each position boolean unique mappability array
    for i, is_unique in enumerate(unique_mappability):
        # If the k-mer is unique for this length
        if is_unique:
            # Add 1 to the multiread mappability array for this position
            # and all positions up to the k-mer length
            multiread_mappability[i:i + kmer_length] += 1.0

    multiread_mappability /= kmer_length
    return multiread_mappability


def write_single_read_bed(kmer_length: int,
                          unique_count_filename: Path,
                          verbose: bool):

    # Get the chromosome name from the unique length filename
    # NB: Assume the chromosome name is the the entire string preceding the
    # ".unique*" part of the unique_count_filename (may contain periods)
    file_basename = unique_count_filename.name
    chr_name = \
        file_basename[:file_basename.find(CHROMOSOME_FILENAME_DELIMITER)]

    # single_read_mappability = np.where(unique_kmer_lengths == kmer_length, 1, 0)

    # NB: The single-read mappability is defined for the entire sequence where
    # a uniquely mappable k-mer would cover. So if a k-mer is uniquely mappable
    # starting at position i, then the single read mappability would be 1 for
    # all positions i to i + kmer_length - 1
    # It follows that the multi-read mappability covers the same positions as
    # the single-read, so any non-zero value would be considered single-read
    # mappable
    multiread_mappability = create_multiread_mappability_from_unique_file(
                            unique_count_filename,
                            kmer_length)
    # NB: Score is only either 0 or 1
    single_read_mappability = np.where(multiread_mappability > 0.0, 1, 0)

    current_start = 0
    current_end = 0
    current_value = single_read_mappability[current_start]

    # TODO: Write out 0-values as well
    # TODO: Add trackline options
    # bed_file.write("track name=\"<type 'type'> k24\"description=\"Single-read mappability with k24-merscolor=%s\n")
    with open("single_read_mappability.bed", "w") as bed_file:
        # For each element in the single read mappability array
        for i, value in enumerate(single_read_mappability):
            # If we have single read mappability scores
            if value > 0:  # value == 1
                # If we are at the start of a new interval
                if current_start == current_end:
                    current_start = i
                    current_end = i + 1
                else:
                    # Otherwise continue merging
                    current_end += 1
            # Otherwise we stopped finding single read mappability
            else:  # value == 0.0
                # Stop merging
                # Write out previous interval if it exists
                # TODO: Move formatting string to constant
                # chr_name, start, end, k-mer length, value
                # BED_FILE_LINE_FORMAT = "{}\t{}\t{}\tk{}\t{}\t.\n"
                if (current_end - current_start) > 0:
                    bed_file.write(BED_FILE_LINE_FORMAT.format(chr_name,
                                                               current_start,
                                                               current_end,
                                                               kmer_length, 1))
                    current_start = current_end

        # Write out the remaining interval if it exists
        if (current_end - current_start) > 0:
            bed_file.write(BED_FILE_LINE_FORMAT.format(chr_name,
                                                       current_start,
                                                       current_end,
                                                       kmer_length, 1))
            current_start = current_end


def main(args):

    kmer_length = args.kmer_length
    unique_count_filename = Path(args.unique_count_file)
    verbose = args.verbose

    write_single_read_bed(kmer_length,
                          Path(unique_count_filename),
                          verbose)
