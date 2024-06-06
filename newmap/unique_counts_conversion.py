from pathlib import Path
import sys
from typing import TextIO

from newmap.util import verbose_print

import numpy as np
import numpy.typing as npt


DEFAULT_KMER_SIZE = 24
CHROMOSOME_FILENAME_DELIMITER = ".unique"

# chr_name, start, end, k-mer length, value
BED_FILE_LINE_FORMAT = "{}\t{}\t{}\tk{}\t{}\t.\n"
WIG_FIXED_STEP_DECLARATION_FORMAT = "fixedStep chrom={} start={} step=1 span=1\n"

STDOUT_FILENAME = "-"


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


def write_single_read_bed(bed_file: TextIO,
                          kmer_length: int,
                          multi_read_mappability: npt.NDArray[np.float64],
                          chr_name: str):
    # NB: Score is only either 0 or 1
    single_read_mappability = np.where(multi_read_mappability > 0.0, 1, 0)

    current_start = 0  # NB: Always 0-based, always start of current interval
    current_position = 0 # NB: Always 0-based, current position in iteration
    current_value = single_read_mappability[current_start]

    # For every value in the single-read mappability array
    # NB: current_position is effectively the current 0-based end position
    for current_position, value in enumerate(single_read_mappability):
        # If the current value is different from the previous
        if value != current_value:
            # We are on the start of a new interval
            # Write out the previous interval
            # NB: End coordinate is 1-based, so the previous end coordinate
            # is the current 0-based position
            previous_end = current_position
            bed_file.write(BED_FILE_LINE_FORMAT.format(chr_name,
                                                       current_start,
                                                       previous_end,
                                                       kmer_length,
                                                       current_value))

            # Update new interval values
            current_start = current_position
            current_value = value

    # Write out the remaining interval if it exists
    if (current_position - current_start) > 0:
        bed_file.write(BED_FILE_LINE_FORMAT.format(chr_name,
                                                   current_start,
                                                   # NB: Convert to 1-based end
                                                   current_position + 1,
                                                   kmer_length,
                                                   current_value))


def write_multi_read_wig(wig_file: TextIO,
                         multi_read_mappability: npt.NDArray[np.float64],
                         chr_name: str):

    # Write out the fixedStep declaration
    wig_file.write(WIG_FIXED_STEP_DECLARATION_FORMAT.format(chr_name, 1))
    for value in multi_read_mappability:
        wig_file.write("{}\n".format(value))


def main(args):
    unique_count_filename = Path(args.unique_count_file)
    kmer_length = args.kmer_length
    single_read_bed_filename = args.single_read_bed_file
    multi_read_wig_filename = args.multi_read_wig_file
    verbose = args.verbose

    # Error if both single-read and multi-read output files are standard output
    if (single_read_bed_filename == STDOUT_FILENAME and
        multi_read_wig_filename == STDOUT_FILENAME):
        raise ValueError("Cannot output both single-read and multi-read files "
                         "to standard output")
    # Error if neither single-read nor multi-read output files are specified
    elif (not single_read_bed_filename and
          not multi_read_wig_filename):
        raise ValueError("Must specify at least one output file")

    # Get the chromosome name from the unique length filename
    # NB: Assume the chromosome name is the the entire string preceding the
    # ".unique*" part of the unique_count_filename (may contain periods)
    file_basename = unique_count_filename.name
    chr_name = \
        file_basename[:file_basename.find(CHROMOSOME_FILENAME_DELIMITER)]

    # NB: The single-read mappability is defined for the entire sequence where
    # a uniquely mappable k-mer would cover. So if a k-mer is uniquely mappable
    # starting at position i, then the single read mappability would be 1 for
    # all positions i to i + kmer_length - 1
    # It follows that the multi-read mappability covers the same positions as
    # the single-read, so any non-zero value would be considered single-read
    # mappable
    verbose_print(verbose,
                  "Calculating mappability regions from minimum unique k-mer "
                  "lengths in file: {}".format(unique_count_filename))

    multi_read_mappability = create_multiread_mappability_from_unique_file(
                             unique_count_filename,
                             kmer_length)

    verbose_print(verbose, "Chromosome size:")
    verbose_print(verbose,
                  "{}\t{}".format(chr_name, multi_read_mappability.shape[0]))

    if single_read_bed_filename:
        verbose_print(verbose, "Writing out single-read mappability regions")

        if single_read_bed_filename == STDOUT_FILENAME:
            write_single_read_bed(sys.stdout,
                                  kmer_length,
                                  multi_read_mappability,
                                  chr_name)
        else:
            with open(single_read_bed_filename, "w") as single_read_bed_file:
                write_single_read_bed(single_read_bed_file,
                                      kmer_length,
                                      multi_read_mappability,
                                      chr_name)

    if multi_read_wig_filename:
        verbose_print(verbose, "Writing out multi-read mappability regions")

        if multi_read_wig_filename == STDOUT_FILENAME:
            write_multi_read_wig(sys.stdout,
                                 multi_read_mappability,
                                 chr_name)
        else:
            with open(multi_read_wig_filename, "w") as multi_read_wig_file:
                write_multi_read_wig(multi_read_wig_file,
                                     multi_read_mappability,
                                     chr_name)
