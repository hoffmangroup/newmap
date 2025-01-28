from pathlib import Path
import sys
from typing import BinaryIO, Union

from newmap.util import verbose_print

import numpy as np
import numpy.typing as npt


DEFAULT_KMER_SIZE = 24
CHROMOSOME_FILENAME_DELIMITER = ".unique"

# chr_name, start, end, k-mer length, value
BED_FILE_LINE_FORMAT = "{}\t{}\t{}\tk{}\t{}\t.\n"
WIG_FIXED_STEP_DECLARATION_FORMAT = \
    "fixedStep chrom={} start={} step=1 span=1\n"

STDOUT_FILENAME = "-"


def create_multiread_mappability_from_unique_file(
     unique_lengths_filename: Path,
     kmer_length: int,
     data_type: Union[np.uint8, np.uint16, np.uint32]):

    # Read the unique k-mer lengths from the unique length file
    unique_kmer_lengths = np.fromfile(str(unique_lengths_filename),
                                      dtype=data_type)
    # Create a boolean array of all values that are not 0 and less than or
    # equal to the k-mer length
    unique_mappability = np.logical_and(unique_kmer_lengths <= kmer_length,
                                        unique_kmer_lengths != 0)

    # Create a zero-array of length of the unique length chromosome file
    multiread_mappability = np.zeros(len(unique_mappability), dtype=np.float64)

    # For every unique k-mer for this length
    unique_mappability_indicies = unique_mappability.nonzero()[0]
    for i in unique_mappability_indicies:
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
    current_position = 0  # NB: Always 0-based, current position in iteration
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


def write_multi_read_wig(wig_file: BinaryIO,
                         multi_read_mappability: npt.NDArray[np.float64],
                         chr_name: str):

    # Write out the fixedStep declaration
    wig_file.write(WIG_FIXED_STEP_DECLARATION_FORMAT
                   .format(chr_name, 1)
                   .encode())
    wig_file.write(b''.join(multi_read_mappability.astype(bytes) + b'\n'))
    wig_file.flush()


def write_mappability_files(unique_count_filenames: list[Path],
                            kmer_length: int,
                            single_read_bed_filename: str,  # Might be stdout
                            multi_read_wig_filename: str,  # Might be stdout
                            verbose: bool):

    # Error if both single-read and multi-read output files are standard output
    if (single_read_bed_filename == STDOUT_FILENAME and
       multi_read_wig_filename == STDOUT_FILENAME):
        raise ValueError("Cannot output both single-read and multi-read files "
                         "to standard output")
    # Error if neither single-read nor multi-read output files are specified
    elif (not single_read_bed_filename and
          not multi_read_wig_filename):
        raise ValueError("Must specify at least one output file")

    # For every unique length file specified
    for unique_count_filename in unique_count_filenames:
        # Get the chromosome name from the unique length filename
        # NB: Assume the chromosome name is the the entire string preceding the
        # ".unique*" part of the unique_count_filename (may contain periods)
        file_basename = unique_count_filename.name
        chr_name = \
            file_basename[:file_basename.find(CHROMOSOME_FILENAME_DELIMITER)]

        # Get the data type from the unique length filename suffix
        data_type_string = unique_count_filename.suffix

        if data_type_string == ".uint8":
            data_type = np.uint8
        elif data_type_string == ".uint16":
            data_type = np.uint16
        elif data_type_string == ".uint32":
            data_type = np.uint32
        else:
            raise ValueError(f"Unknown extension on unique length file: "
                             f"\"{data_type_string}\"")

        # NB: The single-read mappability is defined for the entire sequence
        # where a uniquely mappable k-mer would cover. So if a k-mer is
        # uniquely mappable starting at position i, then the single read
        # mappability would be 1 for all positions i to i + kmer_length - 1
        # It follows that the multi-read mappability covers the same positions
        # as the single-read, so any non-zero value would be considered
        # single-read mappable
        verbose_print(verbose, f"Calculating mappability regions from minimum "
                               f"unique k-mer lengths in file: "
                               f"{unique_count_filename}")

        multi_read_mappability = create_multiread_mappability_from_unique_file(
                                 unique_count_filename,
                                 kmer_length,
                                 data_type)  # type: ignore

        verbose_print(verbose, "Chromosome size:")
        verbose_print(verbose,
                      "{}\t{}".format(chr_name,
                                      multi_read_mappability.shape[0]))

        if single_read_bed_filename:
            verbose_print(verbose, f"Appending single-read mappability "
                                   f"regions to {single_read_bed_filename}")

            if single_read_bed_filename == STDOUT_FILENAME:
                write_single_read_bed(sys.stdout,
                                      kmer_length,
                                      multi_read_mappability,
                                      chr_name)
            else:
                with open(single_read_bed_filename, "a") as \
                          single_read_bed_file:
                    write_single_read_bed(single_read_bed_file,
                                          kmer_length,
                                          multi_read_mappability,
                                          chr_name)

        if multi_read_wig_filename:
            verbose_print(verbose, f"Appending multi-read mappability regions"
                                   f" to {multi_read_wig_filename}")

            if multi_read_wig_filename == STDOUT_FILENAME:
                write_multi_read_wig(sys.stdout.buffer,
                                     multi_read_mappability,
                                     chr_name)
            else:
                with open(multi_read_wig_filename, "ab") as \
                          multi_read_wig_file:
                    write_multi_read_wig(multi_read_wig_file,
                                         multi_read_mappability,
                                         chr_name)


def main(args):
    unique_count_filenames = [Path(filename) for filename in
                              args.unique_count_files]
    kmer_length = args.kmer_length
    single_read_bed_filename = args.single_read_bed_file
    multi_read_wig_filename = args.multi_read_wig_file
    verbose = args.verbose

    write_mappability_files(unique_count_filenames,
                            kmer_length,
                            single_read_bed_filename,
                            multi_read_wig_filename,
                            verbose)
