from pathlib import Path
import sys
from typing import BinaryIO, Union

from newmap.util import DEFAULT_MAPPABILITY_READ_LENGTH, verbose_print

import numpy as np
import numpy.typing as npt


STDOUT_FILENAME = "-"
CHROMOSOME_FILENAME_DELIMITER = ".unique"

# chr_name, start, end, k-mer length, value
BED_FILE_LINE_FORMAT = "{}\t{}\t{}\tk{}\t{}\t.\n"
WIG_FIXED_STEP_DECLARATION_FORMAT = \
    "fixedStep chrom={} start={} step=1 span=1\n"

MULTIREAD_MAPPABILITY_TYPE = np.float64


def create_multiread_mappability_from_unique_file(
     unique_lengths_filename: Path,
     kmer_length: int,
     data_type: Union[np.uint8, np.uint16, np.uint32]):

    # Read the unique k-mer lengths from the unique length file
    unique_kmer_lengths = np.fromfile(str(unique_lengths_filename),
                                      dtype=data_type)

    # NB: This function aims to use numpy's cumulative sum function (cumsum)
    # in order to sum overlapping intervals (because iterating and summing over
    # each k length interval can be very slow)
    # The way this is done in principle is to sum over an array where there are
    # all zeros by default, and for every interval there is a 1 to mark the
    # start and a -1 to mark the end
    # Generally, instead of simply assigning a 1 or -1, these values must
    # added to the current value to account for intervals where starts and ends
    # can overlap

    # NB: Can't just use the default unsigned types, because there may be
    # negative values after subtracting out the end intervals
    # Set a "1" at every position where a interval would start
    interval_marks = (
        (unique_kmer_lengths <= kmer_length) &
        (unique_kmer_lengths != 0)
    ).astype(MULTIREAD_MAPPABILITY_TYPE)

    # And then kmer_length + 1 places away, subtract "1"
    interval_marks[kmer_length:] -= interval_marks[:-kmer_length]

    # Take the cumulative sum to get a count of overlapping intervals
    multiread_mappability = np.cumsum(interval_marks)
    # Get the fraction of overlapping intervals by the length k
    # NB: Cannot sum the fractions, must divide at the end. Would not be the
    # exact same results as Umap had
    multiread_mappability /= kmer_length

    return multiread_mappability


def write_single_read_bed(bed_file: BinaryIO,
                          kmer_length: int,
                          multi_read_mappability: npt.NDArray[np.float64],
                          chr_name: str):
    # NB: Score is only either 0 or 1
    single_read_mappability = np.where(multi_read_mappability > 0.0, 1, 0)

    # Create a mask marking "True" when were are on a new interval
    change_mask = np.concatenate(
         ([True],  # First value is assumed to be a new set
          single_read_mappability[1:] !=  # Shift and compare
          single_read_mappability[:-1]))

    # Get indices where changes in values occur (from 1 to 0 or vice-versa)
    change_start_indices = np.nonzero(change_mask)[0]

    # Get the length of each interval (difference in change indices)
    interval_lengths = np.diff(
        np.append(change_start_indices, single_read_mappability.size)
    )

    # Write out entire set of formatted byte strings
    # BED_FILE_LINE_FORMAT = "{}\t{}\t{}\tk{}\t{}\t.\n"
    bed_file.write(b''.join(
        [f"{chr_name}\t"
         f"{index}\t"
         f"{index+length}\t"
         f"k{kmer_length}\t"
         f"{single_read_mappability[index]}"
         "\t.\n".encode()
         for index, length in zip(change_start_indices, interval_lengths)]
    ))


def write_multi_read_wig(wig_file: BinaryIO,
                         multi_read_mappability: npt.NDArray[np.float64],
                         chr_name: str,
                         decimal_places: int = 2):

    # Write out the fixedStep declaration
    wig_file.write(WIG_FIXED_STEP_DECLARATION_FORMAT
                   .format(chr_name, 1)
                   .encode())

    for mappability_chunk in np.nditer(multi_read_mappability,
                                       flags=['external_loop', 'buffered']):
        wig_file.write(
            b'\n'.join(float_format(value, decimal_places)  # type: ignore
                       for value in mappability_chunk) + b'\n'
        )


def float_format(value: float,
                 decimal_places: int) -> bytes:
    # Special case for 0, ignore trailing decimal places
    if value == 0.0:
        return b"0.0"

    # Otherwise print the specified number of decimal places
    return f"{value:.{decimal_places}f}".encode()


def safe_remove(filename: str):
    if (filename and
       filename != STDOUT_FILENAME and
       Path(filename).exists()):
        Path(filename).unlink()


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

    # Delete any existing mappability files if they exist
    safe_remove(single_read_bed_filename)
    safe_remove(multi_read_wig_filename)

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
            filename = single_read_bed_filename
            if filename == STDOUT_FILENAME:
                filename = "standard output"

            verbose_print(verbose, f"Appending single-read mappability "
                                   f"regions to {filename}")

            if single_read_bed_filename == STDOUT_FILENAME:
                write_single_read_bed(sys.stdout.buffer,
                                      kmer_length,
                                      multi_read_mappability,
                                      chr_name)
            else:
                with open(single_read_bed_filename, "ab") as \
                          single_read_bed_file:
                    write_single_read_bed(single_read_bed_file,
                                          kmer_length,
                                          multi_read_mappability,
                                          chr_name)

        if multi_read_wig_filename:
            filename = multi_read_wig_filename
            if filename == STDOUT_FILENAME:
                filename = "standard output"

            verbose_print(verbose, f"Appending multi-read mappability regions"
                                   f" to {filename}")

            # Calculate the number of decimal places
            decimal_places = int(np.ceil(np.log10(kmer_length)))

            if multi_read_wig_filename == STDOUT_FILENAME:
                write_multi_read_wig(sys.stdout.buffer,
                                     multi_read_mappability,
                                     chr_name,
                                     decimal_places)
            else:
                with open(multi_read_wig_filename, "ab") as \
                          multi_read_wig_file:
                    write_multi_read_wig(multi_read_wig_file,
                                         multi_read_mappability,
                                         chr_name,
                                         decimal_places)


def main(args):
    unique_count_filenames = [Path(filename) for filename in
                              args.unique_count_files]
    kmer_length = args.read_length
    single_read_bed_filename = args.single_read
    multi_read_wig_filename = args.multi_read
    verbose = args.verbose

    # Check to see if the read/k-mer length is valid or specified
    # If not
    if not kmer_length.isdigit():
        # Assume the read/k-mer length on the command line is a unique filename
        unique_count_filenames.insert(0, Path(kmer_length))
        # Set the k-mer length to the default value
        kmer_length = DEFAULT_MAPPABILITY_READ_LENGTH
    else:
        # Convert the k-mer length to an integer
        kmer_length = int(kmer_length)

    # If neither single-read nor multi-read output files are specified
    # Default to single-read standard output
    if (not single_read_bed_filename and
       not multi_read_wig_filename):
        single_read_bed_filename = STDOUT_FILENAME

    write_mappability_files(unique_count_filenames,
                            kmer_length,
                            single_read_bed_filename,
                            multi_read_wig_filename,
                            verbose)
