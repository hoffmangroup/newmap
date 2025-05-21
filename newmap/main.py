from argparse import ArgumentParser
from importlib.metadata import version
import sys

from newmap import generate_index, unique_counts, unique_counts_conversion
from newmap.unique_counts_conversion import STDOUT_FILENAME

# Will throw PackageNotFoundError if package is not installed
__version__ = version("newmap")

# Defaults for FM-index generation
DEFAULT_COMPRESSION_RATIO = 8
DEFAULT_SEED_LENGTH = 12

# Defaults for minimum kmer length counting
# TODO: Check where/if these are used
DEFAULT_KMER_BATCH_SIZE = 1000000
DEFAULT_THREAD_COUNT = 1
DEFAULT_KMER_SEARCH_RANGE = "20:200"

INDEX_EXTENSION = "awfmi"
FASTA_FILE_METAVAR = "fasta_file"

INDEX_SUBCOMMAND = "index"
UNIQUE_LENGTHS_SUBCOMMAND = "search"
GENERATE_MAPPABILITY_SUBCOMMAND = "track"


def parse_subcommands():
    parser = ArgumentParser(
        description="Newmap: A tool for generating mappability "
                    "data for a reference sequence")

    parser.add_argument(
        "--version",
        action="version",
        version=f"%(prog)s {__version__}")

    subparsers = parser.add_subparsers(
        title="To generate mappability data, the following subcommands must "
              "be run in order",
        metavar="",
        required=True)

    # Create a subparser for the "generate-index" command
    generate_index_parser = subparsers.add_parser(
                            INDEX_SUBCOMMAND,
                            help="Create an FM index from a "
                                 "sequence to generate mappability data for.")
    generate_index_parser.set_defaults(func=generate_index.main)

    # TODO: Consider changing to -i and -o for input and output
    generate_index_parser.add_argument(
        "fasta_file",
        help="Filename of reference sequence in Fasta format")

    generate_index_parser.add_argument(
        "--awfmi-base", "-i",
        metavar="BASENAME",
        help="Basename of the index file to write. The default is the "
        f"basename of the input fasta file with the {INDEX_EXTENSION} "
        "extension.")

    fm_index_paramater_group = generate_index_parser.add_argument_group(
        "Indexing parameters",
        "Parameters for the FM-index generation")

    fm_index_paramater_group.add_argument(
        "--compression-ratio", "-c",
        type=int,
        default=DEFAULT_COMPRESSION_RATIO,
        metavar="RATIO",
        help="Compression ratio for the suffix array to be sampled. "
        "Larger ratios reduce file size and increase the average number of "
        "operations per query. "
        "Default is {}.".format(DEFAULT_COMPRESSION_RATIO))

    fm_index_paramater_group.add_argument(
        "--seed-length", "-s",
        type=int,
        default=DEFAULT_SEED_LENGTH,
        metavar="LENGTH",
        help="Length of k-mers to memoize in a lookup table to speed up "
        "searches. Each value increase multiplies memory usage of the index "
        "by 4. "
        "Default is {}.".format(DEFAULT_SEED_LENGTH))

    # Create a subparser for the "search" command
    unique_length_parser = subparsers.add_parser(
                            UNIQUE_LENGTHS_SUBCOMMAND,
                            help="Finds the shortest unique sequence length "
                                 "at each position in the input fasta file. "
                                 "Saves the results to a binary file "
                                 "containing an array with the shortest "
                                 "lengths found within the search range.")

    unique_length_parser.set_defaults(func=unique_counts.main)

    unique_length_parser.add_argument(
        "fasta_file",
        metavar=FASTA_FILE_METAVAR,
        help="Filename of (gzipped) fasta file for kmer generation")

    unique_length_parser.add_argument(
        "index_file",
        nargs="?",
        help="Filename of reference index file to count occurances in. "
             f"Defaults to the basename of the {FASTA_FILE_METAVAR} with "
             f"the {INDEX_EXTENSION} extension.")

    unique_length_output_parameter_group = \
        unique_length_parser.add_argument_group(
            "Output parameters")

    unique_length_output_parameter_group.add_argument(
        "--search-range", "-r",
        metavar="RANGE",
        default=DEFAULT_KMER_SEARCH_RANGE,
        help="Search set of sequence lengths to determine uniqueness. "
             "Use a comma separated list of increasing lengths "
             "or a full inclusive set of lengths separated by a colon. "
             "Examples: 20,24,30 or 20:30. "
             f"Default is {DEFAULT_KMER_SEARCH_RANGE}.")

    unique_length_output_parameter_group.add_argument(
        "--output-directory", "-o",
        metavar="DIR",
        default=".",
        help="Directory to write the binary files containing the 'unique' "
             "lengths to. Default is the current working directory.")

    unique_length_output_parameter_group.add_argument(
        "--include-sequences", "-i",
        metavar="IDS",
        help="A comma separated list of sequence IDs to select from the given "
             "fasta file. Default is to use all sequences when not specified. "
             "Cannot be used with --exclude-sequences.")

    unique_length_output_parameter_group.add_argument(
        "--exclude-sequences", "-x",
        metavar="IDS",
        help="A comma separated list of sequence IDs to exclude from the "
             "given fasta file. Default is to use all sequences when not "
             "specified. Cannot be used with --include-sequences.")

    unique_length_output_parameter_group.add_argument(
        "--verbose", "-v",
        action="store_true",
        help="Print additional information to standard error",)

    unique_length_performance_parameter_group = \
        unique_length_parser.add_argument_group(
            "Performance parameters")

    unique_length_performance_parameter_group.add_argument(
        "--initial-search-length", "-l",
        type=int,
        metavar="LENGTH",
        default=0,
        help="Specify the initial search length. Only valid "
             "when the search range is a continuous range separated by a "
             "colon. Defaults to the midpoint of the range."
    )

    unique_length_performance_parameter_group.add_argument(
        "--kmer-batch-size", "-s",
        default=DEFAULT_KMER_BATCH_SIZE,
        metavar="SIZE",
        type=int,
        help="Maximum number of k-mers to batch per reference sequence from "
             "input fasta file. "
             "Use to control memory usage. "
             "Default is {}".format(DEFAULT_KMER_BATCH_SIZE))

    unique_length_performance_parameter_group.add_argument(
        "--num-threads", "-t",
        default=DEFAULT_THREAD_COUNT,
        metavar="NUM",
        type=int,
        help="Number of threads to parallelize kmer counting. "
             "Default is {}".format(DEFAULT_THREAD_COUNT))

    # Create a subparser for the "generate-mappability" command
    generate_mappability_parser = subparsers.add_parser(
      GENERATE_MAPPABILITY_SUBCOMMAND,
      help="Converts binary array files of unique sequence lengths to "
           "mappability file track(s) for a specific read length")

    generate_mappability_parser.set_defaults(
        func=unique_counts_conversion.main)

    generate_mappability_parser.add_argument(
        "kmer_length",
        type=int,
        help="Kmer length for mappability file output.")

    generate_mappability_parser.add_argument(
        "unique_count_files",
        nargs="+",  # NB: One or more unique files
        help="One or more unique count files to convert to mappability "
             "file(s)")

    # Add (non-positional) arguments for single-read bed file output
    generate_mappability_parser.add_argument(
        "--single-read-bed-file", "-s",
        help="Filename for single-read mappability BED file output. Use '{}' "
             "for standard output.".format(STDOUT_FILENAME))

    # Add (non-positional) arguments for multi-read wiggle file output
    generate_mappability_parser.add_argument(
        "--multi-read-wig-file", "-m",
        help="Filename for multi-read mappability WIG file output. Use '{}' "
             "for standard output.".format(STDOUT_FILENAME))

    generate_mappability_parser.add_argument(
        "--verbose", "-v",
        action="store_true",
        help="Print additional information to standard error",)

    # If there are no arguments, print the help message
    if len(sys.argv) == 1:
        parser.print_help()
    # Otherwise
    else:
        # Parse the arguments
        args = parser.parse_args()
        # Call the function associated with the subcommand
        args.func(args)


if __name__ == "__main__":
    parse_subcommands()
