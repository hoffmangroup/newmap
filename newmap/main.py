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
DEFAULT_MINIMUM_KMER_LENGTH = 20
DEFAULT_MAXIMUM_KMER_LENGTH = 200

# Defaults for mappability output
DEFAULT_KMER_SIZE = 24

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
        help="Basename of the index file to write. The default is the "
        "basename of the input fasta file with the awfmi extension.")

    fm_index_paramater_group = generate_index_parser.add_argument_group(
        "Indexing parameters",
        "Parameters for the FM-index generation")

    fm_index_paramater_group.add_argument(
        "--compression-ratio", "-c",
        type=int,
        default=DEFAULT_COMPRESSION_RATIO,
        help="Compression ratio for the suffix array to be sampled. "
        "Larger ratios reduce file size and increase the average number of "
        "operations per query. "
        "Default is {}.".format(DEFAULT_COMPRESSION_RATIO))

    fm_index_paramater_group.add_argument(
        "--seed-length", "-s",
        type=int,
        default=DEFAULT_SEED_LENGTH,
        help="Length of k-mers to memoize in a lookup table to speed up "
        "searches. Each value increase multiplies memory usage of the index "
        "by 4. "
        "Default is {}.".format(DEFAULT_SEED_LENGTH))

    # Create a subparser for the "search" command
    unique_length_parser = subparsers.add_parser(
                            UNIQUE_LENGTHS_SUBCOMMAND,
                            help="Finds the shortest unique sequence length "
                                 "at each position in the input fasta file. "
                                 "Saves the results to a binary array.")

    unique_length_parser.set_defaults(func=unique_counts.main)

    unique_length_parser.add_argument(
        "kmer_lengths",
        help="Specify k-mer lengths to find unique k-mers. "
             "Use a comma separated list of increasing lengths "
             "or a full inclusive set of lengths separated by a colon. "
             "Example: 20,24,30 or 20:30.")

    unique_length_parser.add_argument(
        "index_file",
        help="Filename of reference index file to count occurances in")

    unique_length_parser.add_argument(
        "fasta_file",
        help="Filename of (gzipped) fasta file for kmer generation")

    unique_length_parser.add_argument(
        "--initial-search-length", "-l",
        type=int,
        default=0,
        help="Specify the initial search length for unique k-mers. Only valid "
             "when the search range is a continuous range separated by a "
             "colon."
    )

    unique_length_parser.add_argument(
        "--include-sequences", "-i",
        help="A comma separated list of sequence IDs to select from the given "
             "fasta file. Default is to use all sequences when not specified. "
             "Cannot be used with --exclude-sequences.")

    unique_length_parser.add_argument(
        "--exclude-sequences", "-x",
        help="A comma separated list of sequence IDs to exclude from the "
             "given fasta file. Default is to use all sequences when not "
             "specified. Cannot be used with --include-sequences.")

    unique_length_parser.add_argument(
        "--kmer-batch-size", "-s",
        default=DEFAULT_KMER_BATCH_SIZE,
        type=int,
        help="Maximum number of k-mers to batch per reference sequence from "
             "input fasta file. "
             "Use to control memory usage. "
             "Default is {}".format(DEFAULT_KMER_BATCH_SIZE))

    unique_length_parser.add_argument(
        "--thread-count", "-t",
        default=DEFAULT_THREAD_COUNT,
        type=int,
        help="Number of threads to parallelize kmer counting. "
             "Default is {}".format(DEFAULT_THREAD_COUNT))

    unique_length_parser.add_argument(
        "--verbose",
        action="store_true",
        help="Print additional information to standard error",)

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
        "--verbose",
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
