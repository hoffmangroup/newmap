from argparse import ArgumentParser

from newmap import generate_index, unique_counts, unique_counts_conversion

# Defaults for FM-index generation
DEFAULT_INDEX_NAME = "index.awfmi"
DEFAULT_SUFFIX_ARRAY_COMPRESSION_RATIO = 8
DEFAULT_KMER_LENGTH_IN_SEED_TABLE = 12

# Defaults for minimum kmer length counting
DEFAULT_KMER_BATCH_SIZE = 1000000
DEFAULT_THREAD_COUNT = 1
DEFAULT_MINIMUM_KMER_LENGTH = 20
DEFAULT_MAXIMUM_KMER_LENGTH = 200

# Defaults for mappability output
DEFAULT_KMER_SIZE = 24


def parse_subcommands():
    parser = ArgumentParser(
        description="Newmap: A tool for generating mappability "
                    "data for a reference sequences")

    subparsers = parser.add_subparsers(
        title="Newmap subcommands",
        description="Valid subcommands",
        help="see --help on each command for additional information",
        required=True)

    # Create a subparser for the "generate-index" command
    generate_index_parser = subparsers.add_parser(
                            "generate-index",
                            description="Generate a FM index from a fasta "
                                        "file")
    generate_index_parser.set_defaults(func=generate_index.main)

    # TODO: Consider changing to -i and -o for input and output
    generate_index_parser.add_argument(
        "fasta_file",
        help="Filename of input fasta file")

    generate_index_parser.add_argument(
        "--index-file", "-i",
        default=DEFAULT_INDEX_NAME,
        help="Filename to output the reference index file for kmer counting. "
        "Default is {}".format(DEFAULT_INDEX_NAME))

    fm_index_paramater_group = generate_index_parser.add_argument_group(
        "FM-index parameters",
        "Parameters for the FM-index generation")

    fm_index_paramater_group.add_argument(
        "--suffix-array-compression-ratio", "-c",
        type=int,
        default=DEFAULT_SUFFIX_ARRAY_COMPRESSION_RATIO,
        help="Compression ratio for the suffix array to be sampled. "
        "Larger ratios reduces file size and increase average number of "
        "operations per query. "
        "Default is {}.".format(DEFAULT_SUFFIX_ARRAY_COMPRESSION_RATIO))

    fm_index_paramater_group.add_argument(
        "--kmer-length-in-seed-table", "-k",
        type=int,
        default=DEFAULT_KMER_LENGTH_IN_SEED_TABLE,
        help="Length of kmers to memoize in a lookup table to speed up "
        "searches. Each value increase multiplies memory usage by 4. "
        "Default is {}.".format(DEFAULT_KMER_LENGTH_IN_SEED_TABLE))

    # Create a subparser for the "unique-lengths" command
    unique_length_parser = subparsers.add_parser(
                            "unique-lengths",
                            description="Creates a binary file with the "
                                        "unique minimum k-mer length each "
                                        "sequence position from the range of "
                                        "k-mer lengths given")

    unique_length_parser.set_defaults(func=unique_counts.main)

    unique_length_parser.add_argument(
        "kmer_lengths",
        help="Specify k-mer lengths to find unique kmers. "
             "Use a comma seperated list of increasing lengths "
             "or a full inclusive set of lengths seperated by a colon. "
             "Example: 20,24,30 or 20-30.")

    unique_length_parser.add_argument(
        "index_file",
        help="Filename of reference index file to count occurances in")

    unique_length_parser.add_argument(
        "fasta_file",
        help="Filename of (gzipped) fasta file for kmer generation")

    unique_length_parser.add_argument(
        "--kmer-batch-size", "-s",
        default=DEFAULT_KMER_BATCH_SIZE,
        type=int,
        help="Maximum number of kmers to batch per reference sequence from "
             "given fasta file. "
             "Use to control memory usage. "
             "Default is {}" .format(DEFAULT_KMER_BATCH_SIZE))

    unique_length_parser.add_argument(
        "--thread-count", "-t",
        default=DEFAULT_THREAD_COUNT,
        type=int,
        help="Number of threads to parallelize kmer counting. "
             "Default is {}" .format(DEFAULT_THREAD_COUNT))

    unique_length_parser.add_argument(
        "--verbose",
        action="store_true",
        help="Print additional information to standard error",)

    # Create a subparser for the "generate-mappability" command
    generate_mappability_parser = subparsers.add_parser(
      "generate-mappability",
      description="Converts minimum unique kmer length files to mappability "
                  "file output")

    generate_mappability_parser.set_defaults(
        func=unique_counts_conversion.main)

    generate_mappability_parser.add_argument(
        "unique_count_file",
        help="Unique count file to convert to bed file")

    generate_mappability_parser.add_argument(
        "--kmer-length", "-k",
        default=DEFAULT_KMER_SIZE,
        type=int,
        help="Kmer length for mappability file output. Default is {}".format(
            DEFAULT_KMER_SIZE))

    # Add (non-positional) arguments for single-read bed file output
    generate_mappability_parser.add_argument(
        "--single-read-bed-file", "-s",
        help="Filename for single-read mappability BED file output. Use '-' "
             "for standard output.")

    # Add (non-positional) arguments for multi-read wiggle file output
    generate_mappability_parser.add_argument(
        "--multi-read-wig-file", "-m",
        help="Filename for multi-read mappability WIG file output. Use '-' "
             "for standard output.")

    generate_mappability_parser.add_argument(
        "--verbose",
        action="store_true",
        help="Print additional information to standard error",)

    # Parse the arguments
    args = parser.parse_args()

    # Call the function associated with the subcommand
    args.func(args)


if __name__ == "__main__":
    parse_subcommands()
