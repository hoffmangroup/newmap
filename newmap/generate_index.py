from argparse import ArgumentParser
import sys

from newmap._c_newmap_generate_index import generate_fm_index

DEFAULT_INDEX_NAME = "index.awfmi"
DEFAULT_SUFFIX_ARRAY_COMPRESSION_RATIO = 8
DEFAULT_KMER_LENGTH_IN_SEED_TABLE = 12


def main():
    args = parse_cmdline()
    fasta_filename = args.fasta_file
    index_filename = args.index_file
    suffix_array_compression_ratio = args.suffix_array_compression_ratio
    kmer_length_in_seed_table = args.kmer_length_in_seed_table

    generate_fm_index(fasta_filename,
                      index_filename,
                      suffix_array_compression_ratio,
                      kmer_length_in_seed_table)


def parse_cmdline():
    parser = ArgumentParser(description="Generate a FM index from a fasta file")

    # TODO: Consider changing to -i and -o for input and output
    parser.add_argument(
        "--fasta-file", "-f",
        required=True,
        help="Filename of input fasta file")

    parser.add_argument(
        "--index-file", "-i",
        default=DEFAULT_INDEX_NAME,
        help="Filename of reference index file for kmer counting. "
        "Default is {}".format(DEFAULT_INDEX_NAME))

    fm_index_paramater_group = parser.add_argument_group(
        "FM index parameters",
        "Parameters for the FM index generation")

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


    return parser.parse_args()


if __name__ == "__main__":
    sys.exit(main())
