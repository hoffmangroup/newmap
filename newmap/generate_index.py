from pathlib import Path

from newmap._c_newmap_generate_index import generate_fm_index
from newmap.util import INDEX_EXTENSION


def main(args):
    fasta_filename = args.fasta_file
    index_filename = args.output

    # If no index filename was specified
    if not index_filename:
        # Use the basename of the fasta file and cwd
        index_filename = Path(fasta_filename).stem + \
                              "." + INDEX_EXTENSION

    suffix_array_compression_ratio = args.compression_ratio
    kmer_length_in_seed_table = args.seed_length

    generate_fm_index(fasta_filename,
                      index_filename,
                      suffix_array_compression_ratio,
                      kmer_length_in_seed_table)
