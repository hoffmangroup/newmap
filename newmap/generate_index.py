from pathlib import Path

from newmap._c_newmap_generate_index import generate_fm_index

INDEX_EXTENSION = ".awfmi"


def main(args):
    fasta_filename = args.fasta_file
    index_basename = args.awfmi_base

    # If a basename is provided, use it to create the index filename
    if index_basename:
        index_filename = Path(index_basename).with_suffix(INDEX_EXTENSION)
    # Otherwise use the basename of the fasta file
    else:
        # NB: Doesn't preserve the original directory structure of fasta file
        fasta_basename = Path(fasta_filename).stem
        index_filename = Path(fasta_basename).with_suffix(INDEX_EXTENSION)

    # Convert to python string for C module
    index_filename = str(index_filename)

    suffix_array_compression_ratio = args.compression_ratio
    kmer_length_in_seed_table = args.seed_length

    generate_fm_index(fasta_filename,
                      index_filename,
                      suffix_array_compression_ratio,
                      kmer_length_in_seed_table)
