from newmap._c_newmap_generate_index import generate_fm_index


def main(args):
    fasta_filename = args.fasta_file
    index_filename = args.index_file
    suffix_array_compression_ratio = args.suffix_array_compression_ratio
    kmer_length_in_seed_table = args.kmer_length_in_seed_table

    generate_fm_index(fasta_filename,
                      index_filename,
                      suffix_array_compression_ratio,
                      kmer_length_in_seed_table)
