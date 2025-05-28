from newmap._c_newmap_generate_index import generate_fm_index


def main(args):
    fasta_filename = args.fasta_file
    index_filename = args.output

    # Convert to python string for C module
    index_filename = str(index_filename)

    suffix_array_compression_ratio = args.compression_ratio
    kmer_length_in_seed_table = args.seed_length

    generate_fm_index(fasta_filename,
                      index_filename,
                      suffix_array_compression_ratio,
                      kmer_length_in_seed_table)
