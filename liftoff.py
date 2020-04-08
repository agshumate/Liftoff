#!/usr/bin/env python3

import map_main_annotation as map_main
from write_new_gff import write_new_gff
import argparse
import map_missing
import map_unplaced
import map_copies
import cleanup_intermediate_files as cif


def parse_args():
    parser = argparse.ArgumentParser(description='Lift genes')
    parser.add_argument('-g', required=True, metavar="gff_or_gtf", help="annotation file to lift over in gff or gtf format" )
    parser.add_argument('-t', required=True, metavar = "target_fasta", help="target fasta genome to lift genes to")
    parser.add_argument('-r', required=True, metavar= "reference_fasta", help = "reference fasta genome to lift genes from")
    parser.add_argument('-target_chroms', required=False, metavar = "target_chroms", help="file with name of chromosomes to be lifted to")
    parser.add_argument('-ref_chroms', required=False, metavar = "reference_chroms", help="file with name of chrosomosomes to be lifted from")
    parser.add_argument('-p', required=False, metavar ="num_processess", help= "processes", default=1)
    parser.add_argument('-o', required=False, metavar= "out_file", help="output file", default = 'stdout')
    parser.add_argument('-w', required=False, metavar = "blast_word_size", help="word size for blast step", default =50)
    parser.add_argument('-unplaced', required=False, metavar="unplaced_sequence", help="file with unplaced sequence names to map onto main assembly")
    parser.add_argument('-copy_num',  action='store_true', required=False, help="look for additional copies of genes in the annotation")
    args = parser.parse_args()
    return args

def parse_chrm_files(original_chrms_file, new_chrms_file):
    return [line.rstrip() for line in open(original_chrms_file, 'r').readlines()], [line.rstrip() for line in open(new_chrms_file, 'r').readlines()]



def main():
   #parse args
    args = parse_args()
    gff = args.g
    target_fasta = args.t
    reference_fasta = args.r
    processes =int(args.p)
    output = args.o
    word_size = int(args.w)
    old_chroms_file = args.ref_chroms
    new_chroms_file = args.target_chroms

    #read chroms
    if old_chroms_file is not None and new_chroms_file is not None:
        old_chrms, new_chrms = parse_chrm_files(old_chroms_file, new_chroms_file)
    else:
        old_chrms=['all']
        new_chrms=['all']

    #lift genes
    final_feature_list, unmapped_genes, gene_db =map_main.map_main_annotation(gff, target_fasta, reference_fasta, old_chrms, new_chrms, processes, word_size)
    if len(unmapped_genes) > 0 and new_chrms != ['all']:
        features_with_missing, unmapped_genes = map_missing.map_unmapped_genes_agaisnt_all(unmapped_genes, target_fasta, reference_fasta, processes, word_size, final_feature_list, gene_db,  gff + "_db", new_chrms)
    else:
        features_with_missing = final_feature_list
    if args.unplaced is not None:
        old_chroms = parse_chrm_files(args.unplaced, new_chroms_file)[0]
        features_with_unplaced, unmapped_genes = map_unplaced.map_unplaced_seqs(target_fasta, reference_fasta, processes, word_size, features_with_missing, gene_db, old_chroms, gff + "_db", new_chrms)
    else:
        features_with_unplaced = features_with_missing

    if args.copy_num is not False:
        old_chroms = ['all']
        features_with_copy_nums, unmapped_genes = map_copies.find_extra_copies(target_fasta, reference_fasta, processes, word_size, features_with_unplaced, gene_db, old_chroms, gff + "_db", new_chrms)
    else:
        features_with_copy_nums = features_with_unplaced
    write_new_gff(features_with_copy_nums, output, gene_db)
    cif.cleanup_intermediate_files()


main()




