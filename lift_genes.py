from blast_genes import blast_all_genes
from get_gene_seqs import get_all_gene_seqs
from make_blast_dbs import build_all_databases
from multiprocessing import Pool
from functools import partial
import gffutils
from process_blast_alignments import process_alignments
from build_coordinate_map import build_coordinate_map
from merge_lifted_features import  merge_lifted_features
from find_best_mapping import find_best_mapping
import numpy as np


def extact_and_align_genes(target_fasta, reference_fasta, old_chrms, new_chrms, processes, word_size, gene_db, search_type, build_db):
    if build_db:
        build_all_databases(target_fasta, processes)
    get_all_gene_seqs(gene_db, reference_fasta, old_chrms)
    all_records = blast_all_genes(old_chrms, new_chrms, processes, word_size, search_type)
    return all_records


def build_gene_database(gff):
    gene_db = gffutils.create_db(gff, gff + "_db", merge_strategy="create_unique", force=True)
    return gene_db


def lift_all_genes(processes, db_name, all_records, exclude_criteria, threshold, weight_threshold):
    pool = Pool(processes)
    func = partial(lift_genes, db_name, exclude_criteria, threshold, weight_threshold)
    feature_list = {}
    all_unmapped = []
    blast_records_array = np.array(list(all_records.values()))
    blast_records_split = np.array_split(blast_records_array, processes)
    for result in pool.imap_unordered(func, blast_records_split):
        feature_list.update(result[0])
        all_unmapped.extend(result[1])
    pool.close()
    pool.join()
    return feature_list, all_unmapped



def lift_genes(gene_db_name,  exclude_criteria, threshold, weight_threshold, blast_records):
    gene_db = gffutils.FeatureDB(gene_db_name)
    features = {}
    unmapped_genes = []
    query_num = 0
    for blast_record in blast_records:
        query_num += 1
        #print(query_num)
        new_parent_name = blast_record.query
        #print(new_parent_name)
        copy_tag_len = len(new_parent_name.split("_")[-1])
        original_parent_name = new_parent_name[:-copy_tag_len-1]
        if new_parent_name in exclude_criteria:
            criteria = exclude_criteria[new_parent_name]
        else:
            criteria = []
        parent = gene_db[original_parent_name]
        exon_alignments = process_alignments(blast_record, gene_db, parent)
        if len(exon_alignments) > 0:
            coordinate_map, alignment_scores = build_coordinate_map(exon_alignments, parent, gene_db)
            mapped_exons, shortest_path_weight = find_best_mapping(coordinate_map, alignment_scores, gene_db, parent,
                                                                    exon_alignments, criteria, weight_threshold)
            features[new_parent_name] = merge_lifted_features(mapped_exons, shortest_path_weight, gene_db, parent, unmapped_genes, threshold, new_parent_name)
        else:
            unmapped_genes.append(parent)
    return features, unmapped_genes