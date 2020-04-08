import numpy as np


def get_seq(fasta, genes, index):
    gene = genes[index]
    return ">"+ gene.id + "\n" + gene.sequence(fasta, use_strand=True) + "\n"


def get_parent_features(gene_db, chrm_name):
    parent_features = []
    if chrm_name != 'all':
        all_feature_list = list(gene_db.region(seqid=chrm_name))
    else:
        all_feature_list = list(gene_db.all_features())
    feature_ids = set([feature.id for feature in all_feature_list])
    for feature in all_feature_list:
        if ("Parent" not in feature.attributes or feature["Parent"][0] not in feature_ids) and feature.featuretype!="region":
            parent_features.append(feature)
    return parent_features


def get_all_seqs_for_chrm(gene_db,fasta, chrm_name, processes):
    parent_features = np.array(get_parent_features(gene_db, chrm_name))
    file_index = 0
    gene_files = []
    if chrm_name == "all":
        parent_features_split = np.array_split(parent_features, processes)
    else:
        parent_features_split=[parent_features]
    for parents_list in parent_features_split:
        fasta_out = open(chrm_name + "_" + str(file_index) + '_genes.fa', 'w')
        gene_files.append(chrm_name + "_" + str(file_index) + '_genes.fa')
        file_index += 1
        for feature in parents_list:
            fasta_out.write(">"+ feature.id + "\n" + feature.sequence(fasta, use_strand=True) + "\n")
    return gene_files



def get_all_gene_seqs(gene_db, fasta, old_chroms, new_chroms, processes):
    blast_file_pairs = {}
    for i in range (len(old_chroms)):
        gene_files = get_all_seqs_for_chrm(gene_db, fasta, old_chroms[i], processes)
        for file in gene_files:
            blast_file_pairs[file] = new_chroms[i]
    return blast_file_pairs
