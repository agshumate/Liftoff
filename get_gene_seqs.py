
def get_seq(fasta, genes, index):
    gene = genes[index]
    return ">"+ gene.id + "\n" + gene.sequence(fasta, use_strand=True) + "\n"


def get_parent_features(gene_db, chrm_name):
    parent_features = []
    if chrm_name != 'all':
        all_feature_list = gene_db.region(seqid=chrm_name)
    else:
        all_feature_list = gene_db.all_features()
    for feature in all_feature_list:
        if len(list(gene_db.parents(feature))) == 0 and feature.featuretype!="region":
            parent_features.append(feature)
    return parent_features


def get_all_seqs_for_chrm(gene_db,fasta, chrm_name):
    fasta_out = open(chrm_name + '_genes.fa', 'w')
    parent_features = get_parent_features(gene_db, chrm_name)
    for feature in parent_features:
        fasta_out.write(">"+ feature.id + "\n" + feature.sequence(fasta, use_strand=True) + "\n")



def get_all_gene_seqs(gene_db, fasta, old_chroms):
    for chrm_name in old_chroms:
        get_all_seqs_for_chrm(gene_db, fasta, chrm_name)