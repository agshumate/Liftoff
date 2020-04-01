
def get_seq(fasta, genes, index):
    gene = genes[index]
    return ">"+ gene.id + "\n" + gene.sequence(fasta, use_strand=True) + "\n"


def get_all_seqs_for_chrm(gene_db, fasta, chrm_name):
    fasta_out = open(chrm_name + '_genes.fa', 'w')
    if chrm_name != "all":
        genes = list(gene_db.region(seqid = chrm_name, featuretype="gene"))
    else:
        genes = list(gene_db.features_of_type(featuretype="gene"))
    for gene in genes:
        fasta_out.write(">"+ gene.id + "\n" + gene.sequence(fasta, use_strand=True) + "\n")



def get_all_gene_seqs(gene_db, fasta, old_chroms):
    for chrm_name in old_chroms:
        get_all_seqs_for_chrm(gene_db, fasta, chrm_name)