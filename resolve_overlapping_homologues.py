from process_blast_alignments import count_overlap
from check_homologues import check_homologues, get_gene_list
import lift_genes as lg


def add_to_exclusion_criteria(remap_genes, remap_gene, exclude_criteria):
    exclude_start, exclude_end, exclude_seq, exclude_strand = remap_genes[remap_gene][0], \
                                                              remap_genes[remap_gene][1], \
                                                              remap_genes[remap_gene][2], \
                                                              remap_genes[remap_gene][3]
    if remap_gene in exclude_criteria:
        exclude_criteria[remap_gene].append((exclude_start, exclude_end,
                                             exclude_seq, exclude_strand))
    else:
        exclude_criteria[remap_gene] = [(exclude_start, exclude_end,
                                         exclude_seq, exclude_strand)]


def resolve_overlapping_homologues(all_records, feature_list, remap_genes,  unmapped_genes, db_name, processes, threshold, weight_threshold):
    exclude_criteria = {}
    while len(remap_genes) > 0:
        blast_records_array = {}
        for remap_gene in remap_genes:
            del feature_list[remap_gene]
            blast_records_array[remap_gene] = all_records[remap_gene]
            add_to_exclusion_criteria(remap_genes, remap_gene, exclude_criteria)
        feature_list_remapped, unmapped_genes_remapped = lg.lift_all_genes(processes, db_name, blast_records_array, exclude_criteria, threshold, weight_threshold)
        feature_list.update(feature_list_remapped)
        unmapped_genes.extend(unmapped_genes_remapped)
        clean_exclude_criteria(feature_list, exclude_criteria)
        remap_genes = check_homologues(feature_list_remapped, feature_list,  db_name, processes)
    return feature_list


def clean_exclude_criteria(lifted_features, exclude_criteria):
    gene_list = get_gene_list(lifted_features)
    gene_dict = {}
    for gene in gene_list:
        if gene.seqid in gene_dict:
            gene_dict[gene.seqid].append(gene)
        else:
            gene_dict[gene.seqid] = [gene]
    for gene_name in exclude_criteria:
        for criteria in exclude_criteria[gene_name]:
            overlapping_features = find_overlapping_features(criteria, gene_dict[criteria[2]])
            if overlapping_features is False:
                exclude_criteria[gene_name].remove(criteria)


def find_overlapping_features(exclude_criteria,  gene_list):
    for gene in gene_list:
        if gene.strand == exclude_criteria[3] and gene.seqid == exclude_criteria[2]:
            if count_overlap(gene.start, gene.end, exclude_criteria[0], exclude_criteria[1]) > 0:
                return True
    return False
