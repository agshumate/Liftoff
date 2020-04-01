import lift_genes as lg
from resolve_overlapping_homologues import resolve_overlapping_homologues
from check_homologues import check_homologues

def find_extra_copies(target_fasta, reference_fasta, processes, word_size, final_feature_list, gene_db, old_chroms, full_db_name):
    for feature in final_feature_list:
        for value in final_feature_list[feature]:
            if value.featuretype == "gene":
                value.score = -1
    all_extra_copies = lg.extact_and_align_genes(target_fasta, reference_fasta, old_chroms, ['all'], processes, word_size,
                                                 gene_db, "copies", False)
    copies_lifted, unmapped_genes = lg.lift_all_genes(processes, full_db_name, all_extra_copies, {}, 1.0, 0.0)
    final_feature_list.update(copies_lifted)
    remap_genes = check_homologues(copies_lifted, final_feature_list, full_db_name, processes)
    final_features_with_copies = resolve_overlapping_homologues(all_extra_copies, final_feature_list,
                                                                  remap_genes,  unmapped_genes, full_db_name,
                                                                  processes, 1.0, 0.0)
    return final_features_with_copies, unmapped_genes
