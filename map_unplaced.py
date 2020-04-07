import lift_genes as lg
from check_homologues import check_homologues
from resolve_overlapping_homologues import resolve_overlapping_homologues
import gffutils

def map_unplaced_seqs(target_fasta, reference_fasta, processes, word_size, final_feature_list, gene_db, old_chroms, full_db_name):
    for feature in final_feature_list:
        for value in final_feature_list[feature]:
            if "Parent" not in value.attributes:
                value.score = -1
    unplaced_features = []
    for seq in old_chroms:
        for feature in list(gene_db.region(seqid=seq)):
            unplaced_features.append(feature)
    unplaced_db = gffutils.create_db(unplaced_features, "unplaced_db", force=True, merge_strategy="create_unique")
    all_unplaced_records = lg.extact_and_align_genes(target_fasta, reference_fasta, old_chroms, ['all'], processes, word_size,
                                                 unplaced_db, "unplaced", False)
    unplaced_lifted, unmapped_genes = lg.lift_all_genes(processes, "unplaced_db", all_unplaced_records, {}, 0.5, None, gene_db)
    final_feature_list.update(unplaced_lifted)
    remap_genes = check_homologues(unplaced_lifted, final_feature_list, full_db_name, processes)
    final_features_with_unplaced = resolve_overlapping_homologues(all_unplaced_records, final_feature_list,
                                                                 remap_genes, unmapped_genes, full_db_name,
                                                                 processes, 0.5, None, gene_db)
    return final_features_with_unplaced, unmapped_genes