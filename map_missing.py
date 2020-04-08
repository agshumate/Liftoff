import lift_genes as lg
from check_homologues import check_homologues
from resolve_overlapping_homologues import resolve_overlapping_homologues
import gffutils

def map_unmapped_genes_agaisnt_all(unmapped_genes, target_fasta, reference_fasta, processes, word_size, final_feature_list, gene_db, full_db_name):
    for feature in final_feature_list:
        for value in final_feature_list[feature]:
            if "Parent" not in value.attributes:
                value.score = -1
    for gene in unmapped_genes:
        for tran in gene_db.children(gene.id, level=1):
            unmapped_genes.append(tran)
    missing_db = gffutils.create_db(unmapped_genes, "missing_db", force=True, merge_strategy="create_unique")
    print("extracting and aligning missing")
    all_missing_records = lg.extact_and_align_genes(target_fasta, reference_fasta, ['all'], ['all'], processes, word_size,
                                                 missing_db, "missing", False)
    print("lifting missing")
    missing_lifted, unmapped_genes = lg.lift_all_genes(processes, "missing_db", all_missing_records, {}, 0.0, None, gene_db)
    final_feature_list.update(missing_lifted)
    print("checking homologues")
    remap_genes = check_homologues(missing_lifted, final_feature_list, full_db_name, processes)
    print("resolving homologues")
    final_features_with_missing = resolve_overlapping_homologues(all_missing_records, final_feature_list, remap_genes, unmapped_genes, full_db_name, processes, 0.0, None, gene_db)
    return final_features_with_missing, unmapped_genes
