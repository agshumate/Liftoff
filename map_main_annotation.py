import lift_genes as lg
from check_homologues import check_homologues
from resolve_overlapping_homologues import resolve_overlapping_homologues


def map_main_annotation(gff, target_fasta, reference_fasta, old_chrms, new_chrms, processes, word_size):
    print("building db")
    gene_db = lg.build_gene_database(gff)
    print("extracting and aligning")
    all_records = lg.extact_and_align_genes(target_fasta, reference_fasta, old_chrms, new_chrms, processes, word_size, gene_db,
                                         "chrm_by_chrm", True)
    print("lifting genes")
    feature_list, unmapped_genes = lg.lift_all_genes(processes, gff + "_db", all_records, {}, 0.5, None, gene_db)
    print("checking homologues")
    remap_genes = check_homologues(feature_list, feature_list, gff + "_db", processes)
    print("resolving homologues")
    final_feature_list = resolve_overlapping_homologues(all_records, feature_list, remap_genes, unmapped_genes, gff + "_db",
                                                        processes, 0.5, None, gene_db)
    return final_feature_list, unmapped_genes, gene_db