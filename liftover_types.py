import lift_features
import fix_overlapping_features
import extract_features
import align_features
import liftoff_utils


def lift_original_annotation(gff, target_fasta, reference_fasta, ref_chroms, target_chroms, processes, db,
                             lifted_feature_list, unmapped_features, infer_transcripts, infer_genes, cov_threshold, seq_threshold,
                             minimap2_path, inter_files):
    if target_chroms[0] == target_fasta:
        cov_threshold, seq_threshold = 0,0
    parent_dict, children_dict, intermediate_dict, feature_db, original_parent_order = extract_features.extract_features_to_lift(
        gff, db, ref_chroms, reference_fasta, processes, infer_transcripts, infer_genes, inter_files)
    aligned_segments = align_features.align_features_to_target(ref_chroms, target_chroms, processes, target_fasta,
                                                               parent_dict, children_dict, "chrm_by_chrm",
                                                               unmapped_features, reference_fasta, minimap2_path, inter_files, True)

    print("lifting features")
    lift_features.lift_all_features(aligned_segments, {}, cov_threshold, feature_db, parent_dict, children_dict,
                                    intermediate_dict, unmapped_features, lifted_feature_list, seq_threshold)
    print("fix homologs")
    fix_overlapping_features.fix_incorrectly_overlapping_features(lifted_feature_list, lifted_feature_list, parent_dict,
                                                                  aligned_segments, unmapped_features,
                                                                  cov_threshold, intermediate_dict, children_dict,
                                                                  feature_db, original_parent_order, seq_threshold, "chrm_by_chrm")
    return feature_db, parent_dict, intermediate_dict, children_dict, original_parent_order


def map_unmapped_genes_agaisnt_all(unmapped_features, target_fasta, reference_fasta, ref_chroms, target_chroms,
                                       processes, lifted_feature_list, feature_db, parent_dict, intermediate_dict,
                                       children_dict, parent_order,minimap2_path,inter_files):
    liftoff_utils.clear_scores(lifted_feature_list, parent_dict)
    unmapped_dict = {}
    for feature in unmapped_features:
        unmapped_dict[feature.id]=feature
    extract_features.get_gene_sequences(unmapped_dict, ref_chroms, reference_fasta, processes, inter_files)
    unmapped_features = []
    aligned_segments=align_features.align_features_to_target(ref_chroms, target_chroms, processes, target_fasta,
                                                               unmapped_dict, children_dict, "missing", unmapped_features, reference_fasta,
                                                             minimap2_path,inter_files, True)
    print("lifting features")
    lift_features.lift_all_features(aligned_segments, {}, 0.0, feature_db, unmapped_dict, children_dict,
                                    intermediate_dict, unmapped_features, lifted_feature_list, 0.0)
    fix_overlapping_features.fix_incorrectly_overlapping_features(lifted_feature_list, lifted_feature_list, parent_dict,
                                                                  aligned_segments, unmapped_features, 0.0,
                                                                  intermediate_dict, children_dict, feature_db,
                                                                   parent_order, 0.0, "missing")
    return unmapped_features


def map_unplaced_genes(unmapped_features, target_fasta, reference_fasta, ref_chroms, target_chroms,
                                       processes, lifted_feature_list, feature_db, parent_dict, intermediate_dict,
                                       children_dict, parent_order,minimap2_path,inter_files):
    liftoff_utils.clear_scores(lifted_feature_list, parent_dict)
    unplaced_dict = {}
    for feature_name in parent_dict:
        feature = parent_dict[feature_name]
        if feature.seqid in ref_chroms:
            unplaced_dict[feature.id] = feature
    extract_features.get_gene_sequences(unplaced_dict, ref_chroms, reference_fasta, processes, inter_files)
    aligned_segments=align_features.align_features_to_target(ref_chroms, target_chroms, processes, target_fasta,
                                                               unplaced_dict, children_dict, "unplaced", unmapped_features,
                                                             reference_fasta,minimap2_path,inter_files, True)
    print("lifting features")
    lift_features.lift_all_features(aligned_segments, {}, 0.0, feature_db, unplaced_dict, children_dict,
                                    intermediate_dict, unmapped_features, lifted_feature_list, 0.0)

    fix_overlapping_features.fix_incorrectly_overlapping_features(lifted_feature_list, lifted_feature_list, parent_dict,
                                                                  aligned_segments, unmapped_features, 0.0,
                                                                  intermediate_dict, children_dict, feature_db, parent_order, 0.0, "unplaced")



def map_extra_copies(target_fasta, reference_fasta, ref_chroms, target_chroms, processes,
                             lifted_feature_list, parent_dict, children_dict, feature_db, intermediate_dict, parent_order,
                     seq_threshold,minimap2_path,inter_files, remap):
    liftoff_utils.clear_scores(lifted_feature_list, parent_dict)
    unmapped_features = []
    extract_features.get_gene_sequences(parent_dict, ref_chroms, reference_fasta, processes, inter_files)
    aligned_segments=align_features.align_features_to_target(ref_chroms, target_chroms, processes, target_fasta,
                                                               parent_dict, children_dict, "copies", unmapped_features, reference_fasta,
                                                             minimap2_path,inter_files, remap
                                                             )

    print("lifting features")
    lift_features.lift_all_features(aligned_segments, {}, 0.0, feature_db, parent_dict, children_dict,
                                    intermediate_dict, unmapped_features, lifted_feature_list, seq_threshold)
    fix_overlapping_features.fix_incorrectly_overlapping_features(lifted_feature_list, lifted_feature_list, parent_dict,
                                                                  aligned_segments, unmapped_features, 0.0,
                                                                  intermediate_dict, children_dict, feature_db, parent_order, seq_threshold, "copies")

