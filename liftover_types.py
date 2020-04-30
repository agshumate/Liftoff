import lift_features
import fix_overlapping_features
import extract_features
import align_features
import liftoff_utils


def lift_original_annotation(gff, target_fasta, reference_fasta, ref_chroms, target_chroms, processes, db,
                             lifted_feature_list, unmapped_features, infer_transcripts):
    if target_chroms[0] != target_fasta:
        threshold = 0.5
    else:
        threshold = 0
    parent_dict, children_dict, intermediate_dict, feature_db, original_parent_order = extract_features.extract_features_to_lift(
        gff, db, ref_chroms, reference_fasta, processes, infer_transcripts)
    aligned_segments = align_features.align_features_to_target(ref_chroms, target_chroms, processes, target_fasta,
                                                               parent_dict, children_dict, "chrm_by_chrm",
                                                               unmapped_features)
    print("lifting features")
    lift_features.lift_all_features(aligned_segments, {}, threshold, feature_db, parent_dict, children_dict,
                                    intermediate_dict, unmapped_features, lifted_feature_list)
    fix_overlapping_features.fix_incorrectly_overlapping_features(lifted_feature_list, lifted_feature_list, parent_dict,
                                                                  aligned_segments, unmapped_features,
                                                                  threshold, intermediate_dict, children_dict,
                                                                  feature_db, original_parent_order)
    return feature_db, parent_dict, intermediate_dict, children_dict, original_parent_order


def map_unmapped_genes_agaisnt_all(unmapped_features, target_fasta, reference_fasta, ref_chroms, target_chroms,
                                       processes, lifted_feature_list, feature_db, parent_dict, intermediate_dict,
                                       children_dict, parent_order):
    liftoff_utils.clear_scores(lifted_feature_list, parent_dict)
    unmapped_dict = {}
    for feature in unmapped_features:
        unmapped_dict[feature.id]=feature
    extract_features.get_gene_sequences(unmapped_dict, ref_chroms, reference_fasta, processes)
    unmapped_features = []
    aligned_segments=align_features.align_features_to_target(ref_chroms, target_chroms, processes, target_fasta,
                                                               unmapped_dict, children_dict, "missing", unmapped_features)
    lift_features.lift_all_features(aligned_segments, {}, 0.0, feature_db, unmapped_dict, children_dict,
                                    intermediate_dict, unmapped_features, lifted_feature_list)
    fix_overlapping_features.fix_incorrectly_overlapping_features(lifted_feature_list, lifted_feature_list, parent_dict,
                                                                  aligned_segments, unmapped_features, 0.0,
                                                                  intermediate_dict, children_dict, feature_db,
                                                                   parent_order)
    return unmapped_features