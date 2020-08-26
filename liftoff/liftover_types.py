from liftoff import fix_overlapping_features, lift_features, liftoff_utils, align_features, extract_features


def lift_original_annotation(ref_chroms, target_chroms, lifted_features_list, args, unmapped_features, parents_to_lift):
    liftover_type = "chrm_by_chrm"
    if target_chroms[0] == args.target and args.exclude_partial == False:
        min_cov, min_seqid = 0.05, 0.05
    else:
        min_cov, min_seqid = args.a, args.s

    feature_hierarchy, feature_db, ref_parent_order = extract_features.extract_features_to_lift(ref_chroms,
                                                                                                liftover_type,
                                                                                                parents_to_lift, args)
    align_and_lift_features(ref_chroms, target_chroms, args, feature_hierarchy, liftover_type, unmapped_features,
                            feature_db,
                            feature_hierarchy.parents, lifted_features_list, ref_parent_order, min_cov, min_seqid,
                            args.overlap)
    return feature_db, feature_hierarchy, ref_parent_order


def align_and_lift_features(ref_chroms, target_chroms, args, feature_hierarchy, liftover_type, unmapped_features,
                            feature_db,
                            features_to_lift, lifted_features_list, ref_parent_order, min_cov, min_seqid, max_overlap):
    aligned_segments = align_features.align_features_to_target(ref_chroms, target_chroms, args, feature_hierarchy,
                                                               liftover_type, unmapped_features)
    print("lifting features")
    feature_locations = None
    lift_features.lift_all_features(aligned_segments, min_cov, feature_db, features_to_lift, feature_hierarchy,
                                    unmapped_features, lifted_features_list, min_seqid, feature_locations, args.d,
                                    ref_parent_order)
    fix_overlapping_features.fix_incorrectly_overlapping_features(lifted_features_list, lifted_features_list,
                                                                  aligned_segments, unmapped_features,
                                                                  min_cov, feature_hierarchy,
                                                                  feature_db, ref_parent_order, min_seqid, args.d,
                                                                   max_overlap)


def map_unmapped_genes_agaisnt_all(unmapped_features, ref_chroms, target_chroms, lifted_features_list, feature_db,
                                   feature_hierarchy, ref_parent_order, args):
    liftoff_utils.clear_scores(lifted_features_list, feature_hierarchy.parents)
    unmapped_dict = get_unmapped_genes(unmapped_features)
    if args.exclude_partial:
        min_cov, min_seqid = args.a, args.s
    else:
        min_cov, min_seqid = 0.05, 0.05
    liftover_type = "unmapped"
    extract_features.get_gene_sequences(unmapped_dict, ref_chroms, args, liftover_type)
    unmapped_features = []
    align_and_lift_features(ref_chroms, target_chroms, args, feature_hierarchy, liftover_type, unmapped_features,
                            feature_db,
                            unmapped_dict, lifted_features_list, ref_parent_order, min_cov, min_seqid, args.overlap)
    return unmapped_features


def get_unmapped_genes(unmapped_features):
    unmapped_dict = {}
    for feature in unmapped_features:
        unmapped_dict[feature.id] = feature
    return unmapped_dict


def map_unplaced_genes(unmapped_features, ref_chroms, target_chroms,
                       lifted_features_list, feature_db, feature_hierarchy, ref_parent_order, args):
    liftoff_utils.clear_scores(lifted_features_list, feature_hierarchy.parents)
    liftover_type = "unplaced"
    unplaced_dict = get_features_from_unplaced_seq(ref_chroms, feature_hierarchy)
    extract_features.get_gene_sequences(unplaced_dict, ref_chroms, args, liftover_type)
    if args.exclude_partial:
        min_cov, min_seqid = args.a, args.s
    else:
        min_cov, min_seqid = 0.05, 0.05
    align_and_lift_features(ref_chroms, target_chroms, args, feature_hierarchy, liftover_type, unmapped_features,
                            feature_db, unplaced_dict, lifted_features_list, ref_parent_order, min_cov, min_seqid,
                            args.overlap)


def get_features_from_unplaced_seq(ref_chroms, feature_hierarchy):
    unplaced_dict = {}
    for feature_name in feature_hierarchy.parents:
        feature = feature_hierarchy.parents[feature_name]
        if feature.seqid in ref_chroms:
            unplaced_dict[feature.id] = feature
    return unplaced_dict


def map_extra_copies(ref_chroms, target_chroms, lifted_features_list, feature_hierarchy, feature_db, ref_parent_order,
                     args):
    liftoff_utils.clear_scores(lifted_features_list, feature_hierarchy.parents)
    unmapped_features = []
    liftover_type = "copies"
    extract_features.get_gene_sequences(feature_hierarchy.parents, ref_chroms, args, liftover_type)
    min_cov, min_seqid = 0, args.sc
    align_and_lift_features(ref_chroms, target_chroms, args, feature_hierarchy, liftover_type, unmapped_features,
                            feature_db,
                            feature_hierarchy.parents, lifted_features_list, ref_parent_order, min_cov, min_seqid,
                            args.overlap)
