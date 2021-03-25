from liftoff import find_best_mapping, liftoff_utils, merge_lifted_features


def lift_all_features(alns, threshold, feature_db,  feature_hierarchy,
                      unmapped_features, lifted_feature_list, seq_id_threshold, feature_locations, args,
                      ref_parent_order):
    features_to_lift = feature_hierarchy.parents
    feature_order = get_feature_order(feature_db)
    alignments = sort_alignments(features_to_lift, alns)
    num_features = 0
    for alignment in alignments:
        num_features += 1
        preivous_feature_start, previous_feature_seq, previous_feature_ref_start = find_neighbor_location(
            feature_hierarchy.parents, alignment, lifted_feature_list, ref_parent_order)

        lifted_features, parent_name = lift_single_feature(threshold, feature_order, features_to_lift,
                                                           feature_hierarchy, preivous_feature_start,
                                                           previous_feature_ref_start,
                                                           previous_feature_seq, unmapped_features,
                                                           alignment, seq_id_threshold, feature_locations,
                                                           lifted_feature_list, args)

        if lifted_features != []:
            lifted_feature_list[parent_name] = lifted_features



def get_feature_order(gene_db):
    feature_types = list(gene_db.featuretypes())
    index = 0
    feature_order = {}
    if 'exon' in feature_types:
        feature_order['exon'] = index
        index += 1
    if 'CDS' in feature_types:
        feature_order['CDS'] = index
        index += 1
    for feature_type in feature_types:
        if feature_type not in feature_order:
            feature_order[feature_type] = index
            index += 1
    return feature_order


def sort_alignments(parent_dict, alignments):
    parent_list = []
    order = 0
    order_dict = {}
    values = list(alignments.values())
    for alignment in alignments:
        parent_list.append(parent_dict[liftoff_utils.convert_id_to_original(alignments[alignment][0].query_name)])
    parent_list.sort(key=lambda x: (x.seqid, x.start))
    for parent in parent_list:
        order_dict[parent.id] = order
        order += 1
    values.sort(key=lambda x: order_dict[liftoff_utils.convert_id_to_original(x[0].query_name)])
    return values


def find_neighbor_location(ref_parents, alignment, lifted_feature_list, ref_parent_order):
    ref_feature = ref_parents[liftoff_utils.convert_id_to_original(alignment[0].query_name)]
    ref_neighbor_name = liftoff_utils.find_nonoverlapping_upstream_neighbor(ref_parent_order, ref_feature.id)
    if ref_neighbor_name is not None:
        ref_neighbor_key = ref_neighbor_name + "_0"
        if ref_neighbor_key in lifted_feature_list:
            previous_feature_start = lifted_feature_list[ref_neighbor_key][0].start
            previous_feature_seq = lifted_feature_list[ref_neighbor_key][0].seqid
            previous_feature_ref_start = ref_parents[ref_neighbor_name].start
            return previous_feature_start, previous_feature_seq, previous_feature_ref_start
    return 0, "", 0


def lift_single_feature(threshold, feature_order, features_to_lift, feature_hierarchy,
                        previous_feature_start, previous_feature_ref_start, previous_gene_seq, unmapped_features,
                        aligned_feature, seq_id_threshold,
                        feature_locations,
                        lifted_features_list, args):
    new_parent_name = aligned_feature[0].query_name
    original_parent_name = liftoff_utils.convert_id_to_original(new_parent_name)
    parent = features_to_lift[original_parent_name]
    if len(aligned_feature) > 0:
        lifted_children, alignment_coverage, seq_id = find_best_mapping.find_best_mapping(aligned_feature,
                                                                                          parent.end - parent.start + 1,
                                                                                          parent,
                                                                                          feature_hierarchy,
                                                                                          previous_feature_start,
                                                                                          previous_feature_ref_start,
                                                                                          previous_gene_seq,
                                                                                          feature_locations,
                                                                                          lifted_features_list,
                                                                                          args)


        lifted_features = merge_lifted_features.merge_lifted_features(lifted_children,
                                                                      parent,
                                                                      unmapped_features, threshold,
                                                                      new_parent_name, feature_order,
                                                                      feature_hierarchy,
                                                                      alignment_coverage, seq_id,
                                                                      seq_id_threshold)



    else:
        unmapped_features.append(parent)
    return lifted_features, aligned_feature[0].query_name
