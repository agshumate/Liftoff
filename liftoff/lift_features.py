from liftoff import find_best_mapping, liftoff_utils, merge_lifted_features



def lift_all_features(alns, threshold, feature_db, features_to_lift, feature_hierarchy,
                      unmapped_features, lifted_feature_list, seq_id_threshold, feature_locations, distance_factor,
                      allow_chrom_split):
    feature_order = get_feature_order(feature_db)
    alignments = sort_alignments(features_to_lift, alns)
    num_features = 0
    previous_ref_chrom = ""
    for alignment in alignments:
        ref_chrom = feature_hierarchy.parents[liftoff_utils.convert_id_to_original(alignment[0].query_name)].seqid
        num_features += 1
        if ref_chrom != previous_ref_chrom:
            previous_gene_seq = ""
            previous_gene_start = 0
        lifted_fragments = lift_single_feature(threshold, feature_order, features_to_lift, feature_hierarchy,
                                                                                previous_gene_start,previous_gene_seq, unmapped_features,
                                                                                alignment, seq_id_threshold, feature_locations,
                                                                                lifted_feature_list, distance_factor,allow_chrom_split)
        for lifted_fragment in lifted_fragments:
            if lifted_fragment[0] != []:
                previous_gene_start = lifted_fragments[0][2]
                previous_gene_seq = lifted_fragments[0][3]
                lifted_feature_list[lifted_fragment[1]] = lifted_fragment[0]
        previous_ref_chrom = ref_chrom


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


def lift_single_feature(threshold, feature_order, features_to_lift, feature_hierarchy,
                        previous_gene_start, previous_gene_seq, unmapped_features, aligned_feature, seq_id_threshold,
                        feature_locations,
                        lifted_features_list, distance_factor, allow_chrom_split):
    new_parent_name = aligned_feature[0].query_name
    original_parent_name = liftoff_utils.convert_id_to_original(new_parent_name)
    parent = features_to_lift[original_parent_name]
    if len(aligned_feature) > 0:
        all_converted_results = find_best_mapping.find_best_mapping(aligned_feature,
                                                                                          parent.end - parent.start + 1,
                                                                                          parent,
                                                                                          feature_hierarchy,
                                                                                          previous_gene_start,
                                                                    previous_gene_seq, feature_locations,
                                                                                          lifted_features_list,
                                                                                          distance_factor, allow_chrom_split)

        lifted_fragments = []
        total_seq_id = sum(result[2] for result in all_converted_results)
        total_coverage = sum(result[1] for result in all_converted_results)
        for i in range (len(all_converted_results)):
            result = all_converted_results[i]
            frag_tag = "_frag" + str(i)
            lifted_children, alignment_coverage, seq_id = result[0], result[1], result[2]
            lifted_features, feature_start, feature_chrm = merge_lifted_features.merge_lifted_features(lifted_children,
                                                                                         parent,
                                                                                         unmapped_features, threshold,
                                                                                         new_parent_name + frag_tag,
                                                                                         feature_order,
                                                                                         feature_hierarchy,
                                                                                         alignment_coverage, seq_id,
                                                                                         seq_id_threshold, i,
                                                                                         total_seq_id, total_coverage)
            lifted_fragments.append([lifted_features, aligned_feature[0].query_name + frag_tag, feature_start, feature_chrm])
    else:
        unmapped_features.append(parent)
        feature_start = 0
    return lifted_fragments


