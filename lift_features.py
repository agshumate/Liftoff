import  merge_lifted_features
import find_best_mapping




def lift_all_features(alns, all_overlapping_features, threshold,  feature_db, parent_dict, children_dict,
                      intermediate_dict, unmapped_features, lifted_feature_list):
    feature_order = get_feature_order(feature_db)
    previous_gene_start = 0
    alignments = sort_alignments(parent_dict, alns)
    for alignment in alignments:
        lifted_features, parent_name, previous_gene_start = lift_features_subset(all_overlapping_features, threshold,
                                                                                 feature_order, parent_dict,
                                                                                 children_dict, intermediate_dict,
                                                                                 previous_gene_start, unmapped_features,
                                                                                 alignment)
        lifted_feature_list[parent_name]=lifted_features


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
            index +=1
    return feature_order


def sort_alignments(parent_dict, alignments):
    parent_list = []
    order = 0
    order_dict = {}
    values = list(alignments.values())
    for alignment in alignments:
       parent_list.append(parent_dict[alignments[alignment][0].original_id])
    parent_list.sort(key=lambda x: (x.seqid, x.start))
    for parent in parent_list:
            order_dict[parent.id] = order
            order+=1
    values.sort(key=lambda x: order_dict[x[0].original_id])
    return values


def lift_features_subset(all_overlapping_features, threshold, feature_order, parent_dict, children_dict,
                             intermediate_dict, previous_gene_start, unmapped_features, aligned_feature):
    new_parent_name = aligned_feature[0].query_name
    copy_tag_len = len(new_parent_name.split("_")[-1])
    original_parent_name = new_parent_name[:-copy_tag_len-1]
    if new_parent_name in all_overlapping_features:
        overlapping_features = all_overlapping_features[new_parent_name]
    else:
        overlapping_features = []
    parent = parent_dict[original_parent_name]
    if len(aligned_feature) > 0:
        lifted_children, shortest_path_weight, alignment_coverage, seq_id = find_best_mapping.find_best_mapping(aligned_feature,
                                                                                    parent.end - parent.start + 1,
                                                                                    parent, overlapping_features,
                                                                                    children_dict, previous_gene_start)
        lifted_feature_list, feature_start = merge_lifted_features.merge_lifted_features(lifted_children,
                                                                                         shortest_path_weight, parent,
                                                                                         unmapped_features, threshold,
                                                                                         new_parent_name, feature_order,
                                                                                         parent_dict, intermediate_dict,alignment_coverage, seq_id)
    else:
        unmapped_features.append(parent)
        feature_start = 0
    return lifted_feature_list, aligned_feature[0].query_name, feature_start
