from liftoff import liftoff_utils, new_feature


def merge_lifted_features(mapped_children, parent, unmapped_features, aln_cov_threshold, copy_id, feature_order,
                          feature_hierarchy, aln_cov, seq_id, seq_id_threshold):
    feature_list, final_features = {}, []
    non_parents = []
    top_target_feature = None
    if len(mapped_children) == 0:
        unmapped_features.append(parent)
        return []
    for child_name in mapped_children:
        child_feature = mapped_children[child_name]
        feature_list[child_feature.id] = child_feature
        if is_top_parent(child_feature, parent) is False:
            non_parents.append((child_feature, child_feature.attributes["Parent"][0]))
        else:
            top_target_feature = child_feature
    while (len(non_parents) != 0):
        non_parents, top_target_feature = create_parents(non_parents, parent, feature_hierarchy, feature_list)
    final_features = process_final_features_list(feature_list, top_target_feature, seq_id, seq_id_threshold, aln_cov,
                                                 aln_cov_threshold,
                                                 unmapped_features, parent, feature_order, copy_id)
    return final_features


def is_top_parent(child, parent):
    return child.id == parent.id


def create_parents(non_parents, top_ref_parent, feature_hierarchy, feature_list):
    added_parent_ids, new_non_parents = [], []
    top_target_feature = None
    for child, parent in non_parents:
        if parent not in added_parent_ids:
            target_parent_feature = make_new_parent(feature_list, parent, feature_hierarchy)
            if is_top_parent(target_parent_feature, top_ref_parent) is False:
                new_non_parents.append((target_parent_feature, target_parent_feature.attributes["Parent"][0]))
            else:
                top_target_feature = target_parent_feature
            added_parent_ids.append(parent)
    return new_non_parents, top_target_feature


def make_new_parent(feature_list, parent, feature_hierarchy):
    children = [feature for feature in feature_list.values() if
                "Parent" in feature.attributes and feature.attributes["Parent"][0] == parent]
    starts, ends = [child.start for child in children], [child.end for child in children]
    ref_parent = get_ref_parent(parent, feature_hierarchy)
    target_parent_feature = new_feature.new_feature(ref_parent.id, ref_parent.featuretype, children[0].seqid,
                                                    'Liftoff', children[0].strand, min(starts), max(ends),
                                                    dict(ref_parent.attributes))
    feature_list[target_parent_feature.id] = target_parent_feature
    return target_parent_feature


def get_ref_parent(parent, feature_hierarchy):
    if parent in feature_hierarchy.parents:
        ref_parent = feature_hierarchy.parents[parent]
    else:
        ref_parent = feature_hierarchy.intermediates[parent]
    return ref_parent


def process_final_features_list(feature_list, top_target_feature, seq_id, seq_id_threshold, aln_cov, aln_cov_threshold,
                                unmapped_features, parent, feature_order, copy_id):
    final_features = [feature for feature in feature_list.values() if feature != top_target_feature]
    final_features.sort(key=lambda x: (x.seqid, x.start))
    final_features.sort(key=lambda x: feature_order[x.featuretype])
    final_features.insert(0, top_target_feature)
    if aln_cov < aln_cov_threshold or seq_id < seq_id_threshold:
        final_features = []
        unmapped_features.append(parent)
    else:
        top_target_feature.score = 1 - seq_id
        top_target_feature.attributes["copy_id"] = [copy_id]
        top_target_feature.attributes["coverage"] = [str(aln_cov)[0:5]]
        top_target_feature.attributes["sequence_ID"] = [str(seq_id)[0:5]]
    return final_features
