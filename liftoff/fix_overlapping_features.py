from liftoff import lift_features, liftoff_utils
import numpy as np
from interlap import InterLap


def fix_incorrectly_overlapping_features(all_lifted_features, features_to_check, all_aligned_segs,
                                         unmapped_features, threshold, feature_hierarchy, feature_db,
                                         ref_parent_order, seq_id_threshold):
    features_to_remap = check_homologues(all_lifted_features, features_to_check, feature_hierarchy.parents,
                                         ref_parent_order, None)
    resolve_overlapping_homologues(all_aligned_segs, all_lifted_features, features_to_remap, unmapped_features,
                                   threshold, feature_hierarchy, feature_db, ref_parent_order, seq_id_threshold)

    # all_feature_list = liftoff_utils.get_parent_list(all_lifted_features, parent_dict)
    # inter = build_interval_list(all_feature_list)
    # for feature in all_feature_list:
    #     overlaps = find_overlaps(feature.start, feature.end, feature.seqid, feature.strand, feature.attributes[
    #     "copy_id"][0], inter)
    #     for overlap in overlaps:
    #         if overlaps_in_ref_annotation(parent_dict[overlap[2][1].id], parent_dict[feature.id]) is False:
    #             print(overlap, feature.id)



def check_homologues(all_lifted_features, lifted_features_to_check, parent_dict, ref_parent_order, inter):
    all_feature_list = liftoff_utils.get_parent_list(all_lifted_features)
    features_to_check_list = liftoff_utils.get_parent_list(lifted_features_to_check)
    target_parent_order = liftoff_utils.find_parent_order(all_feature_list)
    remap_features = {}
    if inter == None:
        inter = build_interval_list(all_feature_list)
    for feature in features_to_check_list:
        overlaps = find_overlaps(feature.start, feature.end, feature.seqid, feature.strand,
                                 feature.attributes["copy_id"][0], inter)
        for overlap in overlaps:
            feature_to_compare = overlap[2][1]
            compare_overlapping_feature(feature_to_compare, parent_dict, feature, remap_features, ref_parent_order,
                                        target_parent_order)
    return remap_features


def build_interval_list(features):
    inter = InterLap()
    feature_coords = [[feature.start, feature.end, [feature.attributes["copy_id"][0], feature]] for feature in
                      features]
    if len(feature_coords) > 0:
        inter.update(feature_coords)
    return inter


def find_overlaps(start, end, chrm, strand, feature_name, intervals):
    all_overlaps = intervals.find((start, end))
    filtered_overlaps = [overlap for overlap in all_overlaps if
                         overlap[2][1].seqid == chrm and overlap[2][1].strand == strand and overlap[2][
                             0] != feature_name]
    return filtered_overlaps


def compare_overlapping_feature(overlapping_feature, ref_parent_list, feature, remap_features, ref_parent_order,
                                target_parent_order):
    ref_feature = ref_parent_list[feature.id]
    ref_overlap_feature = ref_parent_list[overlapping_feature.id]
    if overlaps_in_ref_annotation(ref_feature, ref_overlap_feature) is False:
        feature_to_remap, feature_to_keep = find_feature_to_remap(feature, overlapping_feature,
                                                                  ref_parent_order, target_parent_order, remap_features,
                                                                  ref_feature, ref_overlap_feature)
        remap_features[feature_to_remap.attributes["copy_id"][0]] = [feature_to_keep]


def overlaps_in_ref_annotation(ref_feature1, ref_feature2):
    if ref_feature1.seqid != ref_feature2.seqid:
        return False
    if ref_feature1.strand != ref_feature2.strand:
        return False
    if ref_feature1.id == ref_feature2.id:
        return False
    else:
        return liftoff_utils.count_overlap(ref_feature1.start, ref_feature1.end,
                                           ref_feature2.start,
                                           ref_feature2.end) > 0


def find_feature_to_remap(feature, overlap_feature, ref_parent_order, target_parent_order, remap_features, ref_feature,
                          ref_overlap_feature):
    feature_is_copy = is_copy(feature)
    overlap_feature_is_copy = is_copy(overlap_feature)
    if feature_is_copy and overlap_feature_is_copy is False:
        return feature, overlap_feature
    if overlap_feature_is_copy and feature_is_copy is False:
        return overlap_feature, feature
    if already_in_list(feature, remap_features):
        return feature, overlap_feature
    if already_in_list(overlap_feature, remap_features):
        return overlap_feature, feature
    if has_greater_seq_id(feature, overlap_feature):
        return overlap_feature, feature
    if has_greater_seq_id(overlap_feature, feature):
        return feature, overlap_feature
    if feature_is_copy is False:
        if is_in_order(ref_parent_order, target_parent_order, feature, ref_feature):
            return overlap_feature, feature
    if overlap_feature_is_copy is False:
        if is_in_order(ref_parent_order, target_parent_order, overlap_feature, ref_overlap_feature):
            return feature, overlap_feature
    if is_shorter(feature, overlap_feature):
        return feature, overlap_feature
    if is_shorter(overlap_feature, feature):
        return overlap_feature, feature
    return overlap_feature, feature


def is_copy(feature):
    return feature.attributes["copy_id"][0][-2:] != '_0'


def already_in_list(feature, remap_features):
    return feature.attributes["copy_id"][0] in remap_features


def has_greater_seq_id(feature1, feature2):
    if feature1.score < feature2.score:
        return True
    return False


def is_in_order(ref_parent_order, target_parent_order, feature, ref_feature):
    target_neighbor = find_nonoverlapping_upstream_neighbor(target_parent_order, feature)
    reference_neighbor = find_nonoverlapping_upstream_neighbor(ref_parent_order, ref_feature)
    return target_neighbor == reference_neighbor


def find_nonoverlapping_upstream_neighbor(parent_order, feature):
    feature_location = np.where(parent_order[:, 0] == feature.id)[0][0]
    neighbor_indx = feature_location - 1
    while neighbor_indx > 0:
        neighbor = parent_order[neighbor_indx][1]
        neighbor_id = parent_order[neighbor_indx][0]
        if neighbor.seqid != feature.seqid:
            return None
        if neighbor.end < feature.start:
            return neighbor_id
        neighbor_indx -= 1


def is_shorter(feature1, feature2):
    return (feature1.end - feature1.start) < (feature2.end - feature2.start)


def resolve_overlapping_homologues(all_aligned_segs, lifted_feature_list, features_to_remap, unmapped_features,
                                   threshold, feature_hierarchy, feature_db, ref_parent_order, seq_id_threshold):
    all_overlapping_features = {}
    iter = 0
    max_iter = 10 * len(features_to_remap)
    while len(features_to_remap) > 0:
        print(len(features_to_remap))
        iter += 1
        if iter > max_iter:
            break
        aligned_segs_for_remap = remove_features_and_get_alignments(features_to_remap, lifted_feature_list,
                                                                    all_overlapping_features, all_aligned_segs)
        lift_features.lift_all_features(aligned_segs_for_remap, all_overlapping_features, threshold, feature_db,
                                        feature_hierarchy.parents, feature_hierarchy, unmapped_features,
                                        lifted_feature_list, seq_id_threshold)
        inter = remove_remapped_features_from_overlaps(lifted_feature_list, all_overlapping_features)
        features_to_check = get_successfully_remapped_features(lifted_feature_list, features_to_remap)
        features_to_remap = check_homologues(lifted_feature_list, features_to_check, feature_hierarchy.parents,
                                             ref_parent_order, inter)
    remove_unresolved_features(features_to_remap, feature_hierarchy.parents, lifted_feature_list, unmapped_features)
    return lifted_feature_list


def remove_features_and_get_alignments(features_to_remap, lifted_feature_list, all_overlapping_features,
                                       all_aligned_segs):
    aligned_segs_for_remap = {}
    for feature_to_remap in features_to_remap:
        del lifted_feature_list[feature_to_remap]
        aligned_segs_for_remap[feature_to_remap] = all_aligned_segs[feature_to_remap]
        add_overlapping_feature(features_to_remap, feature_to_remap, all_overlapping_features)
    return aligned_segs_for_remap


def add_overlapping_feature(features_to_remap, feature_to_remap, all_overlapping_features):
    if feature_to_remap in all_overlapping_features:
        all_overlapping_features[feature_to_remap].append(features_to_remap[feature_to_remap][0])
    else:
        all_overlapping_features[feature_to_remap] = [features_to_remap[feature_to_remap][0]]


def remove_remapped_features_from_overlaps(lifted_feature_list, all_overlapping_features):
    parents = liftoff_utils.get_parent_list(lifted_feature_list)
    inter = build_interval_list(parents)
    for feature_name in all_overlapping_features:
        for overlapping_feature in all_overlapping_features[feature_name]:
            updated_overlaps = find_overlaps(overlapping_feature.start - 1, overlapping_feature.end - 1,
                                             overlapping_feature.seqid,
                                             overlapping_feature.strand, feature_name, inter)
            if len(updated_overlaps) == 0:
                all_overlapping_features[feature_name].remove(overlapping_feature)
    return inter


def get_successfully_remapped_features(lifted_feature_list, features_to_remap):
    successfully_remapped_features = {}
    for feature_to_remap in features_to_remap:
        if feature_to_remap in lifted_feature_list:
            successfully_remapped_features[feature_to_remap] = lifted_feature_list[feature_to_remap]
    return successfully_remapped_features


def remove_unresolved_features(features_to_remap, parent_dict, lifted_feature_list, unmapped_features):
    for feature in features_to_remap:
        unmapped_features.append(parent_dict[liftoff_utils.convert_id_to_original(feature)])
        del lifted_feature_list[feature]
