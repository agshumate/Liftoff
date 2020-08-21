from liftoff import lift_features, liftoff_utils
import numpy as np
from interlap import InterLap


def fix_incorrectly_overlapping_features(all_lifted_features, features_to_check, all_aligned_segs,
                                         unmapped_features, threshold, feature_hierarchy, feature_db,
                                         ref_parent_order, seq_id_threshold, distance_factor,
                                         max_overlap):
    features_to_remap, feature_locations = check_homologues(all_lifted_features, features_to_check,
                                                            feature_hierarchy.parents,
                                                            ref_parent_order, max_overlap)
    resolve_overlapping_homologues(all_aligned_segs, all_lifted_features, features_to_remap, unmapped_features,
                                   threshold, feature_hierarchy, feature_db, ref_parent_order, seq_id_threshold,
                                   feature_locations, distance_factor, max_overlap)


def check_homologues(all_lifted_features, lifted_features_to_check, parent_dict, ref_parent_order, max_overlap):
    all_feature_list = liftoff_utils.get_parent_list(all_lifted_features)
    features_to_check_list = liftoff_utils.get_parent_list(lifted_features_to_check)
    target_parent_order = liftoff_utils.find_parent_order(all_feature_list)
    remap_features = set()
    feature_locations = build_interval_list(all_feature_list)
    for feature in features_to_check_list:
        overlaps = liftoff_utils.find_overlaps(feature.start - 1, feature.end - 1, feature.seqid, feature.strand,
                                               feature.attributes["copy_id"][0], feature_locations, parent_dict,
                                               all_lifted_features, max_overlap)
        for overlap in overlaps:
            if overlap[2][0] != feature.attributes["copy_id"][0]:
                feature_to_compare = overlap[2][1]
                compare_overlapping_feature(feature_to_compare, feature, remap_features, ref_parent_order,
                                            target_parent_order)
    return remap_features, feature_locations


def build_interval_list(features):
    inter = InterLap()
    feature_coords = [[feature.start - 1, feature.end - 1, [feature.attributes["copy_id"][0], feature]] for feature in
                      features]
    if len(feature_coords) > 0:
        inter.update(feature_coords)
    return inter


def compare_overlapping_feature(overlapping_feature, feature, remap_features, ref_parent_order,
                                target_parent_order):
    feature_to_remap = find_feature_to_remap(feature, overlapping_feature,
                                             ref_parent_order, target_parent_order, remap_features)

    remap_features.add(feature_to_remap.attributes["copy_id"][0])


def find_feature_to_remap(feature, overlap_feature, ref_parent_order, target_parent_order, remap_features):
    feature_is_copy = is_copy(feature)
    overlap_feature_is_copy = is_copy(overlap_feature)
    if feature_is_copy and overlap_feature_is_copy is False:
        return feature
    if overlap_feature_is_copy and feature_is_copy is False:
        return overlap_feature
    if already_in_list(feature, remap_features):
        return feature
    if already_in_list(overlap_feature, remap_features):
        return overlap_feature
    if has_greater_seq_id(feature, overlap_feature):
        return overlap_feature
    if has_greater_seq_id(overlap_feature, feature):
        return feature
    farthest_neighbor_distance = find_out_of_order_feature(feature_is_copy, ref_parent_order,
                                                           target_parent_order, feature,
                                                           overlap_feature_is_copy, overlap_feature)
    if farthest_neighbor_distance is not None:
        return farthest_neighbor_distance
    if is_shorter(feature, overlap_feature):
        return feature
    if is_shorter(overlap_feature, feature):
        return overlap_feature
    return overlap_feature


def is_copy(feature):
    return feature.attributes["copy_id"][0][-2:] != '_0'


def already_in_list(feature, remap_features):
    return feature.attributes["copy_id"][0] in remap_features


def has_greater_seq_id(feature1, feature2):
    if feature1.score < feature2.score:
        return True
    return False


def find_out_of_order_feature(feature_is_copy, ref_parent_order, target_parent_order, feature,
                              overlap_feature_is_copy, overlap_feature):
    if feature_is_copy is False:
        neighbor_distance_feature = check_order(ref_parent_order, target_parent_order, feature, overlap_feature)
    else:
        neighbor_distance_feature = None
    if overlap_feature_is_copy is False:
        neighbor_distance_overlap_feature = check_order(ref_parent_order, target_parent_order, overlap_feature, feature)
    else:
        neighbor_distance_overlap_feature = None
    if neighbor_distance_feature is not None and neighbor_distance_overlap_feature is not None:
        if neighbor_distance_feature < neighbor_distance_overlap_feature:
            return overlap_feature
        else:
            return feature
    return None


def check_order(ref_parent_order, target_parent_order, feature, overlap_feature):
    ref_neighbor_upstream = liftoff_utils.find_nonoverlapping_upstream_neighbor(ref_parent_order, feature.id)
    if ref_neighbor_upstream is not None:
        ref_distance = find_distance_between_features(ref_neighbor_upstream, target_parent_order, feature,
                                                      overlap_feature)
        return ref_distance
    return None


def find_distance_between_features(ref_neighbor_upstream, target_parent_order, feature, overlap_feature):
    feature_location_arr = np.where(target_parent_order[:, 0] == feature.id)[0]
    neighbor_location_arr = np.where(target_parent_order[:, 0] == ref_neighbor_upstream)[0]
    if len(feature_location_arr) > 0 and len(neighbor_location_arr) > 0:
        feature_location, neighbor_location = feature_location_arr[0], neighbor_location_arr[0]
        features_in_range = target_parent_order[min(feature_location, neighbor_location):max(feature_location,
                                                                                             neighbor_location) + 1]
        distance = np.abs(feature_location - neighbor_location)
        if overlap_feature.id in features_in_range[:, 0]:
            distance -= 1
        return distance
    return None


def is_shorter(feature1, feature2):
    return (feature1.end - feature1.start) < (feature2.end - feature2.start)


def resolve_overlapping_homologues(all_aligned_segs, lifted_feature_list, features_to_remap, unmapped_features,
                                   threshold, feature_hierarchy, feature_db, ref_parent_order, seq_id_threshold,
                                   feature_locations, distance_factor, max_overlap):
    iter = 0
    max_iter = 10 * len(features_to_remap)
    while len(features_to_remap) > 0:
        iter += 1
        if iter > max_iter:
            break
        aligned_segs_for_remap = remove_features_and_get_alignments(features_to_remap, lifted_feature_list,
                                                                    all_aligned_segs)
        lift_features.lift_all_features(aligned_segs_for_remap, threshold, feature_db,
                                        feature_hierarchy.parents, feature_hierarchy, unmapped_features,
                                        lifted_feature_list, seq_id_threshold, feature_locations, distance_factor,
                                        ref_parent_order)
        features_to_check = get_successfully_remapped_features(lifted_feature_list, features_to_remap,
                                                               )
        features_to_remap, feature_locations = check_homologues(lifted_feature_list, features_to_check,
                                                                feature_hierarchy.parents,
                                                                ref_parent_order, max_overlap)
    remove_unresolved_features(features_to_remap, feature_hierarchy.parents, lifted_feature_list, unmapped_features)
    return lifted_feature_list


def remove_features_and_get_alignments(features_to_remap, lifted_feature_list, all_aligned_segs):
    aligned_segs_for_remap = {}
    for feature_to_remap in features_to_remap:
        del lifted_feature_list[feature_to_remap]
        aligned_segs_for_remap[feature_to_remap] = all_aligned_segs[feature_to_remap]
    return aligned_segs_for_remap


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
