from liftoff import lift_features, liftoff_utils
import numpy as np
from interlap import InterLap


def fix_incorrectly_overlapping_features(all_lifted_features, features_to_check, all_aligned_segs,
                                         unmapped_features, threshold, feature_hierarchy, feature_db,
                                         ref_parent_order, seq_id_threshold, distance_factor, allow_chrom_split):
    features_to_remap, feature_locations = check_homologues(all_lifted_features, features_to_check,
                                                            feature_hierarchy.parents,
                                                            ref_parent_order)
    resolve_overlapping_homologues(all_aligned_segs, all_lifted_features, features_to_remap, unmapped_features,
                                   threshold, feature_hierarchy, feature_db, ref_parent_order, seq_id_threshold,
                                   feature_locations, distance_factor, allow_chrom_split)




def check_homologues(all_lifted_features, lifted_features_to_check, parent_dict, ref_parent_order):
    all_feature_list = liftoff_utils.get_parent_list(all_lifted_features)
    features_to_check_list = liftoff_utils.get_parent_list(lifted_features_to_check)
    target_parent_order = liftoff_utils.find_parent_order(all_feature_list)
    remap_features = set()
    feature_locations = build_interval_list(all_feature_list)
    for feature in features_to_check_list:
        overlaps = liftoff_utils.find_overlaps(feature.start - 1, feature.end - 1, feature.seqid, feature.strand,
                                               feature.attributes["copy_id"][0], feature_locations, parent_dict,
                                               all_lifted_features)
        for overlap in overlaps:
            feature_to_compare = overlap[2][1]
            compare_overlapping_feature(feature_to_compare, parent_dict, feature, remap_features, ref_parent_order,
                                        target_parent_order)
    return remap_features, feature_locations


def build_interval_list(features):
    inter = InterLap()
    feature_coords = [[feature.start - 1, feature.end - 1, [feature.attributes["copy_id"][0], feature]] for feature in
                      features]
    if len(feature_coords) > 0:
        inter.update(feature_coords)
    return inter


def compare_overlapping_feature(overlapping_feature, ref_parent_list, feature, remap_features, ref_parent_order,
                                target_parent_order):
    ref_feature = ref_parent_list[feature.id]
    ref_overlap_feature = ref_parent_list[overlapping_feature.id]
    feature_to_remap = find_feature_to_remap(feature, overlapping_feature,
                                             ref_parent_order, target_parent_order, remap_features,
                                             ref_feature, ref_overlap_feature)
    remap_features.add(liftoff_utils.remove_frag_tag(feature_to_remap.attributes["copy_id"][0]))



def find_feature_to_remap(feature, overlap_feature, ref_parent_order, target_parent_order, remap_features, ref_feature,
                          ref_overlap_feature):
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
    farthest_neighbor_distance = find_farthest_neighbor_distance(feature_is_copy, ref_parent_order,
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
    return liftoff_utils.remove_frag_tag(feature.attributes["copy_id"][0])[-2:] != '_0'


def already_in_list(feature, remap_features):
    return feature.attributes["copy_id"][0] in remap_features


def has_greater_seq_id(feature1, feature2):
    if feature1.score < feature2.score:
        return True
    return False


def find_farthest_neighbor_distance(feature_is_copy, ref_parent_order, target_parent_order, feature,
                                    overlap_feature_is_copy, overlap_feature):
    if feature_is_copy is False:
        neighbor_distance_feature = check_order(ref_parent_order, target_parent_order, feature)
    else:
        neighbor_distance_feature = None
    if overlap_feature_is_copy is False:
        neighbor_distance_overlap_feature= check_order(ref_parent_order, target_parent_order, overlap_feature,)
    else:
        neighbor_distance_overlap_feature = None
    if neighbor_distance_feature is not None and neighbor_distance_overlap_feature is not None:
        if neighbor_distance_feature < neighbor_distance_overlap_feature:
            return overlap_feature
        else:
            return feature
    return None


def check_order(ref_parent_order, target_parent_order, feature):
    target_neighbor_upstream = find_nonoverlapping_upstream_neighbor(target_parent_order, feature)
    if target_neighbor_upstream is not None:
        ref_distance = find_distance_between_features(target_neighbor_upstream, ref_parent_order, feature)
        return ref_distance
    return None



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
    return None


def find_distance_between_features(target_neighbor_upstream, ref_parent_order, feature):
    feature_location = np.where(ref_parent_order[:, 0] == feature.id)[0][0]
    neighbor_location = np.where(ref_parent_order[:, 0] == target_neighbor_upstream)[0][0]
    return np.abs(feature_location - neighbor_location)


def is_shorter(feature1, feature2):
    return (feature1.end - feature1.start) < (feature2.end - feature2.start)


def resolve_overlapping_homologues(all_aligned_segs, lifted_feature_list, features_to_remap, unmapped_features,
                                   threshold, feature_hierarchy, feature_db, ref_parent_order, seq_id_threshold,
                                   feature_locations, distance_factor, allow_chrom_split):
    iter = 0
    max_iter = 10 * len(features_to_remap)
    while len(features_to_remap) > 0:
        print(len(features_to_remap))
        iter += 1
        if iter > max_iter:
            break
        aligned_segs_for_remap = remove_features_and_get_alignments(features_to_remap, lifted_feature_list,
                                                                    all_aligned_segs)
        lift_features.lift_all_features(aligned_segs_for_remap, threshold, feature_db,
                                        feature_hierarchy.parents, feature_hierarchy, unmapped_features,
                                        lifted_feature_list, seq_id_threshold, feature_locations, distance_factor,
                                        allow_chrom_split)
        features_to_check = get_successfully_remapped_features(lifted_feature_list, features_to_remap,
                                                               allow_chrom_split)
        features_to_remap, feature_locations = check_homologues(lifted_feature_list, features_to_check,
                                                                feature_hierarchy.parents,
                                                                ref_parent_order)
    remove_unresolved_features(features_to_remap, feature_hierarchy.parents, lifted_feature_list, unmapped_features)
    return lifted_feature_list


def remove_features_and_get_alignments(features_to_remap, lifted_feature_list,
                                       all_aligned_segs):
    aligned_segs_for_remap = {}

    keys_to_remove = [key for key in lifted_feature_list if
                          liftoff_utils.remove_frag_tag(key) in features_to_remap]

    for key in keys_to_remove:
        del lifted_feature_list[key]
        aligned_segs_for_remap[liftoff_utils.remove_frag_tag(key)] = all_aligned_segs[liftoff_utils.remove_frag_tag(key)]
    return aligned_segs_for_remap


def get_successfully_remapped_features(lifted_feature_list, features_to_remap, allow_chrom_split):
    successfully_remapped_features = {}
    for feature_to_remap in features_to_remap:
        if allow_chrom_split:
            remapped_keys = [key for key in lifted_feature_list if liftoff_utils.remove_frag_tag(key) == feature_to_remap]
        else:
            remapped_keys = [key for key in lifted_feature_list if key == feature_to_remap + "_frag0"]
        for key in remapped_keys:
            successfully_remapped_features[key] = lifted_feature_list[key]
    return successfully_remapped_features


def remove_unresolved_features(features_to_remap, parent_dict, lifted_feature_list, unmapped_features):
    for feature in features_to_remap:
        unmapped_features.append(parent_dict[liftoff_utils.convert_id_to_original(feature)])
        del lifted_feature_list[feature]
