from liftoff import lift_features, liftoff_utils
import numpy as np


def fix_incorrectly_overlapping_features(all_lifted_features, features_to_check, parent_dict, all_aligned_segs,
                                         unmapped_features, threshold, intermediate_dict, children_dict, feature_db,
                                         original_parent_order,seq_id_threshold):
    features_to_remap = check_homologues(all_lifted_features, features_to_check, parent_dict, original_parent_order)
    resolve_overlapping_homologues(all_aligned_segs, all_lifted_features, features_to_remap, unmapped_features,
                                threshold,  parent_dict, intermediate_dict,
                                   children_dict, feature_db,original_parent_order, seq_id_threshold)




def check_homologues(all_lifted_features, lifted_features_to_check, parent_dict, original_parent_order):
    all_feature_list = liftoff_utils.get_parent_list(all_lifted_features, parent_dict)
    features_to_check_list = liftoff_utils.get_parent_list(lifted_features_to_check, parent_dict)
    all_feature_list.sort(key = lambda x: (x.seqid, x.start))
    new_parent_order = liftoff_utils.find_parent_order(all_feature_list)
    chrom_index_dict = {}
    feature_index = 0
    for feature in all_feature_list:
        if feature.seqid not in chrom_index_dict:
            chrom_index_dict[feature.seqid] = feature_index
        feature_index += 1
    remap_features = {}
    for feature in features_to_check_list:
        nearby_features = all_feature_list[chrom_index_dict[feature.seqid]:]
        compare_nearby_features(nearby_features, parent_dict, feature, remap_features, original_parent_order,
                                new_parent_order)
    return remap_features




def compare_nearby_features(nearby_features, parent_list, feature, remap_features, original_parent_order, new_parent_order):
    for nearby_feature in nearby_features:
        if nearby_feature.attributes["copy_id"][0] != feature.attributes["copy_id"][0] and feature.strand == nearby_feature.strand:
            if nearby_feature.seqid != feature.seqid or nearby_feature.start > feature.end:
                break
            if nearby_feature.end < feature.start:
                continue
            else:
                original_nearby_feature, original_feature = parent_list[nearby_feature.id], parent_list[feature.id]
                if original_nearby_feature.seqid != original_feature.seqid or original_nearby_feature.strand != original_feature.strand or original_nearby_feature.id == original_feature.id:
                    original_overlap = -1
                else:
                    original_overlap = liftoff_utils.count_overlap(original_feature.start, original_feature.end,
                                                                   original_nearby_feature.start,
                                                                   original_nearby_feature.end)
                lifted_overlap = liftoff_utils.count_overlap(feature.start, feature.end, nearby_feature.start,
                                                             nearby_feature.end)

                if lifted_overlap >0 and original_overlap <= 0:
                    if feature.attributes["copy_id"][0] in remap_features:
                        feature_to_remap, feature_to_keep = feature, nearby_feature
                    elif nearby_feature.attributes["copy_id"][0] in remap_features:
                        feature_to_remap, feature_to_keep = nearby_feature, feature

                    else:
                        feature_to_remap, feature_to_keep = find_feature_to_remap(feature, nearby_feature,
                                                                              original_parent_order, new_parent_order)

                    remap_features[feature_to_remap.attributes["copy_id"][0]]=[feature_to_keep.start - 1,
                                                                                      feature_to_keep.end - 1,
                                                                                      feature_to_keep.seqid,
                                                                                      feature_to_keep.strand]
    return remap_features


def find_neighbors(parent_order, feature, feature_exclude):
    if feature_exclude is not None and feature_exclude.id != feature.id:
        filtered_parent_order = parent_order[parent_order != feature_exclude.id]
    else:
        filtered_parent_order = parent_order
    feature_location = np.where(filtered_parent_order==feature.id)[0][0]
    left_neighbor = max(0, feature_location -1)
    feature_neighbors = filtered_parent_order[left_neighbor].tolist()
    return feature_neighbors

def find_feature_to_remap(feature, overlap_feature, original_parent_order, new_parent_order):
    if feature.score < overlap_feature.score:
        return overlap_feature, feature
    elif overlap_feature.score < feature.score:
        return feature, overlap_feature
    else:
        feature_neighbors = find_neighbors(new_parent_order, feature, overlap_feature)
        overlap_feature_neighbors = find_neighbors(new_parent_order, overlap_feature, feature)
        original_feature_neighbors = find_neighbors(original_parent_order, feature, None)
        original_overlap_feature_neighbors = find_neighbors(original_parent_order, overlap_feature, None)
        if feature_neighbors == original_feature_neighbors and overlap_feature_neighbors != original_overlap_feature_neighbors:
            return overlap_feature, feature
        elif feature_neighbors != original_feature_neighbors and overlap_feature_neighbors == original_overlap_feature_neighbors:
            return feature, overlap_feature

        if overlap_feature.end - overlap_feature.start < feature.end - feature.start:
            return overlap_feature, feature
        elif overlap_feature.end - overlap_feature.start >  feature.end - feature.start:
            return feature, overlap_feature
        else:
            min_feature = min(feature.attributes["copy_id"][0], overlap_feature.attributes["copy_id"][0])
            if min_feature == feature.attributes["copy_id"][0]:
                return feature, overlap_feature
            else:
                return overlap_feature, feature


def resolve_overlapping_homologues(all_aligned_segs, lifted_feature_list, features_to_remap,  unmapped_features,
                                   threshold,  parent_dict, intermediate_dict, children_dict,
                                   feature_db,original_parent_order,seq_id_threshold):
    all_overlapping_features = {}
    starting_remap_feature_num = len(features_to_remap)
    iter = 0
    while len(features_to_remap) > 0:
        iter +=1
        if iter > 10* starting_remap_feature_num:
            break
        features_to_check = {}
        aligned_segs_to_remap = {}
        for feature_to_remap in features_to_remap:
            del lifted_feature_list[feature_to_remap]
            aligned_segs_to_remap[feature_to_remap] = all_aligned_segs[feature_to_remap]
            add_overlapping_feature(features_to_remap, feature_to_remap, all_overlapping_features)
        lift_features.lift_all_features(aligned_segs_to_remap, all_overlapping_features, threshold, feature_db,
                                        parent_dict, children_dict, intermediate_dict, unmapped_features,
                                        lifted_feature_list, seq_id_threshold)
        clean_overlapping_features(lifted_feature_list, all_overlapping_features, parent_dict, features_to_remap, unmapped_features)
        for feature_to_remap in features_to_remap:
            if feature_to_remap in lifted_feature_list:
                features_to_check[feature_to_remap] = lifted_feature_list[feature_to_remap]
        features_to_remap = check_homologues(lifted_feature_list, features_to_check, parent_dict, original_parent_order)
    for feature in features_to_remap:
        unmapped_features.append(parent_dict[liftoff_utils.convert_id_to_original(feature)])
        del lifted_feature_list[feature]
    return lifted_feature_list



def add_overlapping_feature(features_to_remap, feature_to_remap, all_overlapping_features):
    overlap_start, overlap_end, overlap_seq, overlap_strand = features_to_remap[feature_to_remap][0], \
                                                              features_to_remap[feature_to_remap][1], \
                                                              features_to_remap[feature_to_remap][2], \
                                                              features_to_remap[feature_to_remap][3]
    if feature_to_remap in all_overlapping_features:
        all_overlapping_features[feature_to_remap].append((overlap_start, overlap_end,
                                             overlap_seq, overlap_strand))
    else:
        all_overlapping_features[feature_to_remap] = [(overlap_start, overlap_end,
                                             overlap_seq, overlap_strand)]


def clean_overlapping_features(lifted_feature_list, all_overlapping_features, parent_dict, features_to_remap, unmapped_features):
    parent_list = liftoff_utils.get_parent_list(lifted_feature_list, parent_dict)
    feature_by_chrom_dict = {}
    for feature in parent_list:
        if feature.seqid in feature_by_chrom_dict:
            feature_by_chrom_dict[feature.seqid].append(feature)
        else:
            feature_by_chrom_dict[feature.seqid] = [feature]
    for feature_name in all_overlapping_features:
        for overlapping_features in all_overlapping_features[feature_name]:
            if overlapping_features[2] not in feature_by_chrom_dict or  find_overlapping_features(overlapping_features, feature_by_chrom_dict[overlapping_features[2]]) is False:
                all_overlapping_features[feature_name].remove(overlapping_features)
                for feature in unmapped_features:
                    if feature.id == feature_name:
                        unmapped_features.remove(feature)
                        features_to_remap[feature.id]=[(-1, -1, None, None)]


def find_overlapping_features(overlapping_features, feature_list):
    for feature in feature_list:
        if feature.strand == overlapping_features[3] and feature.seqid == overlapping_features[2]:
            if liftoff_utils.count_overlap(feature.start, feature.end, overlapping_features[0], overlapping_features[1]) > 0:
                return True
    return False
