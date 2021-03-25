import numpy as np


def count_overlap(start1, end1, start2, end2):
    overlap = min(end1, end2) - max(start1, start2) +1
    return overlap


def get_relative_child_coord(parent, coord, is_reverse):
    if is_reverse:
        relative_coord = parent.end - coord
    else:
        relative_coord = coord - parent.start
    return relative_coord


def merge_children_intervals(children):
    if len(children) == 0:
        return []
    intervals = [[child.start, child.end] for child in children]
    intervals.sort(key=lambda interval: interval[0])
    merged = [intervals[0]]
    for current in intervals:
        previous = merged[-1]
        if current[0] <= previous[1]:
            previous[1] = max(previous[1], current[1])
        else:
            merged.append(current)
    return merged


def get_parent_list(feature_list):
    parent_list = []
    for feature in feature_list:
        feature_set = feature_list[feature]
        if len(feature_set) > 0:
            parent_list.append(feature_set[0])
    return parent_list


def clear_scores(feature_list, parent_dict):
    for parent in feature_list:
        for feature in feature_list[parent]:
            if feature.id in parent_dict:
                feature.score = -1
    return feature_list


def find_parent_order(parents):
    parents.sort(key=lambda x: (x.seqid, x.start))
    return np.array([[parent.id, parent] for parent in parents])


def convert_id_to_original(id):
    frag_split = id.split("_frag")[0]
    copy_tag_len = len(frag_split.split("_")[-1])
    original_parent_name = frag_split[:-copy_tag_len - 1]
    return original_parent_name


def get_copy_tag(id):
    copy_tag_len = len(id.split("_")[-1])
    copy_tag = id[-copy_tag_len - 1:]
    return copy_tag


def get_strand(aln, parent):
    if aln.is_reverse:
        if parent.strand == "-":
            strand = '+'
        else:
            strand = "-"
    else:
        strand = parent.strand
    return strand


def find_overlaps(start, end, chrm, strand, feature_name, intervals, parent_dict, lifted_features_list,
                  max_overlap_non_copies):
    all_overlaps = intervals.find((start, end))
    incorrect_overlaps = []
    filtered_overlaps = [overlap for overlap in all_overlaps if
                         overlap[2][1].seqid == chrm and overlap[2][1].strand == strand]

    for overlap in filtered_overlaps:
        shortest_feature_length = min(end - start, overlap[1] - overlap[0]) +1
        overlap_amount = count_overlap(start, end, overlap[0], overlap[1])
        ref_feature = parent_dict[convert_id_to_original(feature_name)]
        ref_overlap_feature = parent_dict[convert_id_to_original(overlap[2][0])]
        if get_copy_tag(feature_name) != '_0' or get_copy_tag(overlap[2][0]) != '_0':
            max_overlap = 0
        else:
            max_overlap = max_overlap_non_copies
        if overlaps_in_ref_annotation(ref_feature, ref_overlap_feature) is False and overlap[2][0] in \
                lifted_features_list and overlap_amount / shortest_feature_length > max_overlap:
            incorrect_overlaps.append(overlap)
    return incorrect_overlaps


def overlaps_in_ref_annotation(ref_feature1, ref_feature2):
    if ref_feature1.seqid != ref_feature2.seqid:
        return False
    if ref_feature1.strand != ref_feature2.strand:
        return False
    if ref_feature1.id == ref_feature2.id:
        return False
    else:
        return count_overlap(ref_feature1.start, ref_feature1.end,
                             ref_feature2.start,
                             ref_feature2.end) > 0


def find_nonoverlapping_upstream_neighbor(parent_order, feature_name):
    feature_location = np.where(parent_order[:, 0] == feature_name)[0][0]
    neighbor_indx = feature_location - 1
    neighbor = parent_order[neighbor_indx][1]
    neighbor_id = parent_order[neighbor_indx][0]
    feature = parent_order[neighbor_indx + 1][1]
    if neighbor.seqid != feature.seqid or neighbor_indx < 0:
        return None
    else:
        return neighbor_id
