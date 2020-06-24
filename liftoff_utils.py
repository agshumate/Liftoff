import numpy as np


def count_overlap(start1, end1, start2, end2):
    overlap = min(end1, end2) - max(start1, start2)
    return overlap

def get_relative_child_coord(parent, coord, is_reverse):
    if is_reverse:
        relative_coord = parent.end - coord
    else:
        relative_coord = coord - parent.start
    return relative_coord


def merge_children_intervals(children):
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



def get_parent_list(feature_list, parent_dict):
    parent_list = []
    for feature in feature_list:
        for value in feature_list[feature]:
            if value.id in parent_dict:
                parent_list.append(value)
    return parent_list


def has_child(feature, gene_db):
    for child in gene_db.children(feature.id):
        return True
    return False

def clear_scores(feature_list, parent_dict):
    for parent in feature_list:
        for feature in feature_list[parent]:
            if feature.id in parent_dict:
                feature.score = -1
    return feature_list


def find_parent_order(parents):
    parents.sort(key=lambda x: (x.seqid, x.start))
    return np.array([parent.id for parent in parents])


def convert_id_to_original(id):
    copy_tag_len = len(id.split("_")[-1])
    original_parent_name = id[:-copy_tag_len-1]
    return original_parent_name


def get_copy_tag(id):
    copy_tag_len = len(id.split("_")[-1])
    copy_tag = id[-copy_tag_len-1:]
    return copy_tag


def make_new_feature(copied_feature, new_start, new_end, new_strand, new_seqid):
    copied_feature.seqid = new_seqid
    copied_feature.source="Liftoff"
    copied_feature.start=new_start
    copied_feature.end=new_end
    copied_feature.strand=new_strand
    return copied_feature
