from liftoff import align_features
import numpy as np

def output_variants(alignments_used):
    find_variants(alignments_used)


def find_variants(alignments_used):
    cigar_operations = align_features.get_cigar_operations()
    for feature in alignments_used:
        alignment_list = alignments_used[feature]
        query_length = max([alignment.query_length for alignment in alignment_list])
        expanded_cigar = np.zeros((1,query_length))
        for alignment in alignment_list:
            cigar = alignment.cigar
            feature_index = alignment.query_alignment_start
            for operation, len in cigar:
                for i in range(len):






