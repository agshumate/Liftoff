from multiprocessing import Pool
import math
from functools import partial
import numpy as np
from pyfaidx import Fasta, Faidx
import subprocess
import pysam
import liftoff_utils
import aligned_seg
import copy
import os



def align_features_to_target(ref_chroms, target_chroms, processes, target_fasta, parent_dict, children_dict,
                             search_type, unmapped_features):
    print("aligning features")
    split_target_sequence(target_chroms,target_fasta )
    threads_per_alignment = max(1, math.floor(processes / len(ref_chroms)))
    sam_files = []
    aligned_segments_dict = {}
    pool = Pool(processes)
    func = partial(align_subset, ref_chroms, target_chroms, threads_per_alignment, target_fasta)
    for result in pool.imap_unordered(func, np.arange(0,len(ref_chroms))):
        sam_files.append(result)
    pool.close()
    pool.join()
    for file in sam_files:
        aligned_segments = parse_alignment(file, parent_dict, children_dict, unmapped_features)
        aligned_segments_dict.update(aligned_segments)
    if search_type == "copies":
        aligned_segments = expand_results(aligned_segments_dict)
    else:
        aligned_segments = rename_aligned_segs(aligned_segments_dict)
    return aligned_segments


def split_target_sequence(target_chroms, target_fasta_name):
    Faidx(target_fasta_name)
    target_fasta = Fasta(target_fasta_name, key_function = lambda x: x.split()[0])
    for chrm in target_chroms:
        if chrm != target_fasta_name:
            out=open(chrm+".fa", 'w')
            out.write(">" + chrm + "\n" + str(target_fasta[chrm]))
            out.close()


def align_subset(ref_chroms, target_chroms, threads, target_fasta_name, index):
    features_file = ref_chroms[index] + "_genes.fa"
    if target_chroms[index] == target_fasta_name:
        target_file = target_fasta_name
    else:
        target_file = target_chroms[index] + ".fa"
    out_arg = "-o"+ features_file + "_to_" + target_file
    threads_arg = "-t" + str(threads)
    subprocess.run(['minimap2', out_arg, target_file, features_file, '-a', '--eqx', '-N50', '-p 0.5', threads_arg],
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    #if target_file != target_fasta_name:
        #os.remove(target_file)
    #os.remove(features_file)
    return features_file + "_to_" + target_file


def parse_alignment(file, parent_dict, children_dict, unmapped_features):
    all_aligned_blocks ={}
    sam_file = pysam.AlignmentFile(file,'r')
    sam_file_iter= sam_file.fetch()
    aln_id = 0
    for ref_seq in sam_file_iter:
        if ref_seq.is_unmapped is False:
            aln_id += 1
            aligned_blocks = get_aligned_blocks(ref_seq, aln_id, children_dict, parent_dict)
            if aligned_blocks != []:
                if ref_seq.query_name in all_aligned_blocks:
                    all_aligned_blocks[ref_seq.query_name].extend(aligned_blocks)
                else:
                    all_aligned_blocks[ref_seq.query_name]=aligned_blocks
        else:
            unmapped_features.append(parent_dict[ref_seq.query_name])
    #os.remove(file)
    return all_aligned_blocks


def get_aligned_blocks(alignment, aln_id, children_dict, parent_dict):
    cigar = alignment.cigar
    query_start = alignment.query_alignment_start
    query_end = alignment.query_alignment_end
    reference_block_start = alignment.reference_start
    reference_block_pos = reference_block_start
    if cigar[0][0] == 5:
        query_start += cigar[0][1]
        query_end += cigar[0][1]
    query_block_start = query_start
    query_block_pos = query_block_start
    new_blocks = []
    mismatches = []
    children = children_dict[alignment.query_name]
    parent = parent_dict[alignment.query_name]
    merged_children_coords = liftoff_utils.merge_children_intervals(children)
    for operation, length in cigar:
        if operation == 7 or operation == 8:
            query_block_pos, reference_block_pos = adjust_position(operation, query_block_pos, reference_block_pos,
                                                                   length)
            if operation == 8:
                mismatches.append(query_block_pos)
            if query_block_pos == query_end:
                add_block(query_block_pos, reference_block_pos, aln_id, alignment, query_block_start,
                          reference_block_start, mismatches, new_blocks, merged_children_coords, parent)
                mismatches = []
        elif operation ==1 or operation ==2:
            add_block(query_block_pos, reference_block_pos, aln_id, alignment, query_block_start, reference_block_start,
                      mismatches, new_blocks, merged_children_coords, parent)
            mismatches = []
            query_block_pos, reference_block_pos = adjust_position(operation, query_block_pos, reference_block_pos,
                                                                   length)
            query_block_start = query_block_pos
            reference_block_start = reference_block_pos
    return new_blocks


def adjust_position(operation, query_block_pos, reference_block_pos, length):
    if operation == 7 or operation == 8 or operation ==1:
        query_block_pos += length
    if operation == 7 or operation ==8 or operation ==2:
        reference_block_pos += length
    return query_block_pos, reference_block_pos


def add_block(query_block_pos, reference_block_pos, aln_id, alignment, query_block_start, reference_block_start,
                  mismatches, new_blocks, merged_children_coords, parent):
    query_block_end = query_block_pos - 1
    reference_block_end = reference_block_pos - 1
    new_block = aligned_seg.aligned_seg(aln_id, alignment.query_name, alignment.reference_name, query_block_start, query_block_end,
                                        reference_block_start, reference_block_end, alignment.is_reverse,
                                        np.array(mismatches).astype(int))
    if contains_child(new_block, merged_children_coords, parent):
        new_blocks.append(new_block)



def contains_child(aln, children_coords, parent):
    for child_interval in children_coords:
        relative_start = liftoff_utils.get_relative_child_coord(parent, child_interval[0], aln.is_reverse)
        relative_end = liftoff_utils.get_relative_child_coord(parent, child_interval[1], aln.is_reverse)
        child_start, child_end = min(relative_start, relative_end), max(relative_start, relative_end)
        overlap = liftoff_utils.count_overlap(child_start, child_end, aln.query_block_start, aln.query_block_end)
        if overlap > 0:
            return True
    return False


def expand_results(all_aligned_segs):
    new_aligned_seg_dict = {}
    for ref_seq in all_aligned_segs:
        num_aligned_segs = count_num_segs(all_aligned_segs[ref_seq])
        for i in range(num_aligned_segs):
            seg_id = str(i +2)
            new_aligned_seg_name = ref_seq + "_"+ seg_id
            new_aligned_seg = create_new_aligned_seg(all_aligned_segs[ref_seq], i)
            new_aligned_seg.query = new_aligned_seg_name
            new_aligned_seg_dict[new_aligned_seg] = new_aligned_seg_name
    return new_aligned_seg_dict


def count_num_segs(ref_seq):
    return len(ref_seq)


def create_new_aligned_seg(ref_seq, aln_num):
    new_aligned_seg = copy.deepcopy(ref_seq)
    aln_count = -1
    for aln in new_aligned_seg[:]:
        aln_count += 1
        if aln_count != aln_num:
                new_aligned_seg.remove(aln)
    return new_aligned_seg

def rename_aligned_segs(all_aligned_segs):
    renamed_aligned_segs = {}
    for aln in all_aligned_segs:
        new_name = aln+ "_0"
        renamed_aligned_segs[new_name] = all_aligned_segs[aln]
        for aln in renamed_aligned_segs[new_name]:
            aln.query_name = new_name
    return renamed_aligned_segs
