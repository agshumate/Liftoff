from multiprocessing import Pool
import math
from functools import partial
import numpy as np
from pyfaidx import Fasta, Faidx
import subprocess
import pysam
from liftoff import aligned_seg, liftoff_utils


def align_features_to_target(ref_chroms, target_chroms, processes, target_fasta, parent_dict, children_dict,
                             search_type, unmapped_features, reference_fasta, minimap2_path, inter_files, map):
    print("aligning features")
    genome_size=split_target_sequence(target_chroms,target_fasta, inter_files )
    aligned_segments_dict = {}
    threads_per_alignment = max(1, math.floor(processes / len(ref_chroms)))
    sam_files = []
    pool = Pool(processes)
    func = partial(align_subset, ref_chroms, target_chroms, threads_per_alignment, target_fasta, reference_fasta, minimap2_path, inter_files, map, genome_size, search_type)
    for result in pool.imap_unordered(func, np.arange(0,len(target_chroms))):
        sam_files.append(result)
    pool.close()
    pool.join()
    for file in sam_files:
        aligned_segments = parse_alignment(file, parent_dict, children_dict, unmapped_features, search_type)
        aligned_segments_dict.update(aligned_segments)
    return aligned_segments_dict


def split_target_sequence(target_chroms, target_fasta_name, inter_files):
    Faidx(target_fasta_name)
    genome_size =0
    target_fasta = Fasta(target_fasta_name, key_function = lambda x: x.split()[0])
    for value in target_fasta.values():
        genome_size += len(value)
    for chrm in target_chroms:
        if chrm != target_fasta_name:
            out=open( inter_files + "/" + chrm+".fa", 'w')
            out.write(">" + chrm + "\n" + str(target_fasta[chrm]))
    return genome_size



def align_subset(ref_chroms, target_chroms, threads, target_fasta_name, reference_fasta_name, minimap2_path, inter_files, map, genome_size, liftover_type, index):
    if ref_chroms[index] == reference_fasta_name and (liftover_type == "chrm_by_chrm" or liftover_type == "copies"):
        features_name = 'reference_all'
    elif liftover_type == "unmapped":
        features_name = "unmapped_to_expected_chrom"
    elif liftover_type == "unplaced":
        features_name = "unplaced"
    else:
        features_name = ref_chroms[index]
    features_file =  inter_files+ "/" + features_name + "_genes.fa"
    if liftover_type != "chrm_by_chrm" or target_chroms[0] == target_fasta_name:
        target_file = target_fasta_name
        out_file_target = "target_all"
    else:
        target_file = inter_files + "/" +target_chroms[index] + ".fa"
        out_file_target = target_chroms[index]
    out_arg = inter_files + "/"+ features_name + "_to_" + out_file_target + ".sam"
    threads_arg =  str(threads)
    if map:
        if minimap2_path is None:
            minimap2 = "minimap2"
        else:
            minimap2= minimap2_path
        if genome_size > 4000000000:
            split_prefix = features_name + "_to_" + out_file_target + "_split"
            subprocess.run([minimap2, '-o', out_arg, target_file, features_file, '-a', '--eqx', '-N', '50', '-p',  '0.5', '-t', threads_arg,"--split-prefix", split_prefix],
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        else:
            subprocess.run(
                [minimap2, '-o', out_arg, target_file, features_file, '-a', '--eqx', '-N', '50', '-p', '0.5', '-t',
                 threads_arg])

    return inter_files + "/"+ features_name + "_to_" + out_file_target + ".sam"


def edit_name(search_type, ref_seq, name_dict):
    if search_type != "copies":
        return ref_seq.query_name + "_1"
    else:
        if ref_seq.query_name not in name_dict:
            name_dict[ref_seq.query_name] = 0
        name_dict[ref_seq.query_name] += 1
        return ref_seq.query_name + "_" + str(name_dict[ref_seq.query_name])

def parse_alignment(file, parent_dict, children_dict, unmapped_features, search_type):
    all_aligned_blocks ={}
    sam_file = pysam.AlignmentFile(file,'r',check_sq=False, check_header=False)
    sam_file_iter= sam_file.fetch()
    aln_id = 0
    name_dict = {}
    for ref_seq in sam_file_iter:
        if ref_seq.is_unmapped is False:
            ref_seq.query_name = edit_name(search_type, ref_seq, name_dict)
            aln_id += 1
            aligned_blocks = get_aligned_blocks(ref_seq, aln_id, children_dict, parent_dict, search_type)
            if ref_seq.query_name in all_aligned_blocks:
                all_aligned_blocks[ref_seq.query_name].extend(aligned_blocks)
            else:
                all_aligned_blocks[ref_seq.query_name]=aligned_blocks
        else:
            unmapped_features.append(parent_dict[ref_seq.query_name])
    unaligned_exons = []
    for seq in all_aligned_blocks:
        if all_aligned_blocks[seq]==[]:
            unaligned_exons.append(seq)
            unmapped_features.append(parent_dict[liftoff_utils.convert_id_to_original(seq)])
    for seq in unaligned_exons:
        del all_aligned_blocks[seq]
    return all_aligned_blocks


def get_aligned_blocks(alignment, aln_id, children_dict, parent_dict, search_type):
    cigar = alignment.cigar
    query_start = alignment.query_alignment_start
    query_end = alignment.query_alignment_end
    children = children_dict[liftoff_utils.convert_id_to_original(alignment.query_name)]
    parent = parent_dict[liftoff_utils.convert_id_to_original(alignment.query_name)]
    if parent.end - parent.start +1 != query_end - query_start and search_type == "copies":
        return []
    reference_block_start = alignment.reference_start
    reference_block_pos = reference_block_start
    if cigar[0][0] == 5:
        query_start += cigar[0][1]
        query_end += cigar[0][1]
    query_block_start = query_start
    query_block_pos = query_block_start
    new_blocks = []
    mismatches = []
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


