from multiprocessing import Pool
import math
from functools import partial
import numpy as np
from pyfaidx import Fasta, Faidx
import subprocess
import pysam
from liftoff import aligned_seg, liftoff_utils
from os import path


def align_features_to_target(ref_chroms, target_chroms, args, feature_hierarchy, liftover_type, unmapped_features):
    if args.subcommand == "polish":
        sam_files = [args.dir + "/polish.sam"]
    else:
        target_fasta_dict = split_target_sequence(target_chroms, args.target, args.dir)
        genome_size = get_genome_size(target_fasta_dict)
        threads_per_alignment = max(1, math.floor(int(args.p) / len(ref_chroms)))
        sam_files = []
        pool = Pool(int(args.p))
        print("aligning features")
        func = partial(align_single_chroms, ref_chroms, target_chroms, threads_per_alignment, args, genome_size,
                       liftover_type)
        for result in pool.imap_unordered(func, np.arange(0, len(target_chroms))):
            sam_files.append(result)
        pool.close()
        pool.join()
    return parse_all_sam_files(feature_hierarchy, unmapped_features, liftover_type, sam_files)


def split_target_sequence(target_chroms, target_fasta_name, inter_files):
    Faidx(target_fasta_name)
    target_fasta_dict = Fasta(target_fasta_name, key_function=lambda x: x.split()[0])
    for chrm in target_chroms:
        if chrm != target_fasta_name:
            out = open(inter_files + "/" + chrm + ".fa", 'w')
            out.write(">" + chrm + "\n" + str(target_fasta_dict[chrm]))
    return target_fasta_dict


def get_genome_size(target_fasta_dict):
    genome_size = 0
    for value in target_fasta_dict.values():
        genome_size += len(value)
    return genome_size


def align_single_chroms(ref_chroms, target_chroms, threads, args, genome_size, liftover_type, index):
    max_single_index_size = 4000000000
    features_file, features_name = get_features_file(ref_chroms, args, liftover_type, index)
    target_file, output_file = get_target_file_and_output_file(liftover_type, target_chroms, index, features_name, args)
    threads_arg = str(threads)
    minimap2_path = get_minimap_path(args)
    target_prefix = get_target_prefix_name(target_chroms, index, args, liftover_type)
    if genome_size > max_single_index_size:
        split_prefix = args.dir + "/" + features_name + "_to_" + target_prefix + "_split"
        command = [minimap2_path, '-o', output_file, target_file, features_file] + args.mm2_options.split(" ") + [
            "--split-prefix", split_prefix, '-t', threads_arg]
        subprocess.run(command)
    else:
        minimap2_index = build_minimap2_index(target_file, args, threads_arg, minimap2_path)
        command = [minimap2_path, '-o', output_file, minimap2_index, features_file] + args.mm2_options.split(" ") + [
            '-t', threads_arg]
        subprocess.run(command)
    return output_file


def get_features_file(ref_chroms, args, liftover_type, index):
    if ref_chroms[index] == args.reference and (liftover_type == "chrm_by_chrm" or liftover_type == "copies"):
        features_name = 'reference_all'
    elif liftover_type == "unmapped":
        features_name = "unmapped_to_expected_chrom"
    elif liftover_type == "unplaced":
        features_name = "unplaced"
    else:
        features_name = ref_chroms[index]
    return args.dir + "/" + features_name + "_genes.fa", features_name


def get_target_file_and_output_file(liftover_type, target_chroms, index, features_name, args):
    if liftover_type != "chrm_by_chrm" or target_chroms[0] == args.target:
        target_file = args.target
        out_file_target = "target_all"
    else:
        target_file = args.dir + "/" + target_chroms[index] + ".fa"
        out_file_target = target_chroms[index]
    output_file = args.dir + "/" + features_name + "_to_" + out_file_target + ".sam"
    return target_file, output_file


def get_minimap_path(args):
    if args.m is None:
        minimap2 = "minimap2"
    else:
        minimap2 = args.m
    return minimap2


def get_target_prefix_name(target_chroms, index, args, liftover_type):
    if liftover_type != "chrm_by_chrm" or target_chroms[0] == args.target:
        prefix = "target_all"
    else:
        prefix = target_chroms[index]
    return prefix


def build_minimap2_index(target_file, args, threads, minimap2_path):
    if path.exists(target_file + ".mmi") is False:
        subprocess.run(
            [minimap2_path, '-d', target_file + ".mmi", target_file] + args.mm2_options.split(" ") + ['-t',
             threads ])
    return target_file + ".mmi"


def parse_all_sam_files(feature_hierarchy, unmapped_features, liftover_type, sam_files):
    aligned_segments_dict = {}
    for file in sam_files:
        aligned_segments = parse_alignment(file, feature_hierarchy, unmapped_features, liftover_type)
        aligned_segments_dict.update(aligned_segments)
    return aligned_segments_dict


def parse_alignment(file, feature_hierarchy, unmapped_features, search_type):
    all_aligned_blocks = {}
    sam_file = pysam.AlignmentFile(file, 'r', check_sq=False, check_header=False)
    sam_file_iter = sam_file.fetch()
    aln_id = 0
    name_dict = {}
    align_count_dict = {}
    for ref_seq in sam_file_iter:
        if ref_seq.is_unmapped is False:
            aln_id = add_alignment(ref_seq,  align_count_dict, search_type, name_dict,aln_id, feature_hierarchy,
                                   all_aligned_blocks)
        else:
            unmapped_features.append(feature_hierarchy.parents[ref_seq.query_name])
    remove_alignments_without_children(all_aligned_blocks, unmapped_features, feature_hierarchy)
    return all_aligned_blocks


def add_alignment(ref_seq, align_count_dict, search_type, name_dict, aln_id, feature_hierarchy,
                  all_aligned_blocks):
    ref_seq.query_name = edit_name(search_type, ref_seq, name_dict)
    aln_id += 1
    if ref_seq.query_name in align_count_dict:
        align_count = align_count_dict[ref_seq.query_name] + 1
    else:
        align_count = 0
    align_count_dict[ref_seq.query_name] = align_count
    aligned_blocks = get_aligned_blocks(ref_seq, aln_id, feature_hierarchy, search_type)
    if ref_seq.query_name in all_aligned_blocks:
        all_aligned_blocks[ref_seq.query_name].extend(aligned_blocks)
    else:
        all_aligned_blocks[ref_seq.query_name] = aligned_blocks
    return aln_id


def edit_name(search_type, ref_seq, name_dict):
    if search_type != "copies":
        return ref_seq.query_name + "_0"
    else:
        if ref_seq.query_name not in name_dict:
            name_dict[ref_seq.query_name] = 0
        name_dict[ref_seq.query_name] += 1
        return ref_seq.query_name + "_" + str(name_dict[ref_seq.query_name])


def get_aligned_blocks(alignment, aln_id, feature_hierarchy, search_type):
    cigar_operations = get_cigar_operations()
    cigar = alignment.cigar
    parent = feature_hierarchy.parents[liftoff_utils.convert_id_to_original(alignment.query_name)]
    query_start, query_end = get_query_start_and_end(alignment, cigar, cigar_operations)
    children = feature_hierarchy.children[liftoff_utils.convert_id_to_original(alignment.query_name)]
    end_to_end = is_end_to_end_alignment(parent, query_start, query_end)
    if search_type == "copies" and end_to_end is False:
        return []
    reference_block_start, reference_block_pos = alignment.reference_start, alignment.reference_start
    query_block_start, query_block_pos = query_start, query_start
    new_blocks, mismatches = [], []
    merged_children_coords = liftoff_utils.merge_children_intervals(children)
    for operation, length in cigar:
        if base_is_aligned(operation, cigar_operations):
            query_block_pos, reference_block_pos = add_aligned_base(operation, query_block_pos, reference_block_pos,
                                                                    length, cigar_operations, mismatches)
            if query_block_pos == query_end:
                add_block(query_block_pos, reference_block_pos, aln_id, alignment, query_block_start,
                          reference_block_start, mismatches, new_blocks, merged_children_coords, parent)
                break
        elif is_alignment_gap(operation, cigar_operations):
            add_block(query_block_pos, reference_block_pos, aln_id, alignment, query_block_start, reference_block_start,
                      mismatches, new_blocks, merged_children_coords, parent)
            mismatches, query_block_start, reference_block_start, query_block_pos, reference_block_pos = \
                end_block_at_gap(
                    operation, query_block_pos, reference_block_pos, length, cigar_operations)
    return new_blocks


def get_cigar_operations():
    return {"insertion": 1, "deletion": 2, "hard_clip": 5, "match": 7, "mismatch": 8}


def get_query_start_and_end(alignment, cigar, cigar_operations):
    query_start = alignment.query_alignment_start
    query_end = alignment.query_alignment_end
    if cigar[0][0] == cigar_operations["hard_clip"]:
        query_start += cigar[0][1]
        query_end += cigar[0][1]
    return query_start, query_end


def is_end_to_end_alignment(parent, query_start, query_end):
    return parent.end - parent.start + 1 == query_end - query_start


def base_is_aligned(operation, cigar_operations):
    return operation == cigar_operations["match"] or operation == cigar_operations["mismatch"]


def add_aligned_base(operation, query_block_pos, reference_block_pos, length, cigar_operations, mismatches):
    if operation == cigar_operations["mismatch"]:
        for i in range(query_block_pos, query_block_pos + length):
            mismatches.append(i)
    query_block_pos, reference_block_pos = adjust_position(operation, query_block_pos, reference_block_pos,
                                                           length, cigar_operations)
    return query_block_pos, reference_block_pos


def adjust_position(operation, query_block_pos, reference_block_pos, length, cigar_operations):
    if operation == cigar_operations["match"] or operation == cigar_operations["mismatch"] or operation == \
            cigar_operations["insertion"]:
        query_block_pos += length
    if operation == cigar_operations["match"] or operation == cigar_operations["mismatch"] or operation == \
            cigar_operations["deletion"]:
        reference_block_pos += length
    return query_block_pos, reference_block_pos


def add_block(query_block_pos, reference_block_pos, aln_id, alignment, query_block_start, reference_block_start,
              mismatches, new_blocks, merged_children_coords, parent):
    query_block_end = query_block_pos - 1
    reference_block_end = reference_block_pos - 1
    new_block = aligned_seg.aligned_seg(aln_id, alignment.query_name, alignment.reference_name, query_block_start,
                                        query_block_end,
                                        reference_block_start, reference_block_end, alignment.is_reverse,
                                        np.array(mismatches).astype(int))
    overlapping_children = find_overlapping_children(new_block, merged_children_coords, parent)
    if overlapping_children != []:

        new_blocks.append(new_block)




def find_overlapping_children(aln, children_coords, parent):
    overlapping_children = []
    for child_interval in children_coords:
        relative_start = liftoff_utils.get_relative_child_coord(parent, child_interval[0], aln.is_reverse)
        relative_end = liftoff_utils.get_relative_child_coord(parent, child_interval[1], aln.is_reverse)
        child_start, child_end = min(relative_start, relative_end), max(relative_start, relative_end)
        overlap = liftoff_utils.count_overlap(child_start, child_end, aln.query_block_start, aln.query_block_end)
        if overlap > 0:
            overlapping_children.append(child_start)
            overlapping_children.append(child_end)
    return overlapping_children


def is_alignment_gap(operation, cigar_operations):
    return operation == cigar_operations["insertion"] or operation == cigar_operations["deletion"]


def end_block_at_gap(operation, query_block_pos, reference_block_pos, length, cigar_operations):
    mismatches = []
    query_block_pos, reference_block_pos = adjust_position(operation, query_block_pos, reference_block_pos,
                                                           length, cigar_operations)
    query_block_start = query_block_pos
    reference_block_start = reference_block_pos
    return mismatches, query_block_start, reference_block_start, query_block_pos, reference_block_pos


def remove_alignments_without_children(all_aligned_blocks, unmapped_features, feature_hierarchy):
    features_to_remove = []
    for seq in all_aligned_blocks:
        if all_aligned_blocks[seq] == []:
            features_to_remove.append(seq)
            unmapped_features.append(feature_hierarchy.parents[liftoff_utils.convert_id_to_original(seq)])
    for feature in features_to_remove:
        del all_aligned_blocks[feature]
    return all_aligned_blocks
