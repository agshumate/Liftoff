from collections import defaultdict
import warnings
from Bio.Seq import Seq
from liftoff import liftoff_utils
import parasail
from itertools import groupby
import numpy as np
import re


def polish_annotations(feature_list, ref_faidx, target_faidx, args, feature_heirarchy, target_feature):
    target_sub_features = get_sub_features(feature_list, target_feature)
    ref_sub_features = get_sub_features(feature_heirarchy.children, target_sub_features[0].id)
    ref_gene = feature_heirarchy.parents[target_sub_features[0].id]

    target_sub_features = feature_list[target_feature]
    transcript_list = [feature for feature in target_sub_features if "matches_ref_protein" in feature.attributes]
    good_transcripts = [tran for tran in transcript_list if  tran.attributes["valid_ORF"] == ["True"]]
    if len(transcript_list) != len(good_transcripts):
        output_sam = open(args.dir+"/polish.sam", 'w')
        write_sam_header(target_faidx, output_sam)
        target_gene = target_sub_features[0]
        polish_annotation(ref_gene, target_gene, ref_sub_features, target_sub_features, ref_faidx,
                          target_faidx,output_sam)
        return True
    return False


def get_sub_features(feature_list, feature_name):
    return feature_list[feature_name]


def find_and_check_cds(target_sub_features, ref_sub_features, ref_faidx, target_faidx, feature_list):
    cds_features = get_cds_features(target_sub_features)
    total_cds = 0
    num_good_cds = 0
    if len(cds_features) > 0:
        total_cds, num_good_cds = count_good_cds(cds_features, ref_faidx, target_faidx, ref_sub_features,
                                                 target_sub_features, feature_list)
        target_sub_features[0].attributes["valid_ORFs"] = [str(num_good_cds)]
    return total_cds, num_good_cds


def write_sam_header(target_fa, output_sam):
    for seq in target_fa.keys():
        output_sam.write("@SQ" + "\t" + "SN:" + seq + "\t" + "LN:" + str(len(target_fa[seq])) + "\n")


def get_cds_features(sub_features):
    return [feature for feature in sub_features if feature.featuretype == "CDS"]


def count_good_cds(cds_features, ref_faidx, target_faidx, ref_children, target_sub_features, feature_list):
    grouped_cds = group_cds_by_tran(cds_features)
    good_cds_count = 0
    for cds_group in grouped_cds:
        cds_seq = get_seq(cds_group, target_faidx)
        transcript = [feature for feature in target_sub_features if feature.id == cds_group[0].attributes["Parent"][
            0]][0]
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            protein = cds_seq.translate()
        if matches_ref(cds_group, ref_faidx, ref_children, protein):
            transcript.attributes["matches_ref_protein"] = ['True']
        else:
            transcript.attributes["matches_ref_protein"] = ['False']
        longest_ORF, longest_ORF_coords = get_longest_ORF(cds_seq)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            protein = longest_ORF.translate()
        transcript.attributes["valid_ORF"] = ['False']
        if len(protein) < 3:
            transcript.attributes["partial_ORF"] = ['True']
        elif missing_start(protein):
                transcript.attributes["missing_start_codon"] = ['True']
        elif missing_stop(protein):
            transcript.attributes["missing_stop_codon"] = ['True']
        elif inframe_stop(protein):
            transcript.attributes["inframe_stop_codon"] = ['True']
        else:
            transcript.attributes["valid_ORF"] = ['True']
            good_cds_count +=1
        if longest_ORF != cds_seq:
            adjust_cds_coords(cds_group,longest_ORF_coords, feature_list )
    return len(grouped_cds), good_cds_count


def get_longest_ORF(cds_seq):
    longest_ORF_seq = cds_seq
    longest_ORF = [0, len(cds_seq)]
    ORFs = [[m.start(), m.end()] for m in re.finditer(r'ATG(?:(?!TAA|TAG|TGA)...)*(?:TAA|TAG|TGA)',
                                                      str(cds_seq))]
    ORFs.sort(key=lambda x: x[1] - x[0])
    if len(ORFs) >0 :
        possible_ORF = ORFs[-1]
        if  possible_ORF[1] - possible_ORF[0] != len(cds_seq) and possible_ORF[1] - possible_ORF[0] >= \
                180:
            longest_ORF_seq = cds_seq[longest_ORF[0]:longest_ORF[1]]
            longest_ORF = possible_ORF
    return longest_ORF_seq, longest_ORF


def adjust_cds_coords(cds_features, longest_ORF_coords, feature_list):
    coords = np.hstack(np.array([list(np.arange(cds.start, cds.end +1)) for cds in cds_features], dtype=object))
    for cds in cds_features:
        relative_start = np.where(coords==cds.start)[0][0]
        relative_end = np.where(coords == cds.end)[0][0]
        if relative_end < longest_ORF_coords[0]:
            feature_list.remove([item for item in feature_list if item.id == cds.id][0])
        if relative_start > longest_ORF_coords[1]:
            feature_list.remove([item for item in feature_list if item.id == cds.id][0])
        if longest_ORF_coords[0] >= relative_start and  longest_ORF_coords[0] <= relative_end:
            cds.start = cds.start + (longest_ORF_coords[0] - relative_start)
        if relative_start <= longest_ORF_coords[1] and relative_end >= longest_ORF_coords[1]:
            cds.end = cds.end - (relative_end - longest_ORF_coords[1] +1)


def group_cds_by_tran(cds_features):
    groups = defaultdict(list)
    for cds in cds_features:
        groups[cds.attributes["Parent"][0]].append(cds)
    return groups.values()



def get_seq(cds_group, fasta_index):
    cds_group.sort(key = lambda x: x.start)
    cds_seq = ""
    chrom_seq = fasta_index[cds_group[0].seqid]
    for cds in cds_group:
        sub_seq = chrom_seq[cds.start -1:cds.end]
        cds_seq += sub_seq.seq
    if cds.strand == "-":
        cds_seq = Seq(cds_seq).reverse_complement()
    else:
        cds_seq = Seq(cds_seq)
    return cds_seq.upper()


def missing_start(protein):
    return protein[0] != "M"


def missing_stop(protein):
    return "*" != protein[-1]


def inframe_stop(protein):
    return "*" in protein[:-1]


def matches_ref(cds_group, ref_fa, ref_children, target_protein):
    tran = cds_group[0].attributes["Parent"]
    ref_cds_group= [ref_cds for ref_cds in ref_children if ref_cds.featuretype=="CDS" and ref_cds.attributes[
        "Parent"] ==
                    tran]
    ref_cds_seq = get_seq(ref_cds_group, ref_fa)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        ref_protein = ref_cds_seq.translate()
    return target_protein == ref_protein


def polish_annotation(ref_gene, target_gene, ref_children, target_children, ref_fa, target_fa, output_sam):
    ref_exons = [feature for feature in ref_children if feature.featuretype == "exon"]
    target_exons = [feature for feature in target_children if feature.featuretype == "exon"]
    ref_CDS = [feature for feature in ref_children if feature.featuretype == "CDS"]
    if len(ref_exons) == 0:
        ref_exons = ref_CDS
        target_exons = [feature for feature in target_children if feature.featuretype == "CDS"]
    ref_CDS_intervals = liftoff_utils.merge_children_intervals(ref_CDS)
    splice_sites = add_splice_sites(ref_exons, ref_gene)
    merged_ref_intervals = liftoff_utils.merge_children_intervals(ref_exons)
    exon_group_dict = find_overlapping_exon_groups(merged_ref_intervals, ref_exons)
    matrix = make_scoring_matrix(3)
    for i in range (len(merged_ref_intervals)):
        ref_interval =  merged_ref_intervals[i]
        target_interval = get_target_interval(exon_group_dict[i], target_exons, ref_interval)
        flank = max(0, (ref_interval[1] - ref_interval[0])- (target_interval[1]- target_interval[0])) + 10
        if target_interval != [0,0]:
            is_reverse = ref_gene.strand != target_gene.strand
            ref_seq = get_feature_sequence(ref_interval, ref_fa, ref_gene.seqid,  0 )
            capitalized_ref_seq = cds_and_splice_sites_to_upper(ref_interval, splice_sites, ref_CDS_intervals,
                                                                ref_seq, is_reverse)
            target_seq = get_feature_sequence(target_interval, target_fa, target_gene.seqid, flank)
            alignment = parasail.sg_dx_trace_scan_sat(capitalized_ref_seq, target_seq, 10, 1, matrix)
            write_sam_file(ref_interval, target_interval, ref_gene, target_gene, capitalized_ref_seq, alignment,
                           output_sam, flank)
    remove_splice_sites(ref_exons, ref_gene)


def add_splice_sites(exons, parent):
    splice_sites = []
    for exon in exons:
        if exon.start - 2 >= parent.start:
            exon.start = exon.start - 2
            splice_sites.append([exon.start, exon.start + 1])
        if exon.end + 2 <= parent.end:
            exon.end = exon.end + 2
            splice_sites.append([exon.end - 1, exon.end])
    return splice_sites

def remove_splice_sites(exons, parent):
    for exon in exons:
        if exon.start  != parent.start:
            exon.start = exon.start + 2
        if exon.end  != parent.end:
            exon.end = exon.end - 2



def find_overlapping_exon_groups(merged_ref_intervals, ref_exons):
    exon_group_dict = {}
    for i in range(len(merged_ref_intervals)):
        exon_group_dict[i] = []
        interval = merged_ref_intervals[i]
        for exon in ref_exons:
            if interval[0] <= exon.start <= interval[1]:
                exon_group_dict[i].append(exon.id)
    return exon_group_dict


def make_scoring_matrix(match_reward_coding):
    match_reward_non_coding = 2
    mismatch_penalty = -4
    matrix = parasail.matrix_create("ACGTacgt*", match_reward_non_coding, mismatch_penalty, True)
    matrix[0,4] = match_reward_coding
    matrix[1,5] = match_reward_coding
    matrix[2,6] = match_reward_coding
    matrix[3,7] = match_reward_coding
    matrix[4,0] = match_reward_coding
    matrix[5,1] = match_reward_coding
    matrix[6,2] = match_reward_coding
    matrix[7,3] = match_reward_coding
    return matrix


def get_target_interval(exon_group, target_exons, ref_interval):
    target_group = []
    for ref_exon in exon_group:
        matching_target_exon = [exon for exon in target_exons if exon.id == ref_exon]
        if len(matching_target_exon) !=0:
            target_group.append(matching_target_exon[0])
    if len(target_group) != 0:
        start = min([exon.start for exon in target_group])
        end = max([exon.end for exon in target_group])
    else:
        start = min([exon.start for exon in target_exons]) - (ref_interval[1]-ref_interval[0])
        end = max([exon.end for exon in target_exons]) + (ref_interval[1]-ref_interval[0])
    return [start, end]


def get_feature_sequence(interval, fa, chrm,  flank):
    chrom_seq = fa[chrm]
    flank_start = get_flank_start(interval[0], flank)
    flank_end = get_flank_end(interval[1], flank, len(chrom_seq))
    seq = chrom_seq[flank_start - 1: flank_end]
    final_seq = seq.seq.lower()
    return final_seq


def get_flank_start(start, flank):
    return round(max(1, start - flank))


def get_flank_end(end, flank, chrom_length):
    return round(min(end + flank, chrom_length))


def cds_and_splice_sites_to_upper(ref_interval, splice_sites, cds_intervals, ref_seq, is_reverse):
    cds_and_splice_sites = splice_sites + cds_intervals
    cap_segs = [[start -ref_interval[0], end-ref_interval[0]] for start, end in cds_and_splice_sites if start >=
                ref_interval[0] and end <= ref_interval[1]]
    prev_string = ref_seq
    for start,end in cap_segs:
        new_string = prev_string[:start] + prev_string[start:end+1].upper() + prev_string[end+1:]
        prev_string = new_string
    if is_reverse:
        return str(Seq(prev_string).reverse_complement())
    return prev_string


def write_sam_file(ref_interval, target_interval, ref_gene, target_gene, ref_seq, alignment, output_sam, flank ):
    reference_traceback = alignment.traceback.query
    target_traceback = alignment.traceback.ref
    is_reverse = ref_gene.strand != target_gene.strand
    if is_reverse:
        hard_clip_start = ref_gene.end - ref_interval[1]
    else:
        hard_clip_start = ref_interval[0] - ref_gene.start
    cigar_string, target_start_offset = make_cigar(reference_traceback, target_traceback, hard_clip_start,)
    if is_reverse:
        bit_flag = 16
    else:
        bit_flag = 0
    target_start_with_flank = get_flank_start(target_interval[0], flank)
    sam_fields = np.array([target_gene.id, bit_flag, target_gene.seqid,
                           target_start_with_flank + int(
        target_start_offset),
                           0, cigar_string,"*", 0,0, ref_seq.upper(), "*"])

    output_sam.write('\t'.join(sam_fields.astype(str)) + "\n")


def make_cigar(reference_traceback, target_traceback, hard_clip_start):
    expanded_cigar = ""
    chars = [pos for pos, char in enumerate(reference_traceback) if char != "-"]
    start = chars[0]
    for i in range(start, len(reference_traceback)):
        if reference_traceback[i].upper() == target_traceback[i].upper():
            expanded_cigar+= '='
        elif reference_traceback[i] == "-":
            expanded_cigar += 'D'
        elif reference_traceback[i] != "-" and target_traceback[i]=="-":
            expanded_cigar += 'I'
        elif reference_traceback[i] != "-" and target_traceback[i]!="-":
            expanded_cigar += 'X'
        else:
            print(reference_traceback[i], target_traceback[i])
    cigar = condense_cigar_string(expanded_cigar, hard_clip_start)
    return cigar, start


def condense_cigar_string(expanded_cigar, hard_clip_start):
    cigar_numbers = ([[str(sum(1 for _ in g)), k] for k, g in groupby(expanded_cigar)])
    if cigar_numbers[0][1] == 'D':
        cigar_numbers = cigar_numbers[1:]
    cigar = ''.join([i for sub in cigar_numbers for i in sub])
    if hard_clip_start > 0:
        cigar = str(hard_clip_start) + 'H' + cigar
    return cigar
