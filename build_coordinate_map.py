import numpy as np
import process_blast_alignments as pba

def build_coordinate_map(exon_alignments, gene, gene_db):
    query_length = gene.end - gene.start + 1
    coordinate_map = np.zeros(shape=(query_length, len(exon_alignments)))
    alignment_score_map = np.zeros(shape=(query_length, len(exon_alignments))) + 2
    alignment_score_weights = build_weight_map(query_length, gene, gene_db)
    for i in range(len(exon_alignments)):
        add_hsp_to_map(exon_alignments[i][0], coordinate_map, alignment_score_map, i)
    return coordinate_map, alignment_score_map * alignment_score_weights


def add_hsp_to_map(hsp, query_to_target, alignment_scores, hsp_num):
    query = np.array(list(hsp.query))
    subject = np.array(list(hsp.sbjct))
    if hsp.sbjct_start < hsp.sbjct_end:
        direction = 1
    else:
        direction = -1
    q_start = hsp.query_start - 1
    q_end = hsp.query_end - 1
    s_start = hsp.sbjct_start -1
    query_gaps = np.where(query == '-')
    subject_gaps = np.where(subject == '-')
    mismatches = np.where((query != '-') & (subject != '-') & (query != subject))
    subject_range = np.arange(s_start, s_start +
                              direction * (q_end - q_start + 1), direction)
    query_to_target[q_start: q_end + 1, hsp_num]= subject_range
    alignment_scores[q_start: q_end + 1, hsp_num] = 0
    adjust_coordinates_for_gaps(query_to_target, alignment_scores, subject_gaps, query_gaps, mismatches, hsp_num,
                                q_start, q_end, direction, hsp.query)


def adjust_coordinates_for_gaps(query_to_target, alignment_scores, subject_gaps, query_gaps, mismatches,
                                hsp_num, q_start, q_end, direction, query):
    for s_gap in subject_gaps[0]:
        relative_s_gap = s_gap - query[0:s_gap].count('-')
        query_to_target[q_start + relative_s_gap: q_end+1 , hsp_num] -= direction
        alignment_scores[q_start + relative_s_gap, hsp_num] = 2
    for q_gap in query_gaps[0]:
        relative_q_gap = q_gap - query[0:q_gap].count('-')
        query_to_target[q_start + relative_q_gap: q_end + 1, hsp_num] += direction
    for mismatch in mismatches[0]:
        alignment_scores[q_start + mismatch-query[0:mismatch].count("-"), hsp_num] = 1


def build_weight_map(query_length, gene, gene_db):
    weight_map = np.zeros(shape=(query_length ,1)) + 1
    all_exons, all_cds = pba.get_exons_and_cds(gene_db, gene)
    for exon in all_exons:
        cds_list = pba.find_cds(exon, all_cds)
        for cds in cds_list:
            adjust_cds_weight(cds, gene, query_length, weight_map)
    adjust_start_and_stop_codon_weights(weight_map)
    return weight_map



def adjust_cds_weight(cds, gene, query_length, weight_map):
    relative_start = cds.start - gene.start
    relative_end = cds.end - gene.start
    if gene.strand == "+":
        weight_map[relative_start:relative_end +1] = 2
    else:
        weight_map[(query_length -1) - relative_end: (query_length -1) - relative_start + 1] = 2


def adjust_start_and_stop_codon_weights(weight_map):
    start_codon = np.where(weight_map ==2)[0][0:3]
    stop_codon = np.where(weight_map ==2)[0][-3:]
    for base in start_codon:
        weight_map[base] = 4
    for base in stop_codon:
        weight_map[base]=4

