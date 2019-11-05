from __future__ import print_function
import numpy as np
from Bio.Blast import NCBIXML
from Bio.Blast import Record
from build_new_annotation import add_gene_and_tran_coords
from build_new_annotation import print_gff
from convert_genome_coords import liftover_genes

# Adjust the coordinates of the map for gaps present in the query and the subject. At each alignment position, keep
# track of where the alignment is a gap, match or mismatch
def adjust_coordinates_for_gaps(query_to_target, alignment_types, subject_gaps, query_gaps, mismatches,
                                hsp_num, q_start, q_end, direction, query):
    for s_gap in subject_gaps[0]:
        relative_s_gap = np.sum(query[0:s_gap] != '-')
        query_to_target[q_start + relative_s_gap: q_end+1 , hsp_num] -= direction
        alignment_types[q_start + relative_s_gap, hsp_num] = 2
    for q_gap in query_gaps[0]:
        relative_q_gap = np.sum(query[0:q_gap] != '-')
        query_to_target[q_start + relative_q_gap: q_end + 1, hsp_num] += direction
    for mismatch in mismatches[0]:
        alignment_types[q_start + mismatch-np.sum(query_gaps[0] < mismatch), hsp_num] = 1

# Build a map to convert the coordinates of the gene to the new annotation
def add_to_map(hsp, query_to_target, hsp_num, alignment_types):
    if get_strand(hsp) == "+":
        direction = 1
    else:
        direction = -1
    query = np.array(list(hsp.query))
    subject = np.array(list(hsp.sbjct))
    q_start = hsp.query_start - 1
    q_end = hsp.query_end - 1
    s_start = hsp.sbjct_start -1
    query_gaps = np.where(query == '-')
    subject_gaps = np.where(subject == '-')
    mismatches = np.where((query != '-') & (subject != '-') & (query != subject))
    subject_range = np.arange(s_start, s_start +
                              direction * (q_end - q_start + 1), direction)
    query_to_target[q_start: q_end + 1, hsp_num]= subject_range
    alignment_types[q_start: q_end + 1, hsp_num] = 0
    adjust_coordinates_for_gaps(query_to_target, alignment_types, subject_gaps, query_gaps, mismatches, hsp_num,
                                q_start, q_end, direction, query)


def get_strand(hsp):
    if hsp.sbjct_start > hsp.sbjct_end:
        strand = '-'
    else:
        strand = '+'
    return strand


def has_CDS(gene):
    for tran in gene.Transcripts:
        for exon in tran.Exons:
            if len(exon.CDS) > 0:
                return True
    return False


def get_start_and_end(hsp):
    start = max(hsp.sbjct_start, hsp.sbjct_end) - 1
    end = min(hsp.sbjct_start, hsp.sbjct_end) - 1
    return start, end

# For each high scoring pair in the alignment, find other alignments that are within 1.5 * the gene length
# and on the same strand. Find the number of exons or CDS in the group of alignments, and select the group that
# has the most exons or CDSs. If there are different groups with the same number of exons aligned, break
# the tie by the number of exon bases such that groups with more complete exon/CDS alignments will be chosen
def find_gene_location(alignment, gene_length, gene_name, aligned_list, gff):
    max_num_exons = -1
    max_exon_bps = -1
    best_hsp_in_group = alignment.hsps[0]
    for hsp1 in alignment.hsps:
        hsp1_start, hsp1_end = get_start_and_end(hsp1)
        exons_in_group = {}
        max_exons_in_group = -1
        for hsp2 in alignment.hsps:
            hsp2_start, hsp2_end = get_start_and_end(hsp2)
            hsp_distance = max(hsp1_end, hsp2_end) - min(hsp1_start, hsp2_start)
            if hsp_distance <= 1.5 * gene_length and get_strand(hsp1) == get_strand(hsp2):
                if check_valid_mapping_coords(hsp2, aligned_list, gene_name, gff):
                    num_exons_in_hsp = check_for_exons_or_CDS(hsp2.query_start, hsp2.query_end, gff[gene_name],
                                                              has_CDS(gff[gene_name]), hsp2, exons_in_group)
                    if num_exons_in_hsp < max_exons_in_group:
                        max_exons_in_group = num_exons_in_hsp
                        best_hsp_in_group = hsp2
        num_exons_in_group = len(list(exons_in_group.keys()))
        exons_bps_in_group = sum(list(exons_in_group.values()))
        if num_exons_in_group > max_num_exons or (num_exons_in_group == max_num_exons
                                                  and exons_bps_in_group  > max_exon_bps):
            best_hsp = best_hsp_in_group
            max_num_exons = num_exons_in_group
            max_exon_bps = exons_bps_in_group
    gene_strand = get_strand(best_hsp)
    gene_start, gene_end = get_start_and_end(best_hsp)
    return gene_strand, gene_start, gene_end


# Check if a gene overlaps another gene in the original annotation
def find_containment(gene_name, overlapping_gene, gff):
    gene1_start = gff[gene_name].start
    gene2_start = gff[overlapping_gene].start
    gene1_end = gff[gene_name].end
    gene2_end = gff[overlapping_gene].end
    return min(gene1_end, gene2_end) - max(gene1_start, gene2_start) > 0


# Check if a gene is already mapped this locus. If there is, check the original annotation to determine if those
# genes should overlap
def check_valid_mapping_coords(hsp, aligned_list, gene_name, gff):
    hsp_start, hsp_end = get_start_and_end(hsp)
    hsp_strand = get_strand(hsp)
    for gene in aligned_list:
        if min(hsp_end, gene.end) - max(hsp_start, gene.start) > 0 and hsp_strand == gene.strand:
            overlapping_gene = gene.name
            is_contained = find_containment(gene_name, overlapping_gene, gff)
            return is_contained
    return True

# Check if there are exons in the alignment. Count the number of exon bases covered in the alignment. If the gene has
# at least 1 annotated CDS, count the CDS bases covered instead of exons
def check_for_exons_or_CDS(query_start, query_end, gene, has_CDS, hsp, exons_in_group):
    exon_id = 0
    num_exons_found = 0
    for tran in gene.Transcripts:
        for exon in tran.Exons:
            exon_id += 1
            if has_CDS:
                if len(exon.CDS) > 0:
                    if gene.strand == "+":
                        feature_start = exon.CDS[0].start - gene.start
                        feature_end = exon.CDS[0].end - gene.start
                    else:
                        feature_start = gene.end - exon.CDS[0].end
                        feature_end = gene.end - exon.CDS[0].start
                else:
                    feature_start = 0
                    feature_end = 0
            else:
                if gene.strand == "+":
                    feature_start = exon.start - gene.start
                    feature_end = exon.end-gene.start
                else:
                    feature_start = gene.end - exon.end
                    feature_end = gene.end - exon.start
            overlap = min(query_end - 1, feature_end) - max(query_start - 1, feature_start)
            if overlap > 0:
                num_exons_found += 1
                matches = hsp.match[
                          feature_start - (hsp.query_start - 1): feature_end + 1 - (hsp.query_start - 1)].count('|')
                if exon_id not in exons_in_group.keys():
                    exons_in_group[exon_id] = matches
                else:
                    exons_in_group[exon_id] = max(exons_in_group[exon_id], matches)
    return num_exons_found


def find_full_length_alignment(hsp,gene):
    strand = get_strand(hsp)
    if hsp.query_end - hsp.query_start == gene.end - gene.start and strand == gene.strand:
        return True
    else:
        return False

# Find the high scoring pairs that are in the correct location and have no genes already mapped here
def find_valid_alignments(alignment, gene_strand, gene_start, gene_end, gene_name, gff, aligned_list):
    valid_alignments =[]
    num_hsp = 0
    full_length_alignment_found = False
    length_threshold = 1.5
    for hsp in alignment.hsps:
        num_hsp += 1
        hsp_strand = get_strand(hsp)
        is_valid_mapping = check_valid_mapping_coords(hsp, aligned_list, gene_name, gff)
        is_full_length = find_full_length_alignment(hsp, gff[gene_name])
        hsp_start, hsp_end = get_start_and_end(hsp)
        if is_full_length and is_valid_mapping:
            full_length_alignment_found = True
        distance_between_hsps = max(hsp_end, gene_end) - min(hsp_start, gene_start)
        original_gene_length = gff[gene_name].end - gff[gene_name].start
        if distance_between_hsps <= length_threshold * original_gene_length:
            valid_length = True
            gene_end = max(hsp_end, gene_end)
            gene_start = min(hsp_start, gene_start)
        else:
            valid_length = False
        if (is_full_length or (full_length_alignment_found is False and hsp_strand == gene_strand and valid_length)) \
                and is_valid_mapping:
            valid_alignments.append(hsp)
    return valid_alignments, full_length_alignment_found

# After the general location has been selected, sort the alignments by the number of exon or CDS bps in the alignment
# if there are multiple full length alignments, choose the alignment that has the best overall score for all exons and
# CDSs
def sort_alignments(query_to_target, gene, alignment_types, full_length_alignment_found):
    exon_or_CDS_score = np.zeros(np.shape(query_to_target)[1])
    for i in range(np.shape(query_to_target)[1]):
        alignment_score = 0
        for tran in gene.Transcripts:
            for exon in tran.Exons:
                if has_CDS(gene):
                    if len(exon.CDS) >0:
                        feature = exon.CDS[0]
                    else:
                        break
                else:
                    feature = exon
                if gene.strand == '+':
                    relative_start = feature.start - gene.start
                    relative_end = feature.end - gene.start
                else:
                    relative_start = gene.end - feature.end
                    relative_end = gene.end - feature.start
                alignment_score += np.sum(alignment_types[relative_start:relative_end +1, i])
        exon_or_CDS_score[i] = alignment_score
    best_alignment_order = exon_or_CDS_score.argsort()
    sorted_query_to_target = query_to_target[:, best_alignment_order]
    sorted_alignment_types = alignment_types[:, best_alignment_order]
    if full_length_alignment_found:
        return sorted_query_to_target[:, [0]], sorted_alignment_types[:, [0]]
    else:
        return sorted_query_to_target, sorted_alignment_types


# Get the strand of the mapping location with the best score
def get_final_strand(best_alignment):
    top_alignment_coords = best_alignment[best_alignment != 0]
    if top_alignment_coords[0] < top_alignment_coords[-1]:
        return "+"
    else:
        return "-"


# Map the genes with the best scoring alignments first so a paralogue isn't incorrectly mapped to a location that has a
# better alignment
def find_best_score(record):
    best_score = 0
    for alignment in record.alignments:
        score = alignment.hsps[0].bits
        if score > best_score:
            best_score = score
    return best_score

# Iterate through every gene and every high scoring pair to find the location of the gene on the new assembly
def create_coordinate_conversion(records, aligned_list, gff, attempt, failed_file):
    redo_list = []
    for record in sorted(records, key=lambda blast_record: find_best_score(blast_record), reverse=True):  # type: Record.Blast
        gene_name = record.query[:-3]
        if len(record.alignments) == 0:
            failed_file.write(gene_name + "\n")
        else:
            for alignment in record.alignments: # type: Record.Alignment
                gene_strand, gene_start, gene_end= find_gene_location(alignment, record.query_length, gene_name,
                                                                      aligned_list, gff)
                hsp_num = 0
                valid_alignments, full_length_alignment_found = find_valid_alignments(alignment, gene_strand, gene_start,
                                                                                      gene_end, gene_name,
                                                                                      gff, aligned_list)
                query_to_target = np.zeros(shape=(record.query_length, len(valid_alignments)))
                alignment_types = np.zeros(shape=(record.query_length, len(valid_alignments)))+2
                aligned_chr = alignment.hit_def
                for hsp in valid_alignments:
                    add_to_map(hsp,  query_to_target, hsp_num, alignment_types)
                    hsp_num += 1
                sorted_query_to_target, sorted_alignment_types = sort_alignments(query_to_target, gff[gene_name],
                                                                                 alignment_types,
                                                                                 full_length_alignment_found)
                if hsp_num > 0:
                    gene_strand = get_final_strand(sorted_query_to_target[:,0])
                new_gene = liftover_genes(sorted_query_to_target, gene_strand, aligned_chr, gene_name, gff,
                                          sorted_alignment_types, attempt)
                add_gene_and_tran_coords(new_gene, gff[gene_name])
                if new_gene.status == "incomplete" and attempt < 2:
                    redo_list.append(record)
                else:
                    aligned_list.append(new_gene)
    return redo_list, aligned_list


# Try to find a mapping for the gene three times. Mappings are only selected on the first try if a gene maps completely
# to a location on the same strand as the original annotation. On the second attempt, mappings are allowed if they
# are complete but on the opposite strand. On the third try, the best partial mapping is selected
def attempt_liftover(records, aligned_list, gff, failed_file):
    num_liftover_attempts = 3
    for i in range(num_liftover_attempts):
        redo_list, aligned_list = create_coordinate_conversion(records, aligned_list, gff, i, failed_file)
        records = redo_list
    return aligned_list


def process_blast(blast_output,  gff, output1, output2):
    result_handle = open(blast_output)
    blast_records = NCBIXML.parse(result_handle)
    CDS_records = []
    noncoding_records = []
    aligned_list = []
    f_out = open(output1, 'w')
    failed_file = open(output2, 'w')
    # separate genes with CDS and genes without and prioritize mapping CDS genes
    for record in blast_records:
        if has_CDS(gff[record.query[:-3]]):
            CDS_records.append(record)
        else:
            noncoding_records.append(record)
    aligned_list = attempt_liftover(CDS_records, aligned_list,  gff, failed_file)
    aligned_list = attempt_liftover(noncoding_records, aligned_list, gff, failed_file)
    print_gff(aligned_list, f_out, failed_file)
    f_out.close()
    failed_file.close()



