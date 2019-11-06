import numpy as np

def find_CDS_alignments(exon, alignments, filtered_types):
    alignments_to_score = alignments
    alignment_types = filtered_types
    if len(exon.CDS) != 0:
        start_offset = exon.CDS[0].start - exon.start
        end_offset = exon.end - exon.CDS[-1].end
        if exon.strand == "+":
            CDS_alignments = alignments[start_offset: len(alignments) - end_offset,:]
            CDS_alignment_types = filtered_types[start_offset: len(alignments) - end_offset,:]
        else:
            CDS_alignments = alignments[end_offset: len(alignments) - start_offset,:]
            CDS_alignment_types = filtered_types[end_offset: len(alignments) - start_offset,:]
        alignments_to_score = CDS_alignments
        alignment_types = CDS_alignment_types
    return alignments_to_score, alignment_types


# Find the best alignment for the exon. The best alignment is determined by fewest number of gaps and mismatches. Gaps
# are penalized more than mismatches since they generally introduce a frameshift into the protein sequence
def find_best_alignment(exon, alignments, filtered_types):
    best_alignment = np.array([0])
    min_score = 10000000
    num_hsps = np.shape(alignments)[1]
    alignments_to_score, alignment_types = find_CDS_alignments(exon, alignments, filtered_types)
    for hsp in range(num_hsps):
        if np.sum(alignments_to_score[:,hsp]) != 0:
            score= np.sum(alignment_types[:,hsp] != 0)
            if score < min_score:
                min_score = score
                best_alignment = alignments[:,hsp]
    return best_alignment








