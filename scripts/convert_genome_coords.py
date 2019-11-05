import gff_entry
import numpy as np
import filter

# Use the coordinate conversion map to return all the valid mappings of the interval
def convert_coords(coord_map,alignment_types, gene_start,gene_end, original_start, original_end, gene_strand):
    relative_start = original_start - gene_start
    relative_end = original_end - gene_start
    if gene_strand == '+':
        return coord_map[relative_start: relative_end + 1], alignment_types[relative_start: relative_end + 1]
    else:
        return coord_map[(gene_end-gene_start) - relative_end: (gene_end-gene_start) - relative_start + 1], \
               alignment_types[(gene_end-gene_start) - relative_end: (gene_end-gene_start) - relative_start + 1]

# If this is the first attempt to map a gene, only allow mappings where exons align in the correct order
def filter_alignments(alignments, alignment_types, new_transcript, attempt):
    for i in range(np.shape(alignments)[1]):
        alignment = alignments[:, i]
        valid_alignment = alignment[alignment != 0]
        for exon in new_transcript.Exons:
            if len(valid_alignment) != 0:
                if min(exon.end, np.max(valid_alignment)) - max(exon.start, np.min(valid_alignment)) > 0 or\
                        (np.max(valid_alignment) < exon.start and attempt == 0):
                    alignments[:, i] = 0
                    alignment_types[:, i] = 2
    filtered_alignments = alignments
    filtered_types = alignment_types
    num_complete_alignments = 0
    return filtered_alignments, num_complete_alignments, filtered_types


def adjust_frame(CDS_alignment, lifted_start, original_frame):
    missing_codon_offset = np.where(CDS_alignment == lifted_start)[0][0] % 3
    new_frame = (int(original_frame) + 3 - missing_codon_offset) % 3
    return int(new_frame)

# Find the new start and end coordinate of the CDS based on the best alignment and label it as complete, unaligned, or
# trucated
def convert_CDS_coords(CDS, exon, best_alignment):
    start_offset = CDS.start - exon.start
    end_offset = exon.end - CDS.end
    if exon.strand == "+":
        CDS_alignment = best_alignment[start_offset: len(best_alignment)-end_offset]
    else:
        CDS_alignment = best_alignment[end_offset: len(best_alignment)-start_offset]
    valid_alignments = CDS_alignment[CDS_alignment != 0]
    if len(valid_alignments) == 0:
        lifted_start = -1
        lifted_end = -1
        status = "unaligned"
    else:
        lifted_start = valid_alignments[0]
        lifted_end = valid_alignments[-1]
        if lifted_start == CDS_alignment[0] and lifted_end == CDS_alignment[-1]:
            status = 'complete'
        else:
            status = 'truncated'
            new_frame = adjust_frame(CDS_alignment, lifted_start, CDS.frame)
            CDS.frame = str(new_frame)
    return lifted_start, lifted_end, status


# Find the new start and end coordinate of the exon based on the best alignment and label it as complete, unaligned, or
# trucated
def convert_exon_coords(best_alignment):
    valid_alignments = best_alignment[best_alignment != 0]
    if len(valid_alignments) == 0:
        lifted_start = -1
        lifted_end = -1
        status = "unaligned"
    else:
        lifted_start = valid_alignments[0]
        lifted_end = valid_alignments[-1]
        if lifted_start == best_alignment[0] and lifted_end == best_alignment[-1]:
            status = 'complete'
        else:
            status = 'truncated'
    return lifted_start, lifted_end, status


# Convert the coordinates of each exon in the gene separately. If there is a complete alignment for the gene, use
# that alignment for all exons. If there are multiple complete alignments, chose the alignment that has the fewest
# gaps and mismatches across all exons or CDSs. Otherwise, select the best mapping for each exon individually.
def liftover_genes(query_to_target, strand, chrm , gene_name, gff, alignment_types, attempt):
    gene = gff[gene_name]
    new_gene = gff_entry.Gene(chrm, gene.name, 0, 0, [], strand, "unknown", gene.source)  # type: Gene
    for tran in gene.Transcripts:  # type: Transcript
        new_transcript = gff_entry.Transcript(tran.id, 0, 0, [], strand, "unknown", tran.source) # :type Transcript
        new_gene.Transcripts.append(new_transcript)
        for i in range (len(tran.Exons)):  # type: Exon
            exon = tran.Exons[i]
            exon_alignments, exon_alignment_types = convert_coords(query_to_target, alignment_types, gene.start,
                                                                   gene.end, exon.start, exon.end,
                                                                   gene.strand)
            filtered_alignments, num_complete_alignments, filtered_types = filter_alignments(exon_alignments.copy(),
                                                                                             exon_alignment_types,
                                                                                             new_transcript, attempt)
            if num_complete_alignments != 1:
                best_alignment = filter.find_best_alignment(exon, filtered_alignments, filtered_types)
            else:
                best_alignment = filtered_alignments
            lifted_start_exon, lifted_end_exon, status_exon = convert_exon_coords(best_alignment)
            new_exon = gff_entry.Exon(exon.parent, exon.id, int(min(lifted_start_exon, lifted_end_exon)),
                                      int(max(lifted_start_exon, lifted_end_exon)), status_exon, strand, exon.source, [])
            if len(exon.CDS) != 0:
                for CDS in exon.CDS:
                    lifted_start_CDS, lifted_end_CDS, status_CDS = convert_CDS_coords(CDS, exon, best_alignment)
                    new_CDS = gff_entry.CDS(int(min(lifted_start_CDS, lifted_end_CDS)),
                                            int(max(lifted_start_CDS, lifted_end_CDS)), CDS.frame, status_CDS)
                    new_exon.CDS.append(new_CDS)
            new_transcript.Exons.append(new_exon)
        new_transcript.strand = strand
    new_gene.strand = strand
    return new_gene


