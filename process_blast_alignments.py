def process_alignments(blast_record, gene_db, parent):
    exon_alignments = []
    if len(blast_record.alignments) != 0:
        exon_alignments = find_children_containing_alignments(blast_record.alignments, gene_db, parent)
    return exon_alignments



def find_children_containing_alignments(alignments, gene_db, parent):
    exon_alignments = []
    unique_exons = find_children(gene_db, parent)
    for alignment in alignments:
        target_chrm = alignment.hit_def.split()[0]
        for hsp in alignment.hsps:
            if contains_exon(unique_exons, parent, hsp.query_start, hsp.query_end):
                exon_alignments.append((hsp, target_chrm))
    return exon_alignments

def get_exons_and_cds(gene_db, gene):
    all_exons = gene_db.children(gene, featuretype="exon")
    all_cds = gene_db.children(gene, featuretype="CDS")
    return all_exons, all_cds


def find_children(gene_db, parent):
    lowest_children=[]
    all_children = gene_db.children(parent)
    for child in all_children:
        if len(list(gene_db.children(child))) == 0:
            lowest_children.append(child)

    if len(lowest_children) == 0:
        lowest_children.append(parent)
    lowest_children.sort(key=lambda x: (int(x.start), int(x.end)))
    return lowest_children



def contains_exon(unique_exons, gene, hsp_start, hsp_end):
    for exon in unique_exons:
        exon_start_offset, exon_end_offset = find_offsets(exon, gene)
        if count_overlap(hsp_start-1, hsp_end-1, exon_start_offset, exon_end_offset ) > 0:
            return True
    return False


def count_overlap(start1, end1, start2, end2):
    overlap = min(end1, end2) - max(start1, start2)
    return overlap



def find_offsets(exon, gene):
    if gene.strand == "+":
        start_offset = exon.start - gene.start
        end_offset = exon.end - gene.start
    elif gene.strand == "-":
        start_offset = gene.end - exon.end
        end_offset = gene.end - exon.start
    return start_offset, end_offset






