def process_alignments(blast_record, gene_db, gene):
    exon_alignments = []
    if len(blast_record.alignments) != 0:
        exon_alignments = find_exon_containing_alignments(blast_record.alignments, gene_db, gene)
    return exon_alignments



def find_exon_containing_alignments(alignments, gene_db, gene):
    exon_alignments = []
    unique_exons = find_unique_exons(gene_db, gene)
    for alignment in alignments:
        target_chrm = alignment.hit_def.split()[0]
        for hsp in alignment.hsps:
            if contains_exon(unique_exons,  gene, hsp.query_start, hsp.query_end):
                exon_alignments.append((hsp, target_chrm))
    return exon_alignments

def get_exons_and_cds(gene_db, gene):
    all_exons = gene_db.children(gene, featuretype="exon")
    all_cds = gene_db.children(gene, featuretype="CDS")
    return all_exons, all_cds


def find_unique_exons(gene_db, gene):
    all_exons, all_cds = get_exons_and_cds(gene_db, gene)
    unique_exons_dict = {}
    for exon in all_exons:
        exon_key = str(exon.start) + ":" + str(exon.end)
        if exon_key not in unique_exons_dict:
            unique_exons_dict[exon_key] = exon
        elif exon_has_cds(find_cds(unique_exons_dict[exon_key], all_cds)) is False and exon_has_cds(find_cds(exon, all_cds)) is True:
            unique_exons_dict[exon_key] = exon
    unique_exons = list(unique_exons_dict.values())
    unique_exons.sort(key=lambda x: (int(x.start), int(x.end)))
    return unique_exons



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


def find_cds(exon, all_cds):
    cds_list = []
    for cds in all_cds:
        if cds.start >= exon.start and cds.end <= exon.end:
            cds_list.append(cds)
    return cds_list



def exon_has_cds(cds_list):
    return len(cds_list) > 0