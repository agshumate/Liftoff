import gffutils
from process_blast_alignments import  find_cds, exon_has_cds, find_unique_exons
from find_best_mapping import find_boundary_coords


def build_lifted_gene_model(mapped_exons, shortest_path_weight, gene_db, gene, unmapped_genes, threshold, gene_name):
    original_transcripts = list(gene_db.children(gene, level=1))
    transcript_list = []
    feature_list = []
    for transcript in original_transcripts:
        exon_and_cds_list = []
        add_exons_to_list(transcript, gene_db, mapped_exons, feature_list, exon_and_cds_list)
        if len(exon_and_cds_list) > 0:
            build_new_transcript(exon_and_cds_list, feature_list, transcript, gene, gene_db, transcript_list)
    if len(transcript_list) > 0:
        gene_list = build_new_gene(transcript_list, gene, gene_db, shortest_path_weight, threshold, gene_name)
    else:
        gene_list = []
    if len(gene_list) !=0:

        feature_list.append(gene_list[0])
    else:
        unmapped_genes.append(gene)
    return feature_list


def add_exons_to_list(transcript, gene_db, mapped_exons, feature_list, exon_and_cds_list):
    all_exons = list(gene_db.children(transcript, featuretype="exon"))
    all_cds = list(gene_db.children(transcript, featuretype="CDS"))
    all_exons.sort(key=lambda x: x.start)
    exon_num = 0
    for exon in all_exons:
        exon_num += 1
        exon_feature = build_new_exon(exon, mapped_exons, exon_num, transcript, exon_and_cds_list, feature_list)
        if exon_feature is not None:
            cds_list = find_cds(exon, all_cds)
            for cds in cds_list:
                cds_feature = build_new_cds(cds, exon, mapped_exons, transcript, feature_list, exon_and_cds_list)
            if exon_has_cds(cds_list) and cds_feature is not None:
                    exon_feature.attributes["status"] = cds_feature.attributes["status"]


def get_exon_key(exon):
    return str(exon.start) + ":" + str(exon.end)


def build_new_exon(exon, mapped_exons, exon_num, transcript,  exon_list, feature_list):
    exon_key = get_exon_key(exon)
    if exon_key in mapped_exons:
        lifted_exon_node = mapped_exons[exon_key]
        exon_feature = gffutils.Feature(seqid=lifted_exon_node.chrm, source="Liftoff", featuretype='exon',
                               start=lifted_exon_node.lifted_start + 1, end=lifted_exon_node.lifted_end + 1,
                               strand=lifted_exon_node.strand)
        exon_feature.attributes["Parent"] = transcript.id

        exon_feature.attributes["exon_id"] = [str(exon_num)]
        exon_feature.attributes["status"] = get_exon_status(lifted_exon_node.alignment)
        feature_list.append(exon_feature)
        exon_list.append(exon_feature)
        return exon_feature
    return None


def build_new_cds(cds, exon, mapped_exons, transcript, feature_list, exon_list):
    exon_key = get_exon_key(exon)
    lifted_exon_node = mapped_exons[exon_key]
    cds_alignment = find_cds_alignment(lifted_exon_node.alignment, exon, cds)
    cds_start, cds_end = find_boundary_coords(cds_alignment)
    if cds_start !=0 and cds_end !=0 :
        cds_status = get_exon_status(cds_alignment)
        cds_feature = gffutils.Feature(seqid=lifted_exon_node.chrm, source="Liftoff", featuretype='CDS',
                         start=cds_start +1, end=cds_end +1,
                         strand=lifted_exon_node.strand)
        cds_feature.attributes["Parent"] = transcript.id
        cds_feature.attributes["status"] = cds_status
        feature_list.append(cds_feature)
        exon_list.append(cds_feature)
        return cds_feature
    return None


def find_cds_alignment(alignment, exon, cds):
    start_offset = cds.start - exon.start
    end_offset = exon.end - cds.end
    if exon.strand == "+":
        cds_alignment = alignment[start_offset: len(alignment)-end_offset]
    else:
        cds_alignment = alignment[end_offset: len(alignment)-start_offset]
    return cds_alignment


def get_exon_status(alignment):
    if alignment[0] == 0 or alignment[-1] == 0:
        return ['truncated']
    else:
        return ['complete']


def get_transcript_coverage(original_transcript, transcript_feature, feature_list, gene_db):
    original_length = 0
    lifted_length = 0
    original_exons = gene_db.children(original_transcript, featuretype="exon")
    lifted_exons = [feature for feature in feature_list if feature.featuretype=="exon" and feature.attributes["Parent"][0]==transcript_feature.id]
    for exon in original_exons:
        original_length += (exon.end - exon.start +1)
    for exon in lifted_exons:
        lifted_length += (exon.end - exon.start +1)
    return round(lifted_length/original_length,2)


def build_new_transcript(exon_list, feature_list, transcript, gene, gene_db, transcript_list):
    starts = [exon.start for exon in exon_list]
    ends = [exon.end for exon in exon_list]
    strand = exon_list[0].strand
    chrm = exon_list[0].seqid
    attributes=transcript.attributes
    attributes["Parent"] = gene.id
    attributes["status"] = get_transcript_status(exon_list, transcript, gene_db)
    transcript_feature = gffutils.Feature(id=transcript.id, source="Liftoff", featuretype=transcript.featuretype, start = min(starts), end=max(ends), seqid=chrm, strand=strand, attributes=attributes)
    coverage = get_transcript_coverage(transcript, transcript_feature, feature_list, gene_db)
    transcript_feature.attributes["coverage"] = str(coverage)
    feature_list.append(transcript_feature)
    transcript_list.append(transcript_feature)
    return transcript_list


def get_transcript_status(exon_list, transcript, gene_db):
    if len(list(gene_db.children(transcript.id, featuretype="exon"))) + len(list(gene_db.children(transcript.id, featuretype="CDS")))> len(exon_list):
        return "incomplete"
    for exon in exon_list:
        if "complete" not in exon.attributes["status"]:
            return "incomplete"
    return "complete"


def build_new_gene(transcript_list, gene, gene_db, path_weight, threshold, gene_name):
    starts = [transcript.start for transcript in transcript_list]
    ends = [transcript.end for transcript in transcript_list]
    strand = transcript_list[0].strand
    chrm = transcript_list[0].seqid
    attributes = gene.attributes
    gene_score = calculate_gene_score(path_weight, gene, gene_db)
    gene_feature = gffutils.Feature(id=gene_name, featuretype="gene", source="Liftoff", start=min(starts), end=max(ends), seqid=chrm, strand=strand, attributes=attributes, score=gene_score)
    attributes["coverage"] = str(max([float(transcript.attributes["coverage"][0]) for transcript in transcript_list]))
    attributes["status"] = get_gene_status(transcript_list, gene_db, gene)
    gene_list = []
    if  float(attributes["coverage"][0]) >= threshold:
        gene_list.append(gene_feature)
    return gene_list


def calculate_gene_score(shortest_path_length, gene, gene_db):
    unique_exons = find_unique_exons(gene_db, gene)
    exon_lengths = [exon.end - exon.start for exon in unique_exons]
    total_length = sum(exon_lengths)
    return shortest_path_length/total_length


def get_gene_status(transcript_list, gene_db, gene):
    transcript_statuses = [transcript.attributes["status"][0] for transcript in transcript_list]
    if "complete" not in transcript_statuses:
        return "none"
    if "complete" in transcript_statuses and "incomplete" in transcript_statuses:
        return "partial"
    if "complete" in transcript_statuses and len(list(gene_db.children(gene, level=1))) != len(transcript_list):
        return "partial"
    if "incomplete" not in transcript_statuses:
        return "all"
