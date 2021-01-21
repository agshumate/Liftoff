import argparse
from pyfaidx import Fasta
from liftoff import extract_features, liftoff_utils, new_feature, run_liftoff, write_new_gff
import parasail
from itertools import groupby
import numpy as np
import cProfile, pstats, io
from Bio.Seq import Seq
import gffutils
import ujson as json




def main(arglist=None):
    pr = cProfile.Profile()
    pr.enable()
    args = parse_args(arglist)
    ref_db = extract_features.build_database(args.db1, args.g1, True, True)
    target_db = extract_features.build_database(args.db2, args.g2, True, True)
    ref_fa = Fasta(args.reference_genome)
    target_fa = Fasta(args.target_genome)
    if args.f is not None:
        remap_genes = [gene.rstrip() for gene in open(args.f, 'r').readlines()]
    else:
        remap_genes = []
    reference_parent_dict, target_child_dict, target_parent_dict = find_and_polish_genes(target_db, ref_db, target_fa,
                                                                                   ref_fa, remap_genes)
    liftoff_args = ([args.target_genome, args.reference_genome, '-db', 'polished_features_db', '-subcommand',
                     "polish", "-o", args.o, "-d", '1000000000' ])
    args = run_liftoff.parse_args(liftoff_args)
    polished_annotations = run_liftoff.run_all_liftoff_steps(args)
    print("updating target features")
    update_target_database(polished_annotations, target_db)
    complete_feature_list = update_feature_list(target_db)
    write_new_gff.write_new_gff(complete_feature_list, reference_parent_dict, args, ref_db)
    pr.disable()  # end profiling
    s = io.StringIO()
    sortby = 'cumulative'
    ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
    ps.print_stats()
    print (s.getvalue())


def update_target_database(polished_annotations, target_db):
    to_delete, to_update = [], []
    for ann in polished_annotations:
        gene = polished_annotations[ann][0]
        if gene.attributes["sequence_ID"][0] >= target_db[gene.id].attributes["sequence_ID"][0] and gene.attributes[
            "coverage"][0] >= target_db[gene.id].attributes["coverage"][0]:
            for feature in polished_annotations[ann]:
                to_delete.append(feature.id)
                to_update.append(gffutils.Feature(id=feature.id, source=feature.source,
                                                  featuretype=feature.featuretype,
                                                  seqid=feature.seqid, \
                                                  strand=feature.strand,
                                                  attributes=feature.attributes, start=feature.start,
                                                  end=feature.end))
    target_db.delete(to_delete)
    target_db.update(to_update, merge_strategy="create_unique")



def update_feature_list(db):
    new_feature_list = {}
    c = db.conn.cursor()
    query = "select * from relations join features as a on a.id  = relations.child join features as b on b.id = " \
            "relations.parent "
    results = c.execute(query)
    for result in results:
        result_tup = tuple(result)
        parent_name = result_tup[0]
        child_name = result_tup[1]
        child_feature = result_tup[4:13]
        parent_feature = result_tup[16:25]
        parent_attributes = json.loads(parent_feature[-1])
        child_attributes = json.loads(child_feature[-1])
        if "Parent" not in parent_attributes:
            if parent_attributes["copy_num_ID"][0] not in new_feature_list:
                new_parent = new_feature.new_feature(parent_name, parent_feature[2], parent_feature[0],
                                                     parent_feature[1], parent_feature[6], parent_feature[3],
                                                     parent_feature[4], parent_attributes)
                new_feature_list[new_parent.attributes["copy_num_ID"][0]] = [new_parent]
            new_child = new_feature.new_feature(child_name, child_feature[2], child_feature[0],
                                                     child_feature[1], child_feature[6], child_feature[3],
                                                     child_feature[4], child_attributes)
            new_feature_list[new_parent.attributes["copy_num_ID"][0]].append(new_child)
    return new_feature_list






def parse_args(arglist=None):
    parser = argparse.ArgumentParser(description='Polish lifted over annotations')
    parser.add_argument('reference_genome', help='reference fasta')
    parser.add_argument('target_genome', help='target fasta')
    req_grp1 = parser.add_argument_group('Required input (reference annotation)')
    mxgrp1 = req_grp1.add_mutually_exclusive_group(required=True)

    mxgrp1.add_argument(
        '-g1', metavar='GFF', help='reference annotation file to lift over in GFF or GTF format')
    mxgrp1.add_argument(
        '-db1', metavar='DB', help='name of reference feature database; if not specified, the -ref_g '
                                      'argument must be provided and a database will be built automatically' )

    req_grp2 = parser.add_argument_group('Required input (target annotation)')
    mxgrp2 = req_grp2.add_mutually_exclusive_group(required=True)
    mxgrp2.add_argument(
        '-g2', metavar='GFF', help=' target annotation file to lift over in GFF or GTF format')
    mxgrp2.add_argument(
        '-db2', metavar='DB', help='name of target feature database; if not specified, the -ref_g '
                                      'argument must be provided and a database will be built automatically')
    parser._positionals.title = 'Required input (sequences)'
    parser.add_argument("-f", metavar='TXT', help="text file of gene names to remap ")
    parser.add_argument('-o', default='stdout', metavar='FILE',help='write output to FILE in same format as input; by default, output is written to terminal (stdout)'
    )
    args = parser.parse_args(arglist)
    return args


def find_and_polish_genes(target_db, ref_db, target_fa, ref_fa, genes_to_remap):
    output_sam = open("intermediate_files/polish.sam", 'w')
    write_sam_header(target_fa, output_sam)
    num_genes = 0
    if len(genes_to_remap) == 0:
        f=open("intermediate_files/polished.ids", 'w')
        genes_to_remap = [gene.id for gene in target_db.features_of_type(featuretype="gene") if gene.attributes[
            "sequence_ID"][0] != '1.0']
        [f.write(gene + "\n") for gene in genes_to_remap]
    sub_ref_db_list = []
    ref_child_dict = build_child_dict(ref_db, genes_to_remap, sub_ref_db_list)
    target_child_dict = build_child_dict(target_db, genes_to_remap,[])
    ref_parent_dict = build_parent_dict(ref_db, genes_to_remap, sub_ref_db_list)
    target_parent_dict = build_parent_dict(target_db, genes_to_remap,[])
    for gene in genes_to_remap:
        num_genes +=1
        print(num_genes)
        polish_annotation(ref_fa, target_fa, ref_db, target_db, gene, output_sam, ref_child_dict,
                              target_child_dict, ref_parent_dict, target_parent_dict)
    gffutils.create_db(sub_ref_db_list, "polished_features_db", force=True,merge_strategy="create_unique")
    return ref_parent_dict, target_child_dict, target_parent_dict





def build_child_dict(db, genes_to_remap, db_list):
    c = db.conn.cursor()
    cond = ', '.join('"{0}"'.format(w) for w in genes_to_remap)
    query = "select * from relations join features on features.id  = relations.child where relations.parent IN ({" \
            "})".format(
        cond)
    results = c.execute(query)
    children_dict = {}
    for result in results:
        feature_tup = tuple(result)
        parent = feature_tup[0]
        if parent not in children_dict:
            children_dict[parent] = []
        child_attributes = json.loads(feature_tup[12])
        child_attributes["ID"] = [feature_tup[3]]
        child = new_feature.new_feature(feature_tup[3], feature_tup[6], feature_tup[4], feature_tup[5],
                                            feature_tup[10],
                                            feature_tup[7], feature_tup[8], child_attributes)

        child_feature = gffutils.Feature(id=child.id, source = child.source, featuretype = child.featuretype, \
                                                                                  seqid=child.seqid, \
                                                                                        strand=child.strand,
                                        attributes = child_attributes, start=child.start, end=child.end )
        db_list.append(child_feature)
        children_dict[parent].append(child)
    return children_dict


def build_parent_dict(db, genes_to_remap, db_list):
    c = db.conn.cursor()
    cond = ', '.join('"{0}"'.format(w) for w in genes_to_remap)
    query = "select * from features where features.id IN ({})".format(cond)
    parent_dict = {}
    results = c.execute(query)
    for result in results :
        feature_tup = tuple(result)
        parent = new_feature.new_feature(feature_tup[0], feature_tup[3], feature_tup[1], feature_tup[2],
                                            feature_tup[7],
                                            feature_tup[4], feature_tup[5], json.loads(feature_tup[9]))
        parent_dict[feature_tup[0]] = parent
        db_list.append(gffutils.Feature(id=parent.id, source= parent.source, featuretype= parent.featuretype, \
                                                                                seqid=parent.seqid, \
                                                                                      strand=parent.strand,
                                        attributes = feature_tup[9], start=parent.start, end=parent.end ))
    return parent_dict




def write_sam_header(target_fa, output_sam):
    for seq in target_fa.keys():
        output_sam.write("@SQ" + "\t" + "SN:" + seq + "\t" + "LN:" + str(len(target_fa[seq])) + "\n")

def exons_ok(target_exons, ref_db):
    for target_exon in target_exons:
        ref_exon = ref_db[target_exon.id]
        if ref_exon.end - ref_exon.start != target_exon.end - target_exon.start:
            return False
    return True


def polish_annotation(ref_fa, target_fa, ref_db, target_db, gene_name, output_sam, ref_child_dict, target_child_dict,
                      ref_parent_dict, target_parent_dict):
    ref_gene, target_gene = ref_parent_dict[gene_name], target_parent_dict[gene_name]
    ref_exons = [feature for feature in ref_child_dict[ref_gene.id] if feature.featuretype == "exon"]
    target_exons = [feature for feature in target_child_dict[target_gene.id] if feature.featuretype == "exon"]
    ref_CDS = [feature for feature in ref_child_dict[ref_gene.id] if feature.featuretype == "CDS"]
    ref_CDS_intervals = liftoff_utils.merge_children_intervals(ref_CDS)
    splice_sites = add_splice_sites(ref_exons, ref_gene)
    merged_ref_intervals =liftoff_utils.merge_children_intervals(ref_exons)
    exon_group_dict = find_overlapping_exon_groups(merged_ref_intervals, ref_exons)
    matrix = make_scoring_matrix(3)
    for i in range (len(merged_ref_intervals)):
        ref_interval =  merged_ref_intervals[i]
        target_interval = get_target_interval(exon_group_dict[i], target_exons, target_db)
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






def add_splice_sites(exons, parent):
    splice_sites = []
    for exon in exons:
        if exon.start -2 >= parent.start:
            exon.start = exon.start -2
            splice_sites.append([exon.start, exon.start+1])
        if exon.end +2 <= parent.end:
            exon.end = exon.end +2
            splice_sites.append([exon.end-1, exon.end])
    return splice_sites


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
    matrix = parasail.matrix_create("ACGTacgt", match_reward_non_coding, mismatch_penalty, True)
    matrix[0,4] = match_reward_coding
    matrix[1,5] = match_reward_coding
    matrix[2,6] = match_reward_coding
    matrix[3,7] = match_reward_coding
    matrix[4,0] = match_reward_coding
    matrix[5,1] = match_reward_coding
    matrix[6,2] = match_reward_coding
    matrix[7,3] = match_reward_coding
    return matrix



def get_target_interval(exon_group, target_exons, target_db):
    target_group = []
    for ref_exon in exon_group:
        matching_target_exon = [exon for exon in target_exons if exon.id == ref_exon]
        if len(matching_target_exon) !=0:
            target_group.append(matching_target_exon[0])
    if len(target_group) != 0:
        start = min([exon.start for exon in target_group])
        end = max([exon.end for exon in target_group])
    else:
        start = 0
        end = 0
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



def write_sam_file(ref_interval, target_interval, ref_gene, target_gene, ref_seq, alignment, output_sam, flank ):
    reference_traceback = alignment.traceback.query
    target_traceback = alignment.traceback.ref
    print(reference_traceback, target_traceback)
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
    sam_fields = np.array([target_gene.id, bit_flag, target_gene.seqid, target_start_with_flank + int(target_start_offset),
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


if __name__ == "__main__":
    pr = cProfile.Profile()
    pr.enable()  # start profiling

    main()

    pr.disable()  # end profiling
    s = io.StringIO()
    sortby = 'cumulative'
    ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
    ps.print_stats()
    print (s.getvalue())
