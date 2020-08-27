from liftoff import write_new_gff, liftover_types
import argparse


def main(arglist=None):
    args = parse_args(arglist)
    if args.chroms is not None:
        ref_chroms, target_chroms = parse_chrm_files(args.chroms)
    else:
        ref_chroms = [args.reference]
        target_chroms = [args.target]
    parent_features_to_lift = get_parent_features_to_lift(args.f)
    lifted_feature_list = {}
    unmapped_features = []
    feature_db, feature_hierarchy, ref_parent_order = liftover_types.lift_original_annotation(ref_chroms, target_chroms,
                                                                                              lifted_feature_list, args,
                                                                                              unmapped_features,
                                                                                              parent_features_to_lift)
    unmapped_features = map_unmapped_features(unmapped_features, target_chroms, lifted_feature_list, feature_db,
                                              feature_hierarchy, ref_parent_order, args)
    map_features_from_unplaced_seq(unmapped_features, lifted_feature_list, feature_db, feature_hierarchy,
                                   ref_parent_order, args)
    write_unmapped_features_file(args.u, unmapped_features)
    map_extra_copies(args, lifted_feature_list, feature_hierarchy, feature_db, ref_parent_order)
    write_new_gff.write_new_gff(lifted_feature_list, feature_hierarchy.parents, args)


def parse_args(arglist=None):
    parser = argparse.ArgumentParser(description='Lift features from one genome assembly to another')
    parser.add_argument('target', help='target fasta genome to lift genes to')
    parser.add_argument('reference', help='reference fasta genome to lift genes from')

    refrgrp = parser.add_argument_group('Required input (annotation)')
    mxgrp = refrgrp.add_mutually_exclusive_group(required=True)
    mxgrp.add_argument(
        '-g', metavar='GFF', help='annotation file to lift over in GFF or GTF format'
    )
    mxgrp.add_argument(
        '-db', metavar='DB', help='name of feature database; if not specified, the -g '
        'argument must be provided and a database will be built automatically'
    )

    outgrp = parser.add_argument_group('Output')
    outgrp.add_argument(
        '-o', default='stdout', metavar='FILE',
        help='write output to FILE in GFF3 format; by default, output is written to terminal (stdout)'
    )
    outgrp.add_argument(
        '-u', default='unmapped_features.txt', metavar='FILE',
        help='write unmapped features to FILE; default is "unmapped_features.txt"',
    )
    outgrp.add_argument(
        '-exclude_partial', action='store_true',
        help='write partial mappings below -s and -a threshold to unmapped_features.txt; if true '
        'partial/low sequence identity mappings will be included in the gff file with '
        'partial_mapping=True, low_identity=True in comments'
    )
    outgrp.add_argument(
        '-dir', default='intermediate_files', metavar='DIR',
        help='name of directory to save intermediate fasta and SAM files; default is "intermediate_files"',
    )

    aligngrp = parser.add_argument_group('Alignment filtering')
    aligngrp.add_argument(
        '-a', default=0.5, metavar='A', type=float,
        help='designate a feature mapped only if it aligns with coverage ≥A; by default A=0.5',
    )
    aligngrp.add_argument(
        '-s', default=0.5, metavar='S', type=float,
        help='designate a feature mapped only if its child features (usually exons/CDS) align '
        'with sequence identity ≥S; by default S=0.5'
    )
    aligngrp.add_argument(
        '-n', default=50, metavar='N', type=int,
        help='consider at most N Minimap2 alignments for each feature; by default N=50'
    )
    aligngrp.add_argument(
        '-d', metavar='D', default=2.0, type=float,
        help='distance scaling factor; alignment nodes separated by more than a factor of D in '
        'the target genome will not be connected in the graph; by default D=2.0'
    )
    aligngrp.add_argument(
        '-flank', default=0, metavar='F', type=float, help="amount of flanking sequence to align as a "
        "fraction [0.0-1.0] of gene length. This can improve gene alignment where gene structure  differs between "
                                                           "target and "
                                                         "reference; by default F=0.0")

    parser.add_argument('-V', '--version', help='show program version', action='version', version='v1.4.2')
    parser.add_argument(
        '-p', default=1, type=int, metavar='P', help='use P parallel processes to accelerate alignment; by default P=1'
    )
    parser.add_argument('-m', help='Minimap2 path', metavar='PATH')
    parser.add_argument('-f', metavar='TYPES', help='list of feature types to lift over')
    parser.add_argument(
        '-infer_genes', action='store_true',
        help='use if annotation file only includes transcripts, exon/CDS features'
    )
    parser.add_argument(
        '-infer_transcripts', action='store_true', required=False,
        help='use if annotation file only includes exon/CDS features and does not include transcripts/mRNA'
    )
    parser.add_argument(
        '-chroms', metavar='TXT', help='comma seperated file with corresponding chromosomes in '
        'the reference,target sequences',
    )
    parser.add_argument(
        '-unplaced', metavar='TXT',
        help='text file with name(s) of unplaced sequences to map genes from after genes from '
        'chromosomes in chroms.txt are mapped; default is "unplaced_seq_names.txt"',
    )
    parser.add_argument('-copies', action='store_true', help='look for extra gene copies in the target genome')
    parser.add_argument(
        '-sc', default=1.0, metavar='SC', type=float,
        help='with -copies, minimum sequence identity in exons/CDS for which a gene is considered '
        'a copy; must be greater than -s; default is 1.0',
    )
    parser.add_argument('-overlap', default=0.1, metavar='O', help="maximum fraction [0.0-1.0] of overlap allowed by 2 "
                                                               "features; by default O=0.1", type=float)
    parser._positionals.title = 'Required input (sequences)'
    parser._optionals.title = 'Miscellaneous settings'
    parser._action_groups = [parser._positionals, refrgrp, outgrp, aligngrp, parser._optionals]
    args = parser.parse_args(arglist)
    if (float(args.s) > float(args.sc)):
        parser.error("-sc must be greater than or equal to -s")
    if (args.chroms is None and args.unplaced is not None):
        parser.error("-unplaced must be used with -chroms")
    return args


def parse_chrm_files(chroms_file):
    chroms = open(chroms_file, 'r')
    ref_chroms, target_chroms = [], []
    for line in chroms.readlines():
        ref_and_target_chrom = line.rstrip().split(",")
        ref_chroms.append(ref_and_target_chrom[0])
        if len(ref_and_target_chrom) > 1:
            target_chroms.append(ref_and_target_chrom[1])
    chroms.close()
    return ref_chroms, target_chroms


def get_parent_features_to_lift(feature_types_file):
    feature_types = ["gene"]
    if feature_types_file is not None:
        f = open(feature_types_file)
        for line in f.readlines():
            feature_types.append(line.rstrip())
    return feature_types


def map_unmapped_features(unmapped_features, target_chroms, lifted_feature_list, feature_db, feature_hierarchy,
                          ref_parent_order, args):
    if len(unmapped_features) > 0 and target_chroms[0] != args.target:
        print("mapping unaligned features to whole genome")
        ref_chroms = [args.reference]
        target_chroms = [args.target]
        return liftover_types.map_unmapped_genes_agaisnt_all(unmapped_features, ref_chroms, target_chroms,
                                                             lifted_feature_list, feature_db, feature_hierarchy,
                                                             ref_parent_order, args)
    return unmapped_features


def map_features_from_unplaced_seq(unmapped_features, lifted_feature_list, feature_db, feature_hierarchy,
                                   ref_parent_order, args):
    if args.unplaced is not None and args.chroms is not None:
        print("mapping unplaced genes")
        ref_chroms, target_chroms = parse_chrm_files(args.unplaced)
        target_chroms = [args.target]
        liftover_types.map_unplaced_genes(unmapped_features, ref_chroms, target_chroms,
                                          lifted_feature_list, feature_db, feature_hierarchy, ref_parent_order, args)


def write_unmapped_features_file(out_arg, unmapped_features):
    unmapped_out = open(out_arg, 'w')
    for gene in unmapped_features:
        unmapped_out.write(gene.id + "\n")
    unmapped_out.close()


def map_extra_copies(args, lifted_feature_list, feature_hierarchy, feature_db, ref_parent_order):
    if args.copies:
        print("mapping gene copies")
        ref_chroms = [args.reference]
        target_chroms = [args.target]
        liftover_types.map_extra_copies(ref_chroms, target_chroms, lifted_feature_list, feature_hierarchy, feature_db,
                                        ref_parent_order, args)


if __name__ == "__main__":
    main()
