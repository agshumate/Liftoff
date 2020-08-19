from liftoff import write_new_gff, liftover_types
import argparse


def main():
    args = parse_args()
    if args.chroms is not None:
        ref_chroms, target_chroms = parse_chrm_files(args.chroms)
    else:
        ref_chroms = [args.r]
        target_chroms = [args.t]
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



def parse_args():
    parser = argparse.ArgumentParser(description='Lift features from one genome assembly to another')
    group = parser.add_mutually_exclusive_group(required=True)
    parser.add_argument("-V", "--version", help="show program version", action='version', version="v1.3.0")
    parser.add_argument('-t', required=True, help="target fasta genome to lift genes to", metavar="<target.fasta>")
    parser.add_argument('-r', required=True, help="reference fasta genome to lift genes from",
                        metavar="<reference.fasta>")
    group.add_argument('-g', help="annotation file to lift over in gff or gtf format", metavar="<ref_annotation.gff>")
    parser.add_argument('-chroms', required=False,
                        help="comma seperated file with corresponding chromosomes in the reference,target sequences",
                        metavar="<chroms.txt>")
    parser.add_argument('-p', required=False, help="processes", default=1, metavar=1, type=int)
    parser.add_argument('-o', required=False, help="output file", default='stdout', metavar="<output.gff>")
    group.add_argument('-db',
                       help="name of feature database. If none, -g argument must be provided and a database will be "
                            "built automatically")
    parser.add_argument('-infer_transcripts', action='store_true', required=False,
                        help="use if GTF file only includes exon/CDS features and does not include transcripts/mRNA")
    parser.add_argument('-u', required=False, help="name of file to write unmapped features to",
                        default="unmapped_features.txt", metavar="<unmapped_features.txt>")
    parser.add_argument('-infer_genes', action='store_true', required=False,
                        help="use if GTF file only includes transcripts, exon/CDS features")
    parser.add_argument('-a', required=False, help="minimum alignment coverage to consider a feature mapped [0-1]",
                        default=0.5, metavar=0.5, type=float)
    parser.add_argument('-s', required=False,
                        help="minimum sequence identity in child features (usually exons/CDS) to consider a feature "
                             "mapped [0-1]",
                        default=0.5, metavar=0.5, type=float)
    parser.add_argument('-unplaced', required=False,
                        help="text file with name(s) of unplaced sequences to map genes from after genes from "
                             "chromosomes in chroms.txt are mapped",
                        metavar="<unplaced_seq_names.txt>")
    parser.add_argument('-copies', required=False, action='store_true',
                        help="look for extra gene copies in the target genome")
    parser.add_argument('-sc', required=False,
                        help='with -copies, minimum sequence identity in exons/CDS for which a gene is considered a '
                             'copy. Must be greater than -s',
                        default=1.0, metavar=1.0, type=float)
    parser.add_argument('-m', required=False, help="Minimap2 path", metavar='PATH')
    parser.add_argument('-dir', required=False, help="name of directory to save intermediate fasta and SAM files",
                        default="intermediate_files", metavar="<intermediate_files_dir>")
    parser.add_argument('-n', required=False, default=50, metavar=50,
                        help="max number of Minimap2 alignments to consider for each feature", type=int)
    parser.add_argument('-f', required=False, metavar="feature types", help="list of feature types to lift-over")
    parser.add_argument('-d', required=False, metavar=2, default=2, help="distance scaling factor. Alignment nodes " \
                                                                       "father apart "
                                                              "than this in the target genome will not be connected in "
                                                              "the graph", type=float)
    parser.add_argument('-exclude_partial', default=False, action='store_true',
                        help="write partial mappings below -s and -a threshold to unmapped_features.txt. If true "
                             "partial/low sequence identity mappings will be included in the gff file with "
                             "partial_mapaping=True, low_identity=True in comments")
    parser.add_argument('-flank', default=0, metavar=0, type=float, help="amount of flanking sequence to align as a "
                                                                         "percentage of gene length. This can improve"
                                                                         " gene alignment where gene structure "
                                                                         "differs between target and reference" )
    parser.add_argument('-frag', default =False, action='store_true', help="allow genes to map across contigs for "
                                                                           "fragmented assemblies")
    parser.add_argument('-scaffold', default=False, action='store_true')
    parser.add_argument('-overlap', default=0.1, metavar=0.1, help="maximum fraction of overlap allowed by 2 "
                                                                   "features", type=float)
    args = parser.parse_args()
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
    if len(unmapped_features) > 0 and target_chroms[0] != args.t:
        print("mapping unaligned features to whole genome")
        ref_chroms = [args.r]
        target_chroms = [args.t]
        return liftover_types.map_unmapped_genes_agaisnt_all(unmapped_features, ref_chroms, target_chroms,
                                                             lifted_feature_list, feature_db, feature_hierarchy,
                                                             ref_parent_order, args)
    return unmapped_features


def map_features_from_unplaced_seq(unmapped_features, lifted_feature_list, feature_db, feature_hierarchy,
                                   ref_parent_order, args):
    if args.unplaced is not None and args.chroms is not None:
        print("mapping unplaced genes")
        ref_chroms, target_chroms = parse_chrm_files(args.unplaced)
        target_chroms = [args.t]
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
        ref_chroms = [args.r]
        target_chroms = [args.t]
        liftover_types.map_extra_copies(ref_chroms, target_chroms, lifted_feature_list, feature_hierarchy, feature_db,
                                        ref_parent_order, args)


if __name__ == "__main__":
    main()
