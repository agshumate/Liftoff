from liftoff import write_new_gff, liftover_types
import argparse



def parse_args():
    parser = argparse.ArgumentParser(description='Lift features from one genome assembly to another')
    group = parser.add_mutually_exclusive_group(required=True)
    parser.add_argument("-V", "--version", help="show program version", action='version', version="v1.2.0" )
    parser.add_argument('-t', required=True, help="target fasta genome to lift genes to", metavar="<target.fasta>")
    parser.add_argument('-r', required=True, help = "reference fasta genome to lift genes from", metavar="<reference.fasta>")
    group.add_argument('-g', help="annotation file to lift over in gff or gtf format", metavar="<ref_annotation.gff>")
    parser.add_argument('-chroms', required=False,
                        help="comma seperated file with corresponding chromosomes in the reference,target sequences", metavar="<chroms.txt>")
    parser.add_argument('-p', required=False, help= "processes", default=1, metavar=1)
    parser.add_argument('-o', required=False,  help="output file", default = 'stdout', metavar="<output.gff>")
    group.add_argument('-db',
                       help="name of feature database. If none, -g argument must be provided and a database will be built automatically")
    parser.add_argument('-infer_transcripts', action='store_true', required=False,
                        help="use if GTF file only includes exon/CDS features and does not include transcripts/mRNA")
    parser.add_argument('-u', required=False, help= "name of file to write unmapped features to", default="unmapped_features.txt", metavar="<unmapped_features.txt>")
    parser.add_argument('-infer_genes', action='store_true', required=False,
                        help="use if GTF file only includes transcripts, exon/CDS features")
    parser.add_argument('-a', required=False, help="minimum alignment coverage to consider a feature mapped [0-1]", default=0.5, metavar=0.5)
    parser.add_argument('-s', required=False, help="minimum sequence identity in child features (usually exons/CDS) to consider a feature mapped [0-1]", default=0.5, metavar=0.5)
    parser.add_argument('-unplaced', required=False,
                        help="text file with name(s) of unplaced sequences to map genes from after genes from chromosomes in chroms.txt are mapped",
                        metavar="<unplaced_seq_names.txt>")
    parser.add_argument('-copies', required=False, action='store_true',
                        help="look for extra gene copies in the target genome")
    parser.add_argument('-sc', required=False, help='with -copies, minimum sequence identity in exons/CDS for which a gene is considered a copy. Must be greater than -s', default = 1.0, metavar=1.0 )
    parser.add_argument('-m', required=False, help="Minimap2 path", metavar='PATH')
    parser.add_argument('-dir', required=False, help="name of directory to save intermediate fasta and SAM files", default="intermediate_files", metavar="<intermediate_files_dir>")
    parser.add_argument('-n', required=False, default=50, metavar=50, help="max number of Minimap2 alignments to consider for each feature")
    parser.add_argument('-f', required=False, metavar="feature types", help="list of feature types to lift-over")
    args = parser.parse_args()
    if (float(args.s) > float(args.sc)):
        parser.error("-sc must be greater than or equal to -s")
    if (args.chroms is None and args.unplaced is not None):
        parser.error("-unplaced must be used with -chroms")
    return args


def main():
   #parse args
    args = parse_args()
    gff = args.g
    target_fasta = args.t
    reference_fasta = args.r
    processes =int(args.p)
    output = args.o
    chroms_file = args.chroms
    db = args.db
    infer_transcripts = args.infer_transcripts
    infer_genes = args.infer_genes
    unplaced_seq = args.unplaced
    copies = args.copies
    cov_threshold = float(args.a)
    seq_threshold = float(args.s)
    seq_threshold_copies = float(args.sc)
    minimap2_path = args.m
    inter_files = args.dir
    max_alns = int(args.n)
    feature_types_file=args.f

    #read chroms
    if chroms_file is not None:
        ref_chroms, target_chroms = parse_chrm_files(chroms_file)
    else:
        ref_chroms=[reference_fasta]
        target_chroms= [target_fasta]

    parent_features_to_lift = get_parent_features_to_lift(feature_types_file)
    #lift genes
    lifted_feature_list = {}
    unmapped_features = []


    feature_db, parent_features, intermediate_features, children_features, parent_order = liftover_types.lift_original_annotation(
    gff, target_fasta, reference_fasta, ref_chroms, target_chroms, processes, db, lifted_feature_list,
    unmapped_features, infer_transcripts, infer_genes, cov_threshold, seq_threshold, minimap2_path, inter_files, max_alns, parent_features_to_lift)
    unmapped_out = open(args.u, 'w')
    if len(unmapped_features) > 0 and target_chroms[0] != target_fasta:
       print("mapping unaligned features to whole genome")
       ref_chroms = [reference_fasta]
       target_chroms = [target_fasta]
       unmapped_features = liftover_types.map_unmapped_genes_agaisnt_all(unmapped_features, target_fasta,
                                                                         reference_fasta, ref_chroms, target_chroms,
                                                                         processes, lifted_feature_list, feature_db,
                                                                         parent_features, intermediate_features,
                                                                         children_features, parent_order, minimap2_path,inter_files, max_alns)
    if unplaced_seq is not None and chroms_file is not None:
        print("mapping unplaced genes")
        ref_chroms, target_chroms = parse_chrm_files(unplaced_seq)
        target_chroms = [target_fasta]
        liftover_types.map_unplaced_genes(unmapped_features, target_fasta, reference_fasta, ref_chroms, target_chroms,
                                       processes, lifted_feature_list, feature_db, parent_features, intermediate_features,
                                       children_features, parent_order,minimap2_path,inter_files, max_alns)

    for gene in unmapped_features:
        unmapped_out.write(gene.id + "\n")
    unmapped_out.close()
    if copies:
       print("mapping gene copies")
       ref_chroms = [reference_fasta]
       target_chroms = [target_fasta]
       remap = chroms_file is not None
       liftover_types.map_extra_copies(target_fasta, reference_fasta, ref_chroms, target_chroms, processes,
       lifted_feature_list, parent_features, children_features, feature_db, intermediate_features, parent_order, seq_threshold_copies, minimap2_path,inter_files, remap, max_alns)



    write_new_gff.write_new_gff(lifted_feature_list, output,  parent_features, cov_threshold, seq_threshold)


def parse_chrm_files(chroms_file):
    chroms = open(chroms_file, 'r')
    ref_chroms, target_chroms = [], []
    for line in chroms.readlines():
        ref_and_target_chrom= line.rstrip().split(",")
        ref_chroms.append(ref_and_target_chrom[0])
        if len(ref_and_target_chrom) > 1:
            target_chroms.append(ref_and_target_chrom[1])
    chroms.close()
    return ref_chroms, target_chroms

def get_parent_features_to_lift(feature_types_file):
    if feature_types_file is None:
        return ["gene"]
    else:
        f= open(feature_types_file)
        return [line.rstrip() for line in f.readlines()]

if __name__ == "__main__":
    main()