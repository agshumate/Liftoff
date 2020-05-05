import liftover_types
import write_new_gff
import argparse



def parse_args():
    parser = argparse.ArgumentParser(description='Lift features from one genome assembly to another')
    group = parser.add_mutually_exclusive_group(required=True)
    parser.add_argument('-t', required=True, metavar = "target fasta", help="target fasta genome to lift genes to")
    parser.add_argument('-r', required=True, metavar= "reference fasta", help = "reference fasta genome to lift genes from")
    group.add_argument('-g', metavar="gff or gtf", help="annotation file to lift over in gff or gtf format")
    parser.add_argument('-chroms', required=False, metavar="chromosome names",
                        help="comma seperated file with corresponding chromosomes in the reference,target sequences")
    parser.add_argument('-p', required=False, metavar ="num processess", help= "processes", default=1)
    parser.add_argument('-o', required=False, metavar= "outfile", help="output file", default = 'stdout')
    group.add_argument('-db', metavar="feature database",
                       help="name of feature database. If none, -g argument must be provided and a database will be built automatically")
    parser.add_argument('-infer_transcripts', action='store_true', required=False,
                        help="use if GTF file only includes exon/CDS features")
    args = parser.parse_args()
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

    #read chroms
    if chroms_file is not None:
        ref_chroms, target_chroms = parse_chrm_files(chroms_file)
    else:
        ref_chroms=[reference_fasta]
        target_chroms= [target_fasta]

    #lift genes
    lifted_feature_list = {}
    unmapped_features = []


    feature_db, parent_features, intermediate_features, children_features, parent_order = liftover_types.lift_original_annotation(
    gff, target_fasta, reference_fasta, ref_chroms, target_chroms, processes, db, lifted_feature_list,
    unmapped_features, infer_transcripts)
    unmapped_out = open("unmapped_features", 'w')
    if len(unmapped_features) > 0 and target_chroms[0] != target_fasta:
       ref_chroms = [reference_fasta]
       target_chroms = [target_fasta]
       unmapped_features = liftover_types.map_unmapped_genes_agaisnt_all(unmapped_features, target_fasta,
                                                                         reference_fasta, ref_chroms, target_chroms,
                                                                         processes, lifted_feature_list, feature_db,
                                                                         parent_features, intermediate_features,
                                                                         children_features, parent_order)

    for gene in unmapped_features:
        unmapped_out.write(gene.id + "\n")
    unmapped_out.close()
    write_new_gff.write_new_gff(lifted_feature_list, output,  parent_features)


def parse_chrm_files(chroms_file):
    chroms = open(chroms_file, 'r')
    ref_chroms, target_chroms = [], []
    for line in chroms.readlines():
        ref_and_target_chrom= line.rstrip().split(",")
        ref_chroms.append(ref_and_target_chrom[0])
        target_chroms.append(ref_and_target_chrom[1])
    return ref_chroms, target_chroms


main()