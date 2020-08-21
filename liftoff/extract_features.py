import gffutils
from pyfaidx import Fasta, Faidx
from liftoff import liftoff_utils, feature_hierarchy
import os
import sys


def extract_features_to_lift(ref_chroms, liftover_type, parents_to_lift, args):
    print("extracting features")
    if os.path.exists(args.dir) is False:
        os.mkdir(args.dir)
    feature_db, feature_db_name = create_feature_db_connections(args)
    feature_hierarchy, parent_order = seperate_parents_and_children(feature_db, parents_to_lift)
    get_gene_sequences(feature_hierarchy.parents, ref_chroms, args, liftover_type)
    return feature_hierarchy, feature_db, parent_order


def create_feature_db_connections(args):
    if args.infer_transcripts is True:
        disable_transcripts = False
    else:
        disable_transcripts = True
    if args.infer_genes is True:
        disable_genes = False
    else:
        disable_genes = True
    if args.db is None:
        feature_db = gffutils.create_db(args.g, args.g + "_db", merge_strategy="create_unique", force=True,
                                        disable_infer_transcripts=disable_transcripts,
                                        disable_infer_genes=disable_genes, verbose=True)
        feature_db_name = args.db
    else:
        feature_db = gffutils.FeatureDB(args.db)
        feature_db_name = args.db
    return feature_db, feature_db_name


def seperate_parents_and_children(feature_db, parent_types_to_lift):
    parent_dict, child_dict, intermediate_dict = {}, {}, {}
    feature_types = feature_db.featuretypes()
    parent_types, child_types, intermediate_types = find_feature_type_hierarchy(feature_types, feature_db)
    filtered_parents = [parent for parent in parent_types if parent in parent_types_to_lift]
    for parent_type in filtered_parents:
        for parent in feature_db.features_of_type(featuretype=parent_type):
            add_parent(parent_dict, parent, child_dict)
            add_children(child_types, parent, feature_db, child_dict)
    add_intermediates(intermediate_types, intermediate_dict, feature_db)
    parent_order = liftoff_utils.find_parent_order(
        [parent for parent in list(parent_dict.values()) if parent is not None])
    ref_feature_hierarchy = feature_hierarchy.feature_hierarchy(parent_dict, intermediate_dict, child_dict)
    return ref_feature_hierarchy, parent_order


def find_feature_type_hierarchy(feature_types, feature_db):
    parent_types, intermediate_types, child_types = [], [], []
    for feature_type in feature_types:
        for feature in feature_db.features_of_type(featuretype=feature_type):
            is_child = has_child(feature, feature_db) is False
            is_parent = has_parent(feature, feature_db) is False
            if is_child:
                child_types.append(feature_type)
            if is_parent:
                parent_types.append(feature_type)
            if is_parent is False and is_child is False:
                intermediate_types.append(feature_type)
            break
    return parent_types, child_types, intermediate_types


def has_child(feature, feature_db):
    for child in feature_db.children(feature.id):
        return True
    return False


def has_parent(feature, feature_db):
    for parent in feature_db.parents(feature.id, level=1):
        if feature.id != parent.id:
            return True
    return False


def add_parent(parent_dict, parent, child_dict):
    parent_dict[parent.id] = parent
    child_dict[parent.id] = []


def add_children(child_types, parent, feature_db, child_dict):
    for child_type in child_types:
        for child in feature_db.children(parent, featuretype=child_type):
            child_dict[parent.id].append(child)
            if "Parent" not in child.attributes:
                add_parent_tag(child, feature_db)
    if len(child_dict[parent.id]) == 0:
        child_dict[parent.id].append(parent)


def add_parent_tag(feature, feature_db):
    parent_id = ""
    parents = list(feature_db.parents(feature, level=1))
    if len(parents) > 0:
        parent_id = parents[0].id
    else:
        parents = list(feature_db.parents(feature))
        if len(parents) > 0:
            parent_id = parents[0].id
    feature.attributes["Parent"] = parent_id


def add_intermediates(intermediate_types, intermediate_dict, feature_db):
    for intermediate_type in intermediate_types:
        for intermediate_feature in feature_db.features_of_type(featuretype=intermediate_type):
            intermediate_dict[intermediate_feature.id] = intermediate_feature
            if "Parent" not in intermediate_feature.attributes:
                add_parent_tag(intermediate_feature, feature_db)


def get_gene_sequences(parent_dict, ref_chroms, args, liftover_type):
    fai = Fasta(args.reference)
    if liftover_type == "unplaced":
        open(args.dir + "/unplaced_genes.fa", 'w')
    for chrom in ref_chroms:
        fasta_out = get_fasta_out(chrom, args.reference, liftover_type, args.dir)
        sorted_parents = sorted(list(parent_dict.values()), key=lambda x: x.seqid)
        if len(sorted_parents) == 0:
            sys.exit(
                "GFF does not contain any gene features. Use -f to provide a list of other feature types to lift over.")
        write_gene_sequences_to_file(chrom, args.reference, fai, sorted_parents, fasta_out, args)
        fasta_out.close()


def get_fasta_out(chrom_name, reference_fasta_name, liftover_type, inter_files):
    if chrom_name == reference_fasta_name and (liftover_type == "chrm_by_chrm" or liftover_type == "copies"):
        fasta_out_name = "reference_all"
    elif liftover_type == "unmapped":
        fasta_out_name = "unmapped_to_expected_chrom"
    elif liftover_type == "unplaced":
        fasta_out_name = "unplaced"
    else:
        fasta_out_name = chrom_name
    if liftover_type == "unplaced":
        fasta_out = open(inter_files + "/" + fasta_out_name + "_genes.fa", 'a')
    else:
        fasta_out = open(inter_files + "/" + fasta_out_name + "_genes.fa", 'w')
    return fasta_out


def write_gene_sequences_to_file(chrom_name, reference_fasta_name, reference_fasta_idx, parents, fasta_out, args):
    if chrom_name == reference_fasta_name:
        current_chrom = parents[0].seqid
    else:
        current_chrom = chrom_name
    chrom_seq = reference_fasta_idx[current_chrom]
    for parent in parents:
        if parent.seqid != current_chrom and chrom_name == reference_fasta_name:
            current_chrom = parent.seqid
            chrom_seq = reference_fasta_idx[current_chrom]
        if parent.seqid == chrom_name or chrom_name == reference_fasta_name:
            gene_length = parent.end - parent.start + 1
            parent.start = round(max(1, parent.start - args.flank * gene_length))
            parent.end = round(min(parent.end + args.flank * gene_length, len(chrom_seq) - 1))
            parent_seq = chrom_seq[parent.start - 1: parent.end].seq
            fasta_out.write(">" + parent.id + "\n" + str(parent_seq) + "\n")
