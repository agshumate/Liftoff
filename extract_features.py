import gffutils
from multiprocessing import Pool
from pyfaidx import Fasta, Faidx
from functools import partial
import liftoff_utils


def extract_features_to_lift(g_arg, db_arg, ref_chroms, reference_fasta, processes, infer_transcripts):
    print("extracting features")
    feature_db, feature_db_name = create_feature_db_connections(g_arg, db_arg, infer_transcripts)
    parent_dict, child_dict, intermediate_dict, parent_order= seperate_parents_and_children(feature_db)
    get_gene_sequences(parent_dict, ref_chroms, reference_fasta, processes)
    return parent_dict, child_dict, intermediate_dict, feature_db, parent_order


def create_feature_db_connections(g_arg, db_arg, infer_transcripts):
    if infer_transcripts is True:
        disable = False
    else:
        disable = True
    if db_arg is None:
        feature_db = gffutils.create_db(g_arg, g_arg + "_db", merge_strategy="create_unique", force=True,
                                 disable_infer_transcripts=disable, disable_infer_genes=True, verbose=True)
        feature_db_name = db_arg
    else:
        feature_db = gffutils.FeatureDB(db_arg)
        feature_db_name = db_arg
    return feature_db, feature_db_name


def seperate_parents_and_children(feature_db):
    parent_dict, child_dict, intermediate_dict = {} , {}, {}
    feature_types = feature_db.featuretypes()
    parent_types, child_types, intermediate_types = find_feature_type_order(feature_types, feature_db)
    feature_num = 0
    for parent_type in parent_types:
        for parent in feature_db.features_of_type(featuretype=parent_type):
            feature_num += 1
            parent_dict[parent.id] = parent
            child_dict[parent.id] = []
            for child_type in child_types:
                for child in feature_db.children(parent, featuretype=child_type):
                    child_dict[parent.id].append(child)
                    if "Parent" not in child.attributes:
                        add_parent_tag(child, feature_db)
            if len(child_dict[parent.id]) == 0:
                child_dict[parent.id].append(parent)
    for intermediate_type in intermediate_types:
        for intermediate_feature in feature_db.features_of_type(featuretype=intermediate_type):
            intermediate_dict[intermediate_feature.id] = intermediate_feature
            if "Parent" not in intermediate_feature.attributes:
                add_parent_tag(intermediate_feature, feature_db)
    parent_order = liftoff_utils.find_parent_order(
        [parent for parent in list(parent_dict.values()) if parent is not None])
    return parent_dict, child_dict, intermediate_dict, parent_order


def add_parent_tag(feature, feature_db):
    parent = list(feature_db.parents(feature, level=1))[0]
    feature.attributes["Parent"]= parent.id


def find_feature_type_order(feature_types, feature_db):
    parent_types, intermediate_types, child_types = [],[],[]
    for feature_type in feature_types:
        for feature in feature_db.features_of_type(featuretype=feature_type):
            if has_child(feature, feature_db) is False:
                child_types.append(feature_type)
            elif has_parent(feature, feature_db) is False:
                parent_types.append(feature_type)
            else:
                intermediate_types.append(feature_type)
            break
    return parent_types, child_types, intermediate_types



def has_parent(feature, feature_db):
    for parent in feature_db.parents(feature.id, level=1):
        if feature.id != parent.id:
            return True
    return False


def has_child(feature, feature_db):
    for child in feature_db.children(feature.id):
        return True
    return False



def get_gene_sequences(parent_dict, ref_chroms, reference_fasta_name, processes):
    pool = Pool(processes)
    Faidx(reference_fasta_name)
    func = partial(get_gene_sequences_subset, parent_dict, reference_fasta_name)
    for result in pool.imap_unordered(func, ref_chroms):
        continue
    pool.close()
    pool.join()
    return


def get_gene_sequences_subset(parent_dict, reference_fasta_name, chrom_name):
    reference_fasta = Fasta(reference_fasta_name)
    fasta_out = open(chrom_name + "_genes.fa", 'w')
    sorted_parents = sorted(list(parent_dict.values()), key=lambda x: x.seqid)
    if chrom_name == reference_fasta_name:
        current_chrom = sorted_parents[0].seqid
    else:
        current_chrom = chrom_name
    chrom_seq = reference_fasta[current_chrom]
    for parent in sorted_parents:
        if parent.seqid != current_chrom and chrom_name == reference_fasta_name:
            current_chrom=parent.seqid
            chrom_seq = reference_fasta[current_chrom]
        if parent.seqid == chrom_name or chrom_name == reference_fasta_name:
            parent_seq = chrom_seq[parent.start-1: parent.end].seq
            fasta_out.write(">" + parent.id + "\n" + str(parent_seq) + "\n")
    fasta_out.close()
    return





