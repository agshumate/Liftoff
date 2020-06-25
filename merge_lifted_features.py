import gffutils
import liftoff_utils
import copy


def meets_length_threshold(feature, threshold, gene):
    original_feature = gene
    original_length = original_feature.end - original_feature.start
    lifted_length = feature.end - feature.start
    coverage = lifted_length/float(original_length)
    feature.attributes["coverage"] = coverage
    return coverage >= threshold


def merge_lifted_features(mapped_children, gene, unmapped_genes, threshold,
                          gene_name, feature_order, parent_dict, intermediate_dict, aln_cov, seq_id, seq_id_threshold):
    feature_list, final_features = {}, []
    orphans=[]
    top_feature = None
    if len(mapped_children) == 0:
        unmapped_genes.append(gene)
        return [], 0
    for child_name in mapped_children:
        child_feature = mapped_children[child_name]
        feature_list[child_feature.id]=child_feature
        if child_feature.id != gene.id:
            orphans.append((child_feature, child_feature["Parent"][0]))
        else:
            top_feature = child_feature
    while(len(orphans) != 0):
        orphans, top_feature=create_parents(orphans, parent_dict, feature_list, intermediate_dict, gene)
    final_features = list(feature_list.values())
    final_features.sort(key=lambda x: (x.seqid, x.start))
    final_features.sort(key=lambda x:feature_order[x.featuretype])
    if aln_cov < threshold or seq_id < seq_id_threshold:
        final_features = []
        unmapped_genes.append(gene)
    else:
        top_feature.score = 1- seq_id
        top_feature.attributes["copy_id"] = gene_name
        top_feature.attributes["coverage"] = str(round(aln_cov,5))
        top_feature.attributes["sequence_ID"] = str(round(seq_id,5))
    return final_features, top_feature.start




def create_parents(orphans, parent_dict, feature_list, intermediate_dict, gene):
    added_parent_ids = []
    new_orphans=[]
    top_feature = None
    for orphan, parent in orphans:
        if parent not in added_parent_ids:
            children = [feature for feature in feature_list.values() if
                        "Parent" in feature.attributes and feature.attributes["Parent"][0] == parent]
            starts = [child.start for child in children]
            ends = [child.end for child in children]
            if parent in parent_dict:
                original_parent = parent_dict[parent]
            else:
                original_parent = intermediate_dict[parent]
            parent_feature = liftoff_utils.make_new_feature(copy.deepcopy(original_parent), min(starts), max(ends), children[0].strand, children[0].seqid)
            feature_list[parent_feature.id] = parent_feature
            if parent_feature.id != gene.id:
                new_orphans.append((parent_feature, parent_feature["Parent"][0]))
            else:
                top_feature = parent_feature
            added_parent_ids.append(parent)
    return new_orphans, top_feature









