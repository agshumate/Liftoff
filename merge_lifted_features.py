import gffutils


def meets_length_threshold(feature, threshold, gene):
    original_feature = gene
    original_length = original_feature.end - original_feature.start
    lifted_length = feature.end - feature.start
    coverage = lifted_length/float(original_length)
    feature.attributes["coverage"] = coverage
    return coverage >= threshold


def merge_lifted_features(mapped_children, shortest_path_weight, gene, unmapped_genes, threshold,
                          gene_name, feature_order, parent_dict, intermediate_dict, aln_cov, seq_id):
    feature_list = []
    orphans=[]
    top_feature = None
    if len(mapped_children) == 0:
        unmapped_genes.append(gene)
        return [], 0
    for child_name in mapped_children:
        child_feature = mapped_children[child_name]
        feature_list.append(child_feature)
        if child_feature.id != gene.id:
            orphans.append((child_feature, child_feature["Parent"][0]))
        else:
            top_feature = child_feature
    while(len(orphans) != 0):
        orphans, top_feature=create_parents(orphans, parent_dict, feature_list, intermediate_dict, gene)
    feature_list.sort(key=lambda x: (x.seqid, x.start))
    feature_list.sort(key=lambda x:feature_order[x.featuretype])
    if aln_cov < threshold:
        feature_list = []
        unmapped_genes.append(gene)
    else:
        top_feature.score = shortest_path_weight/(gene.end - gene.start)
        top_feature.attributes["copy_id"] = gene_name
        top_feature.attributes["coverage"] = round(aln_cov,5)
        top_feature.attributes["sequence_ID"] = str(round(seq_id,5))
    return feature_list, top_feature.start



def create_parents(orphans, parent_dict, feature_list, intermediate_dict, gene):
    added_parent_ids = []
    new_orphans=[]
    top_feature = None
    for orphan, parent in orphans:
        if parent not in added_parent_ids:
            children = [feature for feature in feature_list if
                        "Parent" in feature.attributes and feature.attributes["Parent"][0] == parent]
            starts = [child.start for child in children]
            ends = [child.end for child in children]
            if parent in parent_dict:
                original_parent = parent_dict[parent]
            else:
                original_parent = intermediate_dict[parent]
            parent_feature = gffutils.Feature(seqid=children[0].seqid, source="Liftoff",
                                             featuretype=original_parent.featuretype,
                                             start=min(starts), end=max(ends),
                                             strand=children[0].strand, id=original_parent.id,
                                            attributes=original_parent.attributes)
            feature_list.append(parent_feature)
            if parent_feature.id != gene.id:
                new_orphans.append((parent_feature, parent_feature["Parent"][0]))
            else:
                top_feature = parent_feature
            added_parent_ids.append(parent)
    return new_orphans, top_feature









