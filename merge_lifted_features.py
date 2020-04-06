import process_blast_alignments as pba
import gffutils


def meets_length_threshold(feature, threshold, gene_db):
    original_feature = gene_db[feature]
    original_length = original_feature.end - original_feature.start
    lifted_length = feature.end - feature.start
    coverage = lifted_length/float(original_length)
    feature.attributes["coverage"] = coverage
    return coverage >= threshold


def merge_lifted_features(mapped_children, shortest_path_weight, gene_db, gene, unmapped_genes, threshold, gene_name, feature_order):
    feature_list = []
    orphans=[]
    top_feature = None
    if len(mapped_children) == 0:
        unmapped_genes.append(gene)
        return []
    for child_name in mapped_children:
        lifted_child = mapped_children[child_name]
        if lifted_child.chrm != "start" and lifted_child.chrm!="end":
            original_child = gene_db[child_name]
            child_feature = gffutils.Feature(seqid=lifted_child.chrm, source="Liftoff", featuretype=original_child.featuretype,
                                            start=lifted_child.lifted_start + 1, end=lifted_child.lifted_end + 1,
                                            strand=lifted_child.strand, id=original_child.id, attributes=original_child.attributes)
            parent_list = list(gene_db.parents(original_child, level=1))
            feature_list.append(child_feature)
            if len(parent_list) !=0:
                child_feature.attributes["Parent"]=parent_list[0].id
                orphans.append((child_feature, parent_list[0]))
            else:
                top_feature = child_feature
    while(len(orphans) != 0):
        orphans, top_feature=create_parents(orphans, gene_db, feature_list)
    feature_list.sort(key=lambda x: (x.seqid, x.start, feature_order[x.featuretype]))
    if meets_length_threshold(top_feature, threshold, gene_db) is False:
        feature_list = []
        unmapped_genes.append(gene)
    else:
        top_feature.score = shortest_path_weight
        top_feature.attributes["copy_id"] = gene_name
        if "Parent" in top_feature.attributes:
            del top_feature.attributes["Parent"]
    return feature_list




def create_parents(orphans, gene_db, feature_list):
    added_parent_ids = []
    new_orphans=[]
    top_feature = None
    for orphan, parent in orphans:
        if parent.id not in added_parent_ids:
            #print(parent.id)
            children=[feature for feature in feature_list if "Parent" in feature.attributes and feature.attributes["Parent"][0] == parent.id]
            starts = [child.start for child in children]
            ends = [child.end for child in children]
            original_parent = gene_db[parent.id]
            parent_feature = gffutils.Feature(seqid=children[0].seqid, source="Liftoff",
                                             featuretype=original_parent.featuretype,
                                             start=min(starts), end=max(ends),
                                             strand=children[0].strand, id=original_parent.id, attributes=original_parent.attributes)
            new_parent_list = list(gene_db.parents(original_parent, level=1))
            feature_list.append(parent_feature)
            if len(new_parent_list) !=0:
                parent_feature.attributes["Parent"]=new_parent_list[0].id
                new_orphans.append((parent_feature, new_parent_list[0]))
            else:
                top_feature = parent_feature
            added_parent_ids.append(parent.id)
    return new_orphans, top_feature









