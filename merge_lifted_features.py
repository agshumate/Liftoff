import process_blast_alignments as pba
import gffutils


def get_feature_order(gene_db):
    feature_types = list(gene_db.featuretypes())
    index = 0
    feature_order = {}
    if 'exon' in feature_types:
        feature_order['exon'] = index
        index += 1
    if 'CDS' in feature_types:
        feature_order['CDS'] = index
        index += 1
    for feature_type in feature_types:
        if feature_type not in feature_order:
            feature_order[feature_type] = index
            index +=1
    return feature_order



def merge_lifted_features(mapped_children, shortest_path_weight, gene_db, gene, unmapped_genes, threshold, gene_name):
    feature_list = []
    orphans=[]
    for child_name in mapped_children:
        #print(child_name)
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
                child_feature.score=shortest_path_weight
                child_feature.attributes["copy_id"] = gene_name
    while(len(orphans) != 0):
        orphans=create_parents(orphans, gene_db, feature_list, shortest_path_weight, gene_name)
    feature_order = get_feature_order(gene_db)
    feature_list.sort(key=lambda x: (x.seqid, x.start, feature_order[x.featuretype]))
    return feature_list




def create_parents(orphans, gene_db, feature_list, shortest_path_weight, gene_name):
    added_parent_ids = []
    new_orphans=[]
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
                parent_feature.score=shortest_path_weight
                parent_feature.attributes["copy_id"] = gene_name
            added_parent_ids.append(parent.id)
    return new_orphans









