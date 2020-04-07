from get_gene_seqs import get_parent_features
def write_line(line, out_file):
    if out_file == "stdout":
        print(line)
    else:
<<<<<<< HEAD
        out_file.write(str(line))
=======
        out_file.write(line)
>>>>>>> d5e2965468bb302f52c570bda0ab3d92e52c0da1
        out_file.write("\n")





def write_new_gff(lifted_features, out_file, gene_db):
    copy_num_dict ={}
    if out_file != 'stdout':
        f=open(out_file, 'w')
    else:
        f="stdout"
    all_features= []
    for feature in lifted_features:
        all_features.extend(lifted_features[feature])
    all_features.sort(key=lambda x: (x.seqid, x.start))
    parents = [feature for feature in all_features if "Parent" not in feature.attributes]
    for parent in parents:
        child_features = lifted_features[parent.attributes["copy_id"][0]]
        parent.score = "."
        if parent.id in copy_num_dict:
            copy_num_dict[parent.id] +=1
        else:
            copy_num_dict[parent.id] =1

        copy_num=copy_num_dict[parent.id]
        parent.attributes["copy_number"]=str(copy_num)
        if parent.attributes["coverage"][0] < 0.5:
            parent.attributes["partial_mapping"] = "True"
        del parent.attributes["coverage"]
<<<<<<< HEAD
        write_feature([parent], f, child_features)
=======
        write_feature([parent], out_file, child_features)
>>>>>>> d5e2965468bb302f52c570bda0ab3d92e52c0da1





def write_feature(children, outfile, child_features):
    if len(children) == 0:
        return
    else:
        for child in children:
            write_line(child, outfile)
            new_children= [feature for feature in child_features if "Parent" in feature.attributes and feature.attributes["Parent"][0] == child.id]
            write_feature(new_children, outfile, child_features)
