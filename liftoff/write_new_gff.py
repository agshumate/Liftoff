from liftoff import liftoff_utils


def write_line(feature, out_file):
    feature = delete_attributes(feature)
    if out_file == "stdout":
        print(feature)
    else:
        line=make_gff_line(feature)
        out_file.write(line)
        out_file.write("\n")

def make_gff_line(feature):
    attributes_str = ""
    for attr in feature.attributes:
        value_str = ""
        for value in feature.attributes[attr]:
            value_str += value + ","
        attributes_str += (attr +"=" + value_str[:-1] + ";")
    return feature.seqid + "\t" + feature.source + "\t" + feature.featuretype + "\t" + str(feature.start) + \
           "\t" + str(feature.end) + "\t" + "." + "\t" + feature.strand + "\t" + "." + "\t" + attributes_str[:-1]


def delete_attributes(line):
    if "coverage" in line.attributes:
        line.attributes["coverage"] = str(line.attributes["coverage"][0])
    return line


def write_new_gff(lifted_features, out_file,  parent_dict, cov_threshold,seq_threshold ):
    copy_num_dict ={}
    if out_file != 'stdout':
        f=open(out_file, 'w')
    else:
        f="stdout"
    parents = liftoff_utils.get_parent_list(lifted_features, parent_dict)
    parents.sort(key=lambda x: x.id)
    final_parent_list = []
    for parent in parents:
        parent.score = "."
        if parent.id in copy_num_dict:
            copy_num_dict[parent.id] +=1
        else:
            copy_num_dict[parent.id] =0
        copy_num=copy_num_dict[parent.id]
        parent.attributes["extra_copy_number"]=str(copy_num)
        if float(parent.attributes["coverage"][0]) < cov_threshold:
            parent.attributes["partial_mapping"] = "True"
        if float(parent.attributes["sequence_ID"][0]) < seq_threshold:
            parent.attributes["low_identity"] = "True"
        final_parent_list.append(parent)
    final_parent_list.sort(key=lambda x: (x.seqid, x.start))
    for final_parent in final_parent_list:
        child_features = lifted_features[final_parent.attributes["copy_id"][0]]
        parent_child_dict = build_parent_dict(child_features, parent_dict)
        write_feature([final_parent], f, child_features, parent_child_dict)


       
def build_parent_dict(child_features, parent_dict):
    parent_child_dict = {}
    for child in child_features:

        if child.id not in parent_dict:
            if child.attributes["Parent"][0] in parent_child_dict:
                parent_child_dict[child.attributes["Parent"][0]].append(child)
            else :
                parent_child_dict[child.attributes["Parent"][0]] = [child]
    return parent_child_dict



def write_feature(children, outfile, child_features, parent_dict):
    for child in children:
        write_line(child, outfile)
        if child.id in parent_dict:
            new_children= parent_dict[child.id]
            write_feature(new_children, outfile, child_features, parent_dict)
    return
