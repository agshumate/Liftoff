import liftoff_utils

def write_line(line, out_file):
    line = delete_attributes(line)
    if out_file == "stdout":
        print(line)
    else:
        out_file.write(str(line)) 
        out_file.write("\n")



def delete_attributes(line):
    if "coverage" in line.attributes:
        line.attributes["coverage"] = str(line.attributes["coverage"][0])
    return line


def write_new_gff(lifted_features, out_file,  parent_dict ):
    copy_num_dict ={}
    if out_file != 'stdout':
        f=open(out_file, 'w')
    else:
        f="stdout"
    parents = liftoff_utils.get_parent_list(lifted_features, parent_dict)
    for parent in parents:
        child_features = lifted_features[parent.attributes["copy_id"][0]]
        parent_child_dict = build_parent_dict(child_features, parent_dict)
        parent.score = "."
        if parent.id in copy_num_dict:
            copy_num_dict[parent.id] +=1
        else:
            copy_num_dict[parent.id] =1

        copy_num=copy_num_dict[parent.id]
        parent.attributes["copy_number"]=str(copy_num)

        if float(parent.attributes["coverage"][0]) < 0.5:
            parent.attributes["partial_mapping"] = "True"
        write_feature([parent], f, child_features, parent_child_dict)


       
def build_parent_dict(child_features, parent_dict):
    parent_child_dict = {}
    for child in child_features:
        if child.id not in parent_dict:
            if child.attributes["Parent"][0] in parent_child_dict:
                parent_child_dict[child.attributes["Parent"][0] ].append(child)
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
