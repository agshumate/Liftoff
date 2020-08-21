from liftoff import liftoff_utils


def write_new_gff(lifted_features, parents_dict, args):
    if args.o != 'stdout':
        f = open(args.o, 'w')
    else:
        f = "stdout"
    parents = liftoff_utils.get_parent_list(lifted_features)
    parents.sort(key=lambda x: x.id)
    final_parent_list = finalize_parent_features(parents, args)
    final_parent_list.sort(key=lambda x: (x.seqid, x.start))
    for final_parent in final_parent_list:
        child_features = lifted_features[final_parent.attributes["copy_id"][0]]
        parent_child_dict = build_parent_dict(child_features, parents_dict, final_parent)
        write_feature([final_parent], f, child_features, parent_child_dict)


def finalize_parent_features(parents, args):
    final_parent_list = []
    copy_num_dict = {}
    for parent in parents:
        add_to_copy_num_dict(parent, copy_num_dict)
        copy_num = copy_num_dict[parent.id]
        add_attributes(parent, copy_num, args)
        final_parent_list.append(parent)
    return final_parent_list


def add_to_copy_num_dict(parent, copy_num_dict):
    if parent.id in copy_num_dict:
        copy_num_dict[parent.id] += 1
    else:
        copy_num_dict[parent.id] = 0


def add_attributes(parent, copy_num, args):
    parent.score = "."
    parent.attributes["extra_copy_number"] = str(copy_num)
    parent.attributes["copy_num_ID"] = [parent.id + "_" + str(copy_num)]
    if float(parent.attributes["coverage"][0]) < args.a:
        parent.attributes["partial_mapping"] = ["True"]
    if float(parent.attributes["sequence_ID"][0]) < args.s:
        parent.attributes["low_identity"] = ["True"]


def build_parent_dict(child_features, parent_dict, final_parent):
    parent_child_dict = {}
    for child in child_features:
        if child.id not in parent_dict:
            child.attributes["extra_copy_number"] = final_parent.attributes["extra_copy_number"][0]
            if child.attributes["Parent"][0] in parent_child_dict:
                parent_child_dict[child.attributes["Parent"][0]].append(child)
            else:
                parent_child_dict[child.attributes["Parent"][0]] = [child]
    return parent_child_dict


def write_feature(children, outfile, child_features, parent_dict):
    for child in children:
        write_line(child, outfile)
        if child.id in parent_dict:
            new_children = parent_dict[child.id]
            write_feature(new_children, outfile, child_features, parent_dict)
    return


def write_line(feature, out_file):
    line = make_gff_line(feature)
    if out_file == "stdout":
        print(line)
    else:
        out_file.write(line)
        out_file.write("\n")


def make_gff_line(feature):
    attributes_str = ""
    for attr in feature.attributes:
        if attr != "copy_id":
            value_str = ""
            for value in feature.attributes[attr]:
                value_str += value + ","
            attributes_str += (attr + "=" + value_str[:-1] + ";")
    return feature.seqid + "\t" + feature.source + "\t" + feature.featuretype + "\t" + str(feature.start) + \
           "\t" + str(feature.end) + "\t" + "." + "\t" + feature.strand + "\t" + "." + "\t" + attributes_str[:-1]
