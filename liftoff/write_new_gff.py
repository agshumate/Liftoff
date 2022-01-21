from liftoff import liftoff_utils, __version__
import sys

def write_header(f, out_type):
    f.write("# " + " ".join(sys.argv) + "\n")
    f.write("# Liftoff v" + __version__ + "\n")
    if out_type == 'gff3':
        f.write('##gff-version 3' + "\n")


def write_new_gff(lifted_features, args, feature_db):
    if args.o != 'stdout':
        f = open(args.o, 'w')
    else:
        f = "stdout"
    out_type = feature_db.dialect['fmt']
    write_header(f, out_type)
    parents = liftoff_utils.get_parent_list(lifted_features)
    parents.sort(key=lambda x: x.id)
    final_parent_list = finalize_parent_features(parents, args)
    final_parent_list.sort(key=lambda x: (x.seqid, x.start))
    for final_parent in final_parent_list:
        child_features = lifted_features[final_parent.attributes["copy_id"][0]]
        parent_child_dict = build_parent_dict(child_features, final_parent)
        write_feature([final_parent], f, child_features, parent_child_dict, out_type)


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
    keys_to_readd_at_end = ["copy_num_ID", "partial_mapping", "low_identity", "extra_copy_number",
                            "partial_mapping",  "low_identity"]
    if "copy_id" not in parent.attributes:
        parent.attributes["copy_id"] = parent.attributes["copy_num_ID"]
    for key in keys_to_readd_at_end:
        if key in parent.attributes:
            del  parent.attributes[key]
    parent.attributes["extra_copy_number"] = [str(copy_num)]
    parent.attributes["copy_num_ID"] = [parent.id + "_" + str(copy_num)]

    if float(parent.attributes["coverage"][0]) < args.a:
        parent.attributes["partial_mapping"] = ["True"]
    if float(parent.attributes["sequence_ID"][0]) < args.s:
        parent.attributes["low_identity"] = ["True"]


def build_parent_dict(child_features, final_parent):
    parent_child_dict = {}
    for child in child_features:
        if "Parent" in child.attributes:
            child.attributes["extra_copy_number"] = final_parent.attributes["extra_copy_number"]
            if child.attributes["Parent"][0] in parent_child_dict:
                parent_child_dict[child.attributes["Parent"][0]].append(child)
            else:
                parent_child_dict[child.attributes["Parent"][0]] = [child]
    return parent_child_dict


def write_feature(children, outfile, child_features, parent_dict, output_type):
    for child in children:
        write_line(child, outfile, output_type)
        if child.id in parent_dict:
            new_children = parent_dict[child.id]
            write_feature(new_children, outfile, child_features, parent_dict, output_type)
    return


def write_line(feature, out_file, output_type):
    if output_type == 'gff3':
        line = make_gff_line(feature)
    else:
        line = make_gtf_line(feature)
    if out_file == "stdout":
        print(line)
    else:
        out_file.write(line)
        out_file.write("\n")


def make_gff_line(feature):
    edit_copy_ids(feature)
    attributes_str = "ID=" + feature.attributes["ID"][0] + ";" #make ID the first printed attribute
    for attr in feature.attributes:
        if attr != "copy_id":
            value_str = ""
            for value in feature.attributes[attr]:
                    value_str += value + ","
            if attr != "ID":
                attributes_str += (attr + "=" + value_str[:-1] + ";")
    return feature.seqid + "\t" + feature.source + "\t" + feature.featuretype + "\t" + str(feature.start) + \
           "\t" + str(feature.end) + "\t" + "." + "\t" + feature.strand + "\t" + "." + "\t" + attributes_str[:-1]

def edit_copy_ids(feature):
    copy_num = feature.attributes["extra_copy_number"][0]
    if copy_num != '0':
        feature.attributes["ID"] = [feature.attributes["ID"][0]+ "_" + copy_num]
        if "Parent" in feature.attributes:
            feature.attributes["Parent"] = [feature.attributes["Parent"][0] + "_" + copy_num]
        if "gene_id" in feature.attributes:
            feature.attributes["gene_id"] = [feature.attributes["gene_id"][0] + "_" + copy_num]


def make_gtf_line(feature):
    attributes_str = ""
    for attr in feature.attributes:
        if attr != "copy_id":
            if len(feature.attributes[attr]) >0:
                value_str = ""
                for value in feature.attributes[attr]:
                     value_str += value + ","
                attributes_str += (attr + " " + '"' + value_str[:-1] + '"' + "; ")
    return feature.seqid + "\t" + feature.source + "\t" + feature.featuretype + "\t" + str(feature.start) + \
           "\t" + str(feature.end) + "\t" + "." + "\t" + feature.strand + "\t" + "." + "\t" + attributes_str