def write_line(line, out_file):
    if out_file == "stdout":
        print(line)
    else:
        out_file.write(line)
        out_file.write("\n")

def rename_features(lifted_features):
    group_feature_dict = {}
    renamed_features = {}
    for feature in lifted_features:
        original_name = feature[:-len(feature.split("_")[-1]) -1]
        if original_name in group_feature_dict:
            group_feature_dict[original_name].append(lifted_features[feature])
        else:
            group_feature_dict[original_name] = [lifted_features[feature]]
    for gene_name in group_feature_dict:
        for i in range(len(group_feature_dict[gene_name])):
            if i == 0:
                copy_num_tag = ""
            else:
                copy_num_tag = "_" + str(i +1)
            new_gene_name = gene_name + copy_num_tag
            genes = [feature for feature in group_feature_dict[gene_name][i] if feature.featuretype == "gene"]
            new_feature_group = []
            for gene in genes:
                transcripts = [feature for feature in group_feature_dict[gene_name][i] if feature.featuretype!="gene" and feature.attributes["Parent"][0] == gene_name]
                for transcript in transcripts:
                    transcript_num = transcript.id.split(".")[-1]
                    new_tran_name = transcript.id[:-len(transcript.id.split(".")[-1])-1] + copy_num_tag + "." + transcript_num
                    exons = [feature for feature in lifted_features[gene.id] if
                             feature.featuretype == "exon" and feature["Parent"][0] == transcript.id]
                    cdss = [feature for feature in lifted_features[gene.id] if
                            feature.featuretype == "CDS" and feature["Parent"][0] == transcript.id]
                    for exon in exons:
                        exon.attributes["Parent"] = new_tran_name
                        new_feature_group.append(exon)
                    for cds in cdss:
                        cds.attributes["Parent"] = new_tran_name
                        new_feature_group.append(cds)
                    transcript.id = new_tran_name
                    transcript.attributes["ID"] = new_tran_name
                    transcript.attributes["Parent"] = new_gene_name
                    transcript.attributes["Name"] = new_gene_name
                    new_feature_group.append(transcript)
                gene.id = new_gene_name
                gene.attributes["ID"] = new_gene_name
                gene.attributes["Name"] = new_gene_name
                new_feature_group.append(gene)
            renamed_features[new_gene_name] = new_feature_group
    return renamed_features




def write_new_gff(lifted_features, out_file):
    lifted_features = rename_features(lifted_features)
    if out_file != 'stdout':
        f=open(out_file, 'w')
    else:
        f="stdout"
    all_features= []
    for feature in lifted_features:
        all_features.extend(lifted_features[feature])
    genes = [feature for feature in all_features if feature.featuretype == "gene"]
    genes.sort(key=lambda x: (x.seqid, x.start))
    for gene in genes:
        gene.score = "."
        write_line(str(gene), f)
        transcripts = [feature for feature in lifted_features[gene.id] if feature.featuretype!="gene" and feature.attributes["Parent"][0] == gene.id]
        transcripts.sort(key=lambda x: x.start)
        for transcript in transcripts:
            write_line(str(transcript), f)
            exons = [feature for feature in lifted_features[gene.id] if feature.featuretype == "exon" and feature["Parent"][0] == transcript.id]
            exons.sort(key=lambda x: x.start)
            for exon in exons:
                write_line(str(exon), f)
            cdss = [feature for feature in lifted_features[gene.id] if feature.featuretype == "CDS" and feature["Parent"][0] == transcript.id]
            cdss.sort(key=lambda x: x.start)
            for cds in cdss:
                write_line(str(cds), f)
