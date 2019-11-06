import gff_entry as gff

# Take the input file, parse it, and save it in a dictionary of genes


def match_CDS_to_exon(CDS_start, CDS_end, gene_list, log):
    parent_exon = None
    for gene in gene_list.keys():
        for tran in gene_list[gene].Transcripts:
            for exon in tran.Exons:
                if CDS_start >= exon.start and CDS_end <= exon.end:
                    parent_exon = exon
    if parent_exon is None:
        log.write("No parent Exon found for CDS at: " + str(CDS_start) + "-" + str(CDS_end) + "\n")
    else:
        return parent_exon


def add_transcript(gene_list,chrm, source, start, end, strand, desc, log):
    tran_name = (desc.split("ID=")[1]).split(";")[0]
    if "Parent=" in desc:
        parent_name = (desc.split("Parent=")[1]).split(";")[0]
    elif "geneID" in desc:
        parent_name = (desc.split("geneID=")[1]).split(";")[0]
    if parent_name not in gene_list.keys():
        gene_list[parent_name]=gff.Gene(chrm, parent_name,
                                       int(start)-1, int(end) -1, [], strand, "complete", source)
    current_transcript = gff.Transcript(tran_name, int(start) - 1, int(end) - 1, [], strand, "complete", source)
    gene_list[parent_name].Transcripts.append(current_transcript)



def add_exon(gene_list, source, start, end, strand, desc, exon_id, log):
    parent_transcript = None
    if "Parent=" in desc:
        parent_name = (desc.split("Parent=")[1]).split(";")[0]
    elif "transcript_id" in desc:
        parent_name = (desc.split("transcript_id=")[1]).split(";")[0]
    for gene in gene_list.keys():
        for tran in gene_list[gene].Transcripts:
            if tran.id == parent_name:
                parent_transcript = tran
    if parent_transcript is None:
        log.write("No parent transcript for exon: " + parent_name + "\n")
    else:
        current_exon = gff.Exon(parent_name, exon_id,
                                int(start) - 1, int(end) - 1, "complete", strand, source, [])
        parent_transcript.Exons.append(current_exon)


def get_genes_to_lift(file, warnings_file):
    # gff columns
    chrm = 0
    source = 1
    feature = 2
    strand = 6
    frame = 7
    start = 3
    end = 4
    desc = 8
    gene_list = {}
    current_gene = None
    f = open(file, 'r')
    log = open(warnings_file, 'w')
    lines = f.readlines()
    for line in lines:
        if line[0] != "#": #skip header lines
            line = line.rstrip().split("\t")
            # if line[feature] == "gene":
            #     gene_name = (line[desc].split("ID=")[1]).split(";")[0]
            #     current_gene = gff.Gene(line[chrm], gene_name,
            #                             int(line[start])-1, int(line[end])-1, [], line[strand], "complete", line[source])
            #     gene_list[current_gene.name] = current_gene
            if line[feature] == "transcript" or line[feature] == "mRNA":
                exon_id = 0
                add_transcript(gene_list,line[chrm], line[source], line[start], line[end], line[strand],  line[desc], log)
            elif line[feature] == "exon":
                exon_id += 1
                add_exon(gene_list,line[source], line[start], line[end], line[strand], line[desc], exon_id,log)
            elif line[feature] == "CDS":
                parent_exon = match_CDS_to_exon(int(line[start])-1, int(line[end])-1, gene_list, log)
                if parent_exon is not None:
                    parent_exon.CDS.append(gff.CDS(int(line[start])-1, int(line[end])-1, line[frame], "complete"))
    for gene in gene_list.keys():
        gene_list[gene].start = gene_list[gene].Transcripts[0].start
        gene_list[gene].end = gene_list[gene].Transcripts[-1].end
    return gene_list

