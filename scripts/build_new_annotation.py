import numpy as np

# Write the new annotation in gff format
def print_gff(aligned, f, failed_file):
    aligned.sort(key=lambda z: (z.chrm, int(z.start)))
    for gene in aligned:
        chrm = gene.chrm
        if gene.start == -1 and gene.end == -1:
            failed_file.write(gene.name + "\n")
        f.write(chrm + "\t" + gene.source + "\t" + "gene" + "\t" + str(gene.start + 1) + "\t" + str(gene.end + 1)
                + "\t" + "." + "\t" + gene.strand + "\t" + "." + "\t" + "ID="+gene.name + ";" + gene.status + "\n")
        for tran in gene.Transcripts:
            CDS_list = []
            if tran.start != -1 and tran.end != -1:
                f.write(chrm + "\t" + tran.source + "\t" + "transcript" + "\t" + str(tran.start + 1) + "\t" +
                        str(tran.end + 1) + "\t" + "." + "\t" + gene.strand + "\t" + "." + "\t" + "ID="+tran.id + ";" +
                        "Parent=" + gene.name + ";" + tran.status + "\n")
                for exon in tran.Exons:
                    if exon.start != -1 and exon.end != -1 :
                        f.write(chrm + "\t" + exon.source + "\t" + "exon" + "\t" + str(exon.start +1) + "\t" +
                                str(exon.end +1) + "\t" + "." + "\t" + gene.strand + "\t" + "." + "\t" +
                                "Parent=" + tran.id+";" + "exon_id=" + str(exon.id) + ";" + exon.status + "\n")
                        for CDS in exon.CDS:
                            CDS_list.append(CDS)
                for CDS in CDS_list:
                    if CDS.start != -1 and CDS.end != -1 :
                        f.write(chrm + "\t" + exon.source + "\t" + "CDS" + "\t" + str(CDS.start + 1) + "\t"
                                + str(CDS.end +1 ) + "\t" + "." + "\t" + gene.strand + "\t" + CDS.frame +
                                "\t" + "Parent="+tran.id+";" + CDS.status + "\n")

# Aggregate the mappings of each exon into transcripts and genes and label them as complete or incomplete
def add_gene_and_tran_coords(gene, original_gene):
    gene_status = "complete"
    if gene.strand != 'unaligned':
        CDS_statuses = np.array([CDS.status for tran in gene.Transcripts for exon in tran.Exons for CDS in exon.CDS])
        for tran in gene.Transcripts:
            tran_status = "complete"
            tran.Exons.sort(key=lambda x: int(x.start))
            exon_starts = np.array([exon.start for exon in tran.Exons])
            exon_ends = np.array([exon.end for exon in tran.Exons])
            exon_statuses = np.array([exon.status for exon in tran.Exons])
            if len(CDS_statuses) == 0:
                tran_features = exon_statuses
            else:
                tran_features = CDS_statuses
            exon_order = [exon.id for exon in tran.Exons]
            if exon_order != sorted(exon_order) and exon_order != sorted(exon_order, reverse=True):
                tran_status = "incomplete"
                gene_status = "incomplete"
            if "unaligned" in tran_features or "truncated" in tran_features:
                tran_status = "incomplete"
                gene_status = "incomplete"
            valid_exon_starts = exon_starts[exon_starts != -1]
            valid_exon_ends = exon_ends[exon_ends != -1]
            if len(valid_exon_starts) > 0:
                tran.start = np.min(valid_exon_starts)
            else:
                tran.start = -1
            if len(valid_exon_ends) > 0:
                tran.end = np.max(valid_exon_ends)
            else:
                tran.end = -1
            tran.status = tran_status
        gene.Transcripts.sort(key=lambda y: int(y.start))
        tran_starts = np.array([tran.start for tran in gene.Transcripts])
        tran_ends = np.array([tran.end for tran in gene.Transcripts])
        valid_tran_starts = tran_starts[tran_starts != -1]
        valid_tran_ends = tran_ends[tran_ends != -1]
        if len(valid_tran_starts) > 0:
            gene.start = np.min(valid_tran_starts)
        else:
            gene.start = -1
        if len(valid_tran_ends) > 0:
            gene.end = np.max(valid_tran_ends)
        else:
            gene.end = -1
        if gene_status == "complete" and gene.strand != original_gene.strand:
            gene_status = "complete_on_wrong_strand"
        gene.status = gene_status




