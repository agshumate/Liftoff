import Bio.Blast.Applications as bio
from multiprocessing import Pool
from functools import partial
from Bio.Blast import NCBIXML
import math
import copy


def blast_genes_for_chrm(old_chroms, new_chroms, blast_threads, leftover_threads, word_size, search_type, index):
    genes_file = old_chroms[index].split()[0] + "_genes.fa"
    blastdb = new_chroms[index] + "_db"
    if index < leftover_threads:
        blast_threads += 1
    outfile = old_chroms[index].split()[0] + "_to_" + new_chroms[index].split()[0] + ".xml"
    if search_type == "chrm_by_chrm":
        blastn_cline = bio.NcbiblastnCommandline(query=genes_file, db=blastdb, outfmt=5, soft_masking=False, dust='no',

                                                 word_size=word_size, gapopen=3, gapextend=1, xdrop_gap_final=300,
                                                 num_threads=blast_threads, culling_limit =10, out=outfile)
    elif search_type == "missing":
        blastn_cline = bio.NcbiblastnCommandline(query=genes_file, db=blastdb, outfmt=5, soft_masking=False,
                                                 dust='no', gapopen=3, gapextend=1, xdrop_gap_final=300,
                                                 num_threads=blast_threads, culling_limit = 100,
                                                 out=outfile)
    elif search_type == "copies":
        blastn_cline = bio.NcbiblastnCommandline(query=genes_file, db=blastdb, outfmt=5, soft_masking=False,
                                                 dust='no',
                                                 word_size=word_size, gapopen=3, gapextend=1, xdrop_gap_final=300,
                                                 num_threads=blast_threads,
                                                 out=outfile, qcov_hsp_perc=100, culling_limit=100)
    elif search_type == "unplaced":
        blastn_cline = bio.NcbiblastnCommandline(query=genes_file, db=blastdb, outfmt=5, soft_masking=False,
                                                 dust='no',
                                                 word_size=word_size, gapopen=3, gapextend=1, xdrop_gap_final=300,
                                                 num_threads=blast_threads,
                                                 out=outfile, culling_limit = 10)
    blastn_cline()
    return (outfile)



def parse_blast(file):
    all_records = {}
    for blast_record in NCBIXML.parse(open(file, 'r')):
        query_name = blast_record.query
        all_records[query_name] = blast_record
    return all_records


def count_num_records(record):
    count = 0
    for aln in record.alignments:
        count += len(aln.hsps)
    return count


def get_new_record(record, aln_num):
    new_record = copy.deepcopy(record)
    hsp_count = -1
    for aln in new_record.alignments[:]:
        for hsp in aln.hsps[:]:
            hsp_count += 1
            if hsp_count != aln_num:
                aln.hsps.remove(hsp)
    return new_record


def expand_results(all_records):
    new_record_dict = {}
    for record in all_records:
        num_records = count_num_records(all_records[record])
        for i in range(num_records):
            record_id = str(i +2)
            new_record_name = record + "_"+ record_id
            new_record = get_new_record(all_records[record], i)
            new_record.query = new_record_name
            new_record_dict[new_record_name] = new_record
    return new_record_dict


def rename_records(all_records):
    renamed_records = {}
    for record in all_records:
        new_name = record + "_0"
        renamed_records[new_name] = all_records[record]
        renamed_records[new_name].query = new_name
    return renamed_records


def blast_all_genes(old_chroms, new_chroms, threads, word_size, search_type):
    pool = Pool(processes=threads)
    blast_threads = max(1, math.floor(threads / len(old_chroms)))
    leftover_threads = threads % len(old_chroms)
    all_records = {}
    all_files = []
    func = partial(blast_genes_for_chrm, old_chroms, new_chroms, blast_threads, leftover_threads, word_size, search_type)
    for result in pool.imap_unordered(func, range(0, len(old_chroms))):
        all_files.append(result)
    func = partial(parse_blast)
    for result in pool.imap_unordered(func, all_files):
        all_records.update(result)
    pool.close()
    pool.join()
    if search_type == "copies":
        all_records = expand_results(all_records)
    else:
        all_records = rename_records(all_records)
    return all_records
