import Bio.Blast.Applications as bio
from multiprocessing import Pool
from Bio import SeqIO


def build_database(file, name):
    make_blast_db_cline = bio.NcbimakeblastdbCommandline(dbtype="nucl", title=name + "_db",
                                                     out=name + "_db",
                                                     input_file=file)
    make_blast_db_cline()
    return name + "_db"


def get_fasta_for_database(chrm_name):
    db_name = build_database(chrm_name+".fa", chrm_name)
    return db_name

def split_fastas(fasta_dict, new_chrm_list):
    for chrm in new_chrm_list:
        SeqIO.write(fasta_dict[chrm], chrm+".fa", "fasta")


def build_all_databases(new_genome, threads):
    build_database(new_genome, "all")
    missing_fasta = []
    fasta_dict = SeqIO.index(new_genome, "fasta")
    split_fastas(fasta_dict, missing_fasta)
    pool = Pool(processes=threads)
    results_names = []
    for result in pool.imap_unordered(get_fasta_for_database , missing_fasta):
        results_names.append(result)
    pool.close()
    pool.join()




