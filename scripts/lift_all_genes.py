import process_blast_results as pa
import read_gff as rg

def main():
    gff=rg.get_genes_to_lift(snakemake.input[1], snakemake.output[2])
    pa.process_blast(snakemake.input[0], gff, snakemake.output[0], snakemake.output[1])

if __name__ =="__main__":
    main()


