from read_gff import get_genes_to_lift
from build_new_annotation import print_gff

gff_file = get_genes_to_lift(snakemake.input[0], snakemake.output[1])
print_gff(list(gff_file.values()), open(snakemake.output[0],'w'),  open(snakemake.output[2],'w'))
