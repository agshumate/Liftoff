# Liftoff
## Overview
Liftoff is tool to lift over gff3/gtf annotations from one chromosome-level genome assembly to another of the same species. Most lift over tools lift each interval in a gff file over seperately which often results in gene mappings that do not make sense biologically. Liftoff instead is designed to accurately map full gene sequences onto new assemblies. It can map genes onto chromsomal regions that have been rearranged in the new assemblies, as well as resolve the correct mappings of homologous genes. 

### Getting Started

#### Step 1:
Liftoff requires commandline BLAST which can ben installed [here] (https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download). It can also be installed with conda with the following command

```
conda install -c bioconda blast
```



#### Step 2: 
Then clone the workflow into a working directory 
```
git clone https://github.com/agshumate/Liftoff liftoff 
```

#### Step 3:
```
python setup.py install
```

### USAGE
```
usage: liftoff.py [-h] -g gff_or_gtf -t target_fasta -r reference_fasta
                  [-target_chroms target_chroms.txt]
                  [-ref_chroms reference_chroms.txt] [-p num_processess]
                  [-o out_file] [-w blast_word_size]
                  [-unplaced unplaced_sequence] [-copy_num]

Lift genes

optional arguments:
  -h, --help            show this help message and exit
  -g gff_or_gtf         annotation file to lift over in gff or gtf format
  -t target_fasta       target fasta genome to lift genes to
  -r reference_fasta    reference fasta genome to lift genes from
  -target_chroms        file with name of chromosomes to be lifted to
  -ref_chroms           file with name of chromosomes to be lifted from
  -p num_processess     processes
  -o out_file           output gff file
  -w blast_word_size    word size for blast step
  -unplaced unplaced_sequence
                        file with unplaced sequence names. Genes annotated on
                        these sequences will be mapped onto the main assembly
  -copy_num             look for additional copies of genes after the
                        annotation has been lifted over
 
```
The only required inputs are the reference genome sequence(fasta format), the target genome sequence(fasta format) and the reference annotation(gff or gtf format). For chromosome-scale assemblies, a file with a list of chromosomes in the reference genome can be provided in the file reference_chroms.txt file. The corresponding chromosomes in the target assembly should be provided in the file target_chroms.txt. If the chromosomes have the same names in both assemblies, these files will be identical. Providing these files lifts genes over chromosome by chromosome and is generally much faster then mapping all genes to the entire target assembly. 
