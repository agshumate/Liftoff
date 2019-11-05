# Liftoff
## Overview
Liftoff is a [Snakemake](https://snakemake.readthedocs.io/en/stable/) workflow to lift over gff3/gtf annotations from one chromosome-level genome assembly to another. Most lift over tools lift each interval in a gff file over seperately which often results in gene mappings that do not make sense biologically. Liftoff instead is designed to accurately map full gene sequences onto new assemblies. It can map genes onto chromsomal regions that have been rearranged in the new assemblies, as well as resolve the correct mappings of homologous genes.  

### Getting Started

#### Step 1:
Liftoff uses [Snakemake](https://snakemake.readthedocs.io/en/stable/) to manage the workflow and handle all dependencies. The first step to running Liftoff is therefore installing Snakemake through the Python3 version of Miniconda

```
conda install -c bioconda -c conda-forge snakemake
```

More information about installing Snakemake can be found [here](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html)

#### Step 2: 
Then clone the workflow into a working directory 
```
git clone https://github.com/agshumate/Liftoff liftoff 
```

#### Step 3:
Add the original genome fasta file, the gff3/gtf file, and the new genome fasta to this directory 
```
mv my_original_genome.fa liftoff/
mv my_original_genome.gff liftoff/
mv my_new_genome.fa  liftoff/
```

#### Step 4:
Edit 'original_chroms.txt' and 'new_chroms.txt'.  Liftoff performs the gene lift over chromosome by chromosome so it is important to tell Liftoff the names of the chromosomes in the original genome assembly and the names of the corresponding chromosomes in the new assembly. If the chromosomes have the same names in both assemblies, these files will be identical.

#### Step 5:
Edit the config.yaml file with the following information

**gff**: the gtf/gff3 file that will be lifted over to the new assembly  
**output_format**: either 'gtf' of 'gff'  
**original_chroms**: the 'original_chroms.txt' file edited to include the correct chromosome names in the original assembly  
**new_chroms**: the 'new_chroms.txt' file edited to include the correct chromosome names in the original assembly  
**original_fasta**: the genome assembly (fasta format) that the genes will be lifted from   
**new_fasta**: the genome assembly (fasta format) that the genes will be lifted to   
**word_size**: the size of exact match the required by BLAST to align the gene sequence to the chromosome. 50 is recommended for human. A smaller word size will make the BLAST alignment step slower but more sensitive . 

An example config.yaml file for lifting genes over from the human genome assembly version GRCh37 to GRCh38 would be 
```
gff: GRCh37.gff
output_format: gff
original_chroms: original_chroms.txt
new_chroms: new_chroms.txt
original_fasta: GRCh37.fa
new_fasta: GRCh38.fa
word_size: 50
```
### Running Liftoff
Once config.yaml, original_chroms.txt and new_chroms.txt have been edited run the Snakemake workflow with
```
snakemake --use-conda
``` 
Use the -j option to use multiple cores. For example, if you have at least a 24 core machine, you could lift over all human chromsomes in parallel by executing 

```
snakemake --use-conda -j 24
```
Snakemake has many more [options](https://snakemake.readthedocs.io/en/stable/executable.html) you can invoke depending on your machine setup 

### Output
A directory called {chromA}\_to\_{chromB} will be created in the working directory for each chromosome to chromsome liftover containing   
```
{chrom}_genomeB.gff: The lifted over gff or gtf  
{chrom}_failed_genomeB: A list of genes that failed to be lifted over  
warnings.log: A list of problems found in the original gff file such as transcripts or exons without parent features  
```
