# Liftoff

## Overview
Here we introduce Liftoff, an accurate tool that maps annotations in GFF or GTF between assemblies of the same, or closely-related species. Unlike current coordinate lift-over tools which require a pre-generated “chain” file as input, Liftoff is a standalone tool that takes two genome assemblies and a reference annotation as input and outputs an annotation of the target genome. Liftoff uses Minimap2 [(Li, 2018)](https://academic.oup.com/bioinformatics/article/34/18/3094/4994778) to align the gene sequences from a reference genome to the target genome. Rather than aligning whole genomes, aligning only the gene sequences allows genes to be lifted over even if there are many structural differences between the two genomes. For each gene, Liftoff finds the alignments of the exons that maximize sequence identity while preserving the transcript and gene structure.  If two genes incorrectly map to overlapping loci, Liftoff determines which gene is most-likely mis-mapped, and attempts to re-map it. Liftoff can also find additional gene copies present in the target assembly that are not annotated in the reference. 

Liftoff is not limited to genes, transcripts and exons. It can transfer annotations of any feature types present in a GFF file. For any feature or group of features regardless of type, Liftoff aligns the parent feature to the genome, and then converts the coordinates of the child features while preserving the structure of each child and parent. 
<p align="center">
  <img src="https://user-images.githubusercontent.com/29218752/84577010-d0e34680-ad86-11ea-89a2-1638b970dcad.jpg">
</p>

### Getting Started

#### Step 1:
Liftoff requires Python3 and also depends on Minimap2. Minimap2 which can be installed by following instructions [here](https://github.com/lh3/minimap2/releases/tag/v2.17). It can also be installed with conda with the following command

```
conda install -c bioconda minimap2
```
Add minimap2 is in your path after installation or use the -m argument when running Liftoff to specficy a different path 


#### Step 2: 
Then clone the workflow into a working directory 
```
git clone https://github.com/agshumate/Liftoff liftoff 
```

#### Step 3:

```
cd liftoff
python setup.py install
```

### USAGE
```
usage: liftoff [-h] -t <target.fasta> -r <reference.fasta>
                  [-g <ref_annotation.gff>] [-chroms <chroms.txt>] [-p 1]
                  [-o <output.gff>] [-db DB] [-infer_transcripts]
                  [-u <unmapped_features.txt>] [-infer_genes] [-a 0.5]
                  [-s 0.5] [-unplaced <unplaced_seq_names.txt>] [-copies]
                  [-sc 1.0] [-m PATH] [-dir <intermediate_files_dir>]

Lift features from one genome assembly to another

optional arguments:
  -h, --help            show this help message and exit
  -t <target.fasta>     target fasta genome to lift genes to
  -r <reference.fasta>  reference fasta genome to lift genes from
  -g <ref_annotation.gff>
                        annotation file to lift over in gff or gtf format
  -chroms <chroms.txt>  comma seperated file with corresponding chromosomes in
                        the reference,target sequences
  -p 1                  processes
  -o <output.gff>       output file
  -db DB                name of feature database. If none, -g argument must be
                        provided and a database will be built automatically
  -infer_transcripts    use if GTF file only includes exon/CDS features and
                        does not include transcripts/mRNA
  -u <unmapped_features.txt>
                        name of file to write unmapped features to
  -infer_genes          use if GTF file only includes transcripts, exon/CDS
                        features
  -a 0.5                minimum alignment coverage to consider a feature
                        mapped [0-1]
  -s 0.5                minimum sequence identity in child features (usually
                        exons/CDS) to consider a feature mapped [0-1]
  -unplaced <unplaced_seq_names.txt>
                        text file with name(s) of unplaced sequences to map
                        genes from after genes from chromosomes in chroms.txt
                        are mapped
  -copies               look for extra gene copies in the target genome
  -sc 1.0               with -copies, minimum sequence identity in exons/CDS
                        for which a gene is considered a copy. Must be greater
                        than -s
  -m PATH               Minimap2 path
  -dir <intermediate_files_dir>
                        name of directory to save intermediate fasta and SAM
                        files
 
```
The only required inputs are the reference genome sequence(fasta format), the target genome sequence(fasta format) and the reference annotation or feature database. If an annotation file is provided with the -g argument, a feature database will be built automatically and can be used for future lift overs by providing the -db argument. 

By default, a gene (or other parent feature) will only be considered mapped successfully if the alignment coverage and sequence identity in the chlild features (usually exons/CDS) is > 50%. The alignment coverage threshold can be changed with the -a option. Any genes mapping with a coverage below the threshold will be present in the GFF file with the tag partial_mapping=True. The sequence identity threshold can be changed with the -s option. Any genes mapping with a lower sequence identity will be present in the GFF file with the tag low_identity=True. 

Extra gene copies will have the same ID as the reference gene and will be tagged with extra_copy_number={copy_number}

For chromosome-scale assemblies of the same species, performing the lift over chromosome by chromosome is much faster. This is also strongly recommended for repetitive/polyploid genomes. This option can be enabled by providing a  comma seperated file chroms.txt with corresponding chromosome names with the -chroms argument. Each line of the file should follow {ref_chrom_name},{target_chrom_name} for each pair of corresponding chromosomes. For example, a lift over from a Genbank human assembly to a Refseq human assembly would have the following chroms.txt file. 
 ```
chr1,NC_000001.10
chr2,NC_000002.11
chr3,NC_000003.11
chr4,NC_000004.11
chr5,NC_000005.9
chr6,NC_000006.11
chr7,NC_000007.13
chr8,NC_000008.10
chr9,NC_000009.11
chr10,NC_000010.10
chr11,NC_000011.9
chr12,NC_000012.11
chr13,NC_000013.10
chr14,NC_000014.8
chr15,NC_000015.9
chr16,NC_000016.9
chr17,NC_000017.10
chr18,NC_000018.9
chr19,NC_000019.9
chr20,NC_000020.10
chr21,NC_000021.8
chr22,NC_000022.10
chrX,NC_000023.10
chrY,NC_000024.9
```
After the chromosome by chromosome lift over is complete, any genes that did not map will be aligned to the whole genome. 




