# Liftoff
## Overview
Liftoff is tool to lift over gff3/gtf annotations from one genome assembly to another of the same species or closely related species. Most lift over tools lift each interval in a gff file over seperately which often results in gene mappings that do not make sense biologically. Liftoff instead is designed to accurately map full gene sequences onto new assemblies. It can map genes onto chromsomal regions that have been rearranged in the new assemblies, as well as resolve the correct mappings of homologous genes. 

### Getting Started

#### Step 1:
Liftoff requires Minimap2 which can be installed by following instructions [here](https://github.com/lh3/minimap2/releases/tag/v2.17). It can also be installed with conda with the following command

```
conda install -c bioconda minimap2
```
Make sure minimap2 is in your path after installation 


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
usage: liftoff.py [-h] -t target fasta -r reference fasta [-g gff or gtf]
                  [-chroms chromosome names] [-p num processess] [-o outfile]
                  [-db feature database] [-infer_transcripts]

Lift features from one genome assembly to another

optional arguments:
  -h, --help            show this help message and exit
  -t target fasta       target fasta genome to lift genes to
  -r reference fasta    reference fasta genome to lift genes from
  -g gff or gtf         annotation file to lift over in gff or gtf format
  -chroms chromosome names
                        comma seperated file with corresponding chromosomes in
                        the reference,target sequences
  -p num processess     processes
  -o outfile            output file
  -db feature database  name of feature database. If none, -g argument must be
                        provided and a database will be built automatically
  -infer_transcripts    use if GTF file only includes exon/CDS features
 
```
The only required inputs are the reference genome sequence(fasta format), the target genome sequence(fasta format) and the reference annotation or feature database. If an annotation file is provided with the -g argument, a feature database will be built automatically and can be used for future lift overs by providing the -db argument. For chromosome-scale assemblies of the same species, performing the lift over chromosome by chromosome is much faster. This option can be enabled by providing a  comma seperated file chroms.txt with corresponding chromosome names with the -chroms argument. Each line of the file should follow {ref_chrom_name},{target_chrom_name} for each pair of corresponding chromosomes. For example, a lift over from a Genbank human assembly to a Refseq human assembly would have the following chroms.txt file. 
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
After the chromosome by chromosome lift over is complete, any genes that did not map will be aligned agaisnt the whole genome. 
