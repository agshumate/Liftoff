# Liftoff
## Overview
Liftoff is a [Snakemake](https://snakemake.readthedocs.io/en/stable/) workflow to lift over gff3/gtf annotations from one chromosome-level genome assembly to another. Most lift over tools lift each interval in a gff file over seperately which often results in gene mappings that do not make sense biologically. Liftoff instead is designed to accurately map full gene sequences onto new assemblies. It can map genes onto chromsomal regions that have been rearranged in the new assemblies, as well resolve the correct mappings of homologous genes.  
