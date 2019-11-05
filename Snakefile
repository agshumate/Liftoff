configfile: "config.yaml"



NEW_FASTA=config["new_fasta"]
ORIGINAL_FASTA=config["original_fasta"]
GFF=config["gff"]
WORD_SIZE=config["word_size"]
ORIGINAL_CHROMS=open(config["original_chroms"],'r').read().splitlines()
NEW_CHROMS=open(config["new_chroms"], 'r').read().splitlines()
OUTPUT_FORMAT=config["output_format"]


rule all:
	input:
		expand('{original_chrom}_to_{new_chrom}/{new_chrom}_genomeB.' +OUTPUT_FORMAT, zip,original_chrom=ORIGINAL_CHROMS, new_chrom=NEW_CHROMS), 
		expand('{original_chrom}_to_{new_chrom}/{new_chrom}_genomeB_failed', zip ,original_chrom=ORIGINAL_CHROMS, new_chrom=NEW_CHROMS)
		


rule convert_gtf:
	input:
		"{original_chrom}_to_{new_chrom}/{original_chrom}_genomeA."+GFF[:-3]
	output: 
		temp("{original_chrom}_to_{new_chrom}/{original_chrom}_genomeA.gff")
	conda:
		'envs/gffread.yaml'
	
	run:
		if GFF[-3:] == "gtf":
			shell("gffread -O {input} -o {output}")

	
				

rule parse_gff:
	input:
		GFF
	
	output:
		"{original_chrom}_to_{new_chrom}/{original_chrom}_genomeA."+GFF[:-3]
	
	run:
		f=open(input[0],'r')
		lines = f.readlines()
		of=open(output[0], 'w')
		for line in lines:
			values = line.split("\t")
			if values[0] == wildcards.original_chrom:
				of.write(line)
					

rule convert_fasta_to_single_line_fasta:
	input:
		new=NEW_FASTA,
		old=ORIGINAL_FASTA
	output:
		new=temp('single_line_'+ NEW_FASTA),
		old=temp('single_line_'+ ORIGINAL_FASTA)
	shell:
		"""
		awk '/^>/ {{printf("\\n%s\\n",$0);next; }} {{ printf("%s",$0);}}  END {{printf("\\n");}}' < {input.new} > {output.new}
		awk '/^>/ {{printf("\\n%s\\n",$0);next; }} {{ printf("%s",$0);}}  END {{printf("\\n");}}' < {input.old} > {output.old}
		"""

rule split_chrms:
	input:
		old='single_line_' + ORIGINAL_FASTA,
		new='single_line_' + NEW_FASTA

		
	output:
		old=temp('{original_chrom}_to_{new_chrom}/{original_chrom}_genomeA.fa'),
		new=temp('{original_chrom}_to_{new_chrom}/{new_chrom}_genomeB.fa')
	run:
		for i in range(len(input)):	
			lines = open(input[i],'r').readlines()
			outfile=open(output[i],'w')
			for j in range (1,len(lines)):
				line = lines[j].strip()	
				if line == ">" + wildcards[i] or line.split()[0] == ">" + wildcards[i]:
					outfile.write(line + "\n")
					outfile.write(lines[j+1])
			

rule build_blast_dbs:
	input:
		'{original_chrom}_to_{new_chrom}/{new_chrom}_genomeB.fa'
	params:
		prefix='{original_chrom}_to_{new_chrom}/{new_chrom}_genomeB.blastdb'
	output:
		temp('{original_chrom}_to_{new_chrom}/{new_chrom}_genomeB.blastdb.nhr'),	
		temp('{original_chrom}_to_{new_chrom}/{new_chrom}_genomeB.blastdb.nin'),
		temp('{original_chrom}_to_{new_chrom}/{new_chrom}_genomeB.blastdb.nsq')
	
	conda:
		'envs/blast.yaml'

	shell:
		 'makeblastdb -in {input} -dbtype nucl -out {params.prefix}'

rule create_bedfiles:	
	input:
		'{original_chrom}_to_{new_chrom}/{original_chrom}_genomeA.gff'
	output:	
		temp('{original_chrom}_to_{new_chrom}/{original_chrom}_genomeA_genes.bed')
	shell:
		"""awk '{{if($3 == "gene"){{print}}}}' {input} |  awk -F "ID=" '{{print $1"\t"$2}}' | awk -F ";" '{{print $1}}' | awk '{{print $0"\t"$5-$4}}' | sort -nrk10,10| awk '{{print $1"\\t"$4-1"\\t"$5"\\t"$9"\\t"$6"\\t"$7}}' > {output}"""

rule extract_gene_seqs:
	input:
		bed_file='{original_chrom}_to_{new_chrom}/{original_chrom}_genomeA_genes.bed',
		fasta = '{original_chrom}_to_{new_chrom}/{original_chrom}_genomeA.fa'	
	output:
		temp('{original_chrom}_to_{new_chrom}/{original_chrom}_genomeA_genes.fa'),
		temp('{original_chrom}_to_{new_chrom}/{original_chrom}_genomeA.fa.fai')
	conda:
		'envs/bedtools.yaml'
	shell:
		'bedtools getfasta -fi {input.fasta} -bed {input.bed_file} -name -s > {output[0]}'

rule run_blast:
	input: 
		gff='{original_chrom}_to_{new_chrom}/{original_chrom}_genomeA_genes.fa',
		db1='{original_chrom}_to_{new_chrom}/{new_chrom}_genomeB.blastdb.nhr',
		db2='{original_chrom}_to_{new_chrom}/{new_chrom}_genomeB.blastdb.nin',
		db3='{original_chrom}_to_{new_chrom}/{new_chrom}_genomeB.blastdb.nsq'
		
	params:
		blastdb='{original_chrom}_to_{new_chrom}/{new_chrom}_genomeB.blastdb'
	output:
		'{original_chrom}_to_{new_chrom}/{original_chrom}_genomeA_to_{new_chrom}_genomeB.xml'

	shell:
		'blastn -out {output} -outfmt 5 -query {input.gff} -db {params.blastdb} -word_size {WORD_SIZE} -soft_masking false -culling_limit 5 -dust no -gapopen 3 -gapextend 1 -xdrop_gap_final 300'

rule lift_genes:
	input: 
		'{original_chrom}_to_{new_chrom}/{original_chrom}_genomeA_to_{new_chrom}_genomeB.xml', 
		'{original_chrom}_to_{new_chrom}/{original_chrom}_genomeA.gff'
	output:
		'{original_chrom}_to_{new_chrom}/{new_chrom}_genomeB.gff',
		'{original_chrom}_to_{new_chrom}/{new_chrom}_genomeB_failed',
		'{original_chrom}_to_{new_chrom}/warnings.log'
	conda:
		'envs/liftover.yaml'
	script:
		"scripts/lift_all_genes.py"

rule convert_output:
	input:
		'{original_chrom}_to_{new_chrom}/{new_chrom}_genomeB.gff'
	output:
		'{original_chrom}_to_{new_chrom}/{new_chrom}_genomeB.gtf'
	conda:
		'envs/gffread.yaml'
	run:
		if OUTPUT_FORMAT == "gtf":
			shell("gffread {input} -T -o {output}")








		





