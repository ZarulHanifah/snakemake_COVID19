rule bowtiebuild_buildIndex:
	input:
		rules.download_ref_genome_gff.output.fasta
	output:
		directory("results/INDEX")
	conda:
		"../envs/ivar.yaml"
	log:
		"results/log/bowtiebuild_buildIndex/bowtiebuild_buildIndex.log"
	shell:
		"""
		mkdir {output}
		bowtie2-build {input} {output}/idx > {log}
		samtools faidx {input}
		"""

rule bowtie2_align:
	input:
		r1 = rules.trimmomatic_pe.output.r1,
		r2 = rules.trimmomatic_pe.output.r2,
		index = rules.bowtiebuild_buildIndex.output,
		primers_bed = rules.download_primerscheme_reference.output.bed
	output:
		bam1 = "results/ALL_BAM/{sample}.bam",
		bam2 = "results/ALL_BAM/{sample}.primertrim.bam",
		bam3 = "results/ALL_BAM/{sample}.primertrim.sorted.bam"
	conda:
		"../envs/ivar.yaml"
	threads: 8
	log:
		"results/log/bowtie2_align/{sample}.log"
	params:
		prefix = "results/ALL_BAM/{sample}.primertrim"
	shell:
		"""
		bowtie2 -x {input.index}/idx \
			-1 {input.r1} \
			-2 {input.r2} \
			--rg-id {wildcards.sample} \
			-p {threads} \
			--very-sensitive-local 2> {log} | \
		 samtools sort | \
		 samtools view -bS -F 4 -o {output.bam1} 2>> {log}
		
		ivar trim -i {output.bam1} \
			-b {input.primers_bed} \
			-p {params.prefix} &>> {log}

		samtools sort {output.bam2} -o {output.bam3}
		samtools index {output.bam3}
		"""