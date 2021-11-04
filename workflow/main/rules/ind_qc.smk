rule assemblyStats:
	input:
		rules.gen_consensus_fasta.output
	output:
		"results/assemblyStatsOut/{sample}.tsv"
	conda:
		"../envs/assembly-stats.yaml"
	log:
		"results/log/assembly-stats/{sample}.log"
	params:
		sample = "{sample}"
	shell:
		"""
		assembly-stats -s -t {input} | cut -f 2- >> {output}
		"""

rule quast_qc:
	input:
		input_fasta = rules.gen_consensus_fasta.output,
		bam = rules.bowtie2_align.output.bam3,
		ref_fasta = rules.download_ref_genome_gff.output.fasta,
		ref_gff = rules.download_ref_genome_gff.output.gff
	output:
		directory("results/quast/{sample}")
	conda:
		"../envs/quast.yaml"
	shell:
		"""
		quast {input.input_fasta} -r {input.ref_fasta} --features {input.ref_gff} --ref-bam {input.bam} --output-dir {output}
		"""

rule samtools_coverage_qc:
	input:
		rules.bowtie2_align.output.bam3
	output:
		"results/SAMCOV/{sample}.samcov.txt"
	conda:
		"../envs/samcov.yaml"
	params:
		bins = "196"
	shell:
		"""
		samtools coverage {input} -m -w {params.bins} -o {output}
		"""

rule human_align:
	input:
		r1 = rules.trimmomatic_pe.output.r1,
		r2 = rules.trimmomatic_pe.output.r2,
		idx_dir = config["human_index_dir"]
	output:
		"results/human_bowtie2_align/{sample}.bam"
	conda:
		"../envs/ivar.yaml"
	threads: 2
	log:
		"results/log/human_bowtie2_align/{sample}.log"
	params:
		idx_prefix = "GRCh38_noalt_as"
	shell:
		"""
		bowtie2 -x {input.idx_dir}/{params.idx_prefix} -1 {input.r1} -2 {input.r2}  2> {log} |
		 samtools view -bS  - | \
		 samtools sort -O BAM -o {output}
		"""
