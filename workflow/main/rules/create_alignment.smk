rule align_multiple_genomes_to_ref:
	input:
		concat_fasta = rules.assign_pangolin_lineage.output.concat_fasta,
		ref = rules.download_ref_genome_gff.output.fasta,
	output:
		"results/genome_alignment.sam"
	conda:
		"../envs/gofasta.yaml"
	threads: 8
	log:
		"results/log/align_multiple_genomes_to_ref/log.log"
	shell:
		"""
		minimap2 -a {input.ref} \
			{input.concat_fasta} > {output} 2> {log}
		"""

rule create_multialignment:
	input:
		ref = rules.download_ref_genome_gff.output.fasta,
		sam = rules.align_multiple_genomes_to_ref.output
	output:
		"results/genome_alignment.aln.fasta"
	conda:
		"../envs/gofasta.yaml"
	threads: 8
	log:
		"results/log/create_multialignment/log.log"
	shell:
		"""
		gofasta sam toMultiAlign -o {output} \
			-r {input.ref} \
			-s {input.sam} \
			-t {threads} 2> {log}
		"""