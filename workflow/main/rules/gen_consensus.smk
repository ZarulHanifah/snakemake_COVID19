rule gen_consensus_fasta:
	input:
		ref = rules.download_ref_genome_gff.output.fasta,
		bam = rules.bowtie2_align.output.bam3
	output:
		"results/CONSENSUS_FASTA/{sample}.consensus.fa"
	conda:
		"../envs/ivar.yaml"
	params:
		min_base_quality = "20",
		min_depth = "10",
		prefix = "results/CONSENSUS_FASTA/{sample}.consensus"
	shell:
		"""
		samtools mpileup -A -B -Q 0 --reference {input.ref} {input.bam} | \
		 ivar consensus -m {params.min_depth} -q {params.min_base_quality} -p {params.prefix} -n N
        sed -i "s/>Consensus_/>/" {output}
        sed -i "s/\.consensus.*//" {output}
		"""