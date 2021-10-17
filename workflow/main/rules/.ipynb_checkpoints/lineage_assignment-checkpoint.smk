rule assign_pangolin_lineage:
	input:
		input_fasta = expand(rules.gen_consensus_fasta.output, sample = SAMPLE),
		software = rules.download_pangolin.output
	output:
		"results/pangolin_report.tsv"
	conda:
		"../envs/pangolin.yaml"
	params:
		concat_fasta = "results/consensus_sequences.fasta"
	threads: 8
	shell:
		"""
		pangolin --update

		cat {input.input_fasta} |\
		 sed "s/>Consensus_/>/" |\
		 sed "s/\.consensus.*//" > {params.concat_fasta}
		pangolin -t {threads} {params.concat_fasta} --outfile {output}
		"""
