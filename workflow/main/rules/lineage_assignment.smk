rule assign_pangolin_lineage:
	input:
		input_fasta = expand(rules.gen_consensus_fasta.output, sample = samples),
		software = rules.download_pangolin.output
	output:
		concat_fasta = temp("results/consensus_sequences.fasta"),
		report = "results/pangolin_report.tsv"
	conda:
		"../envs/pangolin.yaml"
	params:
	threads: 8
	shell:
		"""
		pangolin --update

		cat {input.input_fasta} |\
		 sed "s/>Consensus_/>/" |\
		 sed "s/\.consensus.*//" > {output.concat_fasta}
		pangolin -t {threads} {output.concat_fasta} --outfile {output.report}
		"""
