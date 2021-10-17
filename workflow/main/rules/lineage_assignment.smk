rule assign_pangolin_lineage:
	input:
		expand(rules.gen_consensus_fasta.output, sample = samples),
	output:
		concat_fasta = temp("results/consensus_sequences.fasta"),
		report = "results/pangolin_report.tsv"
	conda:
		"../envs/pangolin.yaml"
	threads: 8
	shell:
		"""
		pangolin --update

		cat {input} |\
		 sed "s/>Consensus_/>/" |\
		 sed "s/\.consensus.*//" > {output.concat_fasta}
		pangolin -t {threads} {output.concat_fasta} --outfile {output.report}
		"""
