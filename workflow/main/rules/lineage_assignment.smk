rule assign_pangolin_lineage:
	input:
		expand(rules.gen_consensus_fasta.output, sample = samples),
	output:
		concat_fasta = "results/consensus_sequences.fasta",
		report = "results/pangolin_report.tsv"
	conda:
		"../envs/pangolin.yaml"
	threads: 8
	log:
		"results/log/assign_pangolin_lineage/log.log"
	shell:
		"""
		cat {input} |\
		 sed "s/>Consensus_/>/" |\
		 sed "s/\.consensus.*//" > {output.concat_fasta} 2> {log}
		
		pangolin --update
		pangolin -t {threads} \
			{output.concat_fasta} \
			--outfile {output.report} &> {log}
		"""
