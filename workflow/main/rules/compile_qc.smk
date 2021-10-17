rule compile_assemblyStats:
	input:
		expand(rules.assemblyStats.output, sample = samples)
	output:
		"results/assemblyStats_overall.tsv"
	shell:
		"""
		echo " total_length number mean_length longest shortest N_count Gaps N50 N50n N70 N70n N90 N90n" | tr " " "\t" > {output}
		for i in {input}; do
		 sample=$(echo $i | sed "s/.*\///" | sed "s/\.tsv//")
		 echo $sample $(tail -1 $i) | tr " " "\t" >> {output}
		done
		"""

rule compile_samtools_coverage_qc:
	input:
		expand(rules.samtools_coverage_qc.output, sample = samples)
	output:
		"results/samcov_summary.txt"
	params:
		dest_path = "results/SAMCOV"
	shell:
		"""
		rm -rf {output}
		for i in {params.dest_path}/* ; do
		 sample=$(echo $i | sed "s/.*\///" | sed "s/\.samcov\.txt//")
		 echo $sample >> {output}
		 echo >> {output}
		 cat $i | sed "1d" >> {output}
		 echo >> {output}
		done
		"""

rule compile_stats_qc:
	input:
		cov_bam = expand(rules.bowtie2_align.output.bam3, sample = samples),
		human_bam = expand(rules.human_align.log, sample = samples)
	output:
		"results/basic_summary_stats_qc.tsv"
	conda:
		"../envs/htslib.yaml"
	shell:
		"""
		bash src/stats.sh > {output}
		"""
