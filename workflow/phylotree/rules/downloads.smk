rule download_global_tree:
	output:
		global_phylo = "results/downloads/global_phylo.nh",
		global_samples = "results/downloads/global_samples.vcf"
	params:
		global_phylo_link = "https://raw.githubusercontent.com/yatisht/usher/master/test/global_phylo.nh",
		global_samples_link = "https://raw.githubusercontent.com/yatisht/usher/master/test/global_samples.vcf"
	shell:
		"""
		wget -c {param.global_phylo_link} -O {output.global_phylo}
		wget -c {param.global_samples_link} -O {output.global_samples}
		"""