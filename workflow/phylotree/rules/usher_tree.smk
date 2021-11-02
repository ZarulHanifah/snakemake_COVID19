rule convert_global_files_to_mat:
	input:
		global_phylo = rules.download_global_tree.output.global_phylo,
		global_samples = rules.download_global_tree.output.global_samples
	output:
		"results/preprocessing/global_phylo.pb"
	conda:
		"../envs/usher.yaml"
	shell:
		"""
		usher --tree {input.global_phylo} \
			--vcf {input.global_samples} \
			--collapse-tree \
			--save-mutation-annotated-tree {output}
		"""

rule incorporate_sample_vcf_to_tree:
	input:
		input_vcf = config[""]
		global_mat = rules.convert_global_files_to_mat.output
	output:
		""
	conda:
		"../envs/usher.yaml"
	shell:
		"""
		usher --vcf {input.input_vcf} \
			--load-mutation-annotated-tree {input.global_mat} \
			--write-uncondensed-final-tree
		"""


