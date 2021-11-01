rule convert_global_files_to_mat:
	input:
		global_phylo = rules.download_global_tree.output.global_phylo,
		global_samples = rules.download_global_tree.output.global_samples
	output:
		"results/preprocessing/global_phylo.pb"
	conda:
		"../envs/usher.yaml"




usher --tree global_phylo.nh --vcf global_samples.vcf --collapse-tree --save-mutation-annotated-tree global_phylo.pb

usher --vcf new_samples.vcf --load-mutation-annotated-tree global_phylo.pb --write-uncondensed-final-tree

