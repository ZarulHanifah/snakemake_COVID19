rule download_nextera_pe_adaptors_fasta:
	output:
		"results/NexteraPE-PE.fa"
	params:
		trimmomatic_link = "https://github.com/timflutre/trimmomatic",
		nextera_location = "trimmomatic/adapters/NexteraPE-PE.fa"
	shell:
		"""
		git clone {params.trimmomatic_link}
		cp {params.nextera_location} {output}
		rm -rf trimmomatic
		"""

rule download_primerscheme_reference:
	output:
		bed = "results/ARTIC-V3.bed"
	params:
		github_link = "https://github.com/artic-network/artic-ncov2019.git",
		bed_path = "artic-ncov2019/primer_schemes/nCoV-2019/V3/nCoV-2019.scheme.bed",
	shell:
		"""
		git clone {params.github_link}

		awk -v OFS='\t' '$5=60' < {params.bed_path} > {output.bed}
		
		rm -rf artic-ncov2019
		"""

rule download_ref_genome_gff:
	output:
		fasta = "results/Wuhan-Hu.fasta",
		gff = "results/Wuhan-Hu.gff"
	conda:
		"../envs/ncbi.yaml"
	params:
		acc = "MN908947.3"
	shell:
		"""
		ncbi-acc-download -F fasta -o {output.fasta} {params.acc}
		ncbi-acc-download -F gff3 -o {output.gff} {params.acc}
		"""

rule download_pangolin:
	output:
		directory("pangolin")
	params:
		github_link = "https://github.com/cov-lineages/pangolin.git"
	conda:
		"../envs/pangolin.yaml"
	shell:
		"""
		git clone {params.github_link}
		
		cd {output}
		
		pip install .
		"""
