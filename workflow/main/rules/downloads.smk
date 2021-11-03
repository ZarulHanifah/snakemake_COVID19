rule download_adaptors_fasta:
	output:
		temp("results/adapters.fa")
	params:
		link = "https://github.com/BioInfoTools/BBMap.git",
		adapter_file = "BBMap/resources/adapters.fa"
	shell:
		"""
		git clone {params.link}
		cp {params.adapter_file} {output}
		rm -rf $(dirname $(dirname {params.adapter_file}))
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

rule index_fasta:
	input:
		rules.download_ref_genome_gff.output.fasta
	output:
		"results/Wuhan-Hu.fasta.fai"
	conda:
		"../envs/bowtie2.yaml"
	shell:
		"""
		samtools index {input}
		"""