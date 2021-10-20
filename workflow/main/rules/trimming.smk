rule trimmomatic_pe:
	input:
		r1 = os.path.join("../../", config["fastq_dir_path"], "{sample}_R1_.fastq.gz"),
		r2 = os.path.join("../../", config["fastq_dir_path"], "{sample}_R2_.fastq.gz"),
		adapter = rules.download_adaptors_fasta.output
	output:
		r1 = "results/trimmed/{sample}_R1_.fastq.gz",
		r2 = "results/trimmed/{sample}_R2_.fastq.gz"
	conda:
		"../envs/trimmomatic.yaml"
	log:
		"results/log/trimmomatic/{sample}.log"
	params:
		illuminaclip = ":2:30:10",
		slidingwindow = "5:20",
		leading = "10",
		trailing = "10",
		minlen = "100"
	threads: 2
	shell:
		"""
		trimmomatic PE -threads {threads} -phred33 {input.r1} {input.r2} \
		 {output.r1} /dev/null \
		 {output.r2} /dev/null \
		 ILLUMINACLIP:{input.adapter}{params.illuminaclip} \
		 SLIDINGWINDOW:{params.slidingwindow} \
		 LEADING:{params.leading}\
		 TRAILING:{params.trailing}\
		 MINLEN:{params.minlen} &>{log}
		"""
