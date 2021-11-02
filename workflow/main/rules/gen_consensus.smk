rule gen_consensus_fasta:
	input:
		ref = rules.download_ref_genome_gff.output.fasta,
		bam = rules.bowtie2_align.output.bam3
	output:
		"results/CONSENSUS_FASTA/{sample}.consensus.fa"
	conda:
		"../envs/ivar.yaml"
	params:
		min_base_quality = "20",
		min_depth = "10",
		prefix = "results/CONSENSUS_FASTA/{sample}.consensus"
	shell:
		"""
		prefix=$(echo {output} | sed "s/\.fa//")

		samtools mpileup -A -B -Q 0 --reference {input.ref} {input.bam} | \
		 ivar consensus -m {params.min_depth} -q {params.min_base_quality} -p $prefix -n N
        sed -i "s/>Consensus_/>/" {output}
        sed -i "s/\.consensus.*//" {output}

        rm -rf 
		"""

rule ivar_call_variants:
	input:
		ref = rules.download_ref_genome_gff.output.fasta,
		gff = rules.download_ref_genome_gff.output.gff,
		bam = rules.bowtie2_align.output.bam3
	output:
		vcf = "results/variant_calls/{sample}.vcf.gz",
		tsv = "results/variant_calls/{sample}.tsv"
	conda:
		"../envs/ivar.yaml"
	threads: 8
	log:
		"results/log/ivar_call_variants/{sample}.log"
	params:
		min_depth = 0,
		min_baseq = 0,
		ploidy = 1,
		filter_expression = "%QUAL<20 || DP > 10"
	shell:
		"""
		prefix=$(echo {output.tsv} | sed "s/\.tsv//")		

		bcftools mpileup -f {input.ref} \
			--threads {threads} \
			{input.bam} 2> {log} | \
		bcftools call --ploidy {params.ploidy} \
			--threads {threads} -v -m -Oz -o {output.vcf} 2>> {log}

		bcftools index {output.vcf}

		samtools mpileup -aa -A -B -d {params.min_depth} \
			-Q {params.min_baseq} \
			--reference {input.ref} \
			{input.bam} 2>> {log} | \
		ivar variants -p $prefix  2>> {log}
		"""

rule merge_variant_calls:
	input:
		expand(rules.ivar_call_variants.output.vcf, sample = samples)
	output:
		"results/variantcalls.vcf.gz"
	conda:
		"../envs/bowtie2.yaml"
	shell:
		"""
		header=$(for i in {input} ; do basename $i | sed "s/\..*//" ; done | tr "\n" " ")
		headerfile=$(echo {output} | sed "s/\..*/\.header/")
	
		echo $header > $headerfile

		bcftools merge --use-header $headerfile -o {output} \
			{input}
		"""
