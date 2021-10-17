rule mafft_align_genomes:
    input:
        rules.download_ref_genome_gff.output.fasta,
        expand(rules.gen_consensus_fasta.output, sample = samples)
    output:
        "results/mafft_alignment.aln"
    conda:
        "../envs/augur.yaml"
    params:
    	alignment_method = "mafft",
    	tmp_output = "mafft_alignment"
    threads: 8
    shell:
        """
        augur align --sequences {input} --output {params.tmp_output} --nthreads {threads} --method {params.alignment_method}
        mv {params.tmp_output} {output}
        """

rule augur_treeBuilding:
    input:
        rules.mafft_align_genomes.output
    output:
        "results/AUGUR/raw_augur_tree.treefile"
    conda:
        "../envs/augur.yaml"
    params:
        outdir = "results/AUGUR",
        n_init_tree = 2,
        n_iteration = 2,
        loglike_epsilon = 0.05,
        prefix = "results/AUGUR/raw_augur_tree",
        substitution_model = "TESTNEW"
    log:
        "results/log/iqtree.log"
    threads: 8
    shell:
        """
        iqtree -ninit {params.n_init_tree} -n {params.n_iteration} \
         -me {params.loglike_epsilon} -nt {threads} -s {input} \
         -m {params.substitution_model} -pre {params.prefix} -redo > {log}
    	"""
