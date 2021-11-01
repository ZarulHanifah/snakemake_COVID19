#!/bin/bash

echo "Sample Average_coverage Percentage_genome Human_contamination" | tr " " "\t"
for sample in $(ls ../../input_folder/ALL_FASTQ/ | grep "R1" | sed "s/_.*//"); do
	avg_cov=$(samtools depth -r "MN908947.3:1-29903" -a "results/ALL_BAM/""$sample"".primertrim.sorted.bam" | cut -f3 | python -c "import sys; total=sum(int(l) for l in sys.stdin) ; avg=round(total/29903, 2) ; print(avg)")
	min_cov=$(samtools depth -r "MN908947.3:1-29903" -a "results/ALL_BAM/""$sample"".primertrim.sorted.bam" | cut -f3 | sort -n | head -1)
	percentage_genome=$(samtools depth -r "MN908947.3:1-29903" -a "results/ALL_BAM/""$sample"".primertrim.sorted.bam" | awk '$3 >= "10" {print $0}' | wc -l | python -c "import sys; percent= int(sys.stdin.read())/29903*100 ; percent = round(percent, 2) ; percent = str(percent) + '%' ; print(percent)")
	human_contamination=$(cat results/log/human_bowtie2_align/"$sample".log | tail -1 | sed "s/ .*//")
	echo $sample $avg_cov $percentage_genome $human_contamination | tr " " "\t"
done

