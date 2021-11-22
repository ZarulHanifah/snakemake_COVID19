#!/bin/bash

echo "Sample Average_coverage Percentage_genome(%) Human_contamination(%)" | tr " " "\t"

samplelist=$(ls results/CONSENSUS_FASTA/*fa | \
			grep "consensus" | \
			sed "s/.*\///" | \
			sed "s/\..*//")

for sample in $samplelist; do
	bamfile="results/ALL_BAM/""$sample"".primertrim.sorted.bam"
	fastafile="results/ALIGNED_FASTA/""$sample"".aln.fasta"
	humanlogfile="results/log/human_bowtie2_align/"$sample".log"

	avg_cov=$(samtools depth -r "MN908947.3:1-29903" -a $bamfile | \
				cut -f3 | \
				python -c "import sys
total = sum(int(l) for l in sys.stdin)
avg = round(total/29903, 2)
print(avg)")

	min_cov=$(samtools depth -r "MN908947.3:1-29903" -a $bamfile | \
				cut -f3 | \
				sort -n | \
				head -1)

	human_contamination=$(cat $humanlogfile | \
							tail -1 | \
							sed "s/% .*//")
	
	act_percentage_genome=$(echo $fastafile | \
		python -c """
import sys
import re
from Bio.SeqIO.FastaIO import SimpleFastaParser

fastafile = re.sub(\"\n\", \"\", sys.stdin.readlines()[0])
with open(fastafile, \"r\") as f:
    header, seq = list(SimpleFastaParser(f))[0]

    missing = seq.count(\"-\") + seq.count(\"N\")
    coverage = 29903-missing
    coverage = round(coverage/29903*100, 2)
    coverage = f\"{coverage}\"
    print(coverage)
    """
	)
	
	echo $sample $avg_cov $act_percentage_genome $human_contamination | tr " " "\t"
done

