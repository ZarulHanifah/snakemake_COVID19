#!/usr/bin/env python

# Example: python merge_stats.py \
#					basic_stats.tsv \
#					assem_stats.csv \


import sys
import pandas as pd

def merge(basic_table, assem_table):
	basic_s = pd.read_csv(basic_table, sep = "\t", index_col = 0)
	assem_s = pd.read_csv(assem_table, sep = "\t", index_col = 0)
	assem_s = assem_s.loc[:, ["N_count", "Gaps"]]

	basic_s = basic_s.merge(assem_s, left_index = True, right_index = True)
	basic_s.to_csv(sys.stdout, sep = "\t")


if __name__ == "__main__":
	merge(sys.argv[1], sys.argv[2])