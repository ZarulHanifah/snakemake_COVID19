#!/usr/bin/env python

# Example: python merge_df.py \
#					table1.tsv \
#					table2.csv \


import sys
import pandas as pd

def merge(table1, table2):
	df1 = pd.read_csv(table1, sep = "\t", index_col = 0)
	df2 = pd.read_csv(table2, sep = ",", index_col = 0)
	
	df1 = df1.merge(df2, left_index = True, right_index = True)
	df1.to_csv(sys.stdout, sep = "\t")


if __name__ == "__main__":
	merge(sys.argv[1], sys.argv[2])