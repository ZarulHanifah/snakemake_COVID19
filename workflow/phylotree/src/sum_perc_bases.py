import sys
import os
import pandas as pd
import numpy as np

df = pd.read_csv(sys.stdin, sep = "\t",
                 names = ["chrom", "start", "end"])
nsum = (df["end"] - df["start"]).sum()
perc = str(round(nsum/29700*100, 2)) + "%"
print(perc)