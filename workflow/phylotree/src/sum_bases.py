import sys
import os
import pandas as pd
import numpy as np

df = pd.read_csv(sys.stdin, sep = "\t",
                 names = ["chrom", "start", "end"])
print((df["end"] - df["start"]).sum())