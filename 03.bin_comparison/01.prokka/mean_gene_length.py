"""
Script to calculate mean gene length from Prokka output
Dylan Maghini
July 5 2022
"""

import sys
infile = sys.argv[1]
name = sys.argv[2]

gene_length = 0
gene_count = 0
with open(infile, "r") as f:
    f.readline()
    for line in f:
        vals = line.strip().split("\t")
        if vals[1] == "CDS":
            gene_length += int(vals[2])
            gene_count += 1

print(name + "\t" + str(float(gene_length)/float(gene_count)))
