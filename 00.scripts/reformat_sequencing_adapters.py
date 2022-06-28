"""
Reformat the sequencing adapters and contaminants txt file
to be a FASTA format
"""

import sys
infile = sys.argv[1]

with open(infile, "r") as f:
    for line in f:
        line = line.strip()
        if len(line) > 1:
            header=line.split("\t")[0]
            sequence = line.split("\t")[-1]
            print(">" + header + "\n" + sequence)
