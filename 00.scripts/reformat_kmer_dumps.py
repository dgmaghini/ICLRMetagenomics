"""
Script to take jellyfish "dump" files and
convert the fasta format to a TSV format
"""

import sys

infile = sys.argv[1]
outfile = sys.argv[2]

with open(infile, "r") as f:
    with open(outfile, "w") as outf:
        for line in f:
            if line.startswith(">"):
                count=line[1:].strip()
            else:
                kmer = line.strip()
                outf.write(kmer + "\t" + count + "\n")
