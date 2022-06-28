import sys
import math
import gzip

readfile = sys.argv[1]
name = sys.argv[2]

num_reads = 0 
num_bases = 0 
length_dict = {}

for i in range(0, 2000000, 100):
    length_dict[i] = 0

counter = 1
with gzip.open(readfile, "r") as f:
#with open(readfile, "r") as f:
    for line in f:
        if counter == 2:
            num_reads += 1
            length = len(line.strip())
            num_bases += length
            length_dict[math.floor(float(length)/100)*100] += 1
        counter += 1
        if counter == 5: 
            counter = 1

with open(name + "_report.tsv", "w") as outreport:
    outreport.write("Sample\tTotalReads\tTotalBases\n")
    outreport.write(name + "\t" + str(num_reads) + "\t" + str(num_bases) + "\n")

with open(name + "_read_dist.tsv", "w") as outdist:
    outdist.write("Name\tLength\tCount\n")
    for length in length_dict:
        if length_dict[length] != 0:
            outdist.write(name + "\t" + str(length) + "\t" + str(length_dict[length]) + "\n")

 
