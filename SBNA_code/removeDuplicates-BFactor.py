#This code averages the BFactor value for duplicates

import sys
import numpy

file = open(sys.argv[1],'r')
outfile = open(sys.argv[1]+"_duplicatesRemoved",'w')

BFactors = {}

for line in file:
	line = line.strip("\n").split()
	res = line[0]
	if not res in BFactors.keys():
		BFactors[res]=[]
	BFactors[res].append(float(line[1]))

for res in BFactors.keys():
	outfile.write(res+"\t"+str(numpy.mean(BFactors[res]))+"\n")

outfile.close()
