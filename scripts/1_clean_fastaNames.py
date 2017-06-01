#!/usr/bin/pypy
import sys

infileName = sys.argv[1]

outfile = open("cleaned_" + infileName, "w")

infile = open(infileName, "r")
for i in infile:
	if ">" in i:
		i = ">" + i.strip().split(">")[-1::][0] + "\n"
	outfile.write(i)
outfile.close()
infile.close()

