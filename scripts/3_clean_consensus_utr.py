#!/usr/bin/python


from Bio.SeqIO import parse
from Bio.Seq import Seq
import sys

infileName = sys.argv[1]
infile = parse(infileName, "fasta")
Lmin = 100

outfile = open("cleaned_consensus.fas", "w")
for i in infile:
	name = i.id
	seq = str(i.seq)
	newSeq = ""
	for j in seq:
		if j!="N":
			newSeq += j
	if len(newSeq) >= Lmin:
		outfile.write(">{0}\n{1}\n".format(name, newSeq))
outfile.close()

