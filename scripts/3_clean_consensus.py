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
	# multiple of 3
	while(len(seq)%3 != 0):
		seq = seq[:1]
	newSeq = ""
	for j in range(0, len(seq), 3):
		codon = seq[j] + seq[j+1] + seq[j+2]
		if "N" not in codon:
			newSeq += codon
	if len(newSeq) >= Lmin:
		outfile.write(">{0}\n{1}\n".format(name, newSeq))
outfile.close()

