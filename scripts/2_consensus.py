#!/usr/bin/python

from random import sample
from Bio.SeqIO import parse
from Bio.Seq import Seq
import sys

def consensus(x):
	res = ""
	for position in range(len(x[0])):
		bases = {'A':0, 'T':0, 'C':0, 'G':0}
		for sequence in x:
			if sequence[position] in bases:
				bases[sequence[position]] += 1
		if bases['A'] + bases['T'] + bases['C'] + bases['G'] == 0:
			res += 'N'
		else:
			max_count = max(bases.values())
			max_base = []
			for i in ['A', 'T', 'C', 'G']:
				if bases[i] == max_count:
					max_base.append(i)
			if len(max_base) == 1:
				res += max_base[0]
			else:
				res += sample(max_base, 1)[0]
	return(res)

alignement = {}
list_of_contigs = []

infileName = sys.argv[1]
#infileName = "/home/croux/Documents/BUC_v2/data/cleaned_Thymelicus_sylvestris/cleaned_Thymelicus_sylvestris.cds.alr.0.2.fas"

infile = parse(infileName, "fasta")

for i in infile:
	contigName = i.id.split("|")[0]
	if contigName not in alignement:
		list_of_contigs.append(contigName)
		alignement[contigName] = []
	alignement[contigName].append(i.seq)
infile.close()

print("{0} contigs in {1}".format(len(alignement), infileName))

outfile = open("consensus.fas", "w")
for i in list_of_contigs:
	outfile.write(">{0}\n{1}\n".format(i, consensus(alignement[i])))
outfile.close()


