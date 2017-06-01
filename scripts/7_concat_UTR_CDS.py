#!/usr/bin/python

# code to concat 5P_UTR + CDS + 3P_UTR

import sys
from Bio.SeqIO import parse
from Bio.Seq import Seq

species = sys.argv[1]
utrFile = sys.argv[2]
cdsFile = sys.argv[3]
#species = "Abatus_cordatus"
#utrFile = "/home/croux/Documents/BUC_v3/popphyl_contam0.2/utr/Abatus_cordatus/Abatus_cordatus.utr.fas"
#cdsFile = "/home/croux/Documents/BUC_v3/popphyl_contam0.2/cds/Abatus_cordatus/Abatus_cordatus.cds.alr.0.2.fas"

utr_3 = {}
utr_5 = {}
cds = {}

list_of_individuals = []
list_of_contigs = []

# deal with utr 3' and 5'
infile = parse(utrFile, "fasta")
for i in infile:
	tmp = i.id.split("|")
	utr = tmp[0].split("_")[0]
	contig = tmp[0].split("_")[1]
	ind = tmp[2] + "|" + tmp[3]

	if ind not in list_of_individuals:
		list_of_individuals.append(ind)
	
	if contig not in list_of_contigs:
		list_of_contigs.append(contig)
	
	if utr == "5P":
		if contig not in utr_5:
			utr_5[contig] = {}
			utr_5[contig]["L"] = len(i.seq)
		if ind not in utr_5[contig]:
			utr_5[contig][ind] = i.seq
	
	if utr == "3P":
		if contig not in utr_3:
			utr_3[contig] = {}
			utr_3[contig]["L"] = len(i.seq)
		if ind not in utr_3[contig]:
			utr_3[contig][ind] = i.seq
infile.close()

# deal with CDS
infile = parse(cdsFile, "fasta")
for i in infile:
	tmp = i.id.split("|")
	contig = tmp[0]
	ind = tmp[2] + "|" + tmp[3]

	if ind not in list_of_individuals:
		list_of_individuals.append(ind)
	
	if contig not in list_of_contigs:
		list_of_contigs.append(contig)

	if contig not in cds:
		cds[contig] = {}
		cds[contig]["L"] = len(i.seq)
	if ind not in cds[contig]:
		cds[contig][ind] = i.seq
infile.close()

# concatenation:
outfile = open("5P_CDS_3P.fas", "w")

for contig in list_of_contigs:
	for ind in list_of_individuals:
		tmp = ""
 		if contig in utr_5:
			if ind in utr_5[contig]:
				tmp += str(utr_5[contig][ind])
			else:
				tmp += utr_5[contig]["L"] * "N"
		
 		if contig in cds:
			if ind in cds[contig]:
				tmp += str(cds[contig][ind])
			else:
				tmp += cds[contig]["L"] * "N"
		
 		if contig in utr_3:
			if ind in utr_3[contig]:
				tmp += str(utr_3[contig][ind])
			else:
				tmp += utr_3[contig]["L"] * "N"
		
		outfile.write(">{0}|{1}|{2}\n{3}\n".format(contig, species, ind, tmp))
outfile.close()

