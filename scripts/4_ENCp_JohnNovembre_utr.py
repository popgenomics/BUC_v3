#!/usr/bin/python
# script to run on *.orf.new files and producing the "output_ENCprime_JoNo.txt" file
# by launching SeqCount and ENCprime (C codes written by John Novembre)
# ./ENCp_JohnNovembre.py *.orf.new

from Bio.SeqIO import parse
import sys
import os
step = 100

cdsFile = sys.argv[1]
utrFile = sys.argv[2]

list_of_contigs = []

GCcontent = {}
# CDS
alignement_cds = {}
input = parse(cdsFile, "fasta")
for i in input:
	list_of_contigs.append(i.id)
	alignement_cds[i.id] = {}
	alignement_cds[i.id]["name"] = i.id
	alignement_cds[i.id]["seq"] = i.seq
	tmp = str(i.seq)
	pos1 = tmp[0::3]
	pos2 = tmp[1::3]
	pos3 = tmp[2::3]
	GCcontent[i.id] = {}
	GCcontent[i.id]["GCtot"] = (tmp.count("G") + tmp.count("C"))/(len(tmp)*1.0)
	GCcontent[i.id]["GC1"] = (pos1.count("G") + pos1.count("C"))/(len(pos1)*1.0)
	GCcontent[i.id]["GC2"] = (pos2.count("G") + pos2.count("C"))/(len(pos2)*1.0)
	GCcontent[i.id]["GC3"] = (pos3.count("G") + pos3.count("C"))/(len(pos3)*1.0)
input.close()

# UTR
retained_contigs = []
alignement_utr = {}
input = parse(utrFile, "fasta")
for i in input:
	if i.id in list_of_contigs:
		retained_contigs.append(i.id)
		alignement_utr[i.id] = {}
		alignement_utr[i.id]["name"] = i.id
		alignement_utr[i.id]["seq"] = i.seq
		tmp = str(i.seq)
		GCcontent[i.id]["GCutr"] = (tmp.count("G") + tmp.count("C"))/(len(tmp)*1.0)
		GCcontent[i.id]["Lutr"] = len(tmp) 
input.close()

res = {}

cnt = -1
nLocus = len(retained_contigs)
for i in range(0, nLocus, step):
	nLocusTMP = 0
	cnt += 1
	tmp_utr = ""
	tmp_cds = ""
	if(i + step < nLocus):
		limMax = i + step
	else:
		limMax = nLocus
	for j in range(i, limMax, 1):
		nLocusTMP += 1
		locusTMP = retained_contigs[j]
		tmp_utr += ">{0}\n{1}\n".format(alignement_utr[locusTMP]["name"], alignement_utr[locusTMP]["seq"])
		tmp_cds += ">{0}\n{1}\n".format(alignement_cds[locusTMP]["name"], alignement_cds[locusTMP]["seq"])
	
	tmpFile_utr = "tmp_{0}_utr.fas".format(cnt)
	output = open(tmpFile_utr, "w")
	output.write(tmp_utr)
	output.close()
	
	tmpFile_cds = "tmp_{0}_cds.fas".format(cnt)
	output = open(tmpFile_cds, "w")
	output.write(tmp_cds)
	output.close()
	
	command1 = "SeqCount -c {0} {1}".format(tmpFile_cds, nLocusTMP) # codon count 
	command2 = "SeqCount -n {0} {1}".format(tmpFile_utr, nLocusTMP) # nucleotide composition
	command3 = "ENCprime {0}.codcnt {1}.acgtfreq 1 ExResults_tmp_{2} 1 -q".format(tmpFile_cds, tmpFile_utr, cnt)
	testCommand1 = -1
	testCommand2 = -1
	testCommand3 = -1
	testCommand1 = os.system(command1)
	testCommand2 = os.system(command2)
	testCommand3 = os.system(command3)
	if testCommand1 != 0 or testCommand2 != 0 or testCommand3 != 256:
		continue
	
	acgtfile = "{0}.acgtfreq".format(tmpFile_utr)
	input = open(acgtfile, "r")
	for k in input:
		k = k.strip()
		if "Contig" in k:
			k = k.replace(">", "")
			k = k.split(" ")
			contigName = k[0]
			res[contigName] = {}
			res[contigName]["acgtfreq"] = "\t".join(k[1::])
	input.close()
	
	codcntfile = "{0}.codcnt".format(tmpFile_cds)
	input = open(codcntfile, "r")
	for k in input:
		k = k.strip()
		if "Contig" in k:
			k = k.replace(">", " ")
			k = k.split(" ")
			contigName = k[0]
			res[contigName]["codcnt"] = "\t".join(k[1::])
	input.close()
	
	ExResultFile = "ExResults_tmp_{0}".format(cnt)
	input = open(ExResultFile, "r")
	for k in input:
		k = k.strip()
		if "Contig" in k:
			k = k.replace(":", "")
			k = k.split(" ")
			contigName = k[0]
			res[contigName]["encp"] = "\t".join(k[1::])
	input.close()

res2 = ""
for i in res:
	res2 += "{0}\t{1}\t{2}\t{3}\t{4:.5f}\t{5:.5f}\t{6:.5f}\t{7:.5f}\t{8:.5f}\t{9:.5f}\n".format(i, res[i]["encp"], res[i]["acgtfreq"], res[i]["codcnt"], GCcontent[i]['GCtot'], GCcontent[i]['GC1'], GCcontent[i]['GC2'], GCcontent[i]['GC3'], GCcontent[i]['GCutr'], GCcontent[i]['Lutr'])

commandrm = "rm -rf *tmp*"
os.system(commandrm)

header1 = "Nc Ncp ScaledChi SumChi df p B_KM n_codons"
header2 = "fA fC fG fT"
header3 = "TTT TTC TTA TTG TCT TCC TCA TCG TAT TAC TAA TAG TGT TGC TGA TGG CTT CTC CTA CTG CCT CCC CCA CCG CAT CAC CAA CAG CGT CGC CGA CGG ATT ATC ATA ATG ACT ACC ACA ACG AAT AAC AAA AAG AGT AGC AGA AGG GTT GTC GTA GTG GCT GCC GCA GCG GAT GAC GAA GAG GGT GGC GGA GGG"

header1 = "\t".join(header1.split(" "))
header2 = "\t".join(header2.split(" "))
header3 = "\t".join(header3.split(" "))

header = "Contig\t" + header1 + "\t" + header2 + "\t" + header3 + "\tf_GC_tot\tf_GC_1\tf_GC_2\tf_GC_3\tf_GC_utr\tL_utr\n"

res2 = header + res2

outfile = open("output_ENCprime_JoNo.txt", "w")
outfile.write(res2)
outfile.close()

