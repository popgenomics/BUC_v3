#!/usr/bin/python
# script to run on *.orf.new files and producing the "output_ENCprime_JoNo.txt" file
# by launching SeqCount and ENCprime (C codes written by John Novembre)
# ./ENCp_JohnNovembre.py *.orf.new

from Bio.SeqIO import parse
import sys
import os
step = 100

fastaFile = sys.argv[1]

alignement = {}
cnt = -1

GCcontent = {}
input = parse(fastaFile, "fasta")
for i in input:
	cnt += 1
	alignement[cnt] = {}
	alignement[cnt]["name"] = i.id
	alignement[cnt]["seq"] = i.seq
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

res = {}

cnt = -1
nLocus = len(alignement)
for i in range(0, nLocus, step):
	nLocusTMP = 0
	cnt += 1
	tmp = ""
	if(i + step < nLocus):
		limMax = i + step
	else:
		limMax = nLocus
	for j in range(i, limMax, 1):
		nLocusTMP += 1
		tmp += ">{0}\n{1}\n".format(alignement[j]["name"], alignement[j]["seq"])
	tmpFile = "tmp_{0}.fas".format(cnt)
	output = open(tmpFile, "w")
	output.write(tmp)
	output.close()
	command1 = "SeqCount -c {0} {1}".format(tmpFile, nLocusTMP) 
	command2 = "SeqCount -n {0} {1}".format(tmpFile, nLocusTMP)
	command3 = "ENCprime {0}.codcnt {0}.acgtfreq 1 ExResults_{0} 1 -q".format(tmpFile)
	testCommand1 = -1
	testCommand2 = -1
	testCommand3 = -1
	testCommand1 = os.system(command1)
	testCommand2 = os.system(command2)
	testCommand3 = os.system(command3)
	if testCommand1 != 0 or testCommand2 != 0 or testCommand3 != 256:
		continue
	acgtfile = "tmp_{0}.fas.acgtfreq".format(cnt)
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
	codcntfile = "tmp_{0}.fas.codcnt".format(cnt)
	input = open(codcntfile, "r")
	for k in input:
		k = k.strip()
		if "Contig" in k:
			k = k.replace(">", " ")
			k = k.split(" ")
			contigName = k[0]
			res[contigName]["codcnt"] = "\t".join(k[1::])
	input.close()
	ExResultFile = "ExResults_tmp_{0}.fas".format(cnt)
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
	res2 += "{0}\t{1}\t{2}\t{3}\t{4:.5f}\t{5:.5f}\t{6:.5f}\t{7:.5f}\n".format(i, res[i]["encp"], res[i]["acgtfreq"], res[i]["codcnt"], GCcontent[i]['GCtot'], GCcontent[i]['GC1'], GCcontent[i]['GC2'], GCcontent[i]['GC3'])

commandrm = "rm -rf *tmp*.fas*"
os.system(commandrm)

header1 = "Nc Ncp ScaledChi SumChi df p B_KM n_codons"
header2 = "fA fC fG fT"
header3 = "TTT TTC TTA TTG TCT TCC TCA TCG TAT TAC TAA TAG TGT TGC TGA TGG CTT CTC CTA CTG CCT CCC CCA CCG CAT CAC CAA CAG CGT CGC CGA CGG ATT ATC ATA ATG ACT ACC ACA ACG AAT AAC AAA AAG AGT AGC AGA AGG GTT GTC GTA GTG GCT GCC GCA GCG GAT GAC GAA GAG GGT GGC GGA GGG"

header1 = "\t".join(header1.split(" "))
header2 = "\t".join(header2.split(" "))
header3 = "\t".join(header3.split(" "))

header = "Contig\t" + header1 + "\t" + header2 + "\t" + header3 + "\tf_GC_tot\tf_GC_1\tf_GC_2\tf_GC_3\n"

res2 = header + res2

outfile = open("output_ENCprime_JoNo.txt", "w")
outfile.write(res2)
outfile.close()

