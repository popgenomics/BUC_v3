#!/usr/bin/python

import sys

species = sys.argv[1]

infile = open("output_summarized_v1.txt", "r")
header_BUC = infile.readline()
header_BUC = header_BUC.strip()

stats_BUC = infile.readline()
stats_BUC = stats_BUC.strip()
infile.close()


infile = open("/home/croux/Documents/BUC_v3/scripts/life_history_traits.csv", "r")
header_traits = infile.readline()
header_traits = header_traits.strip().replace(",", "\t")

test = 0
for i in infile:
	if species in i:
		stats_traits = i.strip().replace(",", "\t")
		test = 1
infile.close()

if test==0:
	sys.exit("\n\033[1;31m ERROR: {0} was not found in life_history_traits.csv\033[0m\n".format(species))

header = header_BUC + "\t" + header_traits + "\n"
stats = stats_BUC + "\t" + stats_traits + "\n"


outfile = open("output_summarized_v2.txt", "w")
outfile.write(header)
outfile.write(stats)
outfile.close()

