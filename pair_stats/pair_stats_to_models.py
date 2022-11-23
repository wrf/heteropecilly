#!/usr/bin/env python

# pair_stats_to_models.py

'''
pair_stats_to_models.py pair_stats.tab > models.txt
'''

import sys

if len(sys.argv) < 2:
	print >> sys.stderr, __doc__
else:
	for line in open(sys.argv[1],'r'):
		lsplits = line.split("\t")
		if lsplits[0]=="partition":
			continue
		gene = lsplits[1].split("|")[2]
		partition = lsplits[0].split("_")[2]
		print >> sys.stdout, "LGF, {} = {}".format(gene, partition)

