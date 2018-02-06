#!/usr/bin/env python

# sitewise_columns_to_fasta.py

'''sitewise_columns_to_fasta.py last modified 2018-02-02
    convert tabular sitewise log-likelihood to a base16 FASTA string

sitewise_columns_to_fasta.py RAxML_perSiteLLs.matrix1.tab > RAxML_perSiteLLs.matrix1.fasta

    difference in lnL is ln(tree1)-ln(tree2)

    by default, this is in base16 format
    where neutral sites (no strong dlnL for either tree) are 7 or 8
    all values above 8 (meaning 9-f) favor tree 1
    and all values below 7 (meaning 0-6) favor tree 2

    generate the tabular sitewise results using:
sitewise_ll_to_columns.py RAxML_perSiteLLs.matrix1 > RAxML_perSiteLLs.matrix1.tab
'''

import sys
from collections import defaultdict

if len(sys.argv) < 2:
	print >> sys.stderr, __doc__
else:
	fastastring = ""
	for line in open(sys.argv[1],'r'):
		line = line.strip()
		if line:
			lsplits = line.split("\t")
			if lsplits[0] == "site":
				continue
			if lsplits[1] == "const":
				fastastring += "x"
			else:
				dlnl = float(lsplits[1]) - float(lsplits[2])
				absdlnl = abs(dlnl)
				if absdlnl >= 0.5: # meaning above noise threshold
					if dlnl > 0: # favors tree 1
						if absdlnl < 1: # meaning between 0.5 and 1
							adjdlnl = 9
						else: # greater than 1
							adjdlnl = dlnl + 9
					elif dlnl < 0: # favors tree 2
						if absdlnl < 1: # meaning between 0.5 and 1
							adjdlnl = 6
						else: # greater than 1, so add 6.99 so values are less than 6
							adjdlnl = dlnl + 6.99
					# should be nothing else
				else: # meaning below noise threshold, so abs value is 0.5 or less
					adjdlnl = dlnl + 8
				if adjdlnl >= 16:
					adjdlnl = 15
				elif adjdlnl <= 0:
					adjdlnl = 0
				hexdlnl = hex( int(adjdlnl) )[-1]
				#lnl = int( float(lsplits[1]) - float(lsplits[2]) + 5) # was for base 10
				fastastring += str(hexdlnl)
	print >> sys.stdout, ">diffLnL\n{}".format(fastastring)
