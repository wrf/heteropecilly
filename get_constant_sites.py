#!/usr/bin/env python
#
# get_constant_sites.py

'''
get_constant_sites.py simion2017_20CAT_hp_by_site_w_const.tab > c_sites.txt
'''

import sys

if len(sys.argv) < 2:
	print >> sys.stderr, __doc__
else:
	partitions = [] # list of tuples of intervals
	for line in open(sys.argv[1],'r'):
		line = line.strip()
		if line and line[0]!="#":
			lsplits = line.split("\t")
			if lsplits[1]=="10":
				print >> sys.stdout, lsplits[0]
