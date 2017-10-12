#!/usr/bin/env python
#
# combine_heteropecilly.py

'''
combine heteropecilly data, intervals are extracted from files
for instance,
RemPIP90.bor contains the interval for eliminating 90% most heteropecillous 
positions, so creating the alignment of the 10% most homopecillous positions

so the heteropecilly score of intervals in RemPIP10.bor should be max (9 of 9)
while the score of any sites not included in RemPIP90.bor should be 0 of 9

combine_heteropecilly.py -i *.bor > hp_by_site.tab
'''

import sys
import argparse
import time
from itertools import izip,islice
from collections import Counter

def read_site_info(borfiles):
	hpbysite = {}
	for borfile in reversed(borfiles):
		hetpecvalue = 100-int(borfile.split(".")[0][-2:])
		intcount = 0
		sitelist = []
		print >> sys.stderr, "# Reading intervals from {}, with value {}".format(borfile, hetpecvalue), time.asctime()
		for line in open(borfile, 'r'):
			line = line.strip()
			if line and line[0]!="#":
				intervalsplits = line.split(" ")
				if len(intervalsplits)==1:
					intcount = int(intervalsplits[0])
					print >> sys.stderr, "# should contain {} intervals".format(intcount)
				else:
					print >> sys.stderr, "# line contains {} elements".format(len(intervalsplits))
					for val1, val2 in izip(islice(intervalsplits,0,None,2), islice(intervalsplits,1,None,2)):
						sitelist.extend( range( int(val1), int(val2)+1 ) )
					print >> sys.stderr, "# for {} sites".format(len(sitelist))
			for site in sitelist:
				hpbysite[site] = hetpecvalue
		print >> sys.stderr, "# counted {} cumulative sites".format(len(hpbysite)), time.asctime()
		c = Counter(hpbysite.itervalues())
		print >> sys.stderr, c
	return hpbysite

def main(argv, wayout):
	if not len(argv):
		argv.append('-h')
	parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=__doc__)
	parser.add_argument('-f','--fasta', action="store_true", help="print output as a fasta format line for alignments")
	parser.add_argument('-i','--inputs', nargs="*", help="heteropecilly data files")
	parser.add_argument('-l','--length', type=int, help="alignment length, if intervals do not extend to the end")
	args = parser.parse_args(argv)

	allhpbysite = {}
	allhpbysite = read_site_info(args.inputs)

	print >> sys.stderr, "# found information for {} total sites".format(len(allhpbysite)), time.asctime()
	print >> sys.stderr, "# highest numbered site is {}".format(max(allhpbysite.keys())), time.asctime()

	# adjust length for python indexing
	lastsite = args.length if args.length else max(allhpbysite.keys())

	if args.fasta:
		hpstring = "".join( str(allhpbysite.get(i+1,0)/10) for i in range(lastsite) )
		print >> sys.stdout, ">Heteropecilly_score\n{}".format( hpstring )
	else: # otherwise normal tabular format
		for i in range(lastsite):
			print >> sys.stdout, "{}\t{}".format( i+1, allhpbysite.get(i+1,0) )

if __name__ == "__main__":
	main(sys.argv[1:], sys.stdout)
