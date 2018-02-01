#!/usr/bin/env python
#
# sitewise_get_strong_sites.py  created 2018-02-01

'''sitewise_get_strong_sites.py  last modified 2018-02-01
    subset a supermatrix using only strong-sites identified by site-wise RAxML

sitewise_get_strong_sites.py -a matrix1.aln -l RAxML_perSiteLLs.matrix1.tab -o strong_alignment.aln

    matrix can be in alternate formats (use -f), and gzipped

    generate the tabular sitewise results using:
sitewise_ll_to_columns.py RAxML_perSiteLLs.matrix1 > RAxML_perSiteLLs.matrix1.tab
'''

import sys
import argparse
import time
import gzip
from Bio import AlignIO

def make_strong_alignment(fullalignment, alignformat, valsbysite, mindlnl):
	'''read large alignment, and return a new alignment of only strong sites'''
	if fullalignment.rsplit('.',1)[1]=="gz": # autodetect gzip format
		opentype = gzip.open
		print >> sys.stderr, "# reading alignment {} as gzipped".format(fullalignment), time.asctime()
	else: # otherwise assume normal open
		opentype = open
		print >> sys.stderr, "# reading alignment {}".format(fullalignment), time.asctime()
	alignedseqs = AlignIO.read(opentype(fullalignment), alignformat)
	numtaxa = len(alignedseqs)
	allength = alignedseqs.get_alignment_length()

	print >> sys.stderr, "# Alignment contains {} taxa for {} sites, including gaps".format( numtaxa, allength )
	strongcounter = 0
	newalign = alignedseqs[:,0:0] # start with blank alignment

	for i in range(allength): # sites begin at 0 while lnl begins at 1
		if valsbysite[i+1] >= mindlnl:
			strongcounter += 1
			newalign += alignedseqs[:,i:i+1]
	print >> sys.stderr, "# New alignment contains {} strong sites".format( strongcounter )
	return newalign

def read_tabular_ln(lntabular, treelist):
	'''read tabular log-likelihood results and return a dict where keys are position and value is abs dlnL'''
	lndict = {}
	linecounter = 0
	print >> sys.stderr, "# Using columns {} and {} for topologies".format(*treelist)
	print >> sys.stderr, "# Reading log-likelihood by site from {}".format(lntabular), time.asctime()
	for line in open(lntabular,'r'):
		line = line.strip()
		if line and line[0]!="#": # ignore blank and comment lines
			linecounter += 1
			if linecounter < 2:
				continue
			lsplits = line.split('\t')
			pos = int(lsplits[0]) # sites begin at 1
			dlnl = float(lsplits[1]) - float(lsplits[2])
			absdlnl = abs(dlnl)
			lndict[pos] = absdlnl
	print >> sys.stderr, "# Found log-likelihood for {} sites".format( len(lndict) ), time.asctime()
	return lndict

def main(argv, wayout):
	if not len(argv):
		argv.append('-h')
	parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=__doc__)
	parser.add_argument('-a','--alignment', help="supermatrix alignment")
	parser.add_argument('-f','--format', default="fasta", help="alignment format [fasta]")
	parser.add_argument('-l','--log-likelihood', help="tabular log-likelihood data file from RAxML")
	parser.add_argument('-m','--dlnl-minimum', default=0.5, type=float, help="minimum difference in lnL [0.5]")
	parser.add_argument('-o','--output', help="name of output file", required=True)
	parser.add_argument('-t','--trees', default="1,2", help="two columns of trees 1 and 2, as comma-separated ints [1,2]")
	args = parser.parse_args(argv)

	treelist = [int(i) for i in args.trees.split(",")]
	valsbysite = read_tabular_ln(args.log_likelihood, treelist)
	strongalignment = make_strong_alignment(args.alignment, args.format, valsbysite, args.dlnl_minimum)
	AlignIO.write(strongalignment, args.output, args.format)
	print >> sys.stderr, "# New alignment written to {}".format(args.output), time.asctime()

if __name__ == "__main__":
	main(sys.argv[1:], sys.stdout)
