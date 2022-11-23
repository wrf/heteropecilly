#!/usr/bin/env python
#
# calculate_solvent_exposure.py  created 2017-10-09

'''calculate_solvent_exposure.py

calculate_solvent_exposure.py -a

    chain information must first be parsed from uniprot_sprot.dat.gz
parse_swissprot_data.py
'''

import sys
import os
import gzip
import argparse
import time
import re
from glob import glob
from collections import defaultdict
from Bio import AlignIO
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.ResidueDepth import ResidueDepth

def calculate_se(pdbfile, chainlist):
	parser = PDBParser()
	pdbid = pdbfile.split(".")[0].upper()
	structure = parser.get_structure(pdbfile, pdbfile)
	model = structure[0]
	chain = model["A"]
	residue = chain[1]
	atom = residue["CA"]
	hse = HSExposureCB(model)
	for hset in hse.property_list: # hse tuple looks like: (<Residue GLU het=  resseq=1005 icode= >, (3, 19, 0.0))
		residue, segroup = hset # segroup is tuple of 3 values

	# requires installing msms from:
	# http://mgltools.scripps.edu/downloads
	rd = Bio.PDB.ResidueDepth(structure[0])
	for rset in rd.property_list:

	rda = ResidueDepth(model["A"])

def read_prots_w_pdb(pdbtabular):
	'''read tabular PDB summary, return dict of dicts where keys are PDB-ID, prot ID then list of chains'''
	# assume format is:
	# prot name   PDB-ID  prot length  chains      chain coverage  percentage
	# GROA_HUMAN  1MSH    107          A/B=35-106  72              0.673
	chainsbypdb = defaultdict(dict)
	linecounter = 0
	print >> sys.stderr, "# Reading PDB information from {}".format(pdbtabular), time.asctime()
	for line in open(pdbtabular,'r'):
		line = line.strip()
		if line and line[0]!="#": # skip blank and comment lines
			linecounter += 1
			lsplits = line.split("\t")
			protid = lsplits[0]
			pdbid = lsplits[1]
			chaininfo = lsplits[3]
			chains = chaininfo.split('=')[0].split(',')
			chainsbypdb[pdbid][protid] = chains
	print >> sys.stderr, "# Counted {} lines for {} entries".format( linecounter, len(chainsbypdb) ), time.asctime()
	return chainsbypdb

def main(argv, wayout):
	if not len(argv):
		argv.append('-h')
	parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=__doc__)
	parser.add_argument('-a','--alignments', nargs="*", help="alignment files or directory")
	parser.add_argument('-c','--pdb-chains', help="PDB chain information as tabular", required=True)
	parser.add_argument('-f','--format', default="fasta", help="alignment format [fasta]")
	parser.add_argument('-l','--length-cutoff', default=0.5, type=float, help="minimum length coverage for proteins [0.5]")
	parser.add_argument("-p","--pdbs", help="PDB files or directory")
	args = parser.parse_args(argv)



	if os.path.isdir(args.alignments[0]):
		print >> sys.stderr, "# Reading alignments files from directory {}".format(args.alignments[0]), time.asctime()
		globstring = "{}*.aln".format(args.alignments[0])
		alignmentfiles = glob(globstring)
	elif os.path.isfile(args.alignments[0]):
		alignmentfiles = args.alignments
	else:
		raise OSError("ERROR: Unknown alignments, exiting")

	chaininfo = read_prots_w_pdb(args.pdb_chains)

	alignwpdbcounter = 0
	for alignfile in alignmentfiles:
		uniprot_from_alignment(alignfile, pdbentries, protlengthdict, pdblengthdict, args.length_cutoff)









if __name__ == "__main__":
	main(sys.argv[1:], sys.stdout)
