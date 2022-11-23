#!/usr/bin/env python
#
# get_best_structures.py  created 2017-10-09

'''get_best_structures.py  last modified 2018-02-05

get_best_structures.py -u uniprot_sprot.dat.gz -s 9606

    use a Uniprot database summary from parse_swissprot_data.py with -U
get_best_structures.py -U uniprot-9606_prots_w_pdb.tab -a blast_alignments_hp/ -s 9606 -c pdb_hp_commands.sh > human_prots_w_pdb.tab

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

class pdbmatch:
	'''class for storing information about PDB entries'''
	def __init__(self, pdbid, geneid):
		self.id = pdbid
		self.gene = gene
		self.chains = [] # chains should be letters, all same length
		self.span = (1,0)
	def add_chain(self, chainstring):
		csplits = chainstring.split("=")
		self.span = [int(n) for n in csplits[1].split("-")]
	def span_length(self):
		return self.span[1] - self.span[0] + 1

# 1433B_HUMAN	2BQ0	246	A/B=2-239	238	0.967

def parse_uniprot(uniprotdata, speciesfilter):
	if uniprotdata.rsplit('.',1)[1]=="gz": # autodetect gzip format
		opentype = gzip.open
		print >> sys.stderr, "# reading data from {} as gzipped".format(uniprotdata), time.asctime()
	else: # otherwise assume normal open
		opentype = open
		print >> sys.stderr, "# reading data from {}".format(uniprotdata), time.asctime()

	typecounter = defaultdict(int)
	protlengthdict = {} # key is protid, value is length
	pdblengthdict = defaultdict(dict) # key is 4-letter pdb ID, value is dict where key is protID and value is coverage
	pdbspandict = defaultdict(dict) # key is 4-letter pdb ID, value is dict where key is chain and value is tuple
	pdbentries = defaultdict(list)
	skipentry = False
	for line in opentype(uniprotdata,'r'):
		line = line.rstrip()
		if line:
			key = line[0:5].strip()
			if key: # ignorning blank lines which correspond to sequence
				typecounter[key] += 1
				if key=="//": # reset
					skipentry = False
				elif skipentry:
					continue
				if key=="ID":
					idsplits = line[5:].split()
					protid = idsplits[0]
					protlength = int(idsplits[2])
					protlengthdict[protid] = protlength
				elif key=="OX":
					species = line[5:].split(";")[0].split("=")[1]
					if species not in speciesfilter:
						skipentry = True
						continue
				elif key=="DR": # for database references
					drsplits = [dr.strip() for dr in line[5:].split(";")]
					if drsplits[0]=="PDB":
						pdbid = drsplits[1]
						pdbentries[protid].append(pdbid)
						pdbrangestring = drsplits[4].replace(".","")
						chain, structspan = pdbrangestring.split("=")
						try:
							pdbrangelist = [ int(n) for n in re.search("=(\d+)-(\d+)", pdbrangestring).groups() ]
						except AttributeError:
							print >> sys.stderr, "WARNING CANNOT PARSE {}".format(line)
						pdbcoverage = pdbrangelist[1]-pdbrangelist[0]+1
						pdblengthdict[pdbid].update( {protid:pdbcoverage} )
						pdbspandict[pdbid].update( {chain:structspan} )
						#print >> sys.stdout, "{}\t{}\t{}\t{}".format(protid, pdbid, protlength, pdbrange)
# FT   DOMAIN       57    173       DOMON. {ECO:0000255|PROSITE-
# FT                                ProRule:PRU00246}.

	print >> sys.stderr, "# counted data types from {}".format(uniprotdata), time.asctime()
	#for idtype in sorted(typecounter.keys()):
	#	print >> sys.stderr, "{}\t{}".format(idtype, typecounter[idtype])
	print >> sys.stderr, "# counted {} species".format( len(pdbentries) ), time.asctime()
	#for k,v in pdbentries.iteritems():
	#	print >> sys.stderr, "{}\t{}".format(k, len(v) )
	return pdbentries, protlengthdict, pdblengthdict, pdbspandict

def parse_filtered_uniprot(uniprotdata):
	protlengthdict = {} # key is protid, value is length
	pdblengthdict = defaultdict(dict) # key is 4-letter pdb ID, value is dict where key is protID and value is coverage
	pdbspandict = defaultdict(dict) # key is 4-letter pdb ID, value is dict where key is chain and value is tuple
	pdbentries = defaultdict(list)
	print >> sys.stderr, "# reading data from {}".format(uniprotdata), time.asctime()
	# expects data as:
	# MSH2_HUMAN  2O8B       934   A=1-934          934     1.000
	# for:     0     1         2         3            4         5
	# entry       ID  protlength chainSpan  chainLength  chainCov
	for line in open(uniprotdata,'r'):
		line = line.strip()
		if line:
			lsplits = line.split("\t")
			protlengthdict[lsplits[0]] = lsplits[2]
			pdblengthdict[lsplits[1]].update( {lsplits[0]:lsplits[4]} )
			pdbentries[lsplits[0]].append(lsplits[1])
			for chainstspan in lsplits[3].split(","):
				try:
					chains, structspan = chainstspan.strip().split("=")
				except ValueError:
					print >> sys.stderr, line
				pdbspandict[lsplits[0]].update( {chains:structspan} )
	print >> sys.stderr, "# found PDB for {} entries".format( len(pdbentries) ), time.asctime()
	return pdbentries, protlengthdict, pdblengthdict, pdbspandict

def get_swiss_to_trembl(trembldata):
	if trembldata.rsplit('.',1)[1]=="gz": # autodetect gzip format
		opentype = gzip.open
		print >> sys.stderr, "# reading data from {} as gzipped".format(trembldata), time.asctime()
	else: # otherwise assume normal open
		opentype = open
		print >> sys.stderr, "# reading data from {}".format(trembldata), time.asctime()

	typecounter = defaultdict(int)
	pdbentries = defaultdict(list)
	skipentry = False
	for line in opentype(trembldata,'r'):
		line = line.rstrip()
		if line:
			key = line[0:5].strip()

#Created:	2017-09-19
#Organism:	 Homo sapiens
#Taxid:	9606
#Models:	43227
#Structures:	44280
#UniProt Version:	2017_08
#UniProtAc	coordinate_id	provider	from	to	template	qmean	qmean_norm	url
#A0A024R161	59b56ef3c730d488240127a1	swissmodel	111	151	5ukl.1.C	-0.3787332599	0.7345491548	https://swissmodel.expasy.org/repository/uniprot/A0A024R161.pdbfrom=111&to=151&template=5ukl.1.C&provider=swissmodel
def parse_swissmodel_repo(swissmodeldata, trembltoswiss):
	'''read swissmodel data, return a dict'''
	dbtypes = defaultdict(int)
	modelranges = {}
	print >> sys.stderr, "# reading data from {}".format(swissmodeldata), time.asctime()
	for line in opentype(uniprotdata,'r'):
		line = line.rstrip()
		if line and line[0]!="#": # ignore blank and comment lines
			lsplits = line.split("\t")
			accession = lsplits[0]
			database = lsplits[2]
			modelrange = [ int(n) for n in lsplits[3:5] ]

def accessions_to_swiss(accessiontable):
	'''read gene to accession table and return a dict where key is accession and value is the gene'''
	accessiondict = {}
	print >> sys.stderr, "# reading data from {}".format(accessiontable), time.asctime()
	for line in open(accessiontable,'r'):
		line = line.strip()
		if line and line[0]!="#": # ignore blank and comment lines
			lsplits = line.split("\t")
			geneid = lsplits[0]
			for accession in lsplits[1:]:
				accessiondict[accession] = geneid
	print >> sys.stderr, "# counted {} accessions".format( len(accessiondict) ), time.asctime()
	return accessiondict

def uniprot_from_alignment(alignfile, alignformat, pdbentries, protlengthdict, pdblengthdict, pdbspandict, lengthcutoff, commandfile):
	'''read alignment and extract information regarding coverage and swissprot ID'''
	print >> sys.stderr, "# Reading alignment from {}".format( alignfile )
	alignment = AlignIO.read( alignfile, alignformat )
	al_length = alignment.get_alignment_length()
	num_taxa = len(alignment)
	print >> sys.stderr, "# Alignment contains {} taxa for {} sites, including gaps".format( num_taxa, al_length )
	swissid = alignment[1].id.split("|")[2] # assumes format is sp|Q9BUQ8|DDX23_HUMAN
	if pdbentries.get(swissid, False):
		trimmedlength = len( str(alignment[0].seq).replace("-","") )
		trimmedspan = re.search("\w[\w-]+\w", str(alignment[0].seq) ).span()
		spanlength = int(trimmedspan[1]) - int(trimmedspan[0])
		for entry in pdbentries[swissid]:
			pdblength = int(pdblengthdict[entry][swissid])
			if pdblength < 100:
				print >> sys.stderr, "PDB {} LENGTH IS {}, REMOVING".format(entry, pdblength)
				continue
			chains = pdbspandict[swissid].keys()
			structspan = pdbspandict[swissid].values()
			print >> sys.stdout, "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format( swissid, os.path.basename(alignfile), trimmedlength, trimmedspan, spanlength, entry, protlengthdict[swissid], pdblength, chains, structspan )
			if commandfile:
				pdbfile = "PDBDIR/{}.pdb".format(entry.lower())
				pdbwscores = "{}_w_hp.pdb".format(pdbfile.split(".")[0])
				pdbargs = [ "~/git/pdbcolor/pdb_heteropecilly.py", "-a", alignfile, "-p", pdbfile, "-s", swissid , ">", pdbwscores ]
				with open(commandfile,'a') as cf:
					print >> cf, " ".join(pdbargs)

def main(argv, wayout):
	if not len(argv):
		argv.append('-h')
	parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=__doc__)
	parser.add_argument('-a','--alignments', nargs="*", help="alignment files or directory")
	parser.add_argument("-c","--commands", help="name for optional output file of pdb_heteropecilly commands")
	parser.add_argument('-f','--format', default="fasta", help="alignment format [fasta]")
	parser.add_argument('-l','--length-cutoff', default=0.5, type=float, help="minimum length coverage for proteins [0.5]")
	parser.add_argument('-s','--species', nargs="*", help="filter by one or more species")
	parser.add_argument("-u","--uniprot", help="Uniprot data in text format")
	parser.add_argument("-U","--filtered-uniprot", help="tabular information of filtered Uniprot data")
	args = parser.parse_args(argv)

	if args.uniprot:
		pdbentries, protlengthdict, pdblengthdict, pdbspandict = parse_uniprot(args.uniprot, args.species)
	elif args.filtered_uniprot:
		pdbentries, protlengthdict, pdblengthdict, pdbspandict = parse_filtered_uniprot(args.filtered_uniprot)

	if os.path.isdir(args.alignments[0]):
		print >> sys.stderr, "# Reading alignments files from directory {}".format(args.alignments[0]), time.asctime()
		globstring = "{}*.aln".format(args.alignments[0])
		alignmentfiles = glob(globstring)
	elif os.path.isfile(args.alignments[0]):
		alignmentfiles = args.alignments
	else:
		raise OSError("ERROR: Unknown alignments, exiting")

	if args.commands:
		print >> sys.stderr, "# writing pdb_heteropecilly.py commands to {}".format(args.commands)
		with open(args.commands,'w') as cf:
			print >> cf, "#!/bin/bash\n"

	alignwpdbcounter = 0

	headerline = "#swissProtID\talignFile\ttrimmedLength\tspan\tspanLength\tPDBentry\trefProtLength\tstructLength\tchain\tstructSpan"

	print >> sys.stdout, headerline

	for alignfile in alignmentfiles:
		uniprot_from_alignment(alignfile, args.format, pdbentries, protlengthdict, pdblengthdict, pdbspandict, args.length_cutoff, args.commands)









if __name__ == "__main__":
	main(sys.argv[1:], sys.stdout)
