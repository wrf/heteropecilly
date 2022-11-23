#!/usr/bin/env python
#
# parse_swissprot_data.py  created 2017-09-29

'''parse_swissprot_data.py  last modified 2017-11-01

parse_swissprot_data.py -u uniprot_sprot.dat.gz
'''

import sys
import gzip
import argparse
import time
import re
from collections import defaultdict

def parse_uniprot(uniprotdata, haspdb, speciesfilter, writeaccessions):
	if uniprotdata.rsplit('.',1)[1]=="gz": # autodetect gzip format
		opentype = gzip.open
		print >> sys.stderr, "# reading data from {} as gzipped".format(uniprotdata), time.asctime()
	else: # otherwise assume normal open
		opentype = open
		print >> sys.stderr, "# reading data from {}".format(uniprotdata), time.asctime()

	typecounter = defaultdict(int)
	speciescounter = defaultdict(int)
	pdbentries = defaultdict(list)
	genetoaccession = defaultdict(list)
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
				elif key=="AC":
					acsplits = [ac.strip() for ac in line[5:].split(";") if ac.strip()]
					genetoaccession[protid].extend(acsplits)
				elif key=="OX":
					speciesnum = re.search("NCBI_TaxID=(\d+)",line[5:]).group(1)
					#species = line[5:].split(";")[0].split("=")[1]
					speciescounter[speciesnum] += 1
					if speciesnum not in speciesfilter:
						skipentry = True
						continue
				elif key=="DR": # for database references
					drsplits = [dr.strip() for dr in line[5:].split(";")]
					if drsplits[0]=="PDB":
						pdbid = drsplits[1]
						pdbentries[protid].append(pdbid)
						pdbrangestring = drsplits[4].replace(".","")
						if pdbrangestring=="-": # some entries are blank, for some reason
							# DR   PDB; 1GQ5; X-ray; 2.20 A; -.
							continue
						try:
							pdbrangelist = [ int(n) for n in re.search("=(\d+)-(\d+)", pdbrangestring).groups() ]
						except AttributeError:
							print >> sys.stderr, "WARNING CANNOT PARSE {}".format(line)
							continue
						pdbcoverage = pdbrangelist[1]-pdbrangelist[0]+1
						protcoverage = 1.0 * pdbcoverage / protlength
						print >> sys.stdout, "{}\t{}\t{}\t{}\t{}\t{:.3f}".format(protid, pdbid, protlength, pdbrangestring, pdbcoverage, protcoverage)
	print >> sys.stderr, "# reading data from {}".format(uniprotdata), time.asctime()
	for idtype in sorted(typecounter.keys()):
		print >> sys.stderr, "{}\t{}".format(idtype, typecounter[idtype])
	print >> sys.stderr, "# counted {} proteins with PDB entries".format( len(pdbentries) ), time.asctime()
	#for k,v in pdbentries.iteritems():
	#	print >> sys.stderr, "{}\t{}".format(k, len(v) )

	if writeaccessions:
		with open("gene_to_accession.tab",'w') as gt:
			for k,v in genetoaccession.iteritems():
				print >> gt, "{}\t{}".format(k, "\t".join(v) )
	#for k,v in speciescounter.iteritems():
	#	print >> sys.stderr, "{}\t{}".format(k, v )

def main(argv, wayout):
	if not len(argv):
		argv.append('-h')
	parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=__doc__)
	parser.add_argument("-u","--uniprot", help="Uniprot data in text format", required=True)
	parser.add_argument('-f','--format', default="fasta", help="alignment format [fasta]")
	parser.add_argument('-p','--has-pdb', action="store_true", help="require that entry has a PDB structure")
	parser.add_argument('-s','--species', nargs="*", help="filter by one or more species")
	parser.add_argument('-w','--write-accessions', action="store_true", help="generate gene to accession table")
	args = parser.parse_args(argv)

	parse_uniprot(args.uniprot, args.has_pdb, args.species, args.write_accessions)

if __name__ == "__main__":
	main(sys.argv[1:], sys.stdout)
