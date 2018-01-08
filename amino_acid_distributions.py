#!/usr/bin/env python
#
# amino_acid_distributions.py

'''amino_acid_distributions.py  last updated 2017-11-21

amino_acid_distributions.py -a simion2017_alignment.aln > simion2017_aa_counts.tab

'''

import sys
import argparse
import time
import os
from collections import defaultdict,Counter
from Bio import SeqIO
from Bio import AlignIO

def alignment_aa_counts(alignfile, alignformat, aalist):
	print >> sys.stderr, "# Reading alignment from {}".format( alignfile )
	alignment = AlignIO.read( alignfile, alignformat )

	al_length = alignment.get_alignment_length()
	num_taxa = len(alignment)
	print >> sys.stderr, "# Alignment contains {} taxa for {} sites, including gaps".format( num_taxa, al_length )

	countsbyspecies = {} # key is species, value is Counter
	for seqrec in alignment:
		countsbyspecies[seqrec.id] = Counter( str(seqrec.seq) )

	holozoacounter = defaultdict(int)
	constcounter = defaultdict(int)

	print >> sys.stderr, "# Identifying constant sites from {}".format( alignfile ), time.asctime()
	for i in range(al_length):
		alignment_column = alignment[:,i] # all letters per site
		numgaps = alignment_column.count("-")
		nogap_alignment_column = alignment_column.replace("-","").replace("X","") # excluding gaps
		aa_counter = Counter( nogap_alignment_column )
		mostcommonaa = aa_counter.most_common(1)[0][0]

		# to account for sites that are constant in the 90 taxa, but not the 97
		holozoanumbers = [1 , 6 , 19 , 30 , 56 , 70 , 87] # seq number in alignment
		# for Abeoforma, Amoebidium, Capsaspora, Creolimax, Ministeria, Prium, Sphaerofirma
		aligncol_list = list(alignment_column)
		for sn in reversed(holozoanumbers):
			aligncol_list.pop(sn-1) # seq numbers start at 1, not 0, adjust python index
		filtaligncol = "".join(aligncol_list)
		nogap_filtaligncol = filtaligncol.replace("-","").replace("X","") # excluding gaps
		filtaa_counter = Counter( nogap_filtaligncol )

		if len(aa_counter)>1: # meaning site has more than 1 possible AA, so use HP
			if len(filtaa_counter)==1: # meaning filtered site is constant, but unfiltered is not
				holozoacounter[ filtaa_counter.most_common(1)[0][0] ] += 1
		else: # site is constant, mark as 11 so can be colored accordingly
			constcounter[ mostcommonaa ] += 1

	print >> sys.stderr, "# Counted {} constant sites".format( sum(constcounter.values()) ), time.asctime()

	# print number of constant sites for each AA
	print >> sys.stdout, "Species\t{}".format( "\t".join(aalist) )
	for k in countsbyspecies.iterkeys():
		print >> sys.stdout, "{}\t{}".format( k, "\t".join( [ str(countsbyspecies[k].get(x,0)) for x in aalist ]) )
	print >> sys.stdout, "Constants\t{}".format( "\t".join( [ str(constcounter[x]) for x in aalist] ) )
	print >> sys.stdout, "Holozoanconst\t{}".format( "\t".join( [ str(holozoacounter[x]) for x in aalist] ) )

def count_from_fasta(fastafile):
	aacounts = Counter()
	seqcounts = 0
	print >> sys.stderr, "# Reading proteins from {}".format( fastafile ), time.asctime()
	for seqrec in SeqIO.parse( fastafile, "fasta"):
		seqcounts += 1
		aacounts.update( Counter( str(seqrec.seq) ) )
	print >> sys.stderr, "# Counted {} sequences from {}".format( seqcounts, fastafile )
	return aacounts

def main(argv, wayout):
	if not len(argv):
		argv.append('-h')
	parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=__doc__)
	parser.add_argument("-a","--alignment", help="full multiple sequence alignment", required=True)
	parser.add_argument('-f','--format', default="fasta", help="alignment format [fasta]")
	parser.add_argument('-p','--protein-files', nargs="*", help="additional protein files in fasta format")
	args = parser.parse_args(argv)

	AALIST = list("ACDEFGHIKLMNPQRSTVWY-X") # this ignores U selenocysteine
                  # B      J   O    U X Z

	alignment_aa_counts( args.alignment, args.format, AALIST)

	if args.protein_files:
		for protfile in args.protein_files:
			aacounts = count_from_fasta(protfile)
			print >> sys.stdout, "{}\t{}".format( os.path.basename(protfile), "\t".join( [ str(aacounts[x]) for x in AALIST] ) )

if __name__ == "__main__":
	main(sys.argv[1:], sys.stdout)
