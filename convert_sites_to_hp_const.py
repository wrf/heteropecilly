#!/usr/bin/env python
#
# convert_sites_to_hp_const.py

'''convert_sites_to_hp_const.py  last updated 2017-10-30
convert tabular heteropecilly data to include constant alignment sites

convert_sites_to_hp_const.py -p hp_by_site.tab -a full_alignment.aln

    generate hp_by_site.tab using:

combine_heteropecilly.py -i *.bor > hp_by_site.tab
'''

import sys
import argparse
import time
from collections import defaultdict,Counter
from Bio import AlignIO

def read_tabular_hp(hptabular):
	'''read tabular heteropecilly results and return a dict where keys are position and value is HP'''
	hpdict = {}
	print >> sys.stderr, "# Reading heteropecilly by site from {}".format(hptabular), time.asctime()
	for line in open(hptabular,'r'):
		line = line.strip()
		if line and line[0]!="#": # ignore blank and comment lines
			pos, hp = line.split('\t')
			pos = int(pos)
			hpdict[pos] = hp
	print >> sys.stderr, "# Found heteropecilly for {} sites".format( len(hpdict) ), time.asctime()
	return hpdict

def read_raw_hp(rawhptabular):
	'''read raw heteropecilly results and return a dict where keys are position and value is lnPIP'''
	rawhpdict = {}
	print >> sys.stderr, "# Reading raw heteropecilly by site from {}".format(rawhptabular), time.asctime()
	for line in open(rawhptabular,'r'):
		line = line.strip()
		if line and line[0]!="#": # ignore blank and comment lines
			pos, pip, lnpip = line.split('\t')
			pos = int(pos)+1 # adjust position index
			rawhpdict[pos] = lnpip[0:5] # keep lnPIP as string, take 2 or 3 decimal places
	print >> sys.stderr, "# Found heteropecilly for {} sites".format( len(rawhpdict) ), time.asctime()
	return rawhpdict

def full_alignment_hp(alignfile, alignformat, hpgrpdict, rawhpdict, makefasta):
	print >> sys.stderr, "# Reading alignment from {}".format( alignfile )
	alignment = AlignIO.read( alignfile, alignformat )

	al_length = alignment.get_alignment_length()
	num_taxa = len(alignment)

	constcounter = defaultdict(int)
	print >> sys.stderr, "# Alignment contains {} taxa for {} sites, including gaps".format( num_taxa, al_length )

	index_to_hp = {} # key is full alignment position, value is HP or 10 for constant
	gapsbysite = {} # count of gaps, just to check if variable sites tend to be low occupancy

	constsiteoffset = 0 # initially trimmed and untrimmed begin at same site

	if not makefasta: # make header line
		print >> sys.stdout, "#site\tgroup\tlnPIP\tnumTaxa\tgaps\tmostCommon\tfreq\totherCount"

	for i in range(al_length):
		alignment_column = alignment[:,i] # all letters per site
		numgaps = alignment_column.count("-")
		gapsbysite[i+1] = numgaps
		nogap_alignment_column = alignment_column.replace("-","").replace("X","") # excluding gaps
		num_nogap_taxa = len(nogap_alignment_column)
		aa_counter = Counter( nogap_alignment_column )
		mostcommonaa = aa_counter.most_common(1)[0][0]
		lnPIP = "NA" # set default as NA, otherwise reassigned to value from dict

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
				if makefasta: # use lowercase c
					index_to_hp[i+1] = "c"
				else:
					index_to_hp[i+1] = 10
				constcounter[ filtaa_counter.most_common(1)[0][0] ] += 1
				constsiteoffset -= 1 # decrement by 1, since full alignment is now ahead of trimmed
			else:
				# values should be integers from 0 to 9
				index_to_hp[i+1] = int(hpgrpdict.get(i+1+constsiteoffset,0))/10
				lnPIP = rawhpdict.get(i+1+constsiteoffset,0)
		else: # site is constant, mark as 11 so can be colored accordingly
			if makefasta:
				index_to_hp[i+1] = "C"
			else:
				index_to_hp[i+1] = 11

			constcounter[ mostcommonaa ] += 1
			constsiteoffset -= 1 # decrement by 1, since full alignment is now ahead of trimmed
		if not makefasta:
			mostcommoncount = aa_counter.most_common(1)[0][1]
			print >> sys.stdout, "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format( i+1, index_to_hp[i+1], lnPIP, num_nogap_taxa, numgaps, mostcommonaa, mostcommoncount, num_nogap_taxa-mostcommoncount )

	print >> sys.stderr, "# Assigned values for {} sites".format( len(index_to_hp) )
	totalconstsites = sum(constcounter.values())
	print >> sys.stderr, "# Counted {} constant sites".format( totalconstsites )

	# check for sites constant in choanozoa but not holozoa
	if len(hpgrpdict) + totalconstsites != al_length:
		missingdatacount = al_length - totalconstsites - len(hpgrpdict)
		print >> sys.stderr, "# WARNING COULD NOT ACCOUNT FOR {} SITES".format( missingdatacount )

	# print number of constant sites for each AA
	for k,v in constcounter.iteritems():
		print >> sys.stderr, "{},{}".format(k,v)

	return index_to_hp, gapsbysite

def main(argv, wayout):
	if not len(argv):
		argv.append('-h')
	parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=__doc__)
	parser.add_argument("-a","--alignment", help="full multiple sequence alignment", required=True)
	parser.add_argument('-f','--format', default="fasta", help="alignment format [fasta]")
	parser.add_argument('--fasta', action="store_true", help="print output as a fasta format line for alignments")
	parser.add_argument('-p','--heteropecilly', help="tabular heteropecilly data file from combine_heteropecilly.py")
	parser.add_argument('-P','--raw-PIP', help="tabular raw heteropecilly data file")
	args = parser.parse_args(argv)

	hpgroupsbysite = read_tabular_hp(args.heteropecilly)
	lnhpdict = read_raw_hp(args.raw_PIP) if args.raw_PIP else {} # empty dict so default can be 0
	hp_and_const_index, gapcounts = full_alignment_hp( args.alignment, args.format, hpgroupsbysite, lnhpdict, args.fasta )

	if args.fasta:
		hpstring = "".join( str(hp_and_const_index[k]) for k in sorted(hp_and_const_index.keys()) )
		print >> sys.stdout, ">Heteropecilly_score\n{}".format( hpstring )

if __name__ == "__main__":
	main(sys.argv[1:], sys.stdout)
