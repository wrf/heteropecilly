#!/usr/bin/env python
#
# blast_to_align_pairs.py

'''blast_to_align_pairs.py  last modified 2017-10-16

blast_to_align_pairs.py -q query_prots.fasta -s prot_db.fasta -b blastp.tab -d pair_dir -r Homo_sapiens.fasta

    make tabular blast output from:

blastp -query query_prots.fasta -db prot_db.fasta -outfmt 6 -max_target_seqs 1 > blastp.tab

    alignment reference proteins (still aligned) can be extracted with:

split_supermatrix_to_taxa.py -a matrix.phy -p partitions.txt -d taxa -f phylip-relaxed
'''

import sys
import os
import argparse
import time
import subprocess
from Bio import SeqIO,AlignIO
from Bio.Seq import Seq

def read_tabular_hp(hptabular):
	'''read tabular heteropecilly results and return a dict where keys are position and value is HP'''
	hpdict = {}
	print >> sys.stderr, "# Reading heteropecilly by site from {}".format(hptabular), time.asctime()
	for line in open(hptabular,'r'):
		line = line.strip()
		if line and line[0]!="#": # ignore blank and comment lines
			pos, hp = line.split('\t')[0:2] # only take first two columns
			pos = int(pos)
			if hp == "10":
				hp = "c" # redefine constant sites as C instead of 10 for string
			elif hp == "11":
				hp = "C"
			hpdict[pos] = hp
	print >> sys.stderr, "# Found heteropecilly for {} sites".format( len(hpdict) ), time.asctime()
	return hpdict

def run_mafft(MAFFT, rawseqsfile):
	'''generate multiple sequence alignment from fasta and return MSA filename'''
	aln_output = "{}.aln".format(os.path.splitext(rawseqsfile)[0] )
	aligner_args = [MAFFT, "--auto", "--quiet", rawseqsfile]
	print >> sys.stderr, "{}\n{}".format(time.asctime(), " ".join(aligner_args) )
	with open(aln_output, 'w') as msa:
		subprocess.call(aligner_args, stdout=msa)
	print >> sys.stderr, "# alignment of {} completed".format(aln_output), time.asctime()
	if os.path.isfile(aln_output):
		return aln_output
	else:
		raise OSError("Cannot find expected output file {}".format(aln_output) )

def make_heteropecilly_string(aln_file, hpdict, refprotdict):
	'''from partial alignment to full protein, return string of heteropecilly for partial alignment'''
	print >> sys.stderr, "# Reading alignment from {}".format( aln_file )
	alignment = AlignIO.read( aln_file, "fasta" )
	al_length = alignment.get_alignment_length()
	num_taxa = len(alignment)
	print >> sys.stderr, "# Alignment contains {} taxa for {} sites, including gaps".format( num_taxa, al_length )
	rangeindex = [ int(x) for x in os.path.basename(aln_file).split("-")[0:2] ]

	# remove gapped sites from heteropecilly based on aligned and trimmed protein
	hp_no_gaps = ""
	hpoffset = 0
	refid = alignment[0].id # should always be first seq in alignment
	for i,site in enumerate(str(refprotdict[refid].seq)):
		if site!="-": # if site is not gap, keep heteropecilly score
			adjindex = rangeindex[0] + i # + hpoffset
			hp_no_gaps += hpdict[adjindex]
	#	else: # increase offset for each gap
		#	hpoffset += 1

	# assign heteropecilly score to columns without gaps in alignment
	aligned_hpstring = ""
	refoffset = 0
	for i in range(al_length):
		alignment_column = alignment[:,i] # all letters per site
		if alignment_column[0]=="-": # if first character is a gap, assign gap for heteropecilly
			aligned_hpstring += "-"
		else: # otherwise use value from hp_no_gaps
			aligned_hpstring += hp_no_gaps[refoffset]
			refoffset += 1
	return aligned_hpstring

def make_pairs_from_blast(blastfile, querydict, subjectdict, new_aln_dir, mafftbin, hpdict, refdict ):
	'''iterate through blast hits, generate fasta files of each query subject pair and align them'''
	qspairs = {} # key is query, value is subject, assuming only one hit but multiple HSPs
	linecounter = 0
	print >> sys.stderr, "# parsing blast hits from {}".format(blastfile), time.asctime()
	for line in open(blastfile,'r'):
		line = line.strip()
		if line and line[0]!="#": # skip empty and comment lines
			linecounter += 1
			lsplits = line.split("\t")
			query = lsplits[0]
			subject = lsplits[1]
			qspairs[query] = subject
	print >> sys.stderr, "# counted {} blast hits for {} query proteins".format( linecounter, len(qspairs) ), time.asctime()
	filecounter = 0
	for q,s in qspairs.iteritems():
		fasta_filename = "{}-{}.fasta".format(q.split("_")[-1],s.split("|")[-1])
		fasta_output = os.path.join(new_aln_dir,fasta_filename)
		with open(fasta_output,'w') as fo:
			fo.write( querydict[q].format("fasta") )
			if str(subjectdict[s].seq).count("U"): # check for disallowed selenocysteine
				subjectdict[s].seq = Seq(str(subjectdict[s].seq).replace("U","C"))
			fo.write( subjectdict[s].format("fasta") )
		aln_file = run_mafft(mafftbin, fasta_output)
		if aln_file: # check if file was made
			filecounter += 1
			if refdict and hpdict: # if alignment proteins and heteropecilly are given
				hpstring = make_heteropecilly_string(aln_file, hpdict, refdict)
				with open(aln_file, 'a') as af:
					print >> af, ">Heteropecilly_score\n{}".format( hpstring )
	print >> sys.stderr, "# generated {} alignments".format(filecounter), time.asctime()

def main(argv, wayout):
	if not len(argv):
		argv.append('-h')
	parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=__doc__)
	parser.add_argument('-b','--blast', help="tabular blastp results")
	parser.add_argument('-d','--directory', default="./blast_alignments", help="directory for new alignments [autonamed]")
	parser.add_argument('--mafft', default="mafft", help="path to mafft binary [default is in PATH]")
	parser.add_argument('-q','--query', help="fasta file of blast query proteins from supermatrix")
	parser.add_argument('-s','--subject', help="fasta file of blast database proteins; this was the -db file for blastp")
	parser.add_argument('-r','--reference', help="fasta file of proteins from the alignment with gaps, only needed if heteropecilly is computed")
	parser.add_argument('-p','--heteropecilly', help="tabular heteropecilly data file")
	args = parser.parse_args(argv)

	### DIRECTORY FOR NEW OUTPUT
	new_aln_dir = os.path.abspath("{}_{}".format(args.directory, time.strftime("%Y%m%d-%H%M%S") ) )
	if not os.path.exists(new_aln_dir):
		os.mkdir(new_aln_dir)
		print >> sys.stderr, "# Creating directory {}".format(new_aln_dir), time.asctime()
	elif os.path.isdir(new_aln_dir):
		print >> sys.stderr, "# Using directory {}".format(new_aln_dir), time.asctime()

	print >> sys.stderr, "# Reading sequences from {}".format(args.query), time.asctime()
	querydict = SeqIO.to_dict(SeqIO.parse(args.query,"fasta"))
	print >> sys.stderr, "# Counted {} sequences".format( len(querydict) ), time.asctime()
	print >> sys.stderr, "# Reading sequences from {}".format(args.subject), time.asctime()
	subjectdict = SeqIO.to_dict(SeqIO.parse(args.subject,"fasta"))
	print >> sys.stderr, "# Counted {} sequences".format( len(subjectdict) ), time.asctime()
	if args.reference:
		print >> sys.stderr, "# Reading sequences from {}".format(args.reference), time.asctime()
		refdict = SeqIO.to_dict(SeqIO.parse(args.reference,"fasta"))
		print >> sys.stderr, "# Counted {} sequences".format( len(refdict) ), time.asctime()
	else:
		refdict = None

	hpbysite = read_tabular_hp(args.heteropecilly) if args.heteropecilly else None

	make_pairs_from_blast(args.blast, querydict, subjectdict, new_aln_dir, args.mafft, hpbysite, refdict)

if __name__ == "__main__":
	main(sys.argv[1:], sys.stdout)
