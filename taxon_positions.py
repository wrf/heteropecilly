#!/usr/bin/env python
#
# taxa_from_tree.py  v1.0 created 2017-08-21

'''print taxa in order from a nexus-format tree:

taxa_from_tree.py partial_tree.nex
'''

import sys
from itertools import chain
from Bio import Phylo

# trees in nexus format
# #NEXUS
#begin trees;
#	tree tree_1 = [&R] ((AGAP008544-PA:0.02348,AGAP028140-PA:0.0908)[&bs=1.0]:0.34487,(AGAP028059-#PA:0.07144,AGAP008545-PA:0.08111)[&bs=1.0]:0.34047)[&bs=1.0];
#end;

if len(sys.argv) < 2:
	sys.exit(__doc__)
else:
	tree = Phylo.read(sys.argv[1],"nexus")
	target = sys.argv[2]

	Phylo.draw_ascii(tree)

	# list of lists, where each sublist contains the taxa to add for each node
	taxonorder = []
	# assuming the target is taxonA, should print like
	# [['taxonA'], ['taxonB'], ['taxonC', 'taxonD'], ['taxonE', 'taxonF', 'taxonG', 'taxonH']]
	treetaxa = {}
	for clade in chain( reversed(tree.get_path(target)) , tree.root):
		newtaxaperclade = []
		print clade.depths()
		for termclades in clade.get_terminals():
			if not treetaxa.get(termclades.name,False):
				print termclades.name, termclades.branch_length, termclades.depths()
				treetaxa[termclades.name] = True
				newtaxaperclade.append(termclades.name)
		if newtaxaperclade:
			taxonorder.append(newtaxaperclade)
	print taxonorder
#	for clade in tree.get_terminals():
#		print >> sys.stdout, str(clade.name).replace("'","")
