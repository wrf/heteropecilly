# heteropecilly
in depth analysis of heteropecillious sites in supermatrices

## What is heteropecilly? ##
Homopecilly and heteropecilly are concepts proposed by Roure and Philippe (2011)[http://www.ncbi.nlm.nih.gov/pubmed/21235782], suggesting that the substitution matrix for each site may change uniquely in specific clades over evolutionary time. For instance, a position might accept only Asp and Glu in one clade, and only Asp and Asn in another.

This requires binning taxa into clades, and comparing the variation within each site between clades. For this reason, tree topology cannot be tested (such as monophyly/paraphyly of specific groups), as the binning already assumes a topology.

## combine_heteropecilly.py ##
The raw analysis from [Simion 2017](https://github.com/psimion/SuppData_Metazoa_2017) contained two sets of files: a table of the raw/log scores calculated from the supermatrix alignment, and files containing intervals for removing deciles of most heteropecillious sites. For instance, the file `RemPIP90.bor` contains the interval for eliminating 90% most heteropecillous positions, so creating the alignment of the 10% most homopecillous positions. The heteropecilly score of sites found only in `RemPIP10.bor` should be max (9 of 9, the highest decile) while the score of any sites not included in any file (so finally not in `RemPIP90.bor`) should be 0 of 9.

The script has two output options, one is a table of values by site, the other is a fasta file of one line, where numbers 0-9 are assigned for each position.

## convert_sites_to_hp_const.py ##
Appx. 60k constant sites were removed from the original alignment (which was `supermatrix_97sp_401632pos_1719genes.fasta`). This consists of 3447 sites which were constant only in choanozoa and 56642 sites from all species. 7 protist taxa were removed for the heteropecilly analysis.

The script takes the tabular data from the above script, and adds the constant sites. Sites from all species are indicated by `C` while those only from choanozoa are shown with a lowercase `c`. The fasta output can be added as a row in the alignment:

`cat supermatrix_97sp_401632pos_1719genes.fasta simion2017_hp_by_site_w_const.fasta > supermatrix_97sp_401632pos_1719genes_w_hp_const.fasta`
