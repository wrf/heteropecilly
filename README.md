# heteropecilly
in depth analysis of heteropecillious sites in supermatrices

## What is heteropecilly? ##
Homopecilly and heteropecilly are concepts proposed by [Roure and Philippe (2011)](http://www.ncbi.nlm.nih.gov/pubmed/21235782), suggesting that the substitution matrix for each site may change uniquely in specific clades over evolutionary time. For instance, a position might accept only Asp and Glu in one clade, and only Asp and Asn in another.

![heteropecilly_example_v1.png](https://github.com/wrf/heteropecilly/blob/master/heteropecilly_example_v1.png)

This requires binning taxa into groups, and comparing the variation within each site between groups. For this reason, tree topology cannot be tested (such as monophyly/paraphyly of specific groups), as the binning already assumes a topology.

## combine_heteropecilly.py ##
The raw analysis from [Simion 2017](https://github.com/psimion/SuppData_Metazoa_2017) contained two sets of files: a table of the raw/log scores calculated from the supermatrix alignment, and files containing intervals for removing deciles of most heteropecillious sites. For instance, the file `RemPIP90.bor` contains the interval for eliminating 90% most heteropecillous positions, so creating the alignment of the 10% most homopecillous positions. The heteropecilly score of sites found only in `RemPIP10.bor` should be max (9 of 9, the highest decile) while the score of any sites not included in any file (so finally not in `RemPIP90.bor`) should be 0 of 9.

`./combine_heteropecilly.py -i heteropecilly-v2/*bor -f -l 341543 > heteropecilly-v2/hp_by_site.fasta`

The script has two output options, one is a table of values by site, the other is a fasta file of one line (changed with `-f`), where numbers 0-9 are assigned for each position.

![PAUL-90x341543v-C20.lnPIP.png](https://github.com/wrf/heteropecilly/blob/master/PAUL-90x341543v-C20.lnPIP.png)

## convert_sites_to_hp_const.py ##
Appx. 60k constant sites were removed from the original alignment (which was `supermatrix_97sp_401632pos_1719genes.fasta`). This consists of 3447 sites which were constant only in choanozoa and 56642 sites from all species. 7 protist taxa were removed for the heteropecilly analysis.

`./convert_sites_to_hp_const.py -a simion2017_97sp_401632pos_1719genes.fasta -p heteropecilly-v2/hp_by_site.tab --fasta > heteropecilly-v2/hp_by_site_w_const.fasta`

The script takes the tabular data from the above script, and adds the constant sites. Sites from all species are indicated by `C` while those only from choanozoa are shown with a lowercase `c`. The fasta output can be added as a row in the alignment:

`cat supermatrix_97sp_401632pos_1719genes.fasta simion2017_hp_by_site_w_const.fasta > supermatrix_97sp_401632pos_1719genes_w_hp_const.fasta`

## blast_to_align_pairs.py ##
Because of trimming steps, most proteins in a supermatrix do not represent the entire, or even the majority, of the original protein, and the identity of this protein (say the name of a gene) may be unknown. For cases where human was used, the IDs of the human proteins can be extract with `blastp`, as even trimmed proteins will have the top hit to a real human protein with almost 100% identity. Thus, individual alignments of each trimmed protein can be remade with the reference protein.

Human proteins can be extracted from the [SwissProt set](ftp://ftp.expasy.org/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz). Because of the standard naming scheme of Uniprot proteins, the `getAinB.py` script extracts all proteins with the species tag *_HUMAN*, creating a new file of only human proteins.

`getAinB.py _HUMAN uniprot_sprot.fasta -s > human_uniprot.fasta`

Then, generate a file of all human proteins used in the supermatrix. This script `split_supermatrix_to_taxa.py` can be found in the [supermatrix repo](https://github.com/wrf/supermatrix).

`split_supermatrix_to_taxa.py -a simion2017_97sp_401632pos_1719genes.fasta.gz -p simion2017_partitions.txt -d simion2017_taxa`

Because `blastp` does not allow gaps in sequences, and most sequences extracted from the supermatrix will still have a gap, gaps need to be removed. The `degapper.py` script removes gaps and automatically renames the output file.

`degapper.py simion2017_taxa/Homo_sapiens.fasta`

Make the blast protein database with `makeblastdb` and then run `blastp` and report the results as tabular (`-outfmt 6`).

`blastp -query simion2017_taxa/Homo_sapiens.fasta -db human_uniprot.fasta -outfmt 6 -max_target_seqs 1 > simion2017_taxa/hsapiens_vs_uniprot_blastp.tab`

The top hit for each protein should probably be the original reference protein. If this is not the case, then that protein probably should not be used in phylogeny in the first place. Using the blast results, new alignments can be generated for each protein from the supermatrix and the reference protein. Each file is in fasta format and contains three sequences, the protein used in the supermatrix, the reference protein, and a row of the extracted heteropecilly scores from the supermatrix. The folder containing all of the files is automatically named (based on current time).

`blast_to_align_pairs.py -b simion2017_taxa/hsapiens_vs_uniprot_blastp.tab -q simion2017_taxa/Homo_sapiens.fasta.nogaps -s human_uniprot.fasta -r simion2017_taxa/Homo_sapiens.fasta -p hp_by_site_w_const.tab`

This requires `mafft`, though could be modified to run another aligner.

## alignment_pair_stats.py ##
After using the above scripts, summary information about each alignment can be determined and printed into a tabular format.

`./alignment_pair_stats.py -a blast_alignments_20171009-153805/ > simion2017_align_pair_stats.tab`

Each row contains 8 columns:

* partition - partition from the supermatrix
* protID - blast hit, in this case the Uniprot ID
* trimmedLength - number of non-gap characters in this partition
* trimmedPercent - percent relative to the length of the reference protein
* span - first non-gap character to last non-gap character, though may still contain internal gaps
* spanLength - length of the span in amino acids
* spanPercent - percent relative to the length of the reference protein
* refProtLength - length of reference protein

`Homo_sapiens_81364-81595  sp|P42858|HD_HUMAN  232  0.074  (129, 2455)  2326  0.740  3142`

These results can be summarized between multiple supermatrices. It is clear that for many proteins, most of the sequence was trimmed, likely leaving just a single domain. The orthology of these proteins is obviously in question. This may generally be a consequence of the reliance on transcriptomic data, although the Borowiec 2015 study made use of only genomes yet many proteins are still likely incomplete, possibly due to assembly or annotation problems.

![align_pair_stats_simion-whelan-borowiec.png](https://github.com/wrf/heteropecilly/blob/master/align_pair_stats_simion-whelan-borowiec.png)

