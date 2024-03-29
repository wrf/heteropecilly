# heteropecilly
**AS OF NOV 2022: this repo will probably not be updated again and no new data or analysis will be added**

some in depth analysis of heteropecillious sites in supermatrices

## What is heteropecilly? ##
Homopecilly and heteropecilly are concepts proposed by [Roure and Philippe (2011)](http://www.ncbi.nlm.nih.gov/pubmed/21235782), suggesting that the substitution matrix for each site may change uniquely in specific clades over evolutionary time. For instance, a position might accept only Asp and Glu in one clade, and only Asp and Asn in another.

![heteropecilly_example_v1.png](https://github.com/wrf/heteropecilly/blob/master/heteropecilly_example_v1.png)

This requires binning taxa into groups, and comparing the variation within each site between groups. For this reason, tree topology cannot be tested (such as monophyly/paraphyly of specific groups), as the binning already assumes a topology.

## combine_heteropecilly ##
The raw analysis from [Simion 2017](https://github.com/psimion/SuppData_Metazoa_2017) contained two sets of files: a table of the raw/log scores calculated from the supermatrix alignment, and files containing intervals for removing deciles of most heteropecillious sites. For instance, the file `RemPIP90.bor` contains the interval for eliminating 90% most heteropecillous positions, so creating the alignment of the 10% most homopecillous positions. The heteropecilly score of sites found only in `RemPIP10.bor` should be max (9 of 9, the highest decile) while the score of any sites not included in any file (so finally not in `RemPIP90.bor`) should be 0 of 9.

`./combine_heteropecilly.py -i heteropecilly-v2/*bor -f -l 341543 > heteropecilly-v2/hp_by_site.fasta`

The script has two output options, one is a table of values by site, the other is a fasta file of one line (changed with `-f`), where numbers 0-9 are assigned for each position.

![PAUL-90x341543v-C20.lnPIP.png](https://github.com/wrf/heteropecilly/blob/master/PAUL-90x341543v-C20.lnPIP.png)

## convert_sites_to_hp_const ##
Appx. 60k constant sites were removed from the original alignment (which was `supermatrix_97sp_401632pos_1719genes.fasta`). This consists of 3447 sites which were constant only in choanozoa and 56642 sites from all species. 7 protist taxa were removed for the heteropecilly analysis.

`./convert_sites_to_hp_const.py -a simion2017_97sp_401632pos_1719genes.fasta -p heteropecilly-v2/hp_by_site.tab --fasta > heteropecilly-v2/hp_by_site_w_const.fasta`

The script takes the tabular data from the above script, and adds the constant sites. Sites from all species are indicated by `C` while those only from choanozoa are shown with a lowercase `c`. The fasta output can be added as a row in the alignment:

`cat supermatrix_97sp_401632pos_1719genes.fasta simion2017_hp_by_site_w_const.fasta > supermatrix_97sp_401632pos_1719genes_w_hp_const.fasta`

A more complex table including the original lnPIP scores can be generated with another command:

`convert_sites_to_hp_const.py -a simion2017_97sp_401632pos_1719genes.fasta -p heteropecilly-v2/hp_by_site.tab -P heteropecilly-v2/PAUL-90x341543v-C20.proba.tab > heteropecilly-v2/hp_by_site_w_const_20CAT.tab`

Each row contains 8 columns:

* site - number of site, starting with 1
* group - heteropecilly decile (0-9, 9 being most heteropecillious), 10 or 11 if constant
* lnPIP - log score from `PAUL-90x341543v-C20.proba.tab`, NA if 0
* numTaxa - number of taxa with data
* gaps - number of gaps, should be total minus numTaxa
* mostCommon - most common amino acid at that site
* freq - number of times the most common amino acid was found
* otherCount - total count of other amino acids, should be numTaxa minus freq

## sitewise_ll_to_columns ##
Because not all sites provide the same phylogenetic information for all clades (i.e. sites that are constant except in vertebrates would not affect relationships of other taxa), it is necessary to get information about which sites strongly support certain hypotheses of relationships. This can be done in [RAxML](https://sco.h-its.org/exelixis/web/software/raxml/index.html), using the `-f G` option to calculate site-wise likelihoods (the approach used by [Weigert et al 2014](https://academic.oup.com/mbe/article/31/6/1391/1009370) and [Shen et al 2017](https://www.nature.com/articles/s41559-017-0126) ).

`raxmlHPC-PTHREADS-SSE3-8.2.11 -f G -s simion2017_97sp_401632pos_1719genes.phy -m PROTGAMMALG -z tree_97sp_CAT.rooted_combined.tre -n simion2017_97sp_401632pos_1719genes -T 6`

The output of this can be converted to columns, where it can be more easily used by other programs:

`./sitewise_ll_to_columns.py RAxML_perSiteLLs.simion2017_97sp_401632pos_1719genes > RAxML_perSiteLLs.simion2017_97sp_401632pos_1719genes.tab`

While this makes use a different model, this can nonetheless be directly compared to the heteropecilly scores to give some sense of which sites may provide the most information with one model or the other. The two can be plotted together, and colors correspond to the 10 deciles used in the original analysis. Although ML can be calculated for constant sites (which differs between trees, possibly due to gap positions), these sites have no heteropecilly score and are excluded.

![simion2017_hp_vs_absdln_plot.png](https://github.com/wrf/heteropecilly/blob/master/simion2017_hp_vs_absdln_plot.png)

When this is instead clustered as a barplot in the deciles, the pattern becomes more obvious: heteropecillious sites account for substantially more information than homopecillious sites. The 30% most heteropecillious sites account for more phylogenetic information than the rest combined (26598 vs. 24910 dlnL). This includes both sites favoring porifera-sister and ctenophora-sister. Thus, it is not surprising that removal of these sites leads to an alternate topology being favored (coelenterata, see supplement of [Simion et al](http://www.sciencedirect.com/science/article/pii/S0960982217301999) ). Removing sites purely because the model is inadequate is a questionable practice in general, but from these results suggesting that these sites include the most information, it is also probably unwise to remove these sites.

![simion2017_hp_vs_dln_boxplot.png](https://github.com/wrf/heteropecilly/blob/master/simion2017_hp_vs_dln_boxplot.png)

Curiously, the sites that are determined by `RAxML` to have the most information are actually the sites that cause the most problems for the CAT model `phylobayes`. It could very well be that the differing favored topologies occur because the two programs leverage information from totally different sites. If heteropecillious sites violate certain model assumptions of the CAT model (that changes are uniform across a site), then `phylobayes` can best handle sites that change often and neutrally within a biochemical category. Conversely, heteropecillious sites contain more information that can be used by `RAxML`, perhaps because they reflect rare and lineage-specific changes, such as changes between biochemical categories like L to R.

If it were possible to get site-wise probability information from `phylobayes`, then one could directly examine which sites differentially account for the different results between the two programs. Without such information, any argument favoring one or the other program, model or topology is essentially incomplete, and offers very little understanding or predictive power.

## sitewise_get_strong_sites ##
The vast majority of sites do not weigh in on any given node of a tree. This makes sense as not all sites should be expected to be equally informative to all levels of a tree. Fast-evolving sites may be relevant for narrower time scales, genus- or family level, while slower evolving sites would be relevant for phylum level distinctions. Thus, sites critical for resolving a particular node can be extracted using the `sitewise_get_strong_sites.py` script. By default, sites where the absolute value of the difference is greater than 0.5 are kept. For the Simion dataset, this yields 8106 sites out of 344990 potential sites.

`sitewise_get_strong_sites.py -a simion2017_97sp_401632pos_1719genes.aln -l RAxML_perSiteLLs.simion2017_97sp_401632pos_1719genes.tab -o simion2017_97sp_perSiteLLs_strong.aln`

As it may also be useful to identify which sites favor which of the two topologies, a fasta string can be generated where each site is the calculated difference between the two topologies. Because differences span a range of around -8 to +8, this can be converted to a single character in base16. Thus, a score of 8.0 would mean the likelihood is the same for both, all values above 8.0 (meaning 8.1 - f) favor the first tree, while all values below 8.0 (7.9 - 0) favor the second tree.

``sitewise_columns_to_fasta.py RAxML_perSiteLLs.simion2017_97sp_401632pos_1719genes.tab > simion2017_97sp_perSiteLLs.fasta``

The fasta string can be combined with the original alignment, and sites can be extracted as above.

``cat simion2017_97sp_401632pos_1719genes.aln simion2017_97sp_perSiteLLs.fasta > simion2017_97sp_strong_w_lnl.aln``

## sitewise_recode_constant_sites ##
Constant sites should have no difference between lnL of any topology, yet they do (probably due to which species have gaps/missing data). These sites can specifically be recoded (as `const` or `x`) so that they are excluded from some future analyses.

`sitewise_recode_constant_sites.py -a simion2017_97sp_401632pos_1719genes.aln -l RAxML_perSiteLLs.simion2017_97sp_401632pos_1719genes.tab > RAxML_perSiteLLs.simion2017_const_recoded.tab`

The `const_recoded.tab` file can be used for `sitewise_columns_to_fasta.py` or `blast_to_align_pairs.py`.

## blast_to_align_pairs ##
Because of trimming steps, most proteins in a supermatrix do not represent the entire, or even the majority, of the original protein, and the identity of this protein (say the name of a gene) may be unknown. For cases where human was used, the IDs of the human proteins can be extract with `blastp`, as even trimmed proteins will have the top hit to a real human protein with almost 100% identity. Thus, individual alignments of each trimmed protein can be remade with the reference protein.

Because human is being used as the reference set, all human proteins must be extracted from the [SwissProt set](http://www.uniprot.org/downloads). Because of the standard naming scheme of Uniprot proteins, the `getAinB.py` script extracts all proteins with the species tag *_HUMAN*, creating a new file of only human proteins.

`getAinB.py _HUMAN uniprot_sprot.fasta -s > human_uniprot.fasta`

Then, generate a file of all human proteins used in the supermatrix. This script `split_supermatrix_to_taxa.py` can be found in the [supermatrix repo](https://github.com/wrf/supermatrix).

`split_supermatrix_to_taxa.py -a simion2017_97sp_401632pos_1719genes.fasta.gz -p simion2017_partitions.txt -d simion2017_taxa`

Because `blastp` does not allow gaps in sequences, and most sequences extracted from the supermatrix will still have a gap, gaps need to be removed. The `degapper.py` script removes gaps and automatically renames the output file.

`degapper.py simion2017_taxa/Homo_sapiens.fasta`

Make the blast protein database with `makeblastdb` and then run `blastp` and report the results as tabular (`-outfmt 6`).

`blastp -query simion2017_taxa/Homo_sapiens.fasta.nogaps -db human_uniprot.fasta -outfmt 6 -max_target_seqs 1 > simion2017_taxa/hsapiens_vs_uniprot_blastp.tab`

The top hit for each protein should probably be the original reference protein. If this is not the case, then that protein probably should not be used in phylogeny in the first place. Using the blast results, new alignments can be generated for each protein from the supermatrix and the reference protein. Each file is in fasta format and contains three sequences, the protein used in the supermatrix (which has probably been trimmed), the reference protein, and a row of the extracted heteropecilly scores from the supermatrix. The folder containing all of the files is automatically named (based on current time).

`blast_to_align_pairs.py -b simion2017_taxa/hsapiens_vs_uniprot_blastp.tab -q simion2017_taxa/Homo_sapiens.fasta.nogaps -s human_uniprot.fasta -r simion2017_taxa/Homo_sapiens.fasta -p hp_by_site_w_const.tab`

This requires `mafft`, though could be modified to run another aligner.

### for RAxML site-wise likelihoods ###

Just as above, where an additional fasta line is generated for each alignment to show the heteropecilly scores, this can be used to show the difference in log-likelihood from `RAxML` using the `-l` option. This assumes that the relevant trees are the first two columns, and converts the difference into a hexadecimal value at each position, where no difference is 8, 9-f favor topology 1 and 0-7 favor topology 2.

`blast_to_align_pairs.py -b simion2017_taxa/hsapiens_vs_uniprot_blastp.tab -q simion2017_taxa/Homo_sapiens.fasta.nogaps -s human_uniprot.fasta -r simion2017_taxa/Homo_sapiens.fasta -l RAxML_perSiteLLs.simion2017_97sp_401632pos_1719genes.tab`

## alignment_pair_stats ##
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

Remarkably, even though the pipelines in each paper should be recovering single-copy genes (usually starting from genomes to find the orthologs, then adding transcriptomes to fill up the matrix), very few genes are actually shared between the three studies.

![align_pairs_3-way_venn.png](https://github.com/wrf/heteropecilly/blob/master/pair_stats/align_pairs_3-way_venn.png)

## amino_acid_distributions ##
How do the constant sites compare to the rest of the sequences? Are the sequences that were chosen representative of all proteins? Because they were almost all trimmed to some degree, are they even representative of themselves? The frequency of each amino acid and gaps is counted for all taxa, as well as constant sites. Two reference files of the 1499 human proteins used in the alignments as well as all human proteins are also used to compared the frequencies.

`amino_acid_distributions.py -a simion2017_97sp_401632pos_1719genes.fasta -p simion2017_human_uniprot.fasta human_uniprot.fasta > simion2017_aa_counts.tab`

All 20 amino acids are constant in at some point. For several amino acids, the frequency as a constant site is far above what is found on average (notably glycine, also proline, tryptophan and cysteine), consistent with important structural roles for these ones. For others, the occurrence is much lower (L, V, I, S, T), which account for non-polar and small-polar amino acids. Other biochemical groups like R/K and D/E have one that is less frequent (K and E) with one much more frequent (R and D), again possibly due to slightly different biochemical roles (such as R being always charged).

![simion2017_aa_counts.png](https://github.com/wrf/heteropecilly/blob/master/aa_counts/simion2017_aa_counts.png)

Given the normal approach of trimming, it appears that some amino acids are disproportionately removed, whereby the full-length human proteins have far more S, P, and to a lesser extent, Q and E than the trimmed versions in the alignment. As almost 60% of the original amino acids were removed for the Simion set (and 40% in most other datasets), the final alignments ultimately may be biased by the sequence composition of certain domains. [Tan et al 2015](https://academic.oup.com/sysbio/article/64/5/778/1685763) found that alignment trimming generally makes trees worse, as essentially all methods overtrim. Given this, the high average trimming in a number of studies suggests that large amounts of phylogenetic information are being lost by trimming.

## get_best_structures ##
Although most phylogenetics programs treat sites within an alignment as independent entities, biochemists have known for decades that sites within proteins interact. 

A large number of proteins already have crystal structures, or partial structures. This information is stored within the [text-format database for SwissProt](http://www.uniprot.org/downloads), and can be extracted into a summary table with the `parse_swissprot_data.py` script. Option `-s` sets the species filter to `9606`, which is the NCBI ID for *H. sapiens*.

`parse_swissprot_data.py -u uniprot_sprot.dat.gz -s 9606 > uniprot-9606_prots_w_pdb.tab`

This will generate a table for all human proteins that have any structure in PDB. This does not include proteins where a mouse homolog has a crystal structure that can be used as a template for a human protein model. Each row contains 6 columns:

* Uniprot ID, meaning the protein name, not the 6-digit accession number
* PDB ID
* length of reference protein
* string of protein chain and span
* protein coverage by the structure, in residues
* fraction protein coverage by the structure

`EGFR_HUMAN  4LI5  1210  A=696-1020  325  0.269`

Because the crystal set only contains structures for 618 proteins, less than half of the 1499 human proteins in the set have 3D information. However, the heteropecilly information can be mapped on to the existing structures using a series of scripts, just as can be done for any other information about each site.

Then, to sort out which proteins in the supermatrix have a structure, read in the alignment information, both to compare the span of the protein in the supermatrix and the structures, and also to generate a shell script to colorize them in [PyMol](https://pymol.org/2/) with `-c`.

`get_best_structures.py -U ~/db/uniprot-9606_prots_w_pdb.tab -a blast_alignments_20171009-153805/ -s 9606 -c pdb_hp_commands.sh > human_prots_w_pdb.tab`

The shell script will call `pdb_heteropecilly.py` for each PDB file in a directory, returning an error if that file is not there. This will recode the beta factors of each atom for each residue as heteropecilly scores. Heteropecilly color scheme can be visualized within Pymol, using the command in the console `run ~/git/pdbcolor/color_by_heteropecilly.py`

Below is an example from 2o8b.pdb, which is the structure of the mismatch repair protein MSH2/MSH6 heterodimer (see [Warren et al 2007 Structure of the human MutSalpha DNA lesion recognition complex.](https://www.ncbi.nlm.nih.gov/pubmed/?term=17531815)). Heteropecilly scores were only calculated for MSH2, so the other protein is colored in pale gray. Gaps or missing data are dark gray, constant sites are green, and the colors follow the deciles as in the charts above. In this example, the DNA helix is colored yellow, to distinguish it from the protein. Several features are evident. Large sections of the alignment had been removed by trimming, resulting in gaps when compared to the reference protein, and dark gray regions throughout the protein (long helix connecting the "clamp" domain to the ATPase domain). Constant sites form a distinct sector at the bottom of the image, primarily comprising the ATPase domain, probably also involved in the interface with MSH6. Many heteropecillious sites (red) appear to occur on the surface of the protein, perhaps directly interacting with the solvent or other proteins, though the extent of this was not precisely calculated. This may mean that, in general, heteropecillious sites and lineage-specific changes are a reflection of unique interactions *between* proteins.

![2o8b_w_hp.png](https://github.com/wrf/heteropecilly/blob/master/2o8b_w_hp.png)
