# npstat2

## Population genetics tests and estimators for pooled NGS data

#### Sebastian E. Ramos-Onsins, Sara Guirao-Rico, Ahmed Hafez and Luca Ferretti

This code implements some population genetics tests and estimators that can be applied to pooled sequences from Next Generation Sequencing experiments. The statistics are described in the paper "Population genomics from pool sequencing" by L. Ferretti, S.E. Ramos-Onsins and M. Perez-Enciso, Molecular Ecology (2013), DOI: 10.1111/mec.12522.

In this second version, we have implemented analysis multiscaffold, quasi-singletons, population differentation analysis (Fst, when an additional population is included) and we have increased the number of statistics shown, specially for functional regions. In addition, the assignation of functional positions will be updated soon to increase the precission of the analysis on these regions.

## Code under debugging! 

This code is still not completely validated. There are still known problems in the correct calculation of Fst values (a slight overestimation).

## How to compile

	gcc -o npstat2 NPStat-v2.c -lgsl -lgslcblas -lm -lhts
	
## Input format

The main input of the program is in pileup.gz format. It is also mandatory to include a file with the names of the scaffolds, in the same order than in the pipeline.

Other three types of files could be useful:

- Outgroup sequence in FASTA format. The sequence should be aligned with the reference used to align the .bam file. Multiscaffold data is allowed. Scaffolds must be called with the same name than included in the pileup input file. 

- File with a list of positions of filtered SNPs for the chromosome or scaﬀold analyzed. This could be useful to analyze only the SNPs called by some SNP calling software for pools.

- Annotation file in GFF3 format. This allows to perform the McDonald-Kreitman test. Scaffolds must be called with the same name than included in the pileup input file. 

## How to use

Command:

	npstat [flags] file.pileup.gz
	    
Mandatory flags:

    -n : haploid sample size
    -l : window length
    -scaffolds : name of file including scaffolds   
 
 Additional Optional flags:
   
    -outfile : name output file (default ends with extension '.stats.txt')
    -fstpop2 file2.pileup : computes Fst with a second population contained in file2.pileup
    -n2 : sample size of the second population
    -nolowfreq m : filter on minimum allele count mac>m
    -mincov minimum_coverage : filter on minimum coverage (default 4)
    -maxcov maximum_coverage : filter on maximum coverage (default 100)
    -minqual minimum_base_quality : filter on base quality (default 10)
    -outgroup file.fa : outgroup file in FASTA
    -annot file.gff3 : annotation file in GFF3
    -snpfile file.snp : consider SNPs only if present in file.snp

An important option is `-nolowfreq m`. This specifies how many alleles of
low frequency are discarded. The default option is m=2, which means that
alleles appearing in only 2 reads will be discarded. Data at low coverage 
would need lower values (i.e., m=1 for for read depth smaller than 10.
High error rate would need higher values, e.g. m=100 above read depth 100, etc. 
Use m=0 only if the SNPs have already been called by an external SNP caller 
and passed to the program through the option -snpfile.

## Output

We have included additional statistics in relation to the first version. One file is output (by defaut has the same name than the initial input pileup with extension '.stats.txt', unless the outpfile option is used). If a second population is included (option -n2 and -fstpop2), then two additional files are output: a file with the statististics for the second population (same format than first but with the additional extension '\_p2.txt'), and a third file with the differentiation statistics (same name with the extension '\_fst.txt').

The file with the population output contains the next statistics:

	1.scaffold: name of the scaffold.
	2.window: window number.
	3.start: initial position of the analyzed  window.
	4.end: final position of the analyzed  window.
	5.length: number of bases covered in the window.
	6.length_outgroup: number of bases covered and with known outgroup allele.
	7.read_depth: average read depth.
	8.S: number of segregating sites S.
	9.Theta_FL*: Variability estimate considering quasi-singletons from folded SFS.
	10.Watterson: Watterson estimator of theta.
	11.Pi: Tajima’s Pi estimator of heterozygosity.
	12.Tajima_D: Tajima’s D.
	13.unnormFL*test: unnormalized FuLi D* Test (using quasi-singletons). 
	14.var_S: variance of the number of segregating sites.
	15.var_Watterson: variance of the Watterson estimator.
	16.theta_FL: Variability estimate considering quasi-singletons (unfolded SFS). 
	17.thetaH: Fay and Wu's variability estimator.
	18:thetaZengE: ZengE variability estimator.
	19.unnormFLtest: unnormalized FuLi D Test (using quasi-singletons).
	20.unnorm_FayWu_H: unnormalized Fay and Wu’s H Test.
	21.unnormZengEtest: unnormalized Zeng's E Test.
	22.FayWu_H: normalized Fay and Wu’s H.
	23.div: divergence per base (from outgroup). 
	24.nonsyn_S: nonsynonimous polymorphisms.
	25.syn_S: synonimous polymorphisms
	26.nonsyn_div: nonsynonimous divergence
	27.syn_div: synonimous polymorphisms
	28.len_ns: number of nonsynonymous bases covered in the window.
	29.len_out_ns: number of nonsynonymous bases covered and with known outgroup allele.
	30.thetaFL*_ns: quasi-singletons variability estimate from folded SFS (nonsynonymous) 
	31.Watt_ns: Watterson estimator of theta of nonsynonymous positions.
	32.pi_ns: Tajima’s Pi estimator of heterozygosity of nonsynonymous positions.
	33.thetaFL_ns: quasi-singletons variability estimate from unfolded SFS (nonsynonymous)
	34.thetaH_ns: Fay and Wu's variability estimator of nonsynonymous positions.
	35.thetaZE_ns: ZengE variability estimator (nonsynonymous).
	36.div_ns: divergence per nonsynonymous base (from outgroup).
	37.len_syn: number of synonymous bases covered in the window.
	38.len_out_syn: number of synonymous bases covered and with known outgroup allele. 
	39.thetaFL*_syn: quasi-singletons variability estimate from folded SFS (synonymous) 
	40.Watt_syn: Watterson estimator of theta of synonymous positions.
	41.pi_syn: Tajima’s Pi estimator of heterozygosity of synonymous positions.
	42.thetaFL_syn: quasi-singletons variability estimate from unfolded SFS (synonymous)
	43.thetaH_syn: Fay and Wu's variability estimator of synonymous positions
	44.thetaZE_syn: ZengE variability estimator (synonymous).
	45.div_syn: divergence per synonymous base (from outgroup).
	46.alpha: fraction of substitutions fixed by positive selection.
	47.alpha_watt: fraction of substitutions fixed by positive selection considering Watterson variability.
	48.alpha_pi: fraction of substitutions fixed by positive selection considering Tajima's Pi variability.
	49.alpha_H: fraction of substitutions fixed by positive selection considering Fay and Wu's H variability.

All these statistics are computed after filtering for minimum read depth,
qualities and allele count. Note slightly different assumptions in relation to npstat1 version: (i) the consideration of a variant below the -nolowfreq is rejected and not account to consider three mutliple mutations with the outgroup, and (ii) the synonymous and nonsynonymous segregating sites are counted independently to have outgroup or not at that position. The reason is that in this version the alpha estimates are calculated using the variability and divergence levels and not the explicitely the number of mutations.

The HKA test can be obtained by composing data from S (columns 8), Var(S) (column 14) and divergence (column 23). The McDonald-Kreitman
test could be obtained by composing synonymous/nonsynonimous polymorphism/divergence data in a 2 × 2 contingency table. 

Synonymous and non-synonymous positions are calculated considering Nuclear Universal coding and using Nei-Gojobori (1986) method. Only codons with less than 3 mutations were considered.

The file containing the differentiation output contains the next statistics:

	1.scaffold: name of the scaffold.
	2.window: window number.	
	3.start: initial position of the analyzed  window.
	4.end: final position of the analyzed  window.
	5.length: number of bases covered in the window (considering both populations)
	6.nVariants: number of total variants in this window.
	7.pw_diff_12: pairwise differences between population 1 and 2.
	8.Pi_1: Tajima’s Pi estimator of heterozygosity in population 1 for common positions with pop 2.
	9.Pi_2: Tajima’s Pi estimator of heterozygosity in population 2 for common positions with pop 1.
	10.Pi_t: Tajima’s Pi estimator of heterozygosity considering together pop 1 and 2.
	11.Fst: Differentiation statistic (Fst=1-(mean(Pi_1+Pi_2))/Pi_t)
	
## Example

We provide a small example containing with some regions from different scaffolds of Drosophila melanogaster and D. yakuba as the outgroup. the two populations used are coming from the same population. Unzip before running.

	../npstat2 -n 16 -l 10000 -nolowfreq 2 -minqual 18 -outgroup Dyakuba-mel_final_cns.fa -annot dmel-all-r6.12_sorted.gtf -scaffolds scaffold_file.txt -fstpop2 Pool_seq2.mel_mpileup.txt.gz -n2 16 -outfile Pool_seq_npstat2_results Pool_seq1.mel_mpileup.txt.gz
	
A simple R script to plot the results of the output is included. Simply include the name of the npstat output file of interest and the name of the output pdf file as arguments:

	R --vanilla --args [output_npstat2.txt] [output_plots.pdf] < npstat_plot_windows2.R 

In this example you may do:

	R --vanilla --args Pool_seq_npstat2_results Pool_seq_npstat2_results.pdf < ./npstat2_plot_windows.R
	R --vanilla --args Pool_seq_npstat2_results_p2.txt Pool_seq_npstat2_results_p2.pdf < ./npstat2_plot_windows.R
	R --vanilla --args Pool_seq_npstat2_results_fst.txt Pool_seq_npstat2_results_fst.pdf < ./npstat2_plot_windows.R

	
	
