# npstat2

## Population genetics tests and estimators for pooled NGS data

#### Sebastian E. Ramos-Onsins, Sara Guirao-Rico, Ahmed Hafez and Luca Ferretti

This code implements some population genetics tests and estimators that can be applied to pooled sequences from Next Generation Sequencing experiments. The statistics are described in the paper "Population genomics from pool sequencing" by L. Ferretti, S.E. Ramos-Onsins and M. Perez-Enciso, Molecular Ecology (2013), DOI: 10.1111/mec.12522. npstat2 is designed to work with pools containing from few to few hundred individuals and with no more than few hundred reads of mean read depth.

In this second version, we have implemented analysis multiscaffold, quasi-singletons, population differentation analysis (Fst, when an additional population is included) and we have increased the number of statistics shown, specially for functional regions. In addition, functional positions are now estimated using Nei & Gojobori (1986) method. This version allow to read a BED file for analyzing specific windows (*e.g.*, windows collecting a multigene family gene regions or whatever you are interested). We have modified the default value of *nolowfreq* parameter to 2 (eliminate variants with 2 or less reads) to increase accuracy. Note that the parameters *nolowfreq*, *mincov*, *maxcov* and *minqual* are common for both pools.

Remember to filter your mpileup(s) to exclude sites with possible mapping errors or copy-number variants. Mapping quality filtering is not included in the code.

## Comparison with other softwares 

We are validating the software using simulations under different conditions. Fst works acceptably in comparison with other tools. See simulation results included in the repository.

## How to compile

	gcc -o npstat2 NPStat-v2.c -lgsl -lgslcblas -lm -lhts
	
for MacOS, you may need to include the paths to libraries (previously installed using *brew install*). 

Example:
	
	gcc -o npstat NPStat-v2.c -lgsl -lgslcblas -lm -lhts -I/opt/homebrew/Cellar/gsl/2.8/include -I/opt/homebrew/Cellar/htslib/1.22.1/include -L/opt/homebrew/Cellar/gsl/2.8/lib -L/opt/homebrew/Cellar/htslib/1.22.1/lib

	
## Input format

The main input of the program is in pileup.gz format. It is also mandatory to include a file with the names of the scaffolds, in the same order than in the pipeline.

Other three types of files could be useful:

- Outgroup sequence in FASTA format. The sequence should be aligned with the reference used to align the .bam file. Multiscaffold data is allowed. Scaffolds must be called with the same name than included in the pileup input file. 

- File with a list of positions of filtered SNPs for the chromosome or scaﬀold analyzed. This could be useful to analyze only the SNPs called by some SNP calling software for pools.

- Annotation file in GFF3 format. This allows to perform the McDonald-Kreitman test. Scaffolds must be called with the same name than included in the pileup input file. WARNING: The GFF3 file can not work with different alternative splicing forms per gene (that is, include just one of the possible alternative splicing cases) and overlapping gene annotations must be erased from the file.

## How to use

Command:

	npstat2 [flags] file.pileup.gz
		    
Mandatory flags:

    -n : haploid sample size
    -scaffolds : name of file including scaffolds   
    
    And one of these two options:
    -l : window length
    -bedfile: filename of the bedfile (rows should contain scaffold start end)
 
 Additional Optional flags:
   
    -outfile : name output file (default ends with extension '.stats.txt')
    -fstpop2 file2.pileup : computes Fst with a second population contained in file2.pileup
    -n2 : sample size of the second population
    -nolowfreq m : filter on minimum allele count mac>m (default 1) 
    -mincov minimum_coverage : filter on minimum coverage (default 4)
    -maxcov maximum_coverage : filter on maximum coverage (default 100)
    -minqual minimum_base_quality : filter on base quality (default 10)
    -outgroup file.fa : outgroup file in FASTA
    -annot file.gff3 : annotation file in GFF3
    -snpfile file.snp : consider SNPs only if present in file.snp

An important option is `-nolowfreq m`. This specifies how many alleles of
low frequency are discarded. The default option is m=2, which means that
alleles appearing in only 2 reads will be discarded. Data at low coverage 
would need lower values (i.e., m=1 for read depth smaller than 100.
High error rate would need higher values, e.g. m=3 above read depth 1000, etc. 
Use m=0 only if the SNPs have already been called by an external SNP caller 
and passed to the program through the option -snpfile.

## Output

We have included additional statistics in relation to the first version. One file is output (by defaut has the same name than the initial input pileup with extension '.stats.txt', unless the outpfile option is used). If a second population is included (option -n2 and -fstpop2), then two additional files are output: a file with the statististics for the second population (same format than first but with the additional extension '\_p2.txt'), and a third file with the differentiation statistics (same name with the extension '\_fst.txt').

The file containing the population output includes the next statistics:

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
	24.nonsyn_S: nonsynonym polymorphisms.
	25.syn_S: synonym polymorphisms
	26.nonsyn_div: nonsynonym divergence
	27.syn_div: synonym polymorphisms
	28.len_ns: number of nonsynonym bases covered in the window.
	29.len_out_ns: number of nonsynonym bases covered and with known outgroup allele.
	30.thetaFL*_ns: quasi-singletons variability estimate from folded SFS (nonsynonym) 
	31.Watt_ns: Watterson estimator of theta of nonsynonym positions.
	32.pi_ns: Tajima’s Pi estimator of heterozygosity of nonsynonym positions.
	33.thetaFL_ns: quasi-singletons variability estimate from unfolded SFS (nonsynonym)
	34.thetaH_ns: Fay and Wu's variability estimator of nonsynonym positions.
	35.thetaZE_ns: ZengE variability estimator (nonsynonym).
	36.div_ns: divergence per nonsynonym base (from outgroup).
	37.len_syn: number of synonym bases covered in the window.
	38.len_out_syn: number of synonym bases covered and with known outgroup allele. 
	39.thetaFL*_syn: quasi-singletons variability estimate from folded SFS (synonym) 
	40.Watt_syn: Watterson estimator of theta of synonym positions.
	41.pi_syn: Tajima’s Pi estimator of heterozygosity of synonym positions.
	42.thetaFL_syn: quasi-singletons variability estimate from unfolded SFS (synonym)
	43.thetaH_syn: Fay and Wu's variability estimator of synonym positions
	44.thetaZE_syn: ZengE variability estimator (synonym).
	45.div_syn: divergence per synonym base (from outgroup).
	46.alpha: fraction of substitutions fixed by positive selection.
	47.alpha_watt: fraction of substitutions fixed by positive selection considering Watterson variability.
	48.alpha_pi: fraction of substitutions fixed by positive selection considering Tajima's Pi variability.
	49.alpha_H: fraction of substitutions fixed by positive selection considering Fay and Wu's H variability.

All these statistics are computed after filtering for minimum read depth,
qualities and allele count. Note slightly different assumptions in relation to npstat1 version: (i) the consideration of a variant below the -nolowfreq is rejected and not account to consider three mutliple mutations with the outgroup, and (ii) the synonym and nonsynonym segregating sites are counted independently to have outgroup or not at that position. The reason is that in this version the alpha estimates are calculated using the variability and divergence levels and not the explicitely the number of mutations.

The HKA test can be obtained by composing data from S (columns 8), Var(S) (column 14) and divergence (column 23). The McDonald-Kreitman
test could be obtained by composing synonym/nonsynonym polymorphism/divergence data in a 2 × 2 contingency table. 

Synonymous and non-synonymous positions are calculated considering Nuclear Universal coding and using Nei-Gojobori (1986) method. Only codons with less than 3 mutations were considered.

The file containing the differentiation output includes the next statistics:

	1.scaffold: name of the scaffold.
	2.window: window number.	
	3.start: initial position of the analyzed  window.
	4.end: final position of the analyzed  window.
	5.length: number of bases covered in the window (considering both populations)
	6.nVariants: number of total variants in this window.
	7.pw_diff_12: pairwise differences between population 1 and 2.
	8.Pi_1: Tajima’s Pi estimator of heterozygosity/nt in population 1 for common positions with pop 2.
	9.Pi_2: Tajima’s Pi estimator of heterozygosity/nt in population 2 for common positions with pop 1.
	10.Pi_a: pairwise differences between population 1 and 2 divided by the total effective positions.
	11.Pi_t: Tajima’s Pi estimator of heterozygosity/nt considering together pop 1 and 2.
	12.FstT: Differentiation statistic (Fst=1-(mean(Pi_1+Pi_2))/Pi_t)
	13.FstA: Differentiation statistic (Fst=1-(mean(Pi_1+Pi_2))/Pi_a)

By default, Fst is calculated from the estimation of Pi_a to later estimate Pi_total, using the equations developed in Ferretti et al (2013).
	
## Example

We provide a small example containing some regions from different scaffolds of Drosophila melanogaster and D. yakuba as the outgroup. the two populations used are coming from the same population. Unzip before running.

	../npstat2 -n 16 -l 10000 -nolowfreq 1 -minqual 18 -outgroup Dyakuba-mel_final_cns.fa -annot dmel-all-r6.12_sorted.gtf -scaffolds scaffold_file.txt -fstpop2 Pool_seq2.mel_mpileup.txt.gz -n2 16 -outfile Pool_seq_npstat2_results Pool_seq1.mel_mpileup.txt.gz
	
A simple R script to plot the results of the output is included. Simply include the name of the npstat output file of interest and the name of the output pdf file as arguments:

	R --vanilla --args [output_npstat2.txt] [output_plots.pdf] < npstat_plot_windows2.R 

In this example you may do:

	R --vanilla --args Pool_seq_npstat2_results Pool_seq_npstat2_results.pdf < ./npstat2_plot_windows.R
	R --vanilla --args Pool_seq_npstat2_results_p2.txt Pool_seq_npstat2_results_p2.pdf < ./npstat2_plot_windows.R
	R --vanilla --args Pool_seq_npstat2_results_fst.txt Pool_seq_npstat2_results_fst.pdf < ./npstat2_plot_windows.R

	
	
