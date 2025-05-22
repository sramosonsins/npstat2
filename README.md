# npstat2

## Population genetics tests and estimators for pooled NGS data

#### Sebastian E. Ramos-Onsins, Sara Guirao-Rico, Ahmed Hafez and Luca Ferretti

This code implements some population genetics tests and estimators that can be applied to pooled sequences from Next Generation Sequencing experiments. The statistics are described in the paper "Population genomics from pool sequencing" by L. Ferretti, S.E. Ramos-Onsins and M. Perez-Enciso, Molecular Ecology (2013), DOI: 10.1111/mec.12522.

In this second version, we have implemented analysis multiscaffold, population differentation analysis (Fst, if an additional population is included, STILL TESTING) and have increased the number of statistics shown, specially for functional regions. In addition, the assignation of functional positions will be updated soon to increase the precission of the analysis on these regions.

## How to compile

	gcc -o npstat2 NPStat-v2.c -lgsl -lgslcblas -lm
	
## Input format

The main input of the program is in pileup format. 

Other three types of files could be useful:

- Outgroup sequence in FASTA format. The sequence should be aligned with the reference used to align the .bam file. Multiscaffold data is allowed. Scaffolds must be called with the same name than included in the pileup input file. 

- File with a list of positions of filtered SNPs for the chromosome or scaﬀold analyzed. This could be useful to analyze only the SNPs called by some SNP calling software for pools.

- Annotation file in GFF3 format. This allows to perform the McDonald-Kreitman test. Scaffolds must be called with the same name than included in the pileup input file. 

## How to use

Command:

	npstat [options] file.pileup
	
or to read from standard input:
    
    npstat [options] -
    
Options:

    -n samplesize : haploid sample size
    -l windowlength : window length
    -mincov minimum_coverage : filter on minimum coverage (default 4)
    -maxcov maximum_coverage : filter on maximum coverage (default 100)
    -minqual minimum_base_quality : filter on base quality (default 10)
    -nolowfreq m : filter on minimum allele count mac>m
    -outgroup file.fa : outgroup file in FASTA
    -annot file.gff3 : annotation file in GFF3
    -snpfile file.snp : consider SNPs only if present in file.snp
    -scaffolds : name file including scaffolds
    -outfile : name output file (default ends with extension '.stats.txt')
    -fstpop2 file2.pileup : computes Fst with a second population contained in file2.pileup
    -n2 : sample size of the second population

An important option is `-nolowfreq m`. This specifies how many alleles of
low frequency are discarded. The default option is m=1, which means that
alleles appearing in only 1 read will be discarded. Data at high coverage or
high error rate would need higher values, e.g. m=2 above read depth 50,
m=3 above read depth 100, etc. Use m=0 only if the SNPs have already
been called by an external SNP caller and passed to the program through
the option -snpfile.

## Output

We have included additional statistics in relation to the first version. One file is output (by defaut has the same name than the initial input pileup with extension '.stats.txt', unless the outpfile option is used). If a second population is included (option -n2), then two additional files are ouptut: a file with the statististics for the second population (same format than first but with the additional extension '\_p2.txt'), and a third file with the differentiation statistics (same name with the extension '\_fst.txt').

The file with the population output contains the next statistics:

	1.scaffold: name of the scaffold.
	2.window: window number.
	3.length: number of bases covered in the window.
	4.length_outgroup: number of bases covered and with known outgroup allele.
	5.read_depth: average read depth.
	6.S: number of segregating sites S.
	7.Watterson: Watterson estimator of theta.
	8.Pi: Tajima’s Pi estimator of heterozygosity,
	9.Tajima_D: Tajima’s D.
	10.var_S: variance of the number of segregating sites.
	11.var_Watterson: variance of the Watterson estimator.
	12.thetaH: Fay and Wu's variability estimator.
	13.unnorm_FayWu_H: unnormalized Fay and Wu’s H.
	14.FayWu_H: normalized Fay and Wu’s H.
	15.div: divergence per base (from outgroup). 
	16.nonsyn_S: nonsynonimous polymorphisms.
	17.syn_S: synonimous polymorphisms
	18.nonsyn_div: nonsynonimous divergence
	19.syn_div: synonimous polymorphisms
	20.len_ns: number of nonsynonymous bases covered in the window.
	21.len_out_ns: number of nonsynonymous bases covered and with known outgroup allele.
	22.Watt_ns: Watterson estimator of theta of nonsynonymous positions.
	23.pi_ns: Tajima’s Pi estimator of heterozygosity of nonsynonymous positions.
	24.thetaH_ns: Fay and Wu's variability estimator of nonsynonymous positions.
	25.div_ns: divergence per nonsynonymous base (from outgroup).
	26.len_syn: number of synonymous bases covered in the window.
	27.len_out_syn: number of synonymous bases covered and with known outgroup allele. 
	28.Watt_syn: Watterson estimator of theta of synonymous positions.
	29.pi_syn: Tajima’s Pi estimator of heterozygosity of synonymous positions.
	30.thetaH_syn: Fay and Wu's variability estimator of synonymous positions
	31.div_syn: divergence per synonymous base (from outgroup).
	32.alpha: fraction of substitutions fixed by positive selection.
	33.alpha_watt: fraction of substitutions fixed by positive selection considering Watterson variability.
	34.alpha_pi: fraction of substitutions fixed by positive selection considering Tajima's Pi variability.
	35.alpha_H: fraction of substitutions fixed by positive selection considering Fay and Wu's H variability.

All these statistics are computed after filtering for minimum read depth,
qualities and allele count. The HKA test can be obtained by composing data from S (columns 6), Var(S) (column 10) and divergence (column 15). The McDonald-Kreitman
test can be obtained by composing synonimous/nonsynonimous polymorphism/divergence data from columns 16-18 in a 2 × 2 contingency table.

Note that we approximate all aminoacids to be 4-fold degenerate (i.e.
nonsynonimous and synonimous sites actually correspond to the 1st/2nd
base and 3rd base in the codon, respectively). this will be updated to a Nuclear Universal coding code soon.

The file containing the differentiation output contains the next statistics:

	1.scaffold: name of the scaffold.
	2.window: window number.	
	3.length: number of bases covered in the window (considering both populations)
	4.pw_diff_12: pairwise differences between population 1 and 2.
	5.Pi_1: Tajima’s Pi estimator of heterozygosity in population 1.
	6.Pi_2: Tajima’s Pi estimator of heterozygosity in population 1.
	7.Pi_a: Tajima’s Pi estimator of heterozygosity between pop 1 and 2.
	8.Fst: Differentiation statistic
	
## Example

We provide a small example containing with some regions from different scaffolds of Drosophila melanogaster and D. yakuba as the outgroup. the two populations used are coming from the same population. Unzip before running.

	./npstat2 -n 16 -l 10000 -nolowfreq 2 -minqual 18 -outgroup Dyakuba-mel_final_cns.fa -annot dmel-all-r6.12_sorted.gtf -scaffolds scaffold_file.txt -fstpop2 Pool_seq2.mel_mpileup.txt -n2 16 -outfile Pool_seq_npstat2_results.txt Pool_seq1.mel_mpileup.txt
	
A R script plot the results of the output files. Simply include the name of the npstat output file of interest and the name of the output pdf file as arguments:

	R --vanilla --args [output_npstat2.txt] [output_plots.pdf] < npstat_plot_windows2.R 
	
