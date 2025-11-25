/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/
/*--------------------------- NPStat ---------------------------*/
/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/

/* Code to extract pool statistics (theta, neutrality tests)
 from pileup data and a fasta file representing the outgroup sequence */

/* Compile with
 gcc -O3 -o npstat NPStat-vXX.c -lgsl -lgslcblas -lm
 substituting XX with version number
 */

/* Arguments:
 run "npstat" to see the arguments
 */

/* Output stats:
 window number, bases (after filtering), bases with known outgroup allele (after filtering),
 average read depth, number of segregating sites, Watterson theta, Tajima's Pi, Tajima's D, Fay and Wu's H, divergence per base (from outgroup), HKA etc...
 */

/* Include libraries */
#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <assert.h>
#include <time.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_math.h>
#include <htslib/bgzf.h>

#define VERSION_NPSTAT "npstat2 v.20251124\n"

/* Define substitutions */

#define PRINTSNPS 0
#define DEB(x) //x
#define SEB(x) x

#define print_comb 0 /* put 1 if you want to print combinatorial data, 0 otherwise. NOT IMPLEMENTED */
#define maxcheck 10000000 /* NOT IMPLEMENTED */
#define minimum_coverage 2 /* this fixes the minimum acceptable coverage in terms of reads aligned */
#define maximum_coverage 1000 /* this fixes the maximum acceptable coverage in terms of reads aligned */
#define min(a, b) (((a) < (b)) ? (a) : (b))
#define max(a, b) (((a) > (b)) ? (a) : (b))

#define NUM_ASCII 48

/* Declare structures */

struct tests
{
    unsigned long cov;
    double l;
    double l_out;
    double s;
    double num_t;
    double num_p;
    double num_hl;
    double num_hq;
    double num_nu;
    double num_xi;
    double den_t;
    double den_p;
    double den_hl;
    double den_hq;
    double den_nu;
    double den_xi;
};

struct fst_calc
{
    unsigned long la;
    //unsigned long lt;
    double gen_diff;
    double var_diffL;
    double c_s;
    double Lbias;
    double num_p1;
    double num_p1g;
    double num_p2;
    double num_p2g;
    //double num_pt;
    double den_p1;
    double den_p2;
    //double den_pt;
};

struct combinatorial
{
    double *d_t;
    double *d_p;
    double *d_hl;
    double *d_hq;
    double *d_nu;
    double *d_xi;
};

struct combinatorial_fst
{
    double **c_s;
    double *d_pt;
};

struct keep_codon_data {
    int n;
    unsigned long n_ref;
    unsigned long n_alt_allele;
    unsigned long rd;
    unsigned long n_alt[5];
    int ref_base;
    int alt_base;
    int out_base;
};

typedef struct {
    BGZF *fp;
    int pushed_back;
    int push_char;
} BGZFReader;

/* Declare Constants*/

static char tripletsN[64][3] =
{
    {"111"},    {"112"},    {"114"},    {"113"},    {"121"},    {"122"},    {"124"},    {"123"},
    {"141"},    {"142"},    {"144"},    {"143"},    {"131"},    {"132"},    {"134"},    {"133"},
    {"211"},    {"212"},    {"214"},    {"213"},    {"221"},    {"222"},    {"224"},    {"223"},
    {"241"},    {"242"},    {"244"},    {"243"},    {"231"},    {"232"},    {"234"},    {"233"},
    {"411"},    {"412"},    {"414"},    {"413"},    {"421"},    {"422"},    {"424"},    {"423"},
    {"441"},    {"442"},    {"444"},    {"443"},    {"431"},    {"432"},    {"434"},    {"433"},
    {"311"},    {"312"},    {"314"},    {"313"},    {"321"},    {"322"},    {"324"},    {"323"},
    {"341"},    {"342"},    {"344"},    {"343"},    {"331"},    {"332"},    {"334"},    {"333"},
};
/* order: TCGA=1234 */
static char NuclearUniversalCode[64] =
{
    'F', 'F', 'L', 'L',
    'S', 'S', 'S', 'S',
    'Y', 'Y', '*', '*',
    'C', 'C', '*', 'W',
    'L', 'L', 'L', 'L',
    'P', 'P', 'P', 'P',
    'H', 'H', 'Q', 'Q',
    'R', 'R', 'R', 'R',
    'I', 'I', 'I', 'M',
    'T', 'T', 'T', 'T',
    'N', 'N', 'K', 'K',
    'S', 'S', 'R', 'R',
    'V', 'V', 'V', 'V',
    'A', 'A', 'A', 'A',
    'D', 'D', 'E', 'E',
    'G', 'G', 'G', 'G'
};

/* Declare functions */

int generate_pool_covariance_matrix(double *covmat, double *covmatpool, int na, int nb, int n);
int generate_covariance_matrix(double *covmat, int n);

int base_to_num(char base);

int read_line_pileup(BGZFReader *reader /*FILE * bam_file*/, unsigned long min_qual, unsigned long min_mqual, unsigned long * pos_base, unsigned long * n_ref, unsigned long * n_alt_allele, /*unsigned long*/int * n_tot_allele, unsigned long * n_alt, int * ref_base, int * alt_base, char *cchrom, int m_bar);

int extract_outgroup_base(FILE * fasta_out, unsigned long pos, unsigned long oldpos, int fasta_length);

void init_codons(struct keep_codon_data *codon_pop);
int which_Aa(char *codon);
void calculate_NSynSynPos(char *codon, int *nsyn);
void invert_codon(char *codon, char *inv_codon);
void codon_weights_muts(double *cweights, double *dweights, double *cmutnsyn, double *kmutnsyn, double *cmutsyn, double *kmutsyn, struct keep_codon_data *codon_pop, char strand);

void extract_stats(struct tests * test, struct combinatorial * comb, int n0, unsigned long n_ref, unsigned long n_alt_allele, unsigned long rd, unsigned long * n_alt, int ref_base, int alt_base, int out_base, int mb, float weightp, float weightd, double cmut);

void extract_fst(struct fst_calc * fst, struct combinatorial_fst * combfst, int n01, unsigned long n_ref1, unsigned long n_alt_allele1, unsigned long rd1, unsigned long * n_alt_1, int ref_base1, int alt_base1, int n02, unsigned long n_ref2, unsigned long n_alt_allele2, unsigned long rd2, unsigned long * n_alt_2, int ref_base2, int alt_base2, int out_base, int mb, int *snp);

void extract_fstT(struct fst_calc * fst, struct combinatorial_fst * combfst, struct combinatorial *comb, struct combinatorial *comb2, int n01, unsigned long n_ref1, unsigned long n_alt_allele1, unsigned long rd1, unsigned long * n_alt_1, int ref_base1, int alt_base1, int n02, unsigned long n_ref2, unsigned long n_alt_allele2, unsigned long rd2, unsigned long * n_alt_2, int ref_base2, int alt_base2, int out_base, int mb, int mbf, int *snp, int fst_path);

//SNPS
void extract_snps(unsigned long pos, FILE * output_snps, unsigned long * n_alt_1, unsigned long * n_alt_2, int mb);
//GTF
unsigned long extract_gff(FILE *gffinput, char *gchrom, unsigned long *cds_start, unsigned int *phase_cds, char *strand);

//init tests
void init_tests(struct tests *test1, struct tests *test2, struct tests *tests1, struct tests *testn1, struct tests *tests2, struct tests *testn2, struct fst_calc *fst, unsigned long *vec_rd, unsigned long *vec2_rd, double *psyn, double *pnon, double *dsyn, double *dnon, double *psyn2, double *pnon2, double *dsyn2, double *dnon2, unsigned long *div, unsigned long *div2, unsigned long max_cov, int compute_fst);

//calculate vartests
void calculate_vartests(double *var_h, double *var_d, double *var_s, double *var0_s, double *var0_d, double *var0_h, struct tests *test1, unsigned long *vec_rd, double *vec_s, double *vec_p, double *vec_h, double *vec0_s, double *vec0_d, double *vec0_h, double *covmat, int n01, unsigned long min_cov, unsigned long  max_cov, int m_bar);

//calculate stats and tests
void calculate_window_stats(struct tests *test1, double *cov1_val, double *theta1_val, double *pi1_val, double *d1_val, double *thetaH_val, double *h1_val, double *div_val, unsigned long div/**/, double *thetaE_val/**/, double *thetanu_val, double *thetaxi_val, double *e1_val, double *f1_val, double *fo1_val/**/ );

double calculate_combfst_cs(unsigned int rd1, unsigned int rd2, int m_bar, int n01, int n02);

//calculate max
int imax(int i, int j);
float flmax(float i, float j);

int bgzf_reader_getc(BGZFReader *reader);
void bgzf_reader_ungetc(BGZFReader *reader, int c);
int bgzf_getdeline(BGZFReader *reader, char **line, size_t *len, char delim);

void usage(void);

/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/
/*-------------------------- Functions -------------------------*/
/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/

int bgzf_reader_getc(BGZFReader *reader) {
    if (reader->pushed_back) {
        reader->pushed_back = 0;
        return reader->push_char;
    }
    return bgzf_getc(reader->fp);
}

void bgzf_reader_ungetc(BGZFReader *reader, int c) {
    reader->pushed_back = 1;
    reader->push_char = c;
}

int bgzf_getdeline(BGZFReader *reader, char **line, size_t *len, char delim) {
    size_t cap = 256;
    int i = 0;
    int c;
    
    if (*line == NULL) {
        *line = malloc(cap);
        if (*line == NULL) return -1;
    }
    
    while ((c = bgzf_reader_getc(reader)) != -1 && c != delim) {
        if (i + 1 >= cap) {
            cap *= 2;
            *line = realloc(*line, cap);
            if (*line == NULL) return -1;
        }
        (*line)[i++] = c;
    }
    
    if (c == -1 && i == 0) return -1;
    
    (*line)[i] = '\0';
    *len = i;
    return i;
}

//calculate max
int imax(int i, int j) {
    if(i>j) return i;
    else return j;
}
float flmax(float i, float j) {
    if(i>j) return i;
    else return j;
}


int generate_covariance_matrix(double *covmat, int n)
{
    double *a,*b,a2;
    int i,j,n2,ca,cb;
    gsl_matrix *sigma;
    //double s,*a,*b,bb,bz,zz,a2,*xi_null2,*xi_b2;/*a[n+1],b[n]*/
    //int i,j,n2,ca,cb;
    //double theta; /*temporary*/
    //gsl_matrix *c,*invc,*sigma,*ct,*invct;
    /*struct test_r0 *opt;*/
    
    if(n<2) return -10000;
    
    n2=ceil((double)n/2);
    //debug
    //assert(n2*2==n);
    n=(int)n2*2;
    
    a = (double *)calloc((unsigned long int)n+1,sizeof(double));
    b = (double *)calloc((unsigned long int)n,sizeof(double));
    
    a[0]=0;
    a2=0;
    for(i=2;i<=n+1;i++){
        a[i-1]=a[i-2]+1/(double)(i-1);
    }
    for(i=1;i<n;i++){
        a2+=1/(double)(i*i);
        b[i-1]=2*(double)n*(a[n]-a[i-1])/(double)((n-i+1)*(n-i))-2/(double)(n-i);
    }
    b[n-1]=0;
    sigma=gsl_matrix_alloc(n-1,n-1);
    gsl_matrix_set_zero(sigma);
    for (i=1;i<n;i++){
        if (2*i<n) {
            gsl_matrix_set(sigma,i-1,i-1,b[i]);
        } else {
            if (2*i>n) {
                gsl_matrix_set(sigma,i-1,i-1,b[i-1]-1/gsl_pow_2((double)i));
            } else {
                gsl_matrix_set(sigma,i-1,i-1,(a[n-1]-a[i-1])*2/(double)(n-i)-1/gsl_pow_2((double)i));
            }
        }
    }
    for (i=1;i<n;i++){
        for (j=1;j<i;j++){
            if (i+j<n) {
                gsl_matrix_set(sigma,i-1,j-1,(b[i]-b[i-1])/2); gsl_matrix_set(sigma,j-1,i-1,(b[i]-b[i-1])/2);
            } else {
                if (i+j>n) {
                    gsl_matrix_set(sigma,i-1,j-1,(b[j-1]-b[j])/2-1/(double)(i*j));
                    gsl_matrix_set(sigma,j-1,i-1,(b[j-1]-b[j])/2-1/(double)(i*j));
                } else {
                    gsl_matrix_set(sigma,i-1,j-1,(a[n-1]-a[i-1])/(double)(n-i)+(a[n-1]-a[j-1])/(double)(n-j)-(b[i-1]+b[j])/2-1/(double)(i*j));
                    gsl_matrix_set(sigma,j-1,i-1,(a[n-1]-a[i-1])/(double)(n-i)+(a[n-1]-a[j-1])/(double)(n-j)-(b[i-1]+b[j])/2-1/(double)(i*j));
                }
            }
        }
    }
    
    for(ca=1;ca<n;ca++){
        for(cb=1;cb<n;cb++){
            covmat[(cb-1)*(n-1)+ca-1]=gsl_matrix_get(sigma,ca-1,cb-1);
        }
    }
    
    //gsl_ran_poisson_pdf(k,mu)
    
    free(a);
    free(b);
    
    gsl_matrix_free(sigma);
    return 1;
}

int base_to_num(char base)
{
    switch(base)
    {
        case 'A': return 4; break;
        case 'a': return 4; break;
        case 'C': return 2; break;
        case 'c': return 2; break;
        case 'G': return 3; break;
        case 'g': return 3; break;
        case 'T': return 1; break;
        case 't': return 1; break;
        default: return 0;
    };
};

/*--------------------------------------------------------------*/
//int read_line_pileup(FILE * bam_file, unsigned long min_qual, unsigned long min_mqual, unsigned long * pos_base, unsigned long * n_ref, unsigned long * n_alt_allele, unsigned long * n_tot_allele, unsigned long * n_alt, int * ref_base, int * alt_base)

//read_line_pileup(bam_file1, min_qual, min_mqual, &pos_base1, &n_ref1, &n_alt_allele1, &rd1, n_alt_1, &ref_base1, &alt_base1);
//int read_line_pileup(FILE * bam_file, unsigned long min_qual, unsigned long min_mqual, unsigned long * pos_base, unsigned long * n_ref, unsigned long * n_alt_allele, unsigned long * n_tot_allele, unsigned long * n_alt, int * ref_base, int * alt_base) ;

int read_line_pileup(BGZFReader *reader/*FILE * bam_file*/, unsigned long min_qual, unsigned long min_mqual, unsigned long * pos_base, unsigned long * n_ref, unsigned long * n_alt_allele, /*unsigned long*/int * n_tot_allele, unsigned long * n_alt, int * ref_base, int * alt_base, char *cchrom, int m_bar) {
    
    int count_i;
    char *cline, /* *cchrom,*/ *cpileup, *cqual, *cmqual, crefbase;
    unsigned long count_j, n_ins;
    unsigned long nseq;
    char ct1;
    size_t nline;
    
    DEB(printf("entering read routine\n")); //debug
    
    //cchrom=(char *) malloc(100);
    cpileup=(char *) malloc(40);
    cqual=(char *) malloc(20);
    cmqual=(char *) malloc(20);
    cline=(char *) malloc(1);
    
    nline=1;
    cline=(char *)0;
    
    //getdelim(&cline,&nline,9,bam_file);
    //sscanf(cline,"%s\t",cchrom);
    bgzf_getdeline(reader,&cline,&nline,9);
    sscanf(cline,"%lu\t",pos_base);
    if(*pos_base<1) return(0);
    bgzf_getdeline(reader,&cline,&nline,9);
    sscanf(cline,"%c\t",&crefbase);
    bgzf_getdeline(reader,&cline,&nline,9);
    sscanf(cline,"%lu\t",&nseq);
    bgzf_getdeline(reader,&cline,&nline,9);
    if (strlen(cline)>=40) cpileup=(char *)realloc(cpileup,strlen(cline)+1);
    sscanf(cline,"%s\t",cpileup);
    //getdelim(&cline,&nline,9,bam_file);
    bgzf_getdeline(reader,&cline,&nline,10);
    if (strlen(cline)>=20) cqual=(char *)realloc(cqual,strlen(cline)+1);
    //if (strlen(cline)>=20) cmqual=(char *)realloc(cmqual,strlen(cline)+1);
    sscanf(cline,"%s\t",cqual);
    
    ct1=bgzf_reader_getc(reader);
    bgzf_reader_ungetc(reader,ct1);
    if(ct1!=EOF) {
        bgzf_getdeline(reader,&cline,&nline,9);
        sscanf(cline,"%s\t",cchrom);
    }
    else cchrom[0]=EOF;
    
    //sscanf(cline,"%s\t",cmqual);
    //getdelim(&cline,&nline,10,bam_file);
    //if (strlen(cline)>=20) cmqual=(char *)realloc(cmqual,strlen(cline)+1);
    //sscanf(cline,"%s\n",cmqual);
    
    DEB(printf("read data %s\t%lu\t%c\t%lu\t%s\t%s\t%s\n",cchrom,*pos_base,crefbase,nseq,cpileup,cqual,cmqual)); //debug
    
    //for (;(ct!="\n")&&(ct!=EOF);ct=fgetc(bam_file)) {};
    //ct=fgetc(bam_file); printf("%c",ct);
    //ct=fgetc(bam_file); printf("%c",ct);
    //ct=fgetc(bam_file); printf("%c\n",ct);
    
    for(count_i=0;count_i<5;count_i++){
        n_alt[count_i]=0;
    };
    count_j=0;
    *n_ref=0;
    *n_alt_allele=0;
    *n_tot_allele=0;
    
    if((nseq>=minimum_coverage)&&(nseq<=maximum_coverage)) {
        
        for(count_i=0;count_i<strlen(cpileup);count_i++){
            
            switch (cpileup[count_i]) {
                case '^': count_i++; break;
                case '$': break;
                case '*': count_j++; break;
                case '+': {
                    count_i++;
                    for(n_ins=0;isdigit(cpileup[count_i])!=0;count_i++){
                        n_ins=n_ins*10+(cpileup[count_i]-NUM_ASCII); //printf(">"); //debug
                    }; //printf("%lu ",n_ins); //debug
                    for(n_ins=n_ins-1;n_ins>0;n_ins--) {
                        count_i++; //printf("<"); printf("%lu ",n_ins); if(count_i>20) break; //debug
                    }; //printf("%c ",cpileup[count_i]); //debug
                }; break;
                case '-': {
                    count_i++;
                    for(n_ins=0;isdigit(cpileup[count_i])!=0;count_i++){
                        n_ins=n_ins*10+(cpileup[count_i]-NUM_ASCII);
                    };
                    for(n_ins=n_ins-1;n_ins>0;n_ins=n_ins-1) {
                        count_i++;
                    };
                }; break;
                case 'N': count_j++; break;
                case 'n': count_j++; break;
                case '.': {if((cqual[count_j]>=min_qual+33)) (n_alt[0])++; count_j++;}; break;
                case ',': {if((cqual[count_j]>=min_qual+33)) (n_alt[0])++; count_j++;}; break;
                case 'A': {if((cqual[count_j]>=min_qual+33)) (n_alt[4])++; count_j++;}; break;
                case 'C': {if((cqual[count_j]>=min_qual+33)) (n_alt[2])++; count_j++;}; break;
                case 'G': {if((cqual[count_j]>=min_qual+33)) (n_alt[3])++; count_j++;}; break;
                case 'T': {if((cqual[count_j]>=min_qual+33)) (n_alt[1])++; count_j++;}; break;
                case 'a': {if((cqual[count_j]>=min_qual+33)) (n_alt[4])++; count_j++;}; break;
                case 'c': {if((cqual[count_j]>=min_qual+33)) (n_alt[2])++; count_j++;}; break;
                case 'g': {if((cqual[count_j]>=min_qual+33)) (n_alt[3])++; count_j++;}; break;
                case 't': {if((cqual[count_j]>=min_qual+33)) (n_alt[1])++; count_j++;}; break;
            };
            
            
        };
        
        *ref_base=base_to_num(crefbase);
        n_alt[*ref_base]=n_alt[0];
        n_alt[0]=0;

        *n_alt_allele = 0;
        *n_ref=*alt_base=0; //n_alt[0];
        for(count_i=1;count_i<5;count_i++){
            if (*n_ref < n_alt[count_i]) {
                *n_ref = n_alt[count_i];
                *ref_base = count_i;
            };
        };
        for(count_i=1;count_i<5;count_i++){
            if ((*n_alt_allele < n_alt[count_i]) && (count_i!=*ref_base)) {
                *n_alt_allele = n_alt[count_i];
                *alt_base = count_i;
            };
        };
        /*
        if(*n_ref <= m_bar) {
            n_alt[*ref_base] = 0;
            *n_ref = 0;
            *ref_base = 0;
        }
        if(*n_alt_allele <= m_bar) {
            n_alt[*alt_base] = 0;
            *n_alt_allele = 0;
            *alt_base = 0;
        }
        */
        *n_tot_allele = (int)(*n_ref + *n_alt_allele);
    };
    
    DEB(printf("using data %lu\t%lu\t%lu\t%lu\t%lu\t%lu\t%lu\t%lu\t%u\t%u\n", *n_ref, *n_alt_allele, *n_tot_allele, n_alt[0], n_alt[1],n_alt[2],n_alt[3],n_alt[4], *ref_base, *alt_base)); //debug
    
    DEB(printf("exit read routine\n")); //debug
    free(cpileup);
    free(cqual);
    free(cmqual);
    free(cline);
    //free(cchrom);
    
    return(1);
};

int extract_outgroup_base(FILE * fasta_out, unsigned long pos, unsigned long oldpos, int fasta_length)
{
    unsigned long diff;
    ldiv_t pos_newlines, oldpos_newlines;
    pos_newlines=ldiv(pos-1,(unsigned long)fasta_length);
    oldpos_newlines=ldiv(oldpos-1,(unsigned long)fasta_length);
    diff=pos-oldpos+pos_newlines.quot-oldpos_newlines.quot;
    fseek(fasta_out,diff-1,SEEK_CUR);
    return base_to_num(fgetc(fasta_out));
};

/* Gipo's input and ms2pileup functions
 void read_line_ms(FILE *bam_file, unsigned long w_s, int n03, long int * n_pos, unsigned long int int_positions[], int array_counts[])
 {
 int i = 0;
 char *line, * word, * pEnd;;
 line=(char *)malloc(10000);
 size_t nline;
 nline = 1;
 // Four lines that won't be used
 for (i=1; i<=4; i++)
 {
 getline(&line, &nline, bam_file);
 }
 // Line n° 5, number of positions
 getline(&line, &nline, bam_file);
 word = strtok(line, " ");
 word = strtok(NULL, " ");
 *n_pos = strtol(word, &pEnd, 10);
 // The sixth line, that with positions
 getline(&line, &nline, bam_file);
 word = strtok(line, " ");
 for (i=0; i<*n_pos; i++)
 {
 word = strtok(NULL, " ");
 int_positions[i] = ceil(w_s * strtof(word, &pEnd));
 if (i>0)
 {
 while (int_positions[i] <= int_positions[i - 1])
 {
 (int_positions[i])++;
 }
 }
 }
 // Variation matrix
 int var_line, column;
 for (var_line=1; var_line<=n03; var_line++)
 {
 getline(&line, &nline, bam_file);
 char character;
 for (column=1; column<=*n_pos; column++)
 {
 character = line[column - 1];
 if (character == '1')
 {
 (array_counts[column - 1])++;
 }
 }
 }
 }
 
 void base_repetition(unsigned long w_s, int n_ind, unsigned long * pos_base, long int * n_pos, unsigned long int int_positions[], int *n_variation, int array_counts[], double input_coverage, double error_rate, const gsl_rng * s, unsigned long * n_ref, unsigned long * n_alt_allele, unsigned long * n_tot_allele, unsigned long * n_alt, int * ref_base, int * alt_base)
 {
 unsigned long int cur_phys_pos= (unsigned long int)(*pos_base)+1;
 unsigned long int pos_with_variation = int_positions[*n_variation - 1];
 int * cur_array_counts = &array_counts[*n_variation - 1];
 //  for (cur_phys_pos=1; cur_phys_pos<=w_s; cur_phys_pos++)
 //    {
 if (cur_phys_pos == pos_with_variation)
 {
 *n_tot_allele = gsl_ran_poisson(s, input_coverage);
 double ones_mean = (double) (*cur_array_counts) / (double) (n_ind);
 *n_alt_allele = gsl_ran_binomial (s, ones_mean, *n_tot_allele);
 *n_ref = *n_tot_allele - *n_alt_allele;
 *ref_base = 1;
 *alt_base = 2;
 n_alt[0] = 0;
 n_alt[1] = *n_ref;
 n_alt[2] = *n_alt_allele;
 n_alt[3] = 0;
 n_alt[4] = 0;
 cur_array_counts++;
 (*n_variation)++;
 }
 else
 {
 *n_tot_allele = gsl_ran_poisson(s, input_coverage);
 double argument = (double) (*n_tot_allele) * error_rate;
 *n_alt_allele = min(gsl_ran_poisson(s, argument), *n_tot_allele);
 *n_ref = *n_tot_allele - *n_alt_allele;
 *ref_base = 1;
 *alt_base = 2;
 n_alt[0] = 0;
 n_alt[1] = *n_ref;
 n_alt[2] = *n_alt_allele;
 n_alt[3] = 0;
 n_alt[4] = 0;
 }
 //    }
 (*pos_base)++;
 }
 */

/*--------------------------------------------------------------*/

/*
 int dmau_extract_outgroup_base(FILE * fasta_out, unsigned long pos, unsigned long oldpos, unsigned long * fasta_length, char ref_base)
 {
 char *line;
 size_t nline;
 char ref;
 char out;
 long fpos;
 
 nline=1000;
 line=(char *)malloc(nline*sizeof(char));
 fpos=ftell(fasta_out);
 getline(&line,&nline,fasta_out);
 sscanf(line,"%lu %s %s\n",fasta_length,&ref,&out);
 
 if(pos<(*fasta_length)){
 fseek(fasta_out, fpos, SEEK_SET);
 return ref_base;
 } else {
 return base_to_num(out);
 }
 };
 */

/*--------------------------------------------------------------*/
void init_codons(struct keep_codon_data *codon_pop) {
    int j;
    int na;
    for(j=0;j<3;j++) {
        codon_pop[j].n = 0;
        codon_pop[j].n_ref = 0;
        codon_pop[j].n_alt_allele = 0;
        codon_pop[j].rd = 0;
        for(na=0;na<5;na++) codon_pop[j].n_alt[na] = 0;
        codon_pop[j].ref_base = 0;
        codon_pop[j].alt_base = 0;
        codon_pop[j].out_base = 0;
    }
    return;
}

int which_Aa(char *codon) {
    int j;
    int Aa;
    
    //assign triplet and Associated Aa
    j=0;
    do {
        if((char)(codon[0])==tripletsN[j][0] &&
           (char)(codon[1])==tripletsN[j][1] &&
           (char)(codon[2])==tripletsN[j][2])
            break;
        j++;
    } while(j<64);
    
    if(j<64)
        Aa = NuclearUniversalCode[j];
    else
        Aa = -1;//unassigned

    return(Aa);
}

void calculate_NSynSynPos(char *codon, int *nsyn) {
    //calculation of Syn/Nsyn for each codon:
    //calculate for each position of the codon what is the AminoAcid given a single change (try all four nt) and look for the number of Syn and NSyn
    //the output is a vector with the number of NSyn (4 - Syn) per position (3 values between 0 and 4).
    //according Nei and Gojobori 1986

    //tripletsN[64][3]
    //NuclearUniversalCode[64]

    int i,j,k;
    int Aa,Aa_sn;
    //int syn[3];//test
    char cod_sn[3];
    char cod_obs[3];

    //Define Nsyn changes
    nsyn[0]=nsyn[1]=nsyn[2]=0;
    //syn[0]= syn[1]= syn[2]=0;//test

    //assign triplet and Associated Aa
    cod_obs[0]=(char)(codon[0]);
    cod_obs[1]=(char)(codon[1]);
    cod_obs[2]=(char)(codon[2]);
    Aa = which_Aa(cod_obs);
    if(Aa == -1) {
        nsyn[0]=nsyn[1]=nsyn[2]=-1;//unassigned
        printf("Error assigning Aa");
        exit(1);
        return;
    }
    
    //look for Nsyn changes (assume Syn+NSyn=3 ex. for A->TCG)
    for(k=0;k<3;k++) {
        //assign codon
        cod_sn[0]=(char)(cod_obs[0]);
        cod_sn[1]=(char)(cod_obs[1]);
        cod_sn[2]=(char)(cod_obs[2]);
        for(i=NUM_ASCII+1;i<=NUM_ASCII+4;i++) {
            cod_sn[k]=(char)(cod_obs[k]);
            //modify one site
            if(cod_sn[k]!=i)
                cod_sn[k]=i;
            else
                continue;
            //look for Aa
            Aa_sn = which_Aa(cod_sn);
            if(Aa_sn == -1) {
                nsyn[0]=nsyn[1]=nsyn[2]=-1;//unassigned
                printf("Error assigning Aa_sn");
                exit(1);
                return;
            }
            else { //count if Nsyn
                nsyn[k] += (Aa != Aa_sn);
                //syn[k] += (Aa == Aa_sn);//test
            }
        }
    }
    return;
}

void invert_codon(char *codon, char *inv_codon) {
    
    int i;
    // order: TCGA=1234
    for(i=0;i<3;i++) {
        switch(codon[i])
        {
            case '1':
                inv_codon[2-i] = '4';
                break;
            case '2':
                inv_codon[2-i] = '3';
                break;
            case '3':
                inv_codon[2-i] = '2';
                break;
            case '4':
                inv_codon[2-i] = '1';
                break;
            default:
                printf("Error assigning inverted codons\n");
                exit(1);
        }
    }
    return;
}

void codon_weights_muts(double *cweights, double *dweights, double *cmutnsyn, double *kmutnsyn, double *cmutsyn, double *kmutsyn, struct keep_codon_data *codon_pop, char strand) {

    char codon_init[2][3];
    char codon[4][3];
    char codon_d[3];
    int  cmut[3];
    int  kmut[3];
    int  Aa[4];

    int ncodons=0;
    int pos,i,j,k;
    int totcmut,totkmut;
    int *nsyn_sites;
    //int *nsyn_muts;

    //for each of the three positions:
    //look for the ref, alt (if exist) and out (if exist) for both pops.
    //If ref p1 and p2 does not exist then return weights ZERO
    //look at polymorphism: how many alleles. Define Codons
    
    /* //testing
    nsyn_sites=(int *)calloc(3*4,sizeof(int));
    codon[0][0]=2;
    codon[0][1]=3;
    codon[0][2]=1;
    calculate_NSynSynPos(codon[0],nsyn_sites);
    free(nsyn_sites);
    return;
    */ //end testing


    //Obtain codons for polymorphism
    for(pos=0;pos<3; pos++) {
        //no consider positions with no information
        if(codon_pop[pos].n_ref==0) {
            //codon discarded if no maf
            cweights[0] = cweights[1] = cweights[2] = -1.;
            dweights[0] = dweights[1] = dweights[2] = -1.;
            return;
        }
        //obtaining the different alleles for each position
        if(codon_pop[pos].n_ref) {
            codon_init[0][pos]=codon_pop[pos].ref_base+NUM_ASCII;
            codon_init[1][pos]=codon_pop[pos].ref_base+NUM_ASCII;
        }
        //n_alt_allele is only defined when n_ref is defined
        if(codon_pop[pos].n_alt_allele)
            codon_init[1][pos]=codon_pop[pos].alt_base+NUM_ASCII;
    }
    if(strand == '-') {
        //invert the sequence: triplet from 3 to 1 and A->T,C->G etc..
        for(i=0;i<2;i++) {
            invert_codon(codon_init[i],codon[i]);
        }
    }
    else {
        //if(strand == '+') do not change
        for(i=0;i<2;i++) {
            for(pos=0;pos<3;pos++) {
                codon[i][pos]=codon_init[i][pos];
            }
        }
    }

    //Once codons defined, look at the difference between the two codons:
    //if no diff then calculate Syn/NSyn according Li and Gojobori 1986.
    //if one diff then calculate AVERAGE Syn/NSyn.
    //if two differences than define 4 codons according the combinations and calculate AVERAGE Syn/NSyn.
    //if more than two differences then return weights ZERO.
    
    //comparison
    cmut[0]=cmut[1]=cmut[2]=0;
    cmutnsyn[0]=cmutnsyn[1]=cmutnsyn[2]=0;
    cmutsyn[0] =cmutsyn[1] =cmutsyn[2] =0;
    for(pos=0;pos<3;pos++) {
        if(codon[0][pos] != codon[1][pos])
            cmut[pos] += 1;
    }
    totcmut = cmut[0]+cmut[1]+cmut[2];

    nsyn_sites=(int *)calloc(3*4,sizeof(int));
    if(totcmut == 0) {
        ncodons = 1;
        calculate_NSynSynPos(codon[0],nsyn_sites+0);
        cweights[0] = nsyn_sites[0]/3.0;
        cweights[1] = nsyn_sites[1]/3.0;
        cweights[2] = nsyn_sites[2]/3.0;
        //calculate mutations
        cmutnsyn[0] = cmutnsyn[1] = cmutnsyn[2] = 0.;
        cmutsyn[0]  = cmutsyn[1]  = cmutsyn[2]  = 0.;
    }
    else if(totcmut == 1) {
        ncodons = 2;
        calculate_NSynSynPos(codon[0],nsyn_sites+0);
        calculate_NSynSynPos(codon[1],nsyn_sites+3);
        cweights[0] = (nsyn_sites[0]+nsyn_sites[3])/6.0;
        cweights[1] = (nsyn_sites[1]+nsyn_sites[4])/6.0;
        cweights[2] = (nsyn_sites[2]+nsyn_sites[5])/6.0;
        //calculate mutations
        Aa[0] = which_Aa(codon[0]);
        Aa[1] = which_Aa(codon[1]);
        for(pos=0;pos<3;pos++) {
            if(cmut[pos]) {
                if(Aa[0]!=Aa[1]) {cmutnsyn[pos] = 1.; cmutsyn[pos] = 0.;}
                else {cmutnsyn[pos] = 0.; cmutsyn[pos] = 1.;}
            }
            else {cmutnsyn[pos] = 0.; cmutsyn[pos] = 0.;}
        }
    }
    else if(totcmut == 2) {
        /*
        // preliminary: codon discarded if more than 1 mutation
        //printf("Codon eliminated-pol\n");
        cweights[0] = cweights[1] = cweights[2] = -1.;
        dweights[0] = dweights[1] = dweights[2] = -1.;
        free(nsyn_sites);
        return;
        *//**/
        ncodons = 4;
        //define four codons
        //Ex: ATC,AGC,ATT and AGT: ATC-AGC-AGT|ATC-ATT-AGT: 2 ways. Same for AGC to ATT.
        //pos2: aTc-aGc[-aGt NO CHANGE] | [NO CHANGE aTc-]aTt-aGt
        //pos2: aGc-aTc[-aTt NO CHANGE] | [NO CHANGE aGc-]aGt-aTt
        //both are equal
        i=0;
        for(pos=0;pos<3;pos++) {
            if(cmut[pos] && i==0) {
                codon[2][pos]=codon[0][pos];
                codon[3][pos]=codon[1][pos];
                i=1;
            }
            else {
                if(cmut[pos] && i==1) {
                    codon[2][pos]=codon[1][pos];
                    codon[3][pos]=codon[0][pos];
                }
                else {
                    if(cmut[pos] == 0) {
                        codon[2][pos]=codon[3][pos]=codon[0][pos];
                    }
                }
            }
        }
        calculate_NSynSynPos(codon[0],nsyn_sites+0);
        calculate_NSynSynPos(codon[1],nsyn_sites+3);
        calculate_NSynSynPos(codon[2],nsyn_sites+6);
        calculate_NSynSynPos(codon[3],nsyn_sites+9);
        cweights[0] = (nsyn_sites[0]+nsyn_sites[3]+nsyn_sites[6]+nsyn_sites[9] )/12.0;
        cweights[1] = (nsyn_sites[1]+nsyn_sites[4]+nsyn_sites[7]+nsyn_sites[10])/12.0;
        cweights[2] = (nsyn_sites[2]+nsyn_sites[5]+nsyn_sites[8]+nsyn_sites[11])/12.0;
        //calculate mutations
        for(i=0;i<4;i++)
            Aa[i] = which_Aa(codon[i]);
        for(pos=0;pos<3;pos++) {
            if(cmut[pos]) {
                //count only 2 ways with 2 comparisons per position
                double countc_n = 0.;
                double countc_s = 0.;
                for(i=0;i<3;i++) {
                    for(j=i+1;j<4;j++) {
                        if((i==0 && j==1) ||
                           (i==2 && j==3))
                            continue; //diff two mutations: excluded
                        if(codon[i][pos]!=codon[j][pos]) {//count only diffs in pos
                            if(Aa[i]!=Aa[j]) countc_n += 1.;
                            else countc_s += 1.;
                        }
                    }
                }
                cmutnsyn[pos] = countc_n / 2.;
                cmutsyn[pos]  = countc_s / 2.;
            }
            else {cmutnsyn[pos] = 0.; cmutsyn[pos] = 0.;}
        }
        /**/
    }
    else if(totcmut > 2) {
        // codon discarded if more than 2 mutations
        cweights[0] = cweights[1] = cweights[2] = -1.;
        dweights[0] = dweights[1] = dweights[2] = -1.;
        free(nsyn_sites);
        return;
    }

    //Now calculate divergence weigths, if necessary.
    //out_base is EQUAL in p1 than in p2.
    //Same procedure:
    //calculate codons for divergence versus the polymorphism
    //weight calculated using average pol. versus div.
    //assume divergence has only one codon with no variants
    //case1: if pol has one codon (ex. AAA) divergence has one (ex.AAA|AAT|ATA|ATT|TTT[null]).
    //       do mean pol. vs div. except for ATT -> calculate all 4 cases (2 more and average)
    //case2: if pol has two codons, do always average pol-div. Only define one codon. in div
    //       ex. pol. AAA|AAT. if div AAA (mean 2pol ve 1div), AAT (same), TTT (null),
    //                                ATA (path is AAA-ATA|AAT, i.e., mean 2pol 1div),
    //                                ATT (path is AAA-AAT-ATT, i.e., mean 2pol 1div),
    //case3: if pol has 4 codons, look if coincide, if not define null (mean 4pol 1div)
    
    //Obtain codon for divergence
    for(pos=0;pos<3; pos++) {
        if(codon_pop[pos].out_base==0) {
            //codon discarded if no maf
            dweights[0] = dweights[1] = dweights[2] = -1.;
            free(nsyn_sites);
            return;
        }
        //obtaining the different alleles for each position
        codon_init[0][pos]=codon_pop[pos].out_base+NUM_ASCII;
    }
    if(strand == '-') {
        //invert the sequence: triplet from 3 to 1 and A->T,C->G etc..
        invert_codon(codon_init[0],codon_d);
    }
    else {
        for(pos=0;pos<3; pos++) {
            codon_d[pos] = codon_init[0][pos];
        }
    }

    //calculate net number of differences vs outgroup
    kmut[0]=kmut[1]=kmut[2]=0;
    kmutnsyn[0]=kmutnsyn[1]=kmutnsyn[2]=0;
    kmutsyn[0] =kmutsyn[1] =kmutsyn[2] =0;
    for(pos=0;pos<3;pos++) {
        for(i=0;i<ncodons;i++) {
            if(codon[i][pos] != codon_d[pos]) {
                kmut[pos] = 1;
                break; /*assign 1 if whatever codon is different with outgroup*/
            }
        }
    }
    totkmut = kmut[0]+kmut[1]+kmut[2];
    //estimate weights
    if(totkmut==totcmut) {
        //differ in polymorphism, not in divergence
        for(pos=0;pos<3;pos++) {
            dweights[pos] = cweights[pos]; //assume equal
            kmutnsyn[pos] = cmutnsyn[pos];
            kmutsyn[pos]  = cmutsyn[pos];
         }
    }
    /*if(totkmut==0) {//redundant
        dweights[0] = cweights[0];
        dweights[1] = cweights[1];
        dweights[2] = cweights[2];
        //calculate differences
        kmutnsyn[0] = kmutnsyn[1] = kmutnsyn[2] = 0.;
        kmutsyn[0]  = kmutsyn[1]  = kmutsyn[2]  = 0.;
    }*/
    else {
        if(totkmut==1) { //only 1 or 2 codons defined}
            //at least one mutation comes from divergence
            calculate_NSynSynPos(codon_d,nsyn_sites+0);
            dweights[0] = (cweights[0] + nsyn_sites[0]/3.0)/2.0;
            dweights[1] = (cweights[1] + nsyn_sites[1]/3.0)/2.0;
            dweights[2] = (cweights[2] + nsyn_sites[2]/3.0)/2.0;
            //calculate differences
            Aa[0] = which_Aa(codon[0]);
            Aa[2] = which_Aa(codon_d);
            if(ncodons==1) Aa[1] = which_Aa(codon[0]);//repeated
            if(ncodons==2) Aa[1] = which_Aa(codon[1]);
            for(pos=0;pos<3;pos++) {
                if(kmut[pos]) {
                    //count 2 combinations
                    double countk_n = 0.;
                    double countk_s = 0.;
                    for(i=0;i<2;i++) {
                        if(Aa[i]!=Aa[2]) countk_n += 1.;
                        else countk_s += 1.;
                    }
                    kmutnsyn[pos] = countk_n / 2.;
                    kmutsyn[pos]  = countk_s / 2.;
                }
                else {kmutnsyn[pos] = 0.; kmutsyn[pos] = 0.;
                }
            }
        }
        else if(totkmut==2) {
            /*
             //preliminar: codon discarded if more than 1 mutations
             //printf("Codon eliminated-div\n");
             dweights[0] = dweights[1] = dweights[2] = -1.;
             free(nsyn_sites);
             return;
             *//**/
            if(ncodons==1) {
                //define two more codons within pol
                //Ex: ATC and AGT: ways are ATC-AGC-AGT|ATC-ATT-AGT: 2 ways
                //pos2: aTc-aGc[-aGt NO CHANGE] | [NO CHANGE aTc-]aTt-aGt
                i=0;
                for(pos=0;pos<3;pos++) {
                    if(kmut[pos] && i==0) {
                        codon[2][pos]=codon[0][pos];
                        codon[3][pos]=codon_d[pos];
                        i=1;
                    }
                    else {
                        if(kmut[pos] && i==1) {
                            codon[2][pos]=codon_d[pos];
                            codon[3][pos]=codon[0][pos];
                        }
                        else {
                            if(kmut[pos] == 0) {
                                codon[2][pos]=codon[3][pos]=codon[0][pos];
                            }
                        }
                    }
                    codon[1][pos]=codon_d[pos]; //replicate codon_d in codon[1] (before empty now two diff with codon[0])
                }
                //codon[1] is equal than codon_d !
                calculate_NSynSynPos(codon[0],nsyn_sites+0);
                calculate_NSynSynPos(codon[1],nsyn_sites+3);
                calculate_NSynSynPos(codon[2],nsyn_sites+6);
                calculate_NSynSynPos(codon[3],nsyn_sites+9);
                dweights[0] = (nsyn_sites[0]+nsyn_sites[3]+nsyn_sites[6]+nsyn_sites[9] )/12.0;
                dweights[1] = (nsyn_sites[1]+nsyn_sites[4]+nsyn_sites[7]+nsyn_sites[10])/12.0;
                dweights[2] = (nsyn_sites[2]+nsyn_sites[5]+nsyn_sites[8]+nsyn_sites[11])/12.0;
                //calculate differences
                Aa[0] = which_Aa(codon[0]);
                Aa[1] = which_Aa(codon[1]);
                Aa[2] = which_Aa(codon[2]);
                Aa[3] = which_Aa(codon[3]);
                for(pos=0;pos<3;pos++) {
                    if(kmut[pos]) {
                        //count only 2 ways with 4 effective combinations (two per pos)
                        double countk_n = 0.;
                        double countk_s = 0.;
                        for(i=0;i<3;i++) {
                            for(j=i+1;j<4;j++) {
                                if((i==0 && j==1) ||
                                   (i==2 && j==3))
                                    continue; //diff two mutations: excluded
                                if(codon[i][pos]!=codon[j][pos]) {//count only diffs in pos
                                    if(Aa[i]!=Aa[j]) countk_n += 1.;
                                    else countk_s += 1.;
                                }
                            }
                        }
                        kmutnsyn[pos] = countk_n / 2.;
                        kmutsyn[pos]  = countk_s / 2.;
                    }
                    else {kmutnsyn[pos] = 0.; kmutsyn[pos] = 0.;}
                }
            }
            if(ncodons==2) { //if ncodons is larger than 2 then totkmut==totcmut is true
                calculate_NSynSynPos(codon_d,nsyn_sites+0);
                dweights[0] = (cweights[0] + nsyn_sites[0]/3.0)/2.0;
                dweights[1] = (cweights[1] + nsyn_sites[1]/3.0)/2.0;
                dweights[2] = (cweights[2] + nsyn_sites[2]/3.0)/2.0;
                //Ex. AAA,AAT and out ATA: only ONE way AAT-AAA-ATA
                //Ex. AAA,AAT and out ATT: only ONE way AAA-AAT-ATT
                //Ex. AAA,ATT and out ATA, or AAA or AAT: NOT possible because totkmut==totcmut is true
                //calculate differences:
                for(pos=0;pos<3;pos++)
                    codon[3][pos]=codon_d[pos];//define codon3 for outgroup
                for(i=0;i<2;i++) {
                    k=0;
                    for(pos=0;pos<3;pos++) {
                        if(codon[i][pos]!=codon[3][pos])
                            k++;
                    }
                    if(k==2)
                        break; //codon[i] vs codon[3] are extremes of the path
                }
                if(i==0) j=1; else j=0; //codon[j] has 1 diff with codon[i] and 1 diff with codon[3]
                Aa[0] = which_Aa(codon[i]);
                Aa[1] = which_Aa(codon[j]);
                Aa[2] = which_Aa(codon_d);
                for(pos=0;pos<3;pos++) {
                    if(kmut[pos]) {
                        //count 1 combination (one for each pos)
                        double countk_n = 0.;
                        double countk_s = 0.;
                        if(codon[i][pos]!=codon[j][pos]) {//count only diffs in pos
                            if(Aa[i]!=Aa[j]) countk_n += 1.;
                            else countk_s += 1.;
                        }
                        if(codon[j][pos]!=codon[3][pos]) {//count only diffs in pos
                            if(Aa[i]!=Aa[j]) countk_n += 1.;
                            else countk_s += 1.;
                        }
                        kmutnsyn[pos] = countk_n;
                        kmutsyn[pos]  = countk_s;
                    }
                    else {kmutnsyn[pos] = 0.; kmutsyn[pos] = 0.;}
                }
            }
        }
        else if(totkmut>2) {
            //codon discarded if more than 2 mutations
            dweights[0] = dweights[1] = dweights[2] = -1.;
            free(nsyn_sites);
            return;
        }
    }
    free(nsyn_sites);
    return;
}


void extract_stats(struct tests * test, struct combinatorial * comb, int n0, unsigned long n_ref, unsigned long n_alt_allele, unsigned long rd, unsigned long * n_alt, int ref_base, int alt_base, int out_base, int mb, float weightp, float weightd, double cmut)
{
    char is_out=0;
    float gamma;
    int i;
    
    if(rd > 2*mb+1) {
        if(weightp<0.) weightp=0.;
        if(weightd<0.) weightd=0.;
        
        if (out_base>0)
        {
            if (ref_base==out_base) { is_out=1; }
            else if (n_alt[out_base]==n_alt_allele) { is_out=2; };
            // else is_out=3;
        };
        
        gamma = fmax(((float)rd/(float)n0),1.0);
        
        if(weightp) {
            test->cov+=rd;
            test->l+= 1.0 * weightp;
            test->den_t+= weightp * comb->d_t[rd-1];
            test->den_p+= weightp * comb->d_p[rd-1];
            test->den_nu+= weightp * comb->d_nu[rd-1];
            if (is_out>0 && weightd)
            {
                test->l_out+=1.0 * weightd;
                test->den_hl+= weightd * comb->d_hl[rd-1];
                test->den_hq+= weightd * comb->d_hq[rd-1];
                test->den_xi+= weightd * comb->d_xi[rd-1];
            };
        }
        // else printf("outgroup bases %u %u %u\n",is_out,out_base,ref_base);//debug
        
        if ((n_ref>mb)&&(n_ref<rd-mb) && weightp)
        {
            double sum_nr=0.;
            for (i=1;i<5;i++) {
                if(n_alt[i]>mb) sum_nr += n_alt[i];
            }

            test->s+=1.0 * cmut;
            test->num_t+=1.0 * cmut;
            test->num_p+=(double)(2*n_alt_allele*n_ref)/ (double)/**/(sum_nr*(sum_nr-1))/**/ /*(rd*(rd-1))*/ * cmut;
            
            if((float)n_alt_allele<=gamma)
                test->num_nu+=1.0 * cmut;
            
            if (is_out==1 && weightd){
                test->num_hl+=(double)(n_alt_allele)/(double)/**/(sum_nr-1)/**/ /*(rd-1)*/ * cmut;
                test->num_hq+=(double)(n_alt_allele*n_alt_allele)/(double)/**/(sum_nr*(sum_nr-1))/**/ /*(rd*(rd-1))*/ * cmut;
                if((float)n_alt_allele<=gamma)
                    test->num_xi+=1.0 * cmut;
            };
            if (is_out==2 && weightd)
            {
                test->num_hl+=(double)(n_ref)/(double)/**/(sum_nr-1)/**/ /*(rd-1)*/ * cmut;
                test->num_hq+=(double)(n_ref*n_ref)/(double) /**/(sum_nr*(sum_nr-1))/**/ /*(rd*(rd-1))*/ * cmut;
                if((float)n_ref<=gamma)
                    test->num_xi+=1.0 * cmut;
            }
        };
    }
    
    DEB(printf("using data %u\t%lu\t%lu\t%lu\t%lu\t%lu\t%lu\t%lu\t%lu\t%u\t%u\t%u\n", n0, n_ref, n_alt_allele, rd, n_alt[0], n_alt[1],n_alt[2],n_alt[3],n_alt[4], ref_base, alt_base, out_base)); //debug
    
    return;
};

/*--------------------------------------------------------------*/

unsigned long extract_pos_snpinput( FILE * snpinput, /*ADD 20240612*/char *schrom){
    char *line;
    //int c,i,j;
    size_t l_line;
    unsigned long pos;
    line = (char *) malloc (100*sizeof(char));
    l_line=100;
    char test; test=fgetc(snpinput);
    if (test!=EOF) {
        ungetc(test, snpinput);
        //getline(&line,&l_line,snpinput);
        //for(i=1;i<=column;i++){
        getline(&line,&l_line,snpinput);
        //};
        //BEGIN MOD 20240612
        sscanf(line,"%s\t%lu",schrom,&pos);
        //sscanf(line,"%lu",&pos);//sscanf(line,"%lu\t",&pos); // or sscanf(line,"%lu ",pos);
        //END MOD 20240612
        //i=0;
        //for(c=1;c<column;c++){
        // for(;(i<l_line)&&((line[i]==9)||(line[i]==32))&&((line[i+1]!=9)&&(line[i+1]!=32));i++){};
        //};
        //j=i;
        //for(;(i<l_line)&&((line[i]!=9)&&(line[i]!=32))&&((line[i+1]==9)||(line[i+1]==32));i++){};
        //getline(&line,&l_line,snpinput);
    } else {
        pos=2000000000;
    };
    free(line);
    //return pos;
    return(pos);
};

/*FST*/
void extract_fst(struct fst_calc * fst, struct combinatorial_fst * combfst, int n01, unsigned long n_ref1, unsigned long n_alt_allele1, unsigned long rd1, unsigned long * n_alt_1, int ref_base1, int alt_base1, int n02, unsigned long n_ref2, unsigned long n_alt_allele2, unsigned long rd2, unsigned long * n_alt_2, int ref_base2, int alt_base2, int out_base, int mb, int *snp)
{
    int i,j;
    *snp = 0;
    
    fst->la+=1;
    for (i=1;i<5;i++)
    {
        for (j=1;j<5;j++)
        {
            if((i!=j && n_alt_1[i]>mb)&&(n_alt_2[j]>mb)) {
                fst->gen_diff+=(double)(n_alt_1[i]*n_alt_2[j])/(double)(rd1*rd2);
                *snp = 1;
            }
        };
    };
    fst->c_s+=combfst->c_s[rd1-1][rd2-1];
};
void extract_fstT(struct fst_calc * fst, struct combinatorial_fst * combfst, struct combinatorial *comb, struct combinatorial *comb2, int n01, unsigned long n_ref1, unsigned long n_alt_allele1, unsigned long rd1, unsigned long * n_alt_1, int ref_base1, int alt_base1, int n02, unsigned long n_ref2, unsigned long n_alt_allele2, unsigned long rd2, unsigned long * n_alt_2, int ref_base2, int alt_base2, int out_base, int mb,int mbf, int *snp, int fst_path)
{
    int i,j;
    //int nal;
    //unsigned long int a[5],at,a2;
    //double local_piT;
    unsigned long sum_nr1=0;
    unsigned long sum_nr2=0;

    //a[1]=a[2]=a[3]=a[4]=0;
    *snp = 0;
    
    if(rd1>2*mb+1 && rd2>2*mb+1) {
    //if(rd1>mb && rd2>mb) {
        for (i=1;i<5;i++) {
            if(n_alt_1[i]>mb) sum_nr1 += n_alt_1[i];
            if(n_alt_2[i]>mb) sum_nr2 += n_alt_2[i];
        }
        //if(rd1!=sum_nr1 || rd2!=sum_nr2)
        //    printf("\nrd1: %lu sunnr1: %lu rd2: %lu sumnr2: %lu",rd1,sum_nr1,rd2,sum_nr2);
        
        //calculate diff for Pia
        for (i=1;i<5;i++)
        {
            for (j=1;j<5;j++)
            {
                if(i!=j && n_alt_1[i]>mb && n_alt_2[j]>mb) {
                    fst->gen_diff+=(double)(n_alt_1[i] * n_alt_2[j])/((double)sum_nr1/*rd1*/ * (double)sum_nr2 /*rd2*/);
                    *snp = 1;
                }
            };
            if((n_alt_1[i]>mb || n_alt_1[i]==0) && (n_alt_2[i]>mb || n_alt_2[i]==0))
                fst->var_diffL+=(pow((double)(n_alt_1[i]/(double)sum_nr1 - n_alt_2[i]/(double)sum_nr2),2))/4.0; //(/2?)
        };

        fst->la+=1;
        
        if(fst_path==0) {
            if(combfst->c_s[rd1-1][rd2-1] == 0) {
                combfst->c_s[rd1-1][rd2-1] = calculate_combfst_cs((unsigned int)rd1,(unsigned int)rd2,mb,n01,n02);
            }
            fst->c_s+=combfst->c_s[rd1-1][rd2-1];
        }
        
    //}
    /**/
    //if(rd1>2*mb+1 && rd2>2*mb+1) {    //calculate piW and piT using the same positions and conditions (rd1+rd2>2*mb+1)
        //number of different alleles
        /*
        nal=0;
        at=0;
        for (i=1;i<5;i++) {
            if(n_alt_1[i] + n_alt_2[i] > mb) {
                a[i] = n_alt_1[i] + n_alt_2[i];
                at+=a[i];
                nal++;
            }
        }
        if(nal>2) //not counted for more than two alleles
            return;
         */
        //if ((rd1>2*mb+1)) {
            double fnp1=0;
            fst->den_p1+= comb->d_p[rd1-1];
            if ((n_ref1>mb) && (n_alt_allele1>mb)) {
                fnp1 += (double)(2 * n_alt_allele1 * n_ref1)/((double)sum_nr1/*rd1*/ * (double)(sum_nr1/*rd1*/-1.0));
                //the factor (rd1/(rd1-1.0) is implicit when 2pq is calculated
                fst->num_p1 += fnp1;
            }
            fst->num_p1g += ((double)n01/(n01-1.0))*fnp1;
        //}
        //if ((rd2>2*mb+1)) {
            double fnp2=0;
            fst->den_p2+= comb2->d_p[rd2-1];
            if ((n_ref2>mb) && (n_alt_allele2>mb)) {
                fnp2 += (double)(2 * n_alt_allele2 * n_ref2)/((double)sum_nr2/*rd2*/ * (double)(sum_nr2/*rd2*/-1.0));
                //the factor (rd2/(rd2-1.0) is implicit when 2pq is calculated
                fst->num_p2 += fnp2;
            }
            fst->num_p2g += ((double)n02/(n02-1.0))*fnp2;
        //}
        /*
        local_piT=0.;
        for (i=1;i<5;i++) {
            a2 = 0;
            for (j=1;j<5;j++) {
                if(i!=j && a[j]>mb)
                    a2 += a[j];
            }
            if(a[i]>mb && a2>mb)
                local_piT += (double)(a[i]*a2);
        }
        fst->num_pt += local_piT/(at*(at-1));
        //fst->den_pt+= combfst->d_pt[rd1+rd2-1];
        fst->lt+=1;
        */
    }
    /**/
};



//SNPS
void extract_snps(unsigned long pos, FILE * output_snps, unsigned long * n_alt_1, unsigned long * n_alt_2, int mb)
{
    int i, c1, c2;
    
    c2=0;
    c1=0;
    for(i=1;i<5;i++) {
        if ((n_alt_1[i]+n_alt_2[i])>(n_alt_1[c1]+n_alt_2[c1])) {
            c1=i;
        };
    };
    if (c1==0) c2=1;
    for(i=c2+1;i<5;i++) {
        if (((n_alt_1[i]+n_alt_2[i])>(n_alt_1[c2]+n_alt_2[c2]))&&(i!=c1)) {
            c2=i;
        };
    };
    if((n_alt_1[c2]+n_alt_2[c2])>1) {
        fprintf(output_snps, "%lu\t%lu\t%lu\t%lu\t%lu\n", pos, n_alt_1[c2], n_alt_2[c2], n_alt_1[c1], n_alt_2[c1]);
    };
};

unsigned long extract_gff(FILE *gffinput, char *gchrom, unsigned long *cds_start, unsigned int *phase_cds, char *strand) {
    
    unsigned long cds_end;
    char feature[256];
    char *line_gff;
    size_t n_line_gff;
    
    n_line_gff=1;
    line_gff=malloc(sizeof(char));
    if(getline(&line_gff,&n_line_gff,gffinput)!=-1){
        while ((line_gff[0]=='#')||(line_gff[0]=='\0')) {
            getline(&line_gff,&n_line_gff,gffinput);
        };
        //sscanf(line_gff,"%*s\t%*s\t%s\t%lu\t%lu\t%*s\t%c\t%u\t",feature,cds_start,&cds_end,strand,phase_cds);
        sscanf(line_gff,"%s\t%*s\t%s\t%lu\t%lu\t%*s\t%c\t%u\t",gchrom,feature,cds_start,&cds_end,strand,phase_cds);
        if((strcmp(feature,"CDS")!=0)&&(strcmp(feature,"cds")!=0)) {
            cds_end=0;
        };
    } else {cds_end=2000000000;};
    free(line_gff);
    
    return(cds_end);
};

//better do this function for a single pop and execute two commands in main
void init_tests(struct tests *test1, struct tests *test2, struct tests *tests1, struct tests *testn1, struct tests *tests2, struct tests *testn2, struct fst_calc *fst, unsigned long *vec_rd, unsigned long *vec2_rd, double *psyn, double *pnon, double *dsyn, double *dnon, double *psyn2, double *pnon2, double *dsyn2, double *dnon2, unsigned long *div, unsigned long *div2, unsigned long max_cov, int compute_fst) {
    
    unsigned long rd;
    
    test1->cov=0;
    test1->l=0;
    test1->l_out=0;
    test1->s=0;
    test1->num_t=0;
    test1->num_p=0;
    test1->num_hl=0;
    test1->num_hq=0;
    test1->num_nu=0;
    test1->num_xi=0;
    test1->den_t=0;
    test1->den_p=0;
    test1->den_hl=0;
    test1->den_hq=0;
    test1->den_nu=0;
    test1->den_xi=0;
    
    *div=0;
    
    for(rd=1;rd<=max_cov; rd++){
        vec_rd[rd-1]=0;
        if(compute_fst) vec2_rd[rd-1]=0;
        // vec_s[rd-1]=0;
        // vec_p[rd-1]=0;
        // vec_h[rd-1]=0;
    };
    *psyn=0;
    *dsyn=0;
    *pnon=0;
    *dnon=0;
    
    tests1->cov=0;
    tests1->l=0;
    tests1->l_out=0;
    tests1->s=0;
    tests1->num_t=0;
    tests1->num_p=0;
    tests1->num_hl=0;
    tests1->num_hq=0;
    tests1->num_nu=0;
    tests1->num_xi=0;
    tests1->den_t=0;
    tests1->den_p=0;
    tests1->den_hl=0;
    tests1->den_hq=0;
    tests1->den_nu=0;
    tests1->den_xi=0;
    
    testn1->cov=0;
    testn1->l=0;
    testn1->l_out=0;
    testn1->s=0;
    testn1->num_t=0;
    testn1->num_p=0;
    testn1->num_hl=0;
    testn1->num_hq=0;
    testn1->num_nu=0;
    testn1->num_xi=0;
    testn1->den_t=0;
    testn1->den_p=0;
    testn1->den_hl=0;
    testn1->den_hq=0;
    testn1->den_nu=0;
    testn1->den_xi=0;
    
    if(compute_fst) {
        *div2=0;
        
        *psyn2=0;
        *dsyn2=0;
        *pnon2=0;
        *dnon2=0;
        
        tests2->cov=0;
        tests2->l=0;
        tests2->l_out=0;
        tests2->s=0;
        tests2->num_t=0;
        tests2->num_p=0;
        tests2->num_hl=0;
        tests2->num_hq=0;
        tests2->num_nu=0;
        tests2->num_xi=0;
        tests2->den_t=0;
        tests2->den_p=0;
        tests2->den_hl=0;
        tests2->den_hq=0;
        tests2->den_nu=0;
        tests2->den_xi=0;
        
        testn2->cov=0;
        testn2->l=0;
        testn2->l_out=0;
        testn2->s=0;
        testn2->num_t=0;
        testn2->num_p=0;
        testn2->num_hl=0;
        testn2->num_hq=0;
        testn2->num_nu=0;
        testn2->num_xi=0;
        testn2->den_t=0;
        testn2->den_p=0;
        testn2->den_hl=0;
        testn2->den_hq=0;
        testn2->den_nu=0;
        testn2->den_xi=0;
        
        test2->cov=0;
        test2->l=0;
        test2->l_out=0;
        test2->s=0;
        test2->num_t=0;
        test2->num_p=0;
        test2->num_hl=0;
        test2->num_hq=0;
        test2->num_nu=0;
        test2->num_xi=0;
        test2->den_t=0;
        test2->den_p=0;
        test2->den_hl=0;
        test2->den_hq=0;
        test2->den_nu=0;
        test2->den_xi=0;
        
        fst->la=0;
        //fst->lt=0;
        fst->gen_diff=0;
        fst->var_diffL=0;
        fst->c_s=0;
        fst->Lbias=0;
        fst->num_p1=0;
        fst->num_p1g=0;
        fst->num_p2=0;
        fst->num_p2g=0;
        //fst->num_pt=0;
        fst->den_p1=0;
        fst->den_p2=0;
        //fst->den_pt=0;
    }
    
}

void calculate_vartests(double *var_h, double *var_d, double *var_s, double *var0_s, double *var0_d, double *var0_h, struct tests *test1, unsigned long *vec_rd, double *vec_s, double *vec_p, double *vec_h, double *vec0_s, double *vec0_d, double *vec0_h, double *covmat, int n01, unsigned long min_cov, unsigned long  max_cov, int m_bar) {
    
    int k,lt;
    unsigned long rd;
    double vk_s[n01-1], vk_d[n01-1], vk_h[n01-1];
    /*double vk_l[n01-1],var_l,var0_l;*/
    /*double vk_n[n01-1],var_n,var0_n;*/
    /*double vk_x[n01-1],var_x,var0_x;*/
    
    *var_s=0;
    *var_d=0;
    *var_h=0;
    *var0_s=0;
    *var0_d=0;
    *var0_h=0;
    
    if ((test1->den_t>0)&&(test1->den_p>0)) {
        if (test1->den_hq>0) {
            for(k=1;k<n01;k++){
                vk_s[k-1]=0;
                vk_d[k-1]=0;
                vk_h[k-1]=0;
                for(rd=max(min_cov,2*m_bar+2);rd<=max_cov; rd++){
                    vk_s[k-1]+=vec_rd[rd-1]*vec_s[(k-1)+(n01-1)*(rd-1)];
                    vk_d[k-1]+=vec_rd[rd-1]*(vec_p[(k-1)+(n01-1)*(rd-1)]/test1->den_p-vec_s[(k-1)+(n01-1)*(rd-1)]/test1->den_t);
                    vk_h[k-1]+=vec_rd[rd-1]*(vec_h[(k-1)+(n01-1)*(rd-1)]/test1->den_hq-vec_p[(k-1)+(n01-1)*(rd-1)]/test1->den_p);
                };
            };
            for(rd=max(min_cov,2*m_bar+2);rd<=max_cov; rd++){
                *var0_s+=vec_rd[rd-1]*vec0_s[rd-1];
                *var0_d+=vec_rd[rd-1]*vec0_d[rd-1];
                *var0_h+=vec_rd[rd-1]*vec0_h[rd-1];
            };
            for(k=1;k<n01;k++){
                for(lt=1;lt<n01;lt++){
                    *var_s+=vk_s[k-1]*vk_s[lt-1]*covmat[(n01-1)*(lt-1)+k-1];
                    *var_d+=vk_d[k-1]*vk_d[lt-1]*covmat[(n01-1)*(lt-1)+k-1];
                    *var_h+=vk_h[k-1]*vk_h[lt-1]*covmat[(n01-1)*(lt-1)+k-1];
                };
            };
            *var0_d=*var0_d/(double)(test1->l*test1->l);
            *var0_h=*var0_h/(double)(test1->l_out*test1->l_out);
        } else {
            for(k=1;k<n01;k++){
                vk_s[k-1]=0;
                vk_d[k-1]=0;
                for(rd=max(min_cov,2*m_bar+2);rd<=max_cov; rd++){
                    vk_s[k-1]+=vec_rd[rd-1]*vec_s[(k-1)+(n01-1)*(rd-1)];
                    vk_d[k-1]+=vec_rd[rd-1]*(vec_p[(k-1)+(n01-1)*(rd-1)]/test1->den_p-vec_s[(k-1)+(n01-1)*(rd-1)]/test1->den_t);
                };
            };
            for(rd=max(min_cov,2*m_bar+2);rd<=max_cov; rd++){
                *var0_s+=vec_rd[rd-1]*vec0_s[rd-1];
                *var0_d+=vec_rd[rd-1]*vec0_d[rd-1];
            };
            for(k=1;k<n01;k++){
                for(lt=1;lt<n01;lt++){
                    *var_s+=vk_s[k-1]*vk_s[lt-1]*covmat[(n01-1)*(lt-1)+k-1];
                    *var_d+=vk_d[k-1]*vk_d[lt-1]*covmat[(n01-1)*(lt-1)+k-1];
                };
            };
            *var0_d=*var0_d/(double)(test1->l*test1->l);
        };
    };
}

void calculate_window_stats(struct tests *test1, double *cov1_val, double *theta1_val, double *pi1_val, double *d1_val, double *thetaH_val, double *h1_val, double *div_val, unsigned long div/**/, double *thetaE_val/**/, double *thetanu_val, double *thetaxi_val, double *e1_val, double *f1_val, double *fo1_val/**/ )
{
    //less = than 100 effective positions are not considered
    if(test1->l>100) { *cov1_val=(double)(test1->cov)/(double)(test1->l); } else { *cov1_val=-1; };
    if(test1->den_t>100) { *theta1_val=test1->num_t/test1->den_t; } else { *theta1_val=-1; };
    if(test1->den_p>100) { *pi1_val=test1->num_p/test1->den_p; } else { *pi1_val=-1; };
    if((test1->den_t>100)&&(test1->den_p>0)) { *d1_val=*pi1_val-*theta1_val; } else { *d1_val=-1; };
    if(test1->den_hq>100) { *thetaH_val=test1->num_hq/test1->den_hq; } else { *thetaH_val=-1; };
    if((test1->den_p>100)&&(test1->den_hq>0)) { *h1_val=*pi1_val-*thetaH_val; } else { *h1_val=-1; };
    if(test1->den_hl>100) { *thetaE_val=test1->num_hl/test1->den_hl; } else { *thetaE_val=0; };
    /**/if(test1->den_nu>100) { *thetanu_val=test1->num_nu/test1->den_nu; } else { *thetanu_val=-1; };/**/
    /**/if(test1->den_xi>100) { *thetaxi_val=test1->num_xi/test1->den_xi; } else { *thetaxi_val=-1; };/**/
    /**/if((test1->den_t>100)&&(test1->den_hl>0)) { *e1_val=*theta1_val-*thetaE_val; } else { *e1_val=-1; };/**/
    /**/if((test1->den_t>100)&&(test1->den_nu>0)) { *f1_val=*theta1_val-*thetanu_val; } else { *f1_val=-1; };/**/
    /**/if((test1->den_t>100)&&(test1->den_xi>0)) { *fo1_val=*theta1_val-*thetaxi_val; } else { *fo1_val=-1; };/**/
    if(test1->l_out>100) { *div_val=(double)(div)/(double)(test1->l_out); } else { *div_val=-1; };
}

double calculate_combfst_cs(unsigned int rd1, unsigned int rd2, int m_bar, int n01, int n02) {
    int i,k;
    double x1,y1,x2,y2;
    int m1t,m2t;
    double combfst_cs=0.;
    for(k=1;k<=n01+n02-1;k++) { //maximum sample sizes for the sum of the two pops
        for(i=0;i<=k;i++) { //x1 and y1 are the freqs k-i and i-freqs given total, x2 are 1 - these freqs
            x1=(double)(k-i)/(double)(n01+n02);
            x2=(double)(i)/(double)(n01+n02);
            y1=1-x1;
            y2=1-x2;
            for(m1t=0;m1t<=m_bar;m1t++){ //frequencies from 0 to m_bar (more reads than in one pop) are not considered
                for(m2t=0;m2t<=m_bar;m2t++){
                    combfst_cs +=
                    gsl_ran_hypergeometric_pdf(k-i,n01,n02,k)*((double)(m1t*(rd2-m2t)+m2t*(rd1-m1t))/(double)(rd1*rd2))*
                    gsl_sf_choose(rd1,m1t)*gsl_sf_choose(rd2,m2t)*
                    (gsl_pow_int(x2,m1t)*gsl_pow_int(y2,rd1-m1t)+gsl_pow_int(x2,rd1-m1t)*gsl_pow_int(y2,m1t))*
                    (gsl_pow_int(x1,m2t)*gsl_pow_int(y1,rd2-m2t)+gsl_pow_int(x1,rd2-m2t)*gsl_pow_int(y1,m2t))
                    /(double)k; //SI-28
                    //if(combfst_cs<0.) {
                    //    printf("combfst.c_s[%d-1][%d-1]=%f\n",rd1,rd2,combfst_cs);
                    //}
                }
            }
        }
    }
    //printf("combfst.c_s[%d-1][%d-1]=%f\n",rd1,rd2,combfst_cs);
    return(combfst_cs);
}

void usage(void) {
    fprintf(stderr,"Command:\n   npstat2 [options] file.pileup.gz\n  Options:\n    -n samplesize : haploid sample size\n    -l windowlength : window length\n    -scaffolds: name of file containing scaffold names\n     -mincov minimum_coverage : filter on minimum coverage (default 4)\n    -maxcov maximum_coverage : filter on maximum coverage (default 100)\n    -minqual minimum_base_quality : filter on base quality (default 10)\n    -nolowfreq m : filter on minimum allele count mac>m (default 2)\n    -outgroup file.fa : outgroup file in FASTA\n    -annot file.gff3 : annotation file in GFF3\n    -snpfile file.snp : consider SNPs only if present in file.snp\n    -outfile : name output file (default ends with extension '.stats.txt')\n    -fstpop2 file2.pileup : computes Fst with a second population contained in file2.pileup.gz\n    -n2 : sample size of the second population\n    -h : show command and flags\n");
    
    //printf("Command:\n   npstat2 [options] file.pileup.gz\n  Options:\n    -n samplesize : haploid sample size\n    -l windowlength : window length\n    -scaffolds: name of file containing scaffold names\n    -mincov minimum_coverage : filter on minimum coverage (default 4)\n    -maxcov maximum_coverage : filter on maximum coverage (default 100)\n    -minqual minimum_base_quality : filter on base quality (default 10)\n    -nolowfreq m : filter on minimum allele count mac>m (default 2)\n    -outgroup file.fa : outgroup file in FASTA\n    -annot file.gff3 : annotation file in GFF3\n    -snpfile file.snp : consider SNPs only if present in file.snp\n    -outfile : name output file (default ends with extension '.stats.txt')\n    -fstpop2 file2.pileup : computes Fst with a second population contained in file2.pileup.gz\n    -n2 : sample size of the second population\n    -h : show command and flags\n");
    return;
}
/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/
/*------------------------ Main program ------------------------*/
/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/

int main(int argc, char *argv[])
{
    // arguments: bam file 1, bam file 2, fasta outgroup file, window size, haploid sample size, minimum coverage, maximum coverage, minimum base quality, minimum mapping quality, low frequency alleles removed
    
    /* Input flags */
    double input_coverage1;
    /* Variables */
    
    int count_i,count_j;
    BGZF/*FILE*/ *bam_file1=0;
    BGZF/*FILE*/ *bam_file2=0;
    FILE *fasta_out=0;
    FILE *scaffold_file=0;
    FILE *bed_file=0;
    //char outgroup_available;
    char ct1;
    char ct2=EOF;
    int lt;
    unsigned long pos, oldpos, pos_base1, pos_base2;
    unsigned long n_ref1, n_alt_allele1;//, n_tot_allele1;
    unsigned long n_alt_1[5];
    unsigned long n_ref2, n_alt_allele2;//, n_tot_allele2;
    unsigned long n_alt_2[5];
    unsigned long max_cov, min_cov, min_qual, min_mqual;
    int m_bar,m_bar_fst;//, m1t, m2t;
    unsigned long window_size;
    int n01, n02;
    //unsigned long rd, rd1, rd2;
    int rd,rd1,rd2;
    FILE *output_stat1=0;
    FILE *output_stat2=0;
    FILE *output_fst=0;
    char *temp_file_name1, *temp_file_name2;
    //unsigned long
    int fasta_length;
    //int ref_base;
    int ref_base1, ref_base2, alt_base1, alt_base2, out_base;
    /* Gipo's variables */
    //long int n_pos1;
    //int array_counts1[2000] = {0};
    //unsigned long int int_positions1[1000] = {0};
    int n_variation1;
    n_variation1=1;
    //double error_rate1 = 0.001;
    
    // Random number generator
    const gsl_rng_type * type;
    gsl_rng * r;
    gsl_rng_env_setup();
    type = gsl_rng_default;
    r = gsl_rng_alloc (type);
    gsl_rng_env_setup();
    
    unsigned long n_window;
    unsigned long div,div2;
    unsigned long snp_pos,start,end;
    int variant;
    
    unsigned long *vec_rd, *vec2_rd;
    double *vec_s, *vec_p, *vec_h, *vec0_s, *vec0_d, *vec0_h;
    double *vec2_s, *vec2_p, *vec2_h, *vec02_s, *vec02_d, *vec02_h;
    double *vec_n, *vec_x, *vec_e, *vec0_n, *vec0_x, *vec0_e;
    double *vec2_n, *vec2_x, *vec2_e, *vec02_n, *vec02_x, *vec02_e;
    double *covmat;
    double *covmat2;
    
    struct combinatorial comb1, comb2;
    struct combinatorial_fst combfst;
    struct tests test1, test2;
    struct tests tests1, testn1;
    struct tests tests2, testn2;
    struct fst_calc fst;
    
    FILE * list_snps=0;
    unsigned long pos_snp=0;
    
    FILE * gff=0;
    char if_gff;
    unsigned long cds_start, cds_end;
    unsigned int phase_cds, frame;
    char strand, *gff_file_name;
    double psyn, pnon, dsyn, dnon;
    double psyn2, pnon2, dsyn2, dnon2;
    
    double thetaw_s,thetat_s,thetah_s,div_s;
    double thetaw_n,thetat_n,thetah_n,div_n;
    
    char *cchrom,*cchrom2,*c2chrom2;
    char *schrom,*schrom2;
    char *gchrom,*gchrom2;
    char *ochrom;
    char *cchrom_next;
    char *c2chrom_next;
    unsigned long posw;
    
    char *scaffold;
    scaffold=(char *)malloc(256*sizeof(char));
    unsigned long *scaffold_out;
    scaffold_out=(unsigned long *)malloc(1*sizeof(unsigned long));
    int nscaf=0,nscafc=0;
    int gamma1,gamma2;
    
    char *cline;
    size_t nline;
    char *oline;
    size_t noline;
    cline=(char *)malloc(1);
    nline=1;
    oline=(char *)malloc(1);
    noline=1;
    
    int r1_ok=0;
    int r2_ok=0;
    struct keep_codon_data codon_p1[3];
    struct keep_codon_data codon_p2[3];
    int codon_full=0;
    int *codon_frame;
    double *cweights;
    double *dweights;
    double *cmutnsyn,*cmutsyn;
    double *kmutnsyn,*kmutsyn;
    codon_frame=(int *) calloc(3,sizeof(int));
    cweights=(double *) calloc(3,sizeof(double));
    dweights=(double *) calloc(3,sizeof(double));
    cmutsyn=(double *) calloc(3,sizeof(double));
    kmutsyn=(double *) calloc(3,sizeof(double));
    cmutnsyn=(double *) calloc(3,sizeof(double));
    kmutnsyn=(double *) calloc(3,sizeof(double));

    n_alt_allele1=0;n_alt_allele2=0;
    strand=0;phase_cds=0;frame=0;
    n_ref1=0;n_ref2=0;rd2=0;
    alt_base1=0; ref_base1=0;
    ref_base2=0; alt_base2=0;
    out_base=0; div=0;div2=0;
    psyn=pnon=dsyn=dnon=0.;
    psyn2=pnon2=dsyn2=dnon2=0.;
    
    /* Main part of the program */
    
    /* Read arguments */
    window_size=0;
    n01=0;
    n02=0;
    rd=rd1=rd2=0;
    min_cov=4;
    max_cov=100;
    min_qual=10;
    min_mqual=10;
    m_bar=1; //default
    int fst_path = 0;

    char from_stdin;
    char outgroup_available;
    char compute_fst=0;
    char ext_snps=0;
    char map_qual_pileup=0;
    char *pileup_file_name=0;
    char *pileup2_file_name=0;
    char *outgroup_file_name=0;
    char *snp_file_name=0;
    int arg_i;
    from_stdin=2;
    outgroup_available=0;
    compute_fst=0;
    ext_snps=0;
    map_qual_pileup=0;
    if_gff=0;
    
    /* BEG <----ADDED 29112022*/
    char *scaffold_filename;
    scaffold_filename=(char *) calloc(1000,sizeof(char));
    char *bed_filename;
    bed_filename=(char *) calloc(1000,sizeof(char));
    char *outfile;
    outfile=(char *) calloc(1000,sizeof(char));
    cchrom=(char *) calloc(100,sizeof(char));
    cchrom2=(char *) calloc(100,sizeof(char));
    c2chrom2=(char *) calloc(100,sizeof(char));
    schrom=(char *) calloc(100,sizeof(char));
    schrom2=(char *) calloc(100,sizeof(char));
    gchrom=(char *) calloc(100,sizeof(char));
    gchrom2=(char *) calloc(100,sizeof(char));
    ochrom=(char *) calloc(100,sizeof(char));
    cchrom_next=(char *) calloc(100,sizeof(char));
    c2chrom_next=(char *) calloc(100,sizeof(char));
    /* END <----ADDED 29112022*/
    
    for(arg_i=1;arg_i<argc-1;arg_i++)
    {
        if (strcmp(argv[arg_i], "-n") == 0) {arg_i++; sscanf(argv[arg_i], "%u", &n01); }
        else if (strcmp(argv[arg_i], "-h") == 0) {arg_i++; usage(); exit(0);}
        else if (strcmp(argv[arg_i], "-l") == 0) {arg_i++; sscanf(argv[arg_i], "%lu", &window_size); }
        else if (strcmp(argv[arg_i], "-bedfile") == 0) {arg_i++; sscanf(argv[arg_i], "%s",bed_filename); }
        else if (strcmp(argv[arg_i], "-cov") == 0) {arg_i++; sscanf(argv[arg_i], "%lf", &input_coverage1); }
        else if (strcmp(argv[arg_i], "-n2") == 0) {arg_i++; sscanf(argv[arg_i], "%u", &n02); }
        else if (strcmp(argv[arg_i], "-mincov") == 0) {arg_i++; sscanf(argv[arg_i], "%lu", &min_cov); }
        else if (strcmp(argv[arg_i], "-maxcov") == 0) {arg_i++; sscanf(argv[arg_i], "%lu", &max_cov); }
        else if (strcmp(argv[arg_i], "-minqual") == 0) {arg_i++; sscanf(argv[arg_i], "%lu", &min_qual); }
        //     else if (strcmp(argv[arg_i], "-mapqual") == 0) {map_qual_pileup=1; }
        //     else if (strcmp(argv[arg_i], "-minmapqual") == 0) {arg_i++; sscanf(argv[arg_i], "%lu", &min_mqual); }
        else if (strcmp(argv[arg_i], "-nolowfreq") == 0) {arg_i++; sscanf(argv[arg_i], "%u", &m_bar); }
        else if (strcmp(argv[arg_i], "-outgroup") == 0) {outgroup_available=1; arg_i++; outgroup_file_name=argv[arg_i];}
        else if (strcmp(argv[arg_i], "-fstpop2") == 0) {compute_fst=1; arg_i++; pileup2_file_name=argv[arg_i];}
        //else if (strcmp(argv[arg_i], "-pileup") == 0) {}
        else if (strcmp(argv[arg_i], "-snpfile") == 0) {ext_snps=1; arg_i++; snp_file_name=argv[arg_i];}
        else if (strcmp(argv[arg_i], "-annot") == 0) {if_gff=1; arg_i++; gff_file_name=argv[arg_i];}
        //ADDED
        else if (strcmp(argv[arg_i], "-outfile") == 0) {arg_i++; sscanf(argv[arg_i], "%s",outfile); }
        else if (strcmp(argv[arg_i], "-scaffolds") == 0) {arg_i++; sscanf(argv[arg_i], "%s",scaffold_filename); }
        else if (strcmp(argv[arg_i], "-simplefst") == 0) {fst_path=0; arg_i++; sscanf(argv[arg_i], "%d", &fst_path);}
        /*     else if (strcmp(argv[arg_i], "-snpcolumns") == 0)
         {
         arg_i++; sscanf(argv[arg_i], "%u", &col_chrom);
         arg_i++; sscanf(argv[arg_i], "%u", &col_pos);
         };
         */
    }
    if (arg_i==argc-1) {
        if (strcmp(argv[arg_i], "-") == 0)
            from_stdin=1;
        else {
            from_stdin=0;
            pileup_file_name=argv[arg_i];}/*stdin option is excluded using gz files*/
    };
    //from_stdin=0; pileup_file_name=argv[arg_i-1];
    
    if ((n01==0)||(window_size==0 && bed_filename[0]==0)||(from_stdin))/*stdin option is excluded using gz files*/
    {
        fprintf(stderr,"Missing values in command line!\n");
        usage();
        return(-1);
    } else{
        printf(VERSION_NPSTAT);
        printf("#command: npstat2 ");
        for(arg_i=1;arg_i<argc-1;arg_i++) {printf(" %s ",argv[arg_i]);}
        if (arg_i==argc-1) {
            if (strcmp(argv[arg_i], "-") == 0) {printf(" stdout ");} else {printf(" %s ",argv[arg_i]);}
        };
        printf("\n");
    };
    //fprintf(stderr,"Missing values in command line!\n Command:\n    NPStat [options] file.pileup\n or to read from standard input:\n    NPStat [options] -\n Options:\n    -n samplesize : haploid sample size\n    -l windowlength : window length\n    -mapqual : pileup includes mapping quality\n    -mincov minimum_coverage : filter on minimum coverage (default 4)\n    -maxcov maximum_coverage : filter on maximum coverage (default 100)\n    -minqual minimum_base_quality : filter on base quality (default 10)\n    -minmapqual minimum_mapping_quality : filter on mapping quality (default 10)\n    -nolowfreq m : filter on minimum allele count mac>m\n    -outgroup file.fa : outgroup file in FASTA\n    -fstpop2 file2.pileup : computes Fst with a second population \n    contained in file2.pileup\n    -n2 : sample size of the second population\n");
    //-snpfile file.snp : consider SNPs only if present in file.snp\n
    //-snpcolumns chrom pos : columns for chromosome and \n
    //position in file.snp (default values 1, 2)\n
    
    //outgroup_available=1;
    /*
     outgroup_available=0;
     if (argc<3) {
     fprintf(stderr, "Usage: bam_read_routine file.bam reference.fasta outgroup.fasta") ;
     return 1;
     }
     else if (argc>3) {
     outgroup_available=1;
     };
     */
    
    /* Initialize and open file */
    //temp_file_name1=malloc(strlen(pileup_file_name)+strlen(pileup2_file_name)+10);
    //temp_file_name2=malloc(strlen(pileup_file_name)+strlen(pileup2_file_name)+10);
    temp_file_name1=malloc(strlen(pileup_file_name)+strlen(pileup_file_name)+10);
    temp_file_name2=malloc(strlen(pileup_file_name)+strlen(pileup_file_name)+10);
    
    //if (from_stdin==0) {
    bam_file1=/*fopen*/bgzf_open(pileup_file_name,"r");
    if (bam_file1 == NULL){
        fprintf(stderr,"Error: the pileup file cannot be opened!\n");
        usage();
        return(-1);
    };
    //} else {
    //    bam_file1=stdin;
    //    pileup_file_name="NPStat_file";
    //};
    if(compute_fst) {
        bam_file2=/*fopen*/bgzf_open(pileup2_file_name,"r");
        if (bam_file2 == NULL){
            fprintf(stderr,"Error: the pileup file2 cannot be opened!\n");
            usage();
            return(-1);
        };
    }
    
    if(ext_snps==1){
        list_snps=fopen(snp_file_name,"r");
        if (list_snps == NULL){
            fprintf(stderr,"Error: the SNP list file cannot be opened!\n");
            usage();
            return(-1);
        };
        if(if_gff==1) {
            fprintf(stderr,"Sorry: the SNP list file is incompatible with flag gff\n");
            usage();
            return(-1);
       }
    };
    
    if(if_gff==1){
        gff=fopen(gff_file_name,"r");
        if (gff == NULL){
            fprintf(stderr,"Error: the GFF3 file cannot be opened!\n");
            usage();
            return(-1);
        };
    };
    if(n01%2 != 0 || n02%2 != 0) {
        fprintf(stderr,"Error: Sorry to say that npstat2 only accepts even haploid sample sizes\n");
        usage();
        return(-1);
    }
    
    if(outfile[0]!='\0') {
        strcpy(temp_file_name1,outfile);
    }
    else {
        strcpy(temp_file_name1,pileup_file_name);
        strcpy(temp_file_name2,".stats.txt");
        strcat(temp_file_name1,temp_file_name2);
    }
    output_stat1=fopen(temp_file_name1,"w");
    if (output_stat1 == NULL){
        fprintf(stderr,"Error: the output file cannot be written!\n");
        usage();
        return(-1);
    };
    
    if(compute_fst) {
        if(outfile[0]!='\0') {
            strcpy(temp_file_name1,outfile);
            strcat(temp_file_name1,"_p2.txt");
        }
        else {
            strcpy(temp_file_name1,pileup2_file_name);
            strcpy(temp_file_name2,",stats.txt");
            strcat(temp_file_name1,temp_file_name2);
        }
        output_stat2=fopen(temp_file_name1,"w");
        if (output_stat2 == NULL){
            fprintf(stderr,"Error: the output2 file cannot be written!\n");
            usage();
            return(-1);
        };
        
        if(outfile[0]!='\0') {
            strcpy(temp_file_name1,outfile);
            strcat(temp_file_name1,"_fst.txt");
        }
        else {
            strcpy(temp_file_name1,"fst_");
            strcpy(temp_file_name2,pileup_file_name);
            strcat(temp_file_name1,temp_file_name2);
            strcpy(temp_file_name2,"_");
            strcat(temp_file_name1,temp_file_name2);
            strcpy(temp_file_name2,pileup2_file_name);
            strcat(temp_file_name1,temp_file_name2);
        }
        output_fst=fopen(temp_file_name1,"w");
        if (output_fst == NULL){
            fprintf(stderr,"Error: the output_fst file cannot be written!\n");
            usage();
            return(-1);
        };
    }
    scaffold_file=fopen(scaffold_filename,"r");
    if (scaffold_file == NULL){
        fprintf(stderr,"Error: the scaffold file cannot be opened!\n");
        usage();
        return(-1);
    };
    if(bed_filename[0]!=0) {
        bed_file=fopen(bed_filename,"r");
        if (bed_file == NULL){
            fprintf(stderr,"Error: the bedfile cannot be opened!\n");
            usage();
            return(-1);
        };
    }
    /**/
    
    //SNPS
    /*
     strcpy(temp_file_name1,"snps_");
     strcpy(temp_file_name2,pileup_file_name);
     strcat(temp_file_name1,temp_file_name2);
     strcpy(temp_file_name2,"_");
     strcat(temp_file_name1,temp_file_name2);
     strcpy(temp_file_name2,pileup2_file_name);
     strcat(temp_file_name1,temp_file_name2);
     #if PRINTSNPS == 1
     FILE *output_snps;
     output_snps=fopen(temp_file_name1,"w");
     #endif
     */
    
    /*
     bam_file=fopen(argv[1],"r");
     fasta_ref=fopen(argv[2],"r");
     if (outgroup_available==1) {
     fasta_out=fopen(argv[3],"r");
     };
     */
    /*for (;iscntrl(fgetc(fasta_ref))==0;) {};*/
    fasta_length=0;
    if (outgroup_available==1)
    {
        if((fasta_out=fopen(outgroup_file_name,"r"))==NULL) {
            fprintf(stderr,"Error: the outgroup file %s cannot be opened!\n",outgroup_file_name);
            usage();
            return(-1);
        }
        /*Keep all the initial row positions having the character '>'*/
        unsigned long pos_out=0;
        char test;
        nscafc=0;
        test=fgetc(fasta_out);
        for (;test!=EOF;test=fgetc(fasta_out)){
            if(test=='>') {
                scaffold_out[nscafc]=pos_out;
                nscafc++;
                scaffold_out=(unsigned long *)realloc(scaffold_out,(nscafc+1)*sizeof(unsigned long));
            }
            pos_out++;
        }
        rewind(fasta_out);
        for (;iscntrl(fgetc(fasta_out))==0;){};//fasta definition row
        // fasta_pos=ftell(fasta_out);
        for (count_i=0;iscntrl(fgetc(fasta_out))==0;count_i++) {}; //length of fasta lines (all must be equal)
        fasta_length=count_i;
        rewind(fasta_out);
        //for (;iscntrl(fgetc(fasta_out))==0;){};
    };
    //for (;ct!="\n";ct=fgetc(fasta_ref)) {};
    //for (;ct!="\n";ct=fgetc(fasta_out)) {};
    
    /* Load combinatorics */
    printf("Initializing combinatorics...\n");
    
    vec_rd=(unsigned long *)malloc(max_cov*sizeof(unsigned long));
    
    vec_s=(double *)malloc(max_cov*(n01-1)*sizeof(double));
    vec_p=(double *)malloc(max_cov*(n01-1)*sizeof(double));
    vec_h=(double *)malloc(max_cov*(n01-1)*sizeof(double));
    vec_n=(double *)malloc(max_cov*(n01-1)*sizeof(double));
    vec_x=(double *)malloc(max_cov*(n01-1)*sizeof(double));
    vec_e=(double *)malloc(max_cov*(n01-1)*sizeof(double));
    vec0_s=(double *)malloc(max_cov*sizeof(double));
    vec0_d=(double *)malloc(max_cov*sizeof(double));
    vec0_h=(double *)malloc(max_cov*sizeof(double));
    vec0_n=(double *)malloc(max_cov*sizeof(double));
    vec0_x=(double *)malloc(max_cov*sizeof(double));
    vec0_e=(double *)malloc(max_cov*sizeof(double));
    
    comb1.d_t=(double *)malloc(max_cov*sizeof(double));
    comb1.d_p=(double *)malloc(max_cov*sizeof(double));
    comb1.d_hl=(double *)malloc(max_cov*sizeof(double));
    comb1.d_hq=(double *)malloc(max_cov*sizeof(double));
    comb1.d_nu=(double *)malloc(max_cov*sizeof(double));
    comb1.d_xi=(double *)malloc(max_cov*sizeof(double));
    
    /**/
    if(compute_fst) {
        vec2_rd=(unsigned long *)malloc(max_cov*sizeof(unsigned long));
        
        vec2_s=(double *)malloc(max_cov*(n02-1)*sizeof(double));
        vec2_p=(double *)malloc(max_cov*(n02-1)*sizeof(double));
        vec2_h=(double *)malloc(max_cov*(n02-1)*sizeof(double));
        vec2_n=(double *)malloc(max_cov*(n02-1)*sizeof(double));
        vec2_x=(double *)malloc(max_cov*(n02-1)*sizeof(double));
        vec2_e=(double *)malloc(max_cov*(n02-1)*sizeof(double));
        vec02_s=(double *)malloc(max_cov*sizeof(double));
        vec02_d=(double *)malloc(max_cov*sizeof(double));
        vec02_h=(double *)malloc(max_cov*sizeof(double));
        vec02_n=(double *)malloc(max_cov*sizeof(double));
        vec02_x=(double *)malloc(max_cov*sizeof(double));
        vec02_e=(double *)malloc(max_cov*sizeof(double));
        
        comb2.d_t=(double *)malloc(max_cov*sizeof(double));
        comb2.d_p=(double *)malloc(max_cov*sizeof(double));
        comb2.d_hl=(double *)malloc(max_cov*sizeof(double));
        comb2.d_hq=(double *)malloc(max_cov*sizeof(double));
        comb2.d_nu=(double *)malloc(max_cov*sizeof(double));
        comb2.d_xi=(double *)malloc(max_cov*sizeof(double));
        
        combfst.c_s=(double **)malloc(max_cov*sizeof(double *));
        for(count_i=0;count_i<max_cov;count_i++) {
            combfst.c_s[count_i]=(double *)malloc(max_cov*sizeof(double));
        };
        combfst.d_pt=(double *)malloc(2*max_cov*sizeof(double));
    }
    /**/
    
    for(count_i=0;count_i<max_cov;count_i++)
    {
        comb1.d_t[count_i]=0;
        comb1.d_p[count_i]=0;
        comb1.d_hl[count_i]=0;
        comb1.d_hq[count_i]=0;
        comb1.d_nu[count_i]=0;
        comb1.d_xi[count_i]=0;
        
        /**/
        if(compute_fst) {
            comb2.d_t[count_i]=0;
            comb2.d_p[count_i]=0;
            comb2.d_hl[count_i]=0;
            comb2.d_hq[count_i]=0;
            comb2.d_nu[count_i]=0;
            comb2.d_xi[count_i]=0;
            
            for(count_j=0;count_j<max_cov;count_j++) {
                combfst.c_s[count_i][count_j]=0;
            };
            combfst.d_pt[count_i]=0;
        }
        /**/
    };
    for(rd=m_bar+2;rd<=max_cov;rd++)
    {
        int k,j;
        comb1.d_p[rd-1]+=((double)n01-1)/(double)n01;/*SI-13 den*/
        comb1.d_hl[rd-1]+=(double)rd*(n01-1)/(double)(n01*(rd-1));/*SI-20 den*/
        comb1.d_hq[rd-1]+=(double)rd*((double)((n01-1)*(rd+1))/(double)(2*rd))/(double)(n01*(rd-1)); /*SI-17 den*/
        gamma1 = fmax(((float)rd/(float)n01),1.0);
        if(compute_fst) {
            comb2.d_p[rd-1]+=((double)n02-1)/(double)n02;/*SI-13 den*/
            comb2.d_hl[rd-1]+=(double)rd*(n02-1)/(double)(n02*(rd-1));/*SI-20 den*/
            comb2.d_hq[rd-1]+=(double)rd*((double)((n02-1)*(rd+1))/(double)(2*rd))/(double)(n02*(rd-1)); /*SI-17 den*/
            gamma2 = fmax(((float)rd/(float)n02),1.0);
        }
        //p1
        for(k=1;k<n01;k++)
        {
            comb1.d_t[rd-1]+=(1-gsl_pow_int((double)k/(double)n01,rd)-gsl_pow_int(1-(double)k/(double)n01,rd))/(double)k;/*SI-11 den ...*/
            comb1.d_hl[rd-1]+=-(double)rd/(double)(n01*(rd-1))*(gsl_pow_int((double)k/(double)n01,rd-1)); /*SI-20*/
            comb1.d_hq[rd-1]+=-(double)rd/(double)(n01*(rd-1))*(gsl_pow_int((double)k/(double)n01,rd-1)); /*SI-17*/
            //(((double)rd-1)/(double)n01)*((1+(double)(rd+1)/gsl_pow_int((double)(rd-1),2))*(double)k/(double)n01-1)*gsl_pow_int((double)k/(double)n01,rd-2);
            for(lt=1;lt<=m_bar;lt++)
            {
                comb1.d_t[rd-1]+=-gsl_sf_choose(rd,lt)*gsl_pow_int((double)k/(double)n01,lt-1)*gsl_pow_int(1-(double)k/(double)n01,rd-lt-1)/(double)n01;/*SI-12*/
                comb1.d_p[rd-1]+=-2*gsl_sf_choose(rd-2,lt-1)*gsl_pow_int((double)k/(double)n01,lt-1)*gsl_pow_int(1-(double)k/(double)n01,rd-lt-1)/(double)n01;  /*SI-14*/
                comb1.d_hq[rd-1]+=-(double)rd/(double)(n01*(rd-1))*(double)lt/(double)(rd)* gsl_sf_choose(rd-1,lt-1)*gsl_pow_int((double)k/(double)n01,lt-1)*gsl_pow_int(1-(double)k/(double)n01,rd-lt) -(double)rd/(double)(n01*(rd-1))*(double)(rd-lt)/(double)(rd)* gsl_sf_choose(rd-1,rd-lt-1)*gsl_pow_int((double)k/(double)n01,rd-lt-1)*gsl_pow_int(1-(double)k/(double)n01,lt);/*SI-21*/
                comb1.d_hl[rd-1]+=-(double)rd/(double)(n01*(rd-1))* gsl_sf_choose(rd-1,lt-1)*gsl_pow_int((double)k/(double)n01,lt-1)*gsl_pow_int(1-(double)k/(double)n01,rd-lt) -(double)rd/(double)(n01*(rd-1))* gsl_sf_choose(rd-1,rd-lt-1)*gsl_pow_int((double)k/(double)n01,rd-lt-1)*gsl_pow_int(1-(double)k/(double)n01,lt);/*SI-21*/
            };
            for(j=m_bar+1;(float)j<=gamma1;j++) {
                comb1.d_nu[rd-1]+=1.0/(double)n01*gsl_sf_choose(rd,j)*gsl_pow_int((double)k/(double)n01,j-1)*gsl_pow_int(1-(double)k/(double)n01,rd-j);/*SI-15*/
                comb1.d_xi[rd-1]+=1.0/(double)n01*gsl_sf_choose(rd,j)*gsl_pow_int((double)k/(double)n01,j-1)*gsl_pow_int(1-(double)k/(double)n01,rd-j);/*SI-16*/
            }
        };
        /*
        for(k=1;k<n02;k++) {
            comb2.d_t[rd-1]+=(1-gsl_pow_int((double)k/(double)n02,rd)-gsl_pow_int(1-(double)k/(double)n02,rd))/(double)k;
            for(lt=1;lt<=m_bar;lt++)
            {
                comb2.d_t[rd-1]+=+gsl_sf_choose(rd,lt)*gsl_pow_int((double)k/(double)n02,lt-1)*gsl_pow_int(1-(double)k/(double)n02,rd-lt-1)/(double)n02;
                comb2.d_p[rd-1]+=+2*gsl_sf_choose(rd-2,lt-1)*gsl_pow_int((double)k/(double)n02,lt-1)*gsl_pow_int(1-(double)k/(double)n02,rd-lt-1);
            };
            comb2.d_hl[rd-1]+=((double)rd/(double)n02)*((1-2/((double)rd-1))*(double)k/(double)n02-1)*gsl_pow_int((double)k/(double)n02,rd-2);
            comb2.d_hq[rd-1]+=(((double)rd-1)/(double)n02)*((1+(double)(rd+1)/gsl_pow_int((double)(rd-1),2))*(double)k/(double)n02-1)*gsl_pow_int((double)k/(double)n02,rd-2);
            comb2.d_t[rd-1]+=(1-gsl_pow_int((double)k/(double)n02,rd-1)*((double)k/(double)n02+(double)rd)-gsl_pow_int(1-(double)k/(double)n02,rd))/(double)k;
            comb2.d_p[rd-1]+=(-2*gsl_pow_int((double)k/(double)n02,rd-2))/(double)n02;
            comb2.d_hl[rd-1]+=((double)rd/(double)n02)*((1-2/((double)rd-1))*(double)k/(double)n02-1)
            *gsl_pow_int((double)k/(double)n02,rd-2);
            comb2.d_hq[rd-1]+=(((double)rd-1)/(double)n02)*((1+(double)(rd+1)/gsl_pow_int((double)(rd-1),2))*(double)k/(double)n02-1)*gsl_pow_int((double)k/(double)n02,rd-2);
        };
        */
        //p2
        if(compute_fst) {
            for(k=1;k<n02;k++)
            {
                comb2.d_t[rd-1]+=(1-gsl_pow_int((double)k/(double)n02,rd)-gsl_pow_int(1-(double)k/(double)n02,rd))/(double)k;
                comb2.d_hl[rd-1]+=-(double)rd/(double)(n02*(rd-1))*(gsl_pow_int((double)k/(double)n02,rd-1)); /*SI-20*/
                comb2.d_hq[rd-1]+=-(double)rd/(double)(n02*(rd-1))*(gsl_pow_int((double)k/(double)n02,rd-1)); /*SI-17*/
                for(lt=1;lt<=m_bar;lt++)
                {
                    comb2.d_t[rd-1]+=-gsl_sf_choose(rd,lt)*gsl_pow_int((double)k/(double)n02,lt-1)*gsl_pow_int(1-(double)k/(double)n02,rd-lt-1)/(double)n02;
                    comb2.d_p[rd-1]+=-2*gsl_sf_choose(rd-2,lt-1)*gsl_pow_int((double)k/(double)n02,lt-1)*gsl_pow_int(1-(double)k/(double)n02,rd-lt-1)/(double)n02;
                    comb2.d_hq[rd-1]+=-(double)rd/(double)(n02*(rd-1))*(double)lt/(double)(rd)* gsl_sf_choose(rd-1,lt-1)*gsl_pow_int((double)k/(double)n02,lt-1)*gsl_pow_int(1-(double)k/(double)n02,rd-lt) -(double)rd/(double)(n02*(rd-1))*(double)(rd-lt)/(double)(rd)* gsl_sf_choose(rd-1,rd-lt-1)*gsl_pow_int((double)k/(double)n02,rd-lt-1)*gsl_pow_int(1-(double)k/(double)n02,lt);/*SI-21*/
                    comb2.d_hl[rd-1]+=-(double)rd/(double)(n02*(rd-1))* gsl_sf_choose(rd-1,lt-1)*gsl_pow_int((double)k/(double)n02,lt-1)*gsl_pow_int(1-(double)k/(double)n02,rd-lt) -(double)rd/(double)(n02*(rd-1))* gsl_sf_choose(rd-1,rd-lt-1)*gsl_pow_int((double)k/(double)n02,rd-lt-1)*gsl_pow_int(1-(double)k/(double)n02,lt);/*SI-21*/
                };
                for(j=m_bar+1;(float)j<=gamma2;j++) {
                    comb2.d_nu[rd-1]+=1.0/(double)n02*gsl_sf_choose(rd,j)*gsl_pow_int((double)k/(double)n02,j-1)*gsl_pow_int(1-(double)k/(double)n02,rd-j);
                    comb2.d_xi[rd-1]+=1.0/(double)n02*gsl_sf_choose(rd,j)*gsl_pow_int((double)k/(double)n02,j-1)*gsl_pow_int(1-(double)k/(double)n02,rd-j);
                }
            };
         }
    };
    
    //Define Fst (...)
    //if( (float)abs(n01-n02) / ((float)min(n01,n02)) > 1.0 )
    //    fst_path = 0; //if the difference between nsam is too large
    //else
    //    fst_path = 1; //only if assume similar nsam and similar nreads in both pops!
    //
    //fst_path = 0; //TESTING
    
    m_bar_fst=m_bar;
    //if(compute_fst) {
        //if( fst_path == 0 ) { //if the difference between nsam is too large
            //Fst from Pia
            /*
            for(rd1=2;rd1<=max_cov;rd1++) {
                for(rd2=2;rd2<=max_cov;rd2++) {
                    int k,l;
                    for(k=1;k<=n01+n02-1;k++) {
                        for(l=0;l<=k;l++) {
                            for(m1t=0;m1t<=m_bar;m1t++){
                                for(m2t=0;m2t<=m_bar;m2t++){
                                    combfst.c_s[rd1-1][rd2-1]+= (gsl_ran_hypergeometric_pdf(l,n01,n02,k)*((double)(m1t*(rd2-m2t)+m2t*(rd1-m1t))/(double)(rd1*rd2))*gsl_sf_choose(rd1,m1t)*gsl_sf_choose(rd2,m2t) * (gsl_pow_int((double)l/(double)(n01+n02),m1t)*gsl_pow_int(1-((double)l/(double)(n01+n02)),rd1-m1t)+gsl_pow_int((double)l/(double)(n01+n02),rd1-m1t)*gsl_pow_int(1-((double)l/(double)(n01+n02)),m1t)) * (gsl_pow_int((double)(k-l)/(double)(n01+n02),m2t)*gsl_pow_int(1-((double)(k-l)/(double)(n01+n02)),rd2-m2t)+gsl_pow_int((double)(k-l)/(double)(n01+n02),rd2-m2t)*gsl_pow_int(1-((double)(k-l)/(double)(n01+n02)),m2t)))/(double)k;;
                                }
                            }
                        }
                    }
                    //printf("combfst.c_s[%d-1][%d-1]=%f\n",rd1,rd2,combfst.c_s[rd1-1][rd2-1]);
                }
            }
            */
            /*
            //Fst from Pia (same than above)
            for(rd1=m_bar+1;rd1<=max_cov;rd1++) { //all possible read depth combinations
                for(rd2=m_bar+1;rd2<=max_cov;rd2++) {
                    int i,k;
                    double x1,y1,x2,y2;
                    for(k=1;k<=n01+n02-1;k++) { //maximum sample sizes for the sum of the two pops
                        for(i=0;i<=k;i++) { //x1 and y1 are the freqs k-i and i-freqs given total, x2 are 1 - these freqs
                            x1=(double)(k-i)/(double)(n01+n02);
                            x2=(double)(i)/(double)(n01+n02);
                            y1=1-x1;
                            y2=1-x2;
                            for(m1t=0;m1t<=m_bar;m1t++){ //frequencies from 0 to m_bar (more reads than in one pop) are not considered
                                for(m2t=0;m2t<=m_bar;m2t++){
                                    combfst.c_s[rd1-1][rd2-1]+=
                                    gsl_ran_hypergeometric_pdf(k-i,n01,n02,k)*((double)(m1t*(rd2-m2t)+m2t*(rd1-m1t))/(double)(rd1*rd2))*
                                    gsl_sf_choose(rd1,m1t)*gsl_sf_choose(rd2,m2t)*
                                    (gsl_pow_int(x2,m1t)*gsl_pow_int(y2,rd1-m1t)+gsl_pow_int(x2,rd1-m1t)*gsl_pow_int(y2,m1t))*
                                    (gsl_pow_int(x1,m2t)*gsl_pow_int(y1,rd2-m2t)+gsl_pow_int(x1,rd2-m2t)*gsl_pow_int(y1,m2t))
                                    /(double)k; //SI-28
                                    //if(combfst.c_s[rd1-1][rd2-1]<0.) {
                                    //    printf("combfst.c_s[%d-1][%d-1]=%f\n",rd1,rd2,combfst.c_s[rd1-1][rd2-1]);
                                    //}
                                }
                            }
                        }
                    }
                    //printf("combfst.c_s[%d-1][%d-1]=%f\n",rd1,rd2,combfst.c_s[rd1-1][rd2-1]);
                }
            }
            */
        //}
        //else {
            //Fst for piT. Only if assume similar nsam and similar nreads in both pops!
        /*
             for(rd=2*m_bar_fst+2;rd<=2*max_cov;rd++)
            {
                int k;
                combfst.d_pt[rd-1]+=(double)(n01+n02-1)/(double)(n01+n02);//SI-13 den
                for(k=1;k<n01+n02;k++) {
                    //p1+p2
                    for(lt=1;lt<=m_bar_fst;lt++) {
                        combfst.d_pt[rd-1]+= -2*gsl_sf_choose(rd-2,lt-1) * gsl_pow_int((double)k/(double)(n01+n02),lt-1) * gsl_pow_int(1-(double)k/(double)(n01+n02),rd-lt-1)/(double)(n01+n02);  //SI-14
                    }
                }
            }
        */
        //}
    //}/**/
    //END
    
    for(rd=1;rd<=max_cov*(n01-1); rd++){
        vec_s[rd-1]=0;
        vec_p[rd-1]=0;
        vec_h[rd-1]=0;
        vec_n[rd-1]=0;
        vec_x[rd-1]=0;
        vec_e[rd-1]=0;
    };
    for(rd=2*m_bar+1;rd<=max_cov; rd++){
        int k,covt;
        vec0_s[rd-1]=0;
        vec0_d[rd-1]=0;
        vec0_h[rd-1]=0;
        vec0_n[rd-1]=0;
        vec0_x[rd-1]=0;
        vec0_e[rd-1]=0;
        for(k=1;k<n01;k++){
            vec_s[(k-1)+(n01-1)*(rd-1)]=(1-gsl_pow_int((double)k/(double)n01,rd)-gsl_pow_int(1-(double)k/(double)n01,rd));
            vec_p[(k-1)+(n01-1)*(rd-1)]=(double)(2*k*(n01-k))/(double)(n01*n01);
            vec_h[(k-1)+(n01-1)*(rd-1)]=(double)(k*k)/(double)(n01*n01)+(double)(k)/(double)(n01*(rd-1))-(double)rd/(double)(rd-1)*gsl_pow_int((double)(k)/(double)(n01),rd);
            for(covt=1;covt<=m_bar;covt++){
                vec_s[(k-1)+(n01-1)*(rd-1)]+=-( gsl_ran_binomial_pdf(covt,	(double)k/(double)n01,rd)+gsl_ran_binomial_pdf(rd-covt,(double)k/(double)n01,rd));
                vec_p[(k-1)+(n01-1)*(rd-1)]+=-( 2*(double)(covt*(rd-covt))/(double)(rd*(rd-1))*(gsl_ran_binomial_pdf(covt,	(double)k/(double)n01,rd)+gsl_ran_binomial_pdf(rd-covt,(double)k/(double)n01,rd)));
                vec_h[(k-1)+(n01-1)*(rd-1)]+=-( (double)(covt*covt)/(double)(rd*(rd-1))*gsl_ran_binomial_pdf(covt,	(double)k/(double)n01,rd)+(double)(rd-covt)*(double)(rd-covt)/(double)(rd*(rd-1))*gsl_ran_binomial_pdf(rd-covt,(double)k/(double)n01,rd));
            };
            for(covt=m_bar+1;covt<=rd-m_bar-1;covt++){
                vec0_s[rd-1]+=gsl_ran_binomial_pdf(covt,(double)k/(double)n01,rd)/k;
                vec0_d[rd-1]+=gsl_pow_int(2*(double)(covt*(rd-covt))/(double)(rd*(rd-1))/comb1.d_p[rd-1]-1/comb1.d_t[rd-1],2)*gsl_ran_binomial_pdf(covt,(double)k/(double)n01,rd)/k;
                vec0_h[rd-1]+=gsl_pow_int((double)(covt*covt)/(double)(rd*(rd-1))/comb1.d_hq[rd-1]-2*(double)(covt*(rd-covt))/(double)(rd*(rd-1))/comb1.d_p[rd-1],2)*gsl_ran_binomial_pdf(covt,(double)k/(double)n01,rd)/k;
            };
        };
    };
    if(compute_fst) {
        for(rd=1;rd<=max_cov*(n02-1); rd++){
            vec2_s[rd-1]=0;
            vec2_p[rd-1]=0;
            vec2_h[rd-1]=0;
            vec2_n[rd-1]=0;
            vec2_x[rd-1]=0;
            vec2_e[rd-1]=0;
        };
        for(rd=2*m_bar+1;rd<=max_cov; rd++){
            int k,covt;
            vec02_s[rd-1]=0;
            vec02_d[rd-1]=0;
            vec02_h[rd-1]=0;
            vec02_x[rd-1]=0;
            vec02_n[rd-1]=0;
            vec02_e[rd-1]=0;
            for(k=1;k<n02;k++){
                vec2_s[(k-1)+(n02-1)*(rd-1)]=(1-gsl_pow_int((double)k/(double)n02,rd)-gsl_pow_int(1-(double)k/(double)n02,rd));
                vec2_p[(k-1)+(n02-1)*(rd-1)]=(double)(2*k*(n02-k))/(double)(n02*n02);
                vec2_h[(k-1)+(n02-1)*(rd-1)]=(double)(k*k)/(double)(n02*n02)+(double)(k)/(double)(n02*(rd-1))-(double)rd/(double)(rd-1)*gsl_pow_int((double)(k)/(double)(n02),rd);
                for(covt=1;covt<=m_bar;covt++){
                    vec2_s[(k-1)+(n02-1)*(rd-1)]+=-( gsl_ran_binomial_pdf(covt,    (double)k/(double)n02,rd)+gsl_ran_binomial_pdf(rd-covt,(double)k/(double)n02,rd));
                    vec2_p[(k-1)+(n02-1)*(rd-1)]+=-( 2*(double)(covt*(rd-covt))/(double)(rd*(rd-1))*(gsl_ran_binomial_pdf(covt,    (double)k/(double)n02,rd)+gsl_ran_binomial_pdf(rd-covt,(double)k/(double)n02,rd)));
                    vec2_h[(k-1)+(n02-1)*(rd-1)]+=-( (double)(covt*covt)/(double)(rd*(rd-1))*gsl_ran_binomial_pdf(covt,    (double)k/(double)n02,rd)+(double)(rd-covt)*(double)(rd-covt)/(double)(rd*(rd-1))*gsl_ran_binomial_pdf(rd-covt,(double)k/(double)n02,rd));
                };
                for(covt=m_bar+1;covt<=rd-m_bar-1;covt++){
                    vec02_s[rd-1]+=gsl_ran_binomial_pdf(covt,(double)k/(double)n02,rd)/k;
                    vec02_d[rd-1]+=gsl_pow_int(2*(double)(covt*(rd-covt))/(double)(rd*(rd-1))/comb2.d_p[rd-1]-1/comb2.d_t[rd-1],2)*gsl_ran_binomial_pdf(covt,(double)k/(double)n02,rd)/k;
                    vec02_h[rd-1]+=gsl_pow_int((double)(covt*covt)/(double)(rd*(rd-1))/comb2.d_hq[rd-1]-2*(double)(covt*(rd-covt))/(double)(rd*(rd-1))/comb2.d_p[rd-1],2)*gsl_ran_binomial_pdf(covt,(double)k/(double)n02,rd)/k;
                };
            };
        };
    }
    
    covmat=(double *)malloc((n01-1)*(n01-1)*sizeof(double));
    generate_covariance_matrix(covmat,n01);
    if(compute_fst) {
        covmat2=(double *)malloc((n02-1)*(n02-1)*sizeof(double));
        generate_covariance_matrix(covmat2,n02);
    }
    
    //read_line_ms(bam_file1, window_size, n01, &n_pos1, int_positions1, array_counts1);
    
    fprintf(output_stat1, "scaffold\twindow\tstart\tend\tlength\tlength_outgroup\tread_depth\tS\ttheta_FL*\tWatterson\tPi\tTajima_D\tunnormFL*test\tvar_S\tvar_Watterson\ttheta_FL\tthetaH\tthetaZE\tunnormFLtest\tunnorm_FayWu_H\tunnormZengEtest\tFayWu_H\tdiv\tnonsyn_S\tsyn_S\tnonsyn_div\tsyn_div\tlen_ns\tlen_out_ns\tthetaFL*_ns\tWatt_ns\tpi_ns\tthetaFL_ns\tthetaH_ns\tthetaZE_ns\tdiv_ns\tlen_syn\tlen_out_syn\tthetaFL*_syn\tWatt_syn\tpi_syn\tthetaFL_syn\tthetaH_syn\tthetaZE_syn\tdiv_syn\talpha\talpha_watt\talpha_pi\talpha_H\n");
    
    if(compute_fst) {
        fprintf(output_stat2, "scaffold\twindow\tstart\tend\tlength\tlength_outgroup\tread_depth\tS\ttheta_FL*\tWatterson\tPi\tTajima_D\tunnormFL*test\tvar_S\tvar_Watterson\ttheta_FL\tthetaH\tthetaZE\tunnormFLtest\tunnorm_FayWu_H\tunnormZengEtest\tFayWu_H\tdiv\tnonsyn_S\tsyn_S\tnonsyn_div\tsyn_div\tlen_ns\tlen_out_ns\tthetaFL*_ns\tWatt_ns\tpi_ns\tthetaFL_ns\tthetaH_ns\tthetaZE_ns\tdiv_ns\tlen_syn\tlen_out_syn\tthetaFL*_syn\tWatt_syn\tpi_syn\tthetaFL_syn\tthetaH_syn\tthetaZE_syn\tdiv_syn\talpha\talpha_watt\talpha_pi\talpha_H\n");

        fprintf(output_fst, "scaffold\twindow\tstart\tend\tlength\tnVariants\tpw_diff_12\tPi_1\tPi_2\tPi_a\tPi_t\tFstT\tFstA\n");
    }
    printf("Computing statistics for the window...\n");
    
    
    /* Initialize variables */
    pos=0;
    posw=0;
    oldpos=0;
    pos_base1=0;
    pos_base2=0;
    n_window=1;
    cds_start=0;
    cds_end=0;
    start = 1;
    snp_pos=0;
    
    /*Read name of scaffolds*/
    char *line_sc;
    size_t n_line_sc=1;
    line_sc=malloc(sizeof(char));
    
    getline(&line_sc,&n_line_sc,scaffold_file); //first scaffold
    sscanf(line_sc,"%s",scaffold);
    
    /* Run across all bases */
    BGZFReader reader1 = {bam_file1, 0, 0};
    BGZFReader reader2 = {bam_file2, 0, 0};
    ct1=bgzf_reader_getc(&reader1);
    bgzf_reader_ungetc(&reader1,ct1);
    if(ct1!=EOF) {
        bgzf_getdeline(&reader1,&cline,&nline,9);
        sscanf(cline,"%s\t",cchrom_next);
        strcpy(cchrom,  scaffold); //define first scaffold
        strcpy(schrom,  scaffold);
        strcpy(schrom2, scaffold);
        strcpy(gchrom,  scaffold);
        strcpy(gchrom2, scaffold);
        
        c2chrom_next[0]=EOF;
        if(compute_fst) {
            ct2=bgzf_reader_getc(&reader2);
            bgzf_reader_ungetc(&reader2,ct2);
            if(ct2!=EOF) {
                bgzf_getdeline(&reader2,&cline,&nline,9);
                sscanf(cline,"%s\t",c2chrom_next);
            }
        }
    }
    if(gff) {
        init_codons(codon_p1);
        if(compute_fst) init_codons(codon_p2);
        codon_full = 0;
        codon_frame[0]=0;
        codon_frame[1]=0;
        codon_frame[2]=0;
    }
    nscaf=0;
    for(pos=1;(ct1!=EOF || ct2!=EOF); pos++)
    {
        DEB(printf("new line\n")); //debug
        posw +=1;
        
        /* update cchrom if scaffold is different than current pos in both pileupa*/
        while(strcmp(cchrom,cchrom_next) && strcmp(cchrom,c2chrom_next)) {
            getline(&line_sc,&n_line_sc,scaffold_file);
            sscanf(line_sc,"%s",scaffold);
            if(*line_sc==0) {
                cchrom[0]=0;
                break;
            }
            else strcpy(cchrom,scaffold); //next scaffold
            if(strcmp(cchrom,cchrom_next)==0)  pos_base1=0;
            if(strcmp(cchrom,c2chrom_next)==0) pos_base2=0;
            posw = 1;
            n_window = 1;
        }
        if(*line_sc==0) {
            break;
        }
        
        // INITIALIZE EVERYTHING HERE!!!!
        //include here the windows from the bed_file (that is, one option is window size and the other option is bed_file)
        if (posw==(n_window-1)*window_size+1)
        {
            DEB(printf("inizialize\n")); //debug
            printf(" %s.%lu\t", cchrom, n_window);
            if ( n_window % 10 == 0 ) printf("\n");
            init_tests(&test1, &test2, &tests1, &testn1, &tests2, &testn2, &fst, vec_rd, vec2_rd,\
                       &psyn, &pnon, &dsyn, &dnon, &psyn2, &pnon2, &dsyn2, &dnon2, &div, &div2,\
                       max_cov, compute_fst);
            start = posw;
            snp_pos = 0;
        };
        /*IDENTIFY SELECTED SNPS IF AVAILABLE*/
        if(ext_snps==1){/*SCAFFOLDS MUST BE IN THE SAME ORDER THAN PILEUP!*/
            if(strcmp(schrom2,cchrom)) { /*a new chromosome in cchrom*/
                while (strcmp(schrom,cchrom) && pos_snp<2e9/*EOF*/) pos_snp=extract_pos_snpinput(list_snps,schrom);
                strcpy(schrom2,schrom);
            }
            while (pos_snp<posw && !strcmp(schrom,cchrom)) pos_snp=extract_pos_snpinput(list_snps,schrom);
            //if the new snp defined is out of the current chromosome, schrom2 is the current chrom in pos_snp
        };
        /*DEFINE CDS FROM GTF IF AVAILABLE4*/
        if(if_gff==1){
            /* READ multiscaffold*/
            /* GTF/GFF3 MUST BE SORTED. DESIRABLE NOT TO INCLUDE ALTERNATIVE SPLICING CDS.
             ONLY FIRST CDS (AND LONGEST END) PER GENE ARE READ*/
            if(strcmp(gchrom2,cchrom)) { /*a new chromosome in cchrom*/
                while (strcmp(gchrom,cchrom) && cds_end!=2e9/*EOF*/) {
                    cds_end = extract_gff(gff, gchrom, &cds_start, &phase_cds, &strand);}
                strcpy(gchrom2,gchrom);
            }
            while (cds_end<posw && !strcmp(gchrom,cchrom)) {
                cds_end = extract_gff(gff, gchrom, &cds_start, &phase_cds, &strand);}
        }
        /*DEFINE OUTGROUP IF AVAILABLE*/
        if(outgroup_available==1) {/*define scaffold*/
            /*READ MULTISCAFFOLD FASTA OUTGROUP*/
            while (strcmp(ochrom,cchrom) && nscaf<nscafc) { /*assign a new chromosome in cchrom*/
                fseek(fasta_out,scaffold_out[nscaf],SEEK_SET);//file pointer located at the next '>' character
                getline(&oline,&noline,fasta_out); //leave file pointer located at the first row of the fasta scaffold
                sscanf(oline,">%s",ochrom);
                nscaf++;
                oldpos=0;
            }
        }
        /*READ MPILEUP FILE(S)*/
        r1_ok = 0;
        if (pos_base1<posw && strcmp(cchrom,cchrom_next)==0)
        {
            DEB(printf("reading new base from file 1\n")); //debug
            strcpy(cchrom2,cchrom_next);
            if(read_line_pileup(&reader1, min_qual, min_mqual, &pos_base1, &n_ref1, &n_alt_allele1, &rd1, n_alt_1, &ref_base1, &alt_base1, cchrom_next,m_bar)==0){
                //printf("Exiting file\n");
                ct1=EOF;
            }
            else r1_ok = 1;
        }

        if(compute_fst) {
            r2_ok = 0;
            if (pos_base2<posw && strcmp(cchrom,c2chrom_next)==0)
            {
                DEB(printf("reading new base from file 2\n")); //debug
                strcpy(c2chrom2,c2chrom_next);
                if(read_line_pileup(&reader2, min_qual, min_mqual, &pos_base2, &n_ref2, &n_alt_allele2, &rd2, n_alt_2, &ref_base2, &alt_base2, c2chrom_next,m_bar)==0) {
                    //printf("Exiting file\n");
                    ct2=EOF;
                }
                else r2_ok = 1;
            }
        }
        if(outgroup_available==1){
            out_base=extract_outgroup_base(fasta_out,posw,oldpos,fasta_length); //1;
            oldpos=posw;
        }
        
        if (if_gff==1){
            if (cds_start<=posw){
                /*if gff and cds start, it is also necessary to know the outgroup positions*/
                if(strand=='+'){
                    frame=((posw-cds_start+3-phase_cds)%3)+1;
                } else {
                    if(strand=='-'){
                        frame=((cds_end+3-phase_cds-posw)%3)+1;
                    } else frame=0;
                };
                //keep codon information until having all three positions. Then calculate variability
                if(frame) {
                    if(r1_ok && rd1>=max(min_cov,2*m_bar+2) && (pos_base1==posw) &&!strcmp(cchrom,cchrom2)) {
                        codon_p1[frame-1].n = n01;
                        codon_p1[frame-1].n_ref = n_ref1;
                        codon_p1[frame-1].n_alt_allele = n_alt_allele1;
                        codon_p1[frame-1].rd = rd1;
                        int na;
                        for(na=0;na<5;na++) codon_p1[frame-1].n_alt[na] = n_alt_1[na];
                        codon_p1[frame-1].ref_base = ref_base1;
                        codon_p1[frame-1].alt_base = alt_base1;
                        codon_p1[frame-1].out_base = out_base;
                    }
                    if(compute_fst && r2_ok && rd2>=max(min_cov,2*m_bar+2) && (pos_base2==posw) &&!strcmp(cchrom,c2chrom2)) {
                        codon_p2[frame-1].n = n02;
                        codon_p2[frame-1].n_ref = n_ref2;
                        codon_p2[frame-1].n_alt_allele = n_alt_allele2;
                        codon_p2[frame-1].rd = rd2;
                        int na;
                        for(na=0;na<5;na++) codon_p2[frame-1].n_alt[na] = n_alt_2[na];
                        codon_p2[frame-1].ref_base = ref_base2;
                        codon_p2[frame-1].alt_base = alt_base2;
                        codon_p2[frame-1].out_base = out_base;
                    }
                    
                    codon_frame[codon_full] = frame;
                    codon_full +=1;
                    //in case having different sequence than 123 positions. Re-init codons and eliminate from analysis
                    //codon_frame must be in order 1,2,3 or 3,2,1
                    if(!((strand=='+' && codon_frame[codon_full-1]==codon_full) ||
                         (strand=='-' && codon_frame[codon_full-1]==3-(codon_full-1)))) {
                        //fprintf(stderr,"Frame broken at %c strand: discard codon at position %lu\n",strand,posw);
                        init_codons(codon_p1);
                        if(compute_fst) init_codons(codon_p2);
                        codon_full = 0;
                        codon_frame[0]=0;
                        codon_frame[1]=0;
                        codon_frame[2]=0;
                    }
                }
            } else frame=0;
            
            /*EXTRACT STATISTICS FOR CODON POSITIONS (SYN/NSYN)*/
            //calculate variability for codons and initialize again
            if(codon_full == 3) {
                //obtain weights for codon positions (syn/nsyn) for p1
                codon_weights_muts(cweights,dweights,cmutnsyn,kmutnsyn,cmutsyn,kmutsyn,codon_p1,strand);
                //extract statistics for codon positions in p1
                int pp;
                for(pp=0;pp<3;pp++) {
                    dsyn += kmutsyn[pp];
                    dnon += kmutnsyn[pp];
                    psyn += cmutsyn[pp];
                    pnon += cmutnsyn[pp];
                    //if(cmutsyn[pp])  printf("S:%lu\n",posw-2+pp);
                    //if(cmutnsyn[pp]) printf("N:%lu\n",posw-2+pp);
                    extract_stats(&tests1, &comb1, codon_p1[pp].n, codon_p1[pp].n_ref, codon_p1[pp].n_alt_allele, codon_p1[pp].rd, codon_p1[pp].n_alt, codon_p1[pp].ref_base, codon_p1[pp].alt_base, codon_p1[pp].out_base, m_bar,1.-cweights[pp],1.-dweights[pp],cmutsyn[pp]);
                    extract_stats(&testn1, &comb1, codon_p1[pp].n, codon_p1[pp].n_ref, codon_p1[pp].n_alt_allele, codon_p1[pp].rd, codon_p1[pp].n_alt, codon_p1[pp].ref_base, codon_p1[pp].alt_base, codon_p1[pp].out_base, m_bar,cweights[pp],dweights[pp],cmutnsyn[pp]);
                }
                
                if(compute_fst) {
                    //obtain weights for codon positions (syn/nsyn) for p1
                    codon_weights_muts(cweights,dweights,cmutnsyn,kmutnsyn,cmutsyn,kmutsyn,codon_p2,strand);
                    //extract statistics for codon positions in p1
                    int pp;
                    for(pp=0;pp<3;pp++) {
                        dsyn2 += kmutsyn[pp];
                        dnon2 += kmutnsyn[pp];
                        psyn2 += cmutsyn[pp];
                        pnon2 += cmutnsyn[pp];
                        extract_stats(&tests2, &comb2, codon_p2[pp].n, codon_p2[pp].n_ref, codon_p2[pp].n_alt_allele, codon_p2[pp].rd, codon_p2[pp].n_alt, codon_p2[pp].ref_base, codon_p2[pp].alt_base, codon_p2[pp].out_base, m_bar,1.-cweights[pp],1.-dweights[pp],cmutsyn[pp]);
                        extract_stats(&testn2, &comb2, codon_p2[pp].n, codon_p2[pp].n_ref, codon_p2[pp].n_alt_allele, codon_p2[pp].rd, codon_p2[pp].n_alt, codon_p2[pp].ref_base, codon_p2[pp].alt_base, codon_p2[pp].out_base, m_bar,cweights[pp],dweights[pp],cmutnsyn[pp]);
                    }

                }
                //re-initialize codon structs
                init_codons(codon_p1);
                if(compute_fst) init_codons(codon_p2);
                codon_full = 0;
                codon_frame[0]=0;
                codon_frame[1]=0;
                codon_frame[2]=0;
            }
        };
        
        /*EXTRACT STATISTICS FOR TOTAL POSITIONS*/
        /*POP1*/
        if ((pos_base1==posw) && (rd1>=min_cov) && (rd1<=max_cov) && !strcmp(cchrom,cchrom2))
        {
            if((ext_snps==1)&&((pos_snp!=posw)||strcmp(schrom,cchrom))) {
                if (n_ref1>=n_alt_allele1){
                    n_ref1+=n_alt_allele1;
                    n_alt_allele1=0;
                } else {
                    n_alt_allele1+=n_ref1;
                    n_ref1=0;
                };
            };
            //printf("out base: %u\n",out_base); //debug
            //printf("sample 1: "); //debug
            if (rd1>=max(min_cov,2*m_bar+2))
            {
                extract_stats(&test1, &comb1, n01, n_ref1, n_alt_allele1, rd1, n_alt_1, ref_base1, alt_base1, out_base, m_bar,1.0,1.0,1.0); //+rd1*(pos_snp!=pos));
                if(out_base!=0){
                    if (((n_alt_allele1>=rd1-m_bar)&&(alt_base1!=out_base))||((n_ref1>=rd1-m_bar)&&(ref_base1!=out_base))){
                        div++;
                    };
                }
                //SNPS
#if PRINTSNPS == 1
                extract_snps(pos, output_snps, n_alt_1, n_alt_2, m_bar);
#endif
               vec_rd[rd1-1]++;
            };
        };
        /*EXTRACT STATISTICS*/
        /*POP2*/
        if(compute_fst) {
            if ((pos_base2==posw)&&(rd2>=min_cov)&&(rd2<=max_cov)&&!strcmp(cchrom,c2chrom2))
            {
                if((ext_snps==1)&&((pos_snp!=posw)||strcmp(schrom,cchrom))) {
                    if (n_ref2>=n_alt_allele2){
                        n_ref2+=n_alt_allele2;
                        n_alt_allele2=0;
                    } else {
                        n_alt_allele2+=n_ref2;
                        n_ref2=0;
                    };
                };
                //printf("out base: %u\n",out_base); //debug
                //printf("sample 2: "); //debug
                if (rd2>=max(min_cov,2*m_bar+2))
                {
                    extract_stats(&test2, &comb2, n02, n_ref2, n_alt_allele2, rd2, n_alt_2, ref_base2, alt_base2, out_base, m_bar,1.0,1.0,1.0); //+rd1*(pos_snp!=pos));
                    if(out_base!=0){
                        if (((n_alt_allele2>=rd2-m_bar)&&(alt_base2!=out_base))||((n_ref2>=rd2-m_bar)&&(ref_base2!=out_base))){
                            div2++;
                        };
                    }
                    vec2_rd[rd2-1]++;
                };
            }
            /*CALCULATE FST if there are reads in both pileups */ /*
            if (((pos_base1==posw)&&(rd1>=min_cov)&&(rd1<=max_cov)&&!strcmp(cchrom, cchrom2)&&(rd1>=max(min_cov,2*m_bar+1))) &&  ((pos_base2==posw)&&(rd2>=min_cov)&&(rd2<=max_cov)&&!strcmp(cchrom,c2chrom2)&&(rd2>=max(min_cov,2*m_bar+1)))) {
                extract_fst(&fst, &combfst, n01, n_ref1, n_alt_allele1, rd1, n_alt_1, ref_base1, alt_base1, n02, n_ref2, n_alt_allele2, rd2, n_alt_2, ref_base2, alt_base2, out_base, m_bar, &variant);
                snp_pos += (unsigned long) variant;
            }*/
            if (((pos_base1==posw)&&(rd1<=max_cov)&&!strcmp(cchrom, cchrom2)/*&&(rd1>=max(min_cov,m_bar+1))*/) &&  ((pos_base2==posw)&&(rd2<=max_cov)&&!strcmp(cchrom,c2chrom2)/*&&(rd2>=max(min_cov,m_bar+1))*/)) {
                extract_fstT(&fst, &combfst, &comb1, &comb2, n01, n_ref1, n_alt_allele1, rd1, n_alt_1, ref_base1, alt_base1, n02, n_ref2, n_alt_allele2, rd2, n_alt_2, ref_base2, alt_base2, out_base, m_bar,m_bar_fst, &variant, fst_path);
                snp_pos += (unsigned long) variant;
            }
        }
        /* PRINT OUTPUT(S) */
        ct1=bgzf_reader_getc(&reader1);
        bgzf_reader_ungetc(&reader1,ct1);
        if(compute_fst) {
            ct2=bgzf_reader_getc(&reader2);
            bgzf_reader_ungetc(&reader2,ct2);
        }
        //include here the windows from the bed_file (that is, one option is window size and the other option is bed_file)
        if (((posw==(n_window*window_size)) || (ct1==EOF && ct2==EOF) || (strcmp(cchrom_next,cchrom) && strcmp(c2chrom_next,cchrom))))
        {
            //calculate window values
            end = posw;
            double pi1t_val,pi2t_val,pit_val,pi1t_val_,pi2t_val_,pit_val_,theta1_val, pi1_val, d1_val, thetaH1_val, h1_val, theta2_val, pi2_val, d2_val, thetaH2_val, h2_val, pia_val, pis_val,pia_val_, pis_val_, fst_val,fst_val2, cov1_val, cov2_val, div_val, var_h, var_d, var_s, var0_s, var0_d, var0_h, thetaE1_val,thetaE2_val/**/, thetanu1_val, thetaxi1_val, e1_val, f1_val, fo1_val, thetanu2_val, thetaxi2_val, e2_val, f2_val, fo2_val/*, vk_s[n01-1], vk_d[n01-1], vk_h[n01-1]*//*,pis_val_ponderated*/;
            DEB(printf("printing output\n")); //debug
            
            /*POP1*/
            calculate_vartests(&var_h, &var_d, &var_s, &var0_s, &var0_d, &var0_h, &test1, vec_rd, vec_s, vec_p, vec_h, vec0_s, vec0_d, vec0_h, covmat, n01, min_cov, max_cov, m_bar);
            /*TOTAL_STATS*/
            calculate_window_stats(&test1,&cov1_val,&theta1_val,&pi1_val,&d1_val,&thetaH1_val,&h1_val,&div_val,div/**/,&thetaE1_val/**/,&thetanu1_val, &thetaxi1_val, &e1_val, &f1_val, &fo1_val/**/);
            
            var0_s=var0_s*theta1_val;
            var0_d=var0_d*theta1_val;
            var0_h=var0_h*theta1_val;
            var_s=var_s*theta1_val*theta1_val;
            var_d=var_d*theta1_val*theta1_val;
            var_h=var_h*theta1_val*theta1_val;
            
            pi1t_val=pi1_val;
            
            //PRINT RESULTS
            fprintf(output_stat1, "%s\t%lu\t%lu\t%lu\t%.0f\t%.0f",cchrom2, n_window, start, end, test1.l, test1.l_out);
            /*printf("%lu\t%lu", test1.l, test1.l_out);*/
            
            if (test1.l>0.) {
                if(thetanu1_val!=-1)
                    fprintf(output_stat1, "\t%f\t%.0f\t%f\t%f\t%f\t%f\t%f\t%.3e\t%.3e", cov1_val, test1.s, thetanu1_val, theta1_val, pi1_val, d1_val/sqrt(var0_d+var_d),f1_val,var0_s+var_s, (var0_s+var_s)/(test1.den_t*test1.den_t));
                else
                    fprintf(output_stat1, "\t%f\t%.0f\tNA\t%f\t%f\t%f\tNA\t%.3e\t%.3e", cov1_val, test1.s, theta1_val, pi1_val, d1_val/sqrt(var0_d+var_d),var0_s+var_s, (var0_s+var_s)/(test1.den_t*test1.den_t));
            } else {
                fprintf(output_stat1, "\tNA\t0\tNA\tNA\tNA\tNA\tNA\tNA\tNA");
            };
            if (test1.l_out>0.) {
                if(thetaxi1_val!=-1)
                    fprintf(output_stat1, "\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f", thetaxi1_val,thetaH1_val, thetaE1_val,fo1_val,h1_val,e1_val, h1_val/sqrt(var0_h+var_h), div_val);
                else
                    fprintf(output_stat1, "\tNA\t%f\t%f\tNA\t%f\t%f\t%f\t%f",thetaH1_val, thetaE1_val,h1_val,e1_val, h1_val/sqrt(var0_h+var_h), div_val);
            } else {
                fprintf(output_stat1, "\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA");
            };
            if (if_gff==1) {
                //re-initialize codon structs
                init_codons(codon_p1);
                if(compute_fst) init_codons(codon_p2);
                codon_full = 0;
                codon_frame[0]=0;
                codon_frame[1]=0;
                codon_frame[2]=0;

                //print results
                fprintf(output_stat1, "\t%.1f\t%.1f\t%.1f\t%.1f", pnon, psyn, dnon, dsyn);
                
                /*NONSYN STATS*/
                calculate_window_stats(&testn1,&cov1_val,&theta1_val,&pi1_val,&d1_val,&thetaH1_val,&h1_val,&div_val,dnon/**/,&thetaE1_val/**/,&thetanu1_val, &thetaxi1_val, &e1_val, &f1_val, &fo1_val/**/);
                
                //PRINT RESULTS
                fprintf(output_stat1, "\t%.1f\t%.1f",testn1.l, testn1.l_out);
                if (testn1.l>0.) {
                    if(thetanu1_val!=-1)
                        fprintf(output_stat1, "\t%f\t%f\t%f", thetanu1_val, theta1_val, pi1_val);
                    else
                        fprintf(output_stat1, "\tNA\t%f\t%f", theta1_val, pi1_val);
                } else {
                    fprintf(output_stat1, "\tNA\tNA\tNA");
                };
                if (testn1.l_out>0.) {
                    if(thetaxi1_val!=-1)
                        fprintf(output_stat1, "\t%f\t%f\t%f\t%f", thetaxi1_val,thetaH1_val,thetaE1_val, div_val);
                    else
                        fprintf(output_stat1, "\tNA\t%f\t%f\t%f",thetaH1_val,thetaE1_val, div_val);
                } else {
                    fprintf(output_stat1, "\tNA\tNA\tNA\tNA");
                };
                thetaw_n=theta1_val;
                thetat_n=pi1_val;
                thetah_n=h1_val+pi1_val;
                div_n=div_val;
                
                /*SYN STATS*/
                calculate_window_stats(&tests1,&cov1_val,&theta1_val,&pi1_val,&d1_val,&thetaH1_val,&h1_val,&div_val,dsyn/**/,&thetaE1_val/**/,&thetanu1_val, &thetaxi1_val, &e1_val, &f1_val, &fo1_val/**/);
                
                fprintf(output_stat1, "\t%.1f\t%.1f",tests1.l, tests1.l_out);
                if (tests1.l>0.) {
                    if(thetanu1_val!=-1)
                        fprintf(output_stat1, "\t%f\t%f\t%f", thetanu1_val, theta1_val, pi1_val);
                    else
                        fprintf(output_stat1, "\tNA\t%f\t%f", theta1_val, pi1_val);
                } else {
                    fprintf(output_stat1, "\tNA\tNA\tNA");
                };
                if (tests1.l_out>0.) {
                    if(thetaxi1_val!=-1)
                        fprintf(output_stat1, "\t%f\t%f\t%f\t%f", thetaxi1_val,thetaH1_val,thetaE1_val, div_val);
                    else
                        fprintf(output_stat1, "\tNA\t%f\t%f\t%f",thetaH1_val,thetaE1_val, div_val);
                } else {
                    fprintf(output_stat1, "\tNA\tNA\tNA\tNA");
                };
                thetaw_s=theta1_val;
                thetat_s=pi1_val;
                thetah_s=h1_val+pi1_val;
                div_s=div_val;
                
                //PRINT RESULTS
                if(psyn*dnon==0){
                    if(dsyn*pnon==0){
                        fprintf(output_stat1, "\tNA\tNA\tNA\tNA");
                    } else {
                        fprintf(output_stat1, "\t-Inf\t-Inf\t-Inf\t-Inf");
                    };
                } else {
                    fprintf(output_stat1, "\t%f", 1.0-(double)(dsyn*pnon)/(double)(psyn*dnon));
                    fprintf(output_stat1, "\t%f", 1.0-(double)((div_s)*(thetaw_n))/(double)((thetaw_s)*(div_n)));
                    fprintf(output_stat1, "\t%f", 1.0-(double)((div_s)*(thetat_n))/(double)((thetat_s)*(div_n)));
                    fprintf(output_stat1, "\t%f", 1.0-(double)((div_s)*(thetah_n))/(double)((thetah_s)*(div_n)));
                };
            } else {
                fprintf(output_stat1, "\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA");
            };
            fprintf(output_stat1, "\n");
            
            if(compute_fst) {
                /*POP2*/
                calculate_vartests(&var_h, &var_d, &var_s, &var0_s, &var0_d, &var0_h, &test2, vec2_rd, vec2_s, vec2_p, vec2_h, vec02_s, vec02_d, vec02_h, covmat2, n02, min_cov, max_cov, m_bar);
                /*TOTAL STATS*/
                calculate_window_stats(&test2,&cov2_val,&theta2_val,&pi2_val,&d2_val,&thetaH2_val,&h2_val,&div_val,div2/**/,&thetaE2_val/**/,&thetanu2_val, &thetaxi2_val, &e2_val, &f2_val, &fo2_val/**/);
                
                pi2t_val=pi2_val;
                //if((test2.den_hl>0)&&(test2.den_hq>0)) { h2_val=test2.num_hq/test2.den_hq-test2.num_hl/test2.den_hl; } else { h2_val=0; };
                
                //Calculate Fst:
                if(fst.den_p1>0 && fst.den_p2>0) {
                    pis_val_ = (fst.num_p1/fst.den_p1 + fst.num_p2/fst.den_p2)/2.0; //12
                } else {pis_val_=-1;}
                if(fst.gen_diff>=0) {
                    pia_val_ = fst.gen_diff/(fst.la -  fst.c_s); //SI-27
                } else {pia_val_=-1;}
                if(pis_val_>=0 && pia_val_>=0) {
                    pit_val_ = 1.0/2.0 * (pis_val_ + pia_val_); //13
                } else {pit_val_=-1;}
                if(pit_val_>0) {
                     fst_val2 = 1.0 - pis_val_ / pit_val_; //11
                } else {
                    fst_val2 = -10000;
                }

                //Calculate Fst: Using new Luca approach?
                if(fst.num_p1g>=0 && fst.num_p2g>=0) {
                    pis_val = (fst.num_p1g + fst.num_p2g)/2.0;
                } else {pis_val=-1;}
                if(fst.gen_diff>=0) {
                    pia_val = fst.gen_diff;
                } else {pia_val=-1;}
                if(pis_val>=0 && pia_val>=0) {
                    pit_val = 1.0/2.0 * (pis_val + pia_val);
                } else {pit_val=-1;}
                if(pit_val>0) {
                    //fst_val = 1 - pis_val / pit_val; //Grenedalf ?
                    fst_val = fst.var_diffL / pit_val;
                } else {
                    fst_val = -10000;
                }
                
                
                var0_s=var0_s*theta2_val;
                var0_d=var0_d*theta2_val;
                var0_h=var0_h*theta2_val;
                var_s=var_s*theta2_val*theta2_val;
                var_d=var_d*theta2_val*theta2_val;
                var_h=var_h*theta2_val*theta2_val;
                
                //PRINT RESULTS
                fprintf(output_stat2, "%s\t%lu\t%lu\t%lu\t%.0f\t%.0f",cchrom2, n_window, start, end, test2.l, test2.l_out);
                if (test2.l>0.) {
                    if(thetanu2_val!=-1)
                        fprintf(output_stat2, "\t%f\t%.0f\t%f\t%f\t%f\t%f\t%f\t%.3e\t%.3e", cov2_val, test2.s, thetanu2_val,theta2_val, pi2_val, d2_val/sqrt(var0_d+var_d), f2_val,var0_s+var_s, (var0_s+var_s)/(test2.den_t*test2.den_t));
                    else
                        fprintf(output_stat2, "\t%f\t%.0f\tNA\t%f\t%f\t%f\tNA\t%.3e\t%.3e", cov2_val, test2.s,theta2_val, pi2_val, d2_val/sqrt(var0_d+var_d),var0_s+var_s, (var0_s+var_s)/(test2.den_t*test2.den_t));
                } else {
                    fprintf(output_stat2, "\tNA\t0\tNA\tNA\tNA\tNA\tNA\tNA\tNA");
                };
                if (test1.l_out>0.) {
                    if(thetaxi2_val!=-1)
                        fprintf(output_stat2, "\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f", thetaxi2_val,thetaH2_val,thetaE2_val, fo2_val,h2_val, e2_val,h2_val/sqrt(var0_h+var_h), div_val);
                    else
                        fprintf(output_stat2, "\tNA\t%f\t%f\tNA\t%f\t%f\t%f\t%f",thetaH2_val,thetaE2_val, h2_val, e2_val,h2_val/sqrt(var0_h+var_h), div_val);
                } else {
                    fprintf(output_stat2, "\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA");
                };
                if (if_gff==1) {
                    fprintf(output_stat2, "\t%.1f\t%.1f\t%.1f\t%.1f", pnon2, psyn2, dnon2, dsyn2);
                    
                    /*NONSYN STATS*/
                    calculate_window_stats(&testn2,&cov2_val,&theta2_val,&pi2_val,&d2_val,&thetaH2_val,&h2_val,&div_val,dnon2/**/,&thetaE2_val/**/,&thetanu2_val, &thetaxi2_val, &e2_val, &f2_val, &fo2_val/**/);
                    
                    //PRINT RESULTS
                    fprintf(output_stat2, "\t%.1f\t%.1f",testn2.l, testn2.l_out);
                    if (testn2.l>0.) {
                        if(thetanu2_val!=-1)
                            fprintf(output_stat2, "\t%f\t%f\t%f", thetanu2_val,theta2_val, pi2_val);
                        else
                            fprintf(output_stat2, "\tNA\t%f\t%f",theta2_val, pi2_val);
                    } else {
                        fprintf(output_stat2, "\tNA\tNA\tNA");
                    };
                    if (testn2.l_out>0.) {
                        if(thetaxi2_val!=-1)
                            fprintf(output_stat2, "\t%f\t%f\t%f\t%f", thetaxi2_val,thetaH2_val, thetaE2_val,div_val);
                        else
                            fprintf(output_stat2, "\tNA\t%f\t%f\t%f",thetaH2_val, thetaE2_val,div_val);
                    } else {
                        fprintf(output_stat2, "\tNA\tNA\tNA\tNA");
                    };
                    thetaw_n=theta2_val;
                    thetat_n=pi2_val;
                    thetah_n=h2_val+pi2_val;
                    div_n=div_val;
                    
                    /*SYN STATS*/
                    calculate_window_stats(&tests2,&cov2_val,&theta2_val,&pi2_val,&d2_val,&thetaH2_val,&h2_val,&div_val,dsyn2/**/,&thetaE2_val/**/,&thetanu2_val, &thetaxi2_val, &e2_val, &f2_val, &fo2_val/**/);
                    
                    //PRINT RESULTS
                    fprintf(output_stat2, "\t%.1f\t%.1f",tests2.l, tests2.l_out);
                    if (tests2.l>0.) {
                        if(thetanu2_val!=-1)
                            fprintf(output_stat2, "\t%f\t%f\t%f", thetanu2_val,theta2_val, pi2_val);
                        else
                            fprintf(output_stat2, "\tNA\t%f\t%f",theta2_val, pi2_val);
                    } else {
                        fprintf(output_stat2, "\tNA\tNA\tNA");
                    };
                    if (tests2.l_out>0.) {
                        if(thetaxi2_val!=-1)
                            fprintf(output_stat2, "\t%f\t%f\t%f\t%f", thetaxi2_val,thetaH2_val,thetaE2_val, div_val);
                        else
                            fprintf(output_stat2, "\tNA\t%f\t%f\t%f",thetaH2_val,thetaE2_val, div_val);
                    } else {
                        fprintf(output_stat2, "\tNA\tNA\tNA\tNA");
                    };
                    thetaw_s=theta2_val;
                    thetat_s=pi2_val;
                    thetah_s=h2_val+pi2_val;
                    div_s=div_val;
                    
                    if(psyn2*dnon2==0){
                        if(dsyn2*pnon2==0){
                            fprintf(output_stat2, "\tNA\tNA\tNA\tNA");
                        } else {
                            fprintf(output_stat2, "\t-Inf\t-Inf\t-Inf\t-Inf");
                        };
                    } else {
                        fprintf(output_stat2, "\t%f", 1.0-(double)(dsyn*pnon)/(double)(psyn*dnon));
                        fprintf(output_stat2, "\t%f", 1.0-(double)((div_s)*(thetaw_n))/(double)((thetaw_s)*(div_n)));
                        fprintf(output_stat2, "\t%f", 1.0-(double)((div_s)*(thetat_n))/(double)((thetat_s)*(div_n)));
                        fprintf(output_stat2, "\t%f", 1.0-(double)((div_s)*(thetah_n))/(double)((thetah_s)*(div_n)));
                    };
                    
                } else {
                    fprintf(output_stat2, "\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA");
                };
                fprintf(output_stat2,"\n");
                
                /*RESULTS FST FILE */
                fprintf(output_fst, "%s\t%lu\t%lu\t%lu",cchrom2,n_window, start, end);
                if(fst_path==1) {//Luca approach? // Grenedalf approach?
                    fprintf(output_fst, "\t%lu",fst.la);
                    if(fst.la>0) {
                        fprintf(output_fst,"\t%lu\t%f",snp_pos,fst.var_diffL);//fst.gen_diff);
                        if(fst.num_p1g>=0.) fprintf(output_fst,"\t%f",fst.num_p1g); else fprintf(output_fst,"\tNA");
                        if(fst.num_p2g>=0.) fprintf(output_fst,"\t%f",fst.num_p2g); else fprintf(output_fst,"\tNA");
                        if(pia_val>=0.)  fprintf(output_fst,"\t%f",pia_val);  else fprintf(output_fst,"\tNA");
                        if(pit_val>=0.)  fprintf(output_fst,"\t%f",pit_val);  else fprintf(output_fst,"\tNA");
                        if(fst_val !=-10000) fprintf(output_fst, "\t%f",fst_val); else fprintf(output_fst,"\tNA");
                        if(pia_val >0.) fprintf(output_fst, "\t%f",1.- pis_val/pia_val); else fprintf(output_fst,"\tNA");
                    }
                    else {
                        fprintf(output_fst, "\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA");
                    }
                }
                else {
                    fprintf(output_fst, "\t%lu",fst.la);
                    if(fst.la>0) { //initial npstat approach
                        fprintf(output_fst,"\t%lu\t%f",snp_pos,fst.gen_diff);
                        if(fst.den_p1>0.) fprintf(output_fst,"\t%f",fst.num_p1/fst.den_p1); else fprintf(output_fst,"\tNA");
                        if(fst.den_p2>0.) fprintf(output_fst,"\t%f",fst.num_p2/fst.den_p2); else fprintf(output_fst,"\tNA");
                        if(pia_val_>=0.)  fprintf(output_fst,"\t%f",pia_val_);  else fprintf(output_fst,"\tNA");
                        if(pit_val_>=0.)  fprintf(output_fst,"\t%f",pit_val_);  else fprintf(output_fst,"\tNA");
                        if(fst_val2 !=-10000) fprintf(output_fst, "\t%f",fst_val2); else fprintf(output_fst,"\tNA");
                        if(pia_val_>0.) fprintf(output_fst, "\t%f",1.- pis_val_/pia_val_); else fprintf(output_fst,"\tNA");
                    }
                    else {
                        fprintf(output_fst, "\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA");
                    }
                }
                fprintf(output_fst,"\n");
            }
            
            if(strcmp(cchrom_next,cchrom) && strcmp(c2chrom_next,cchrom)) {
                printf("\n");
            }
            n_window++;
        };
    };
    /* Close files, free gsl random number generator */
    bgzf_close(bam_file1);
    if(compute_fst) {
        bgzf_close(bam_file2);
    }
    fclose(scaffold_file);
    if (outgroup_available==1) { fclose(fasta_out); };
    if (ext_snps==1) { fclose(list_snps); };
    if (if_gff==1) { fclose(gff); };
    fclose(output_stat1);
    if(compute_fst) {
        fclose(output_stat2);
        fclose(output_fst);
    }
    gsl_rng_free(r);
    free(cline);
    free(oline);
    free(cchrom);free(cchrom2);
    free(cchrom_next); free(c2chrom_next);
    free(schrom);free(schrom2);
    free(gchrom);free(gchrom2);
    free(ochrom);
    free(scaffold_out);
    free(outfile);
    free(scaffold);
    free(scaffold_filename);
    free(bed_filename);
    free(line_sc);
    free(codon_frame);
    free(cweights);
    free(dweights);
    free(cmutnsyn);
    free(kmutnsyn);
    free(cmutsyn);
    free(kmutsyn);

    free(vec_s); free(vec_p); free(vec_h);
    free(vec0_s); free(vec0_d); free(vec0_h);
    free(covmat);
    if(compute_fst) {
        free(vec2_s); free(vec2_p); free(vec2_h);
        free(vec02_s); free(vec02_d); free(vec02_h);
        free(covmat2);
        for(count_i=0;count_i<max_cov;count_i++) {
            free(combfst.c_s[count_i]);
        };
        free(combfst.c_s);
        free(combfst.d_pt);
    }
    
    //SNPS
#if PRINTSNPS == 1
    fclose(output_snps);
#endif
    printf("\n");
    printf("Computation of statistics completed.\n");
}
