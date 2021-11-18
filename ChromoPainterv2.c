#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <ctype.h>
#include <time.h>
#include <string.h>
#include "zlib.h"

#define PI 3.141593

/*****************************************************************
// TAKES A PHASE-STYLE SET OF n PHASED HAPLOTYPES, SUCH THAT THE TOP x HAPLOTYPES ARE POTENTIAL DONORS TO THE BOTTOM n-x HAPLOTYPES (THESE n-x HAPLOTYPES ARE ASSUMED TO BE DIPLOID INDIVIDUALS, WITH EACH INDIVIDUAL'S TWO HAPLOTYPES ON CONSECUTIVE LINES). THEN PAINTS EACH OF THE n-x BOTTOM HAPLOTYPES AS A MOSAIC OF THE x DONOR HAPLOTYPES, WITH OR WITHOUT MAXIMIXING OVER N_e (DEFAULT) OR COPYING PROPORTIONS (IF USER SELECTS) USING E-M, OUTPUTTING AS MANY SAMPLES AS DESIRED (AFTER MAXIMIZATION, IF CHOSEN -- MAXIMIZATION IS WHERE THE DIPLOID ASSUMPTION COMES IN, AS N_E ESTIMATES OR COPYING PROPORTIONS ARE AVERAGED OVER AN INDIVIDUAL'S TWO HAPS AT EACH E-M STEP)

// to compile:  gcc -Wall -o ChromoPainterv2 ChromoPainterv2.c -lm -lz

// to run: "./ChromoPainterv2" with following options:
//                 -g geno.filein (PHASEish-style file; no default, required)
//                 -r recommap.filein (no default, required)
//                 -t file listing id and population labels for each individual (no default, required)
//                 -f f_1 f_2 file listing donor and recipient populations; paint recipient individuals f_1 through f_2 using all donor population haplotypes (use '-f <donorlist.filein> 0 0' to paint all recipient inds) (no default, required)
//                 -i number of EM iterations for estimating parameters (default=0)									
//                 -in maximize over average switch rate parameter using E-M
//                 -ip maximize over copying proportions using E-M
//                 -im maximize over donor population mutation (emission) probabilities using E-M
//                 -iM maximize over global mutation (emission) probability using E-M
//                 -s number of samples per recipient haplotype (default=10)
//                 -n average switch rate parameter constant (default=400000 divided by total number of donor haplotypes included in analysis)
//                 -p specify to use prior copying probabilities in donor list file
//                 -m specify to use donor population mutation (emission) probabilities in population list (-f) file
//                 -M global mutation (emission) probability (default=Li & Stephen's (2003) fixed estimate)
//                 -k specify number of expected chunks to define a 'region' (default=100)
//                 -j specify that individuals are haploid
//                 -u specify that data are unlinked
//                 -a a_1 a_2 paint individuals a_1 through a_2 using every other individual (use '-a 0 0' to paint all inds)
//                 -b print-out zipped file  with suffix '.copyprobsperlocus.out' containing prob each recipient copies each donor at every SNP (note: file can be quite large)
//                 -o outfile-prefix (default = "geno.filein")
//                 --help print help menu


// example: 
//        ./ChromoPainterv2 -g example/BrahuiYorubaSimulationChrom22.haplotypes -r example/BrahuiYorubaSimulationChrom22.recomrates -f example/BrahuiYorubaSimulation.poplist.txt 0 0 -t example/BrahuiYorubaSimulation.idfile.txt -o example/BrahuiYorubaSimulationChrom22.chromopainter.out

*******************************************************************/

int reading(st, format, res)
    char **st, *format;
    void *res;
{
    int i;
    char *rs;
    rs = *st;
    for(i = 0; isspace(rs[i]); i++) ; 
    if (!rs[i]) return 0; 
    for(; !isspace(rs[i]); i++) ;
    if (rs[i]) rs[i++] = 0;  
    if (!sscanf(*st, format, res)) return 0; 
    *st += i;
    return 1;
}


struct data_t {
  int nhaps;
  int ndonorhaps;
  int nsnps;
  double *positions;
  double *lambda;
  double *copy_prob;
  double *copy_probSTART;
  double *MutProb_vec;
  int *pop_vec;
  int *hap_label_vec;
  int **cond_chromosomes;
  int **ind_chromosomes;
  double *copy_prob_new;
  double *copy_prob_newSTART;
  double *MutProb_vec_new;
  double **back_prob;
  int *ndonorhaps_vec;
};
 
struct data_t *ReadData(FILE *fd, int ind_val, int ploidy, int * include_ind_vec, char ** pop_label_vec, int num_donor_pops, char ** donor_pop_vec, int * pop_vec_tot, double * copy_prob_tot, double * copy_probSTART_tot, double * MutProb_vec_tot, int * ndonorhaps_tot, int all_versus_all_ind){
  //int line_max=100000000;
  int line_max;
  struct data_t *dat;
  char * firstline = malloc(1000 * sizeof(char));
  char *step;
  char waste[400];
  int i,j,k;
  int nind, num_cond_haps;
  int cond_hap_count, ind_hap_count, include_hap;
  
  dat=malloc(sizeof(struct data_t));
  /* if dat==NULL etc ... */

  /* Number of haplotypes */
  fgets(firstline,2047,fd);
  if (firstline==NULL) { printf("error with PHASE-style input file\n"); exit(1);}
  sscanf(firstline,"%d",&dat->nhaps);
  nind = (int) dat->nhaps/ploidy;
  if ((dat->nhaps <= 0) || (((int)dat->nhaps)!= dat->nhaps)) { printf("Number of total haplotypes must be an integer value and > 0. Exiting...\n"); exit(1);}

  /* Number of SNPs */
  fgets(firstline,2047,fd);
  if (firstline==NULL) { printf("error with PHASE-style input file\n"); exit(1);}
  sscanf(firstline,"%d",&dat->nsnps);
  if ((dat->nsnps <= 0) || (((int)dat->nsnps)!= dat->nsnps)) { printf("Number of sites must be an integer value and > 0. Exiting...\n"); exit(1);}
  line_max=dat->nsnps*8*2;
  char * line = malloc(line_max * sizeof(char));

  dat->ndonorhaps_vec=malloc(num_donor_pops*sizeof(int));
  for (k=0; k < num_donor_pops; k++)
    dat->ndonorhaps_vec[k]=0;
  num_cond_haps=0;
  for (i=0; i < nind; i++)
    {
      for (k=0; k < num_donor_pops; k++)
	{
	  if (include_ind_vec[i]!=0 && i != ind_val && all_versus_all_ind==0 && strcmp(donor_pop_vec[k],pop_label_vec[i])==0)
	    {
	      num_cond_haps=num_cond_haps+ploidy;
	      dat->ndonorhaps_vec[k]=dat->ndonorhaps_vec[k]+ploidy;
	      break;
	    }
	  if (include_ind_vec[i]!=0 && i != ind_val && all_versus_all_ind==1)
	    {
	      num_cond_haps=num_cond_haps+ploidy;
	      dat->ndonorhaps_vec[k]=dat->ndonorhaps_vec[k]+ploidy;
	      break;
	    }
	}
    }
  dat->ndonorhaps=num_cond_haps;
  dat->positions=malloc(dat->nsnps*sizeof(double));
  dat->lambda=malloc((dat->nsnps-1)*sizeof(double));
  dat->cond_chromosomes=malloc(dat->ndonorhaps*sizeof(int *));
  dat->ind_chromosomes=malloc(ploidy*sizeof(int *));
  dat->copy_prob=malloc(dat->ndonorhaps*sizeof(double));
  dat->copy_probSTART=malloc(dat->ndonorhaps*sizeof(double));
  dat->MutProb_vec=malloc(dat->ndonorhaps*sizeof(double));
  dat->pop_vec=malloc(dat->ndonorhaps*sizeof(int));
  dat->hap_label_vec=malloc(dat->ndonorhaps*sizeof(int));
  for (i=0; i<dat->ndonorhaps; i++) 
    dat->cond_chromosomes[i]=malloc(dat->nsnps*sizeof(int));
  for (i=0; i<ploidy; i++) 
    dat->ind_chromosomes[i]=malloc(dat->nsnps*sizeof(int));
  dat->copy_prob_new = malloc(dat->ndonorhaps * sizeof(double));
  dat->copy_prob_newSTART = malloc(dat->ndonorhaps * sizeof(double));
  dat->MutProb_vec_new = malloc(dat->ndonorhaps * sizeof(double));
  dat->back_prob = malloc(8 * sizeof(double *));
  for (i=0; i < 8; i++)
    dat->back_prob[i] = malloc((dat->ndonorhaps+1) * sizeof(double));

  /* Positions */
  fgets(line,line_max,fd);
  if (line==NULL) { printf("error with PHASE-style input file\n"); exit(1);}
  step=line;
  reading(&step,"%s",waste);
  for (i=0; i<dat->nsnps; i++)
    { 
      reading(&step,"%lf",&dat->positions[i]);
      if (dat->positions[i] < 0) { printf("Basepair positions must be >= 0. Exiting...\n"); exit(1);}
      if (i < (dat->nsnps-1)) dat->lambda[i] = 1.0;
    } 

  cond_hap_count=0;
  ind_hap_count=0;
  for(i = 0; i < dat->nhaps; i++)
    {
      fgets(line,line_max,fd);
      if ((line==NULL)||(strlen(line)!=(dat->nsnps+1))) { printf("error with PHASE-style input file\n"); exit(1);}
      include_hap=0;
      for (k=0; k < num_donor_pops; k++)
	{
	  if (include_ind_vec[((int)floor(i/ploidy))]!=0 && (((int)floor(i/ploidy)) != ind_val) && all_versus_all_ind==0 && strcmp(donor_pop_vec[k],pop_label_vec[((int)floor(i/ploidy))])==0)
	    {
	      include_hap=1;
	      break;
	    }
	  if (include_ind_vec[((int)floor(i/ploidy))]!=0 && (((int)floor(i/ploidy)) != ind_val) && all_versus_all_ind==1)
	    {
	      include_hap=1;
	      break;
	    }
	}
      if (include_hap==1)
	{
	  dat->hap_label_vec[cond_hap_count]=i+1;
	  dat->copy_prob[cond_hap_count]=copy_prob_tot[i];
	  dat->copy_probSTART[cond_hap_count]=copy_probSTART_tot[i];
	  dat->MutProb_vec[cond_hap_count]=MutProb_vec_tot[i];
	  if (i < (ploidy*ind_val)) dat->pop_vec[cond_hap_count]=pop_vec_tot[i];
	  if (i > (ploidy*ind_val)) dat->pop_vec[cond_hap_count]=pop_vec_tot[i]-all_versus_all_ind;
	  for (j=0; j<dat->nsnps;j++) {
	    if(line[j]=='0')
	      dat->cond_chromosomes[cond_hap_count][j]=0;
	    if(line[j]=='1')
	      dat->cond_chromosomes[cond_hap_count][j]=1;
	    if(line[j]=='A')
	      dat->cond_chromosomes[cond_hap_count][j]=2;
	    if(line[j]=='C')
	      dat->cond_chromosomes[cond_hap_count][j]=3;
	    if(line[j]=='G')
	      dat->cond_chromosomes[cond_hap_count][j]=4;
	    if(line[j]=='T')
	      dat->cond_chromosomes[cond_hap_count][j]=5;
	    if(line[j]=='?')
	      dat->cond_chromosomes[cond_hap_count][j]=9;
	    if ((line[j]!='0')&&(line[j]!='1')&&(line[j]!='A')&&(line[j]!='G')&&(line[j]!='C')&&(line[j]!='T')&&(line[j]!='?'))
	      {
		printf("Allele-type invalid for hap%d, snp%d. Exiting...\n",i+1,j+1);
		exit(1);
	      }
	  }
	  cond_hap_count=cond_hap_count+1;
	}
      if (((int)floor(i/ploidy)) == ind_val)
	{
	  for (j=0; j<dat->nsnps;j++) {
	    if(line[j]=='0')
	      dat->ind_chromosomes[ind_hap_count][j]=0;
	    if(line[j]=='1')
	      dat->ind_chromosomes[ind_hap_count][j]=1;
	    if(line[j]=='A')
	      dat->ind_chromosomes[ind_hap_count][j]=2;
	    if(line[j]=='C')
	      dat->ind_chromosomes[ind_hap_count][j]=3;
	    if(line[j]=='G')
	      dat->ind_chromosomes[ind_hap_count][j]=4;
	    if(line[j]=='T')
	      dat->ind_chromosomes[ind_hap_count][j]=5;
	    if(line[j]=='?')
	      dat->ind_chromosomes[ind_hap_count][j]=9;
	    if ((line[j]!='0')&&(line[j]!='1')&&(line[j]!='A')&&(line[j]!='G')&&(line[j]!='C')&&(line[j]!='T')&&(line[j]!='?'))
	      {
		printf("Allele-type invalid for hap%d, snp%d. Exiting...\n",i+1,j+1);
		exit(1);
	      }
	  }
	  ind_hap_count=ind_hap_count+1;
	}
    }
  if (ind_hap_count<ploidy)
    {
      printf("Could not find haplotype(s) of individual %d in genotype input file! Row missing?? Exiting...\n",ind_val+1);
      exit(1);
    }
  if (cond_hap_count!=dat->ndonorhaps)
    {
      printf("Something wrong with genotype input file. Exiting...\n");
      exit(1);
    }

  free(firstline);
  free(line);
  return dat;

}

void DestroyData(struct data_t *dat, int ploidy)
{
  int i;
  for (i=0; i<dat->ndonorhaps; i++) 
    free(dat->cond_chromosomes[i]);
  for (i=0; i<ploidy; i++) 
    free(dat->ind_chromosomes[i]);
  free(dat->cond_chromosomes);
  free(dat->ind_chromosomes);
  free(dat->positions);
  free(dat->lambda);
  free(dat->copy_prob);
  free(dat->copy_probSTART);
  free(dat->MutProb_vec);
  free(dat->pop_vec);
  free(dat->hap_label_vec);
  free(dat->copy_prob_new);
  free(dat->copy_prob_newSTART);
  free(dat->MutProb_vec_new);
  for (i=0; i < 8; i++)
    free(dat->back_prob[i]);
  free(dat->back_prob);
  free(dat->ndonorhaps_vec);
}

double ** sampler(int * newh, int ** existing_h, int *p_Nloci, int *p_Nhaps, int *p_nchr,  double p_rhobar, double * MutProb_vec, int * allelic_type_count_vec, double * lambda, double * pos, int nsampTOT, double * copy_prob, double * copy_probSTART, int * pop_vec, int ndonorpops, double region_size, int run_num, int run_thres, int all_versus_all_ind, int haploid_ind, int unlinked_ind, int ind_val, int * hap_label_vec, int print_file9_ind, FILE *fout, FILE *fout3, FILE *fout9)
{
  double rounding_val = 1.0/10000000.0;  // for regional_counts; c is a bit lame
  //double small_missing_val=0.000000000000001;    // prob used when donor data is missing at a SNP (!!! NOT CURRENTLY IMPLEMENTED !!!)

  int i, j, locus;
  double sum, Theta;
  double prob, total_prob, total_gen_dist;
  double total_prob_from_i_to_i,total_prob_to_i_exclude_i,total_prob_from_i_exclude_i,constant_exclude_i,constant_from_i_to_i,constant_exclude_i_both_sides,total_prob_from_any_to_any_exclude_i;
  double total_regional_chunk_count, total_ind_sum;
  int num_regions;
  double * TransProb = malloc( ((*p_Nloci)-1) * sizeof(double));
  double N_e_new;
  int * sample_state = malloc(*p_Nloci * sizeof(int));
  double ObsStateProb, ObsStateProbPREV;
  double delta;    
             //correction to PAC-A rho_est
  double random_unif, random_unifSWITCH;
  double no_switch_prob;
  double ** Alphamat = malloc(*p_Nhaps * sizeof(double *));
  double * BetavecPREV = malloc(*p_Nhaps * sizeof(double));
  double * BetavecCURRENT = malloc(*p_Nhaps * sizeof(double));
  double Alphasum, Alphasumnew, Betasum, Betasumnew;
  double large_num;
  double * copy_prob_new = malloc(*p_Nhaps * sizeof(double));
  double * copy_prob_newSTART = malloc(*p_Nhaps * sizeof(double));
  double * Alphasumvec = malloc(*p_Nloci * sizeof(double));
  double * expected_transition_prob = malloc((*p_Nloci-1)*sizeof(double));
  double * corrected_chunk_count = malloc(*p_Nhaps * sizeof(double));
  double * regional_chunk_count = malloc(*p_Nhaps * sizeof(double));
  double * expected_chunk_length = malloc(*p_Nhaps * sizeof(double));
  double * expected_differences = malloc(*p_Nhaps * sizeof(double));
  double * regional_chunk_count_sum = malloc(ndonorpops * sizeof(double));
  double * regional_chunk_count_sum_final = malloc(ndonorpops * sizeof(double));
  double * regional_chunk_count_sum_squared_final = malloc(ndonorpops * sizeof(double));
  double * ind_snp_sum_vec = malloc(ndonorpops * sizeof(double));
  double * snp_info_measure=malloc(ndonorpops * sizeof(double));
  double * exp_copy_pop=malloc(ndonorpops * sizeof(double));
  double expected_chunk_length_sum, sum_prob;

  double ** copy_prob_new_mat = malloc(8 * sizeof(double *));
  for (i=0; i < 8; i++)
    copy_prob_new_mat[i] = malloc((*p_Nhaps+1) * sizeof(double));

  for(i=0 ; i< *p_Nhaps ; i++)
    {
      Alphamat[i] = malloc(*p_Nloci * sizeof(double));
    }
  for (i=0; i < ndonorpops; i++)
    {
      regional_chunk_count_sum[i] = 0.0;
      regional_chunk_count_sum_final[i] = 0.0;
      regional_chunk_count_sum_squared_final[i] = 0.0;
      ind_snp_sum_vec[i]=0.0;
      snp_info_measure[i]=0.0;
    }

				// Theta as given in Li and Stephens 
  sum = 0;
  for(i = 1; i < *p_nchr; i++){
    sum = sum + 1.0/i; 
  }
  Theta = 1.0 / sum;

  for (i=0; i < *p_Nhaps; i++) 
    {
      if (MutProb_vec[i]<0) MutProb_vec[i]=0.5 * Theta/(*p_Nhaps + Theta);
      //if (MutProb_vec[i]<rounding_val) MutProb_vec[i]=rounding_val;
    }

    // TransProb[i] is probability of copying mechanism "jumping" between
  //   loci i and i+1 
  delta = 1.0;

  if (unlinked_ind==0 && lambda[0] >= 0) TransProb[0] = 1 - exp(-1 * (pos[1]-pos[0]) * delta * p_rhobar*lambda[0]);
  if (unlinked_ind==1 || lambda[0] < 0) TransProb[0] = 1.0;
  /*
  if (TransProb[0]<rounding_val)
	{
	  printf("Transition prob is too low; will likely cause rounding errors. Exiting...\n");
	  exit(1);
	}
  */

  for(locus = 1; locus < *p_Nloci - 1; locus++)
    {
      delta = 1.0;
      if (unlinked_ind==0 && lambda[locus] >= 0) TransProb[locus] = 1 - exp(-1 * (pos[locus+1]-pos[locus]) * delta * p_rhobar*lambda[locus]);
      if (unlinked_ind==1 || lambda[locus] < 0) TransProb[locus] = 1.0;
      /*
      if (TransProb[locus]<rounding_val)
	{
	  printf("Transition prob is too low; will likely cause rounding errors. Exiting...\n");
	  exit(1);
	}
      */
    }
      /* FORWARDS ALGORITHM: (Rabiner 1989, p.262) */
      /* INITIALIZATION: */
  Alphasum = 0.0;
  for (i=0; i < *p_Nhaps; i++)
    {
      //if (newh[0]!=9 && existing_h[i][0] != -9) ObsStateProb = (1-MutProb_vec[i]) * (newh[0] == existing_h[i][0]) + MutProb_vec[i] * (newh[0] != existing_h[i][0]);
      if (newh[0]!=9) ObsStateProb = (1-MutProb_vec[i]) * (newh[0] == existing_h[i][0]) + MutProb_vec[i] * (newh[0] != existing_h[i][0]);
      if (newh[0]==9) ObsStateProb = 1.0;
      //if (existing_h[i][0]==9) ObsStateProb = small_missing_val;
      /*
      if (ObsStateProb<rounding_val)
	{
	  printf("Mutation (emission) rate is too low; will likely cause rounding errors. Exiting...\n");
	  exit(1);
	}
      */
      
      Alphamat[i][0] = log(copy_probSTART[i]*ObsStateProb);
      Alphasum = Alphasum + exp(Alphamat[i][0])*TransProb[0];
    }
     
       /* INDUCTION: */
  Alphasum = log(Alphasum);
  for (locus=1; locus < *p_Nloci; locus++)
    {
      Alphasumnew = 0.0;
      large_num = -1.0*Alphasum;
      for (i=0; i < *p_Nhaps; i++)
	{
	  if (newh[locus]!=9) ObsStateProb = (1-MutProb_vec[i]) * (newh[locus] == existing_h[i][locus]) + MutProb_vec[i] * (newh[locus] != existing_h[i][locus]);
	  if (newh[locus]==9) ObsStateProb = 1.0;
	  /*
	  if (ObsStateProb<rounding_val)
	    {
	      printf("Mutation (emission) rate is too low; will likely cause rounding errors. Exiting...\n");
	      exit(1);
	    }
	  */

	  Alphamat[i][locus] = log(ObsStateProb*copy_prob[i]*exp(Alphasum+large_num) + ObsStateProb*(1-TransProb[(locus-1)])*exp(Alphamat[i][(locus-1)]+large_num)) - large_num;
	  if (locus < (*p_Nloci - 1)) Alphasumnew = Alphasumnew + exp(Alphamat[i][locus]+large_num)*TransProb[locus];
	  if (locus == (*p_Nloci - 1)) Alphasumnew = Alphasumnew + exp(Alphamat[i][locus]+large_num);
	}
      Alphasum = log(Alphasumnew)-large_num;
    }
  //if (Alphasum == (Alphasum*5))
  if (isnan(Alphasum))
    {
      printf("Negative or NaN likelihood. Could be because emission or transition probabilities are too low??...Exiting...\n");
      //for (i=0; i < *p_Nhaps; i++) printf("%d %lf %lf %lf\n",i,copy_prob[i],log(copy_prob[i]),log(MutProb_vec[i]));
      exit(1);
    }
  fprintf(fout3," %.10lf",Alphasum);

  for (i = 0; i < *p_Nhaps; i++)
    {
      copy_prob_new[i] = 0.0;
      corrected_chunk_count[i] = 0.0;
      expected_chunk_length[i] = 0.0;
      expected_differences[i] = 0.0;
      regional_chunk_count[i] = 0.0;
    }
  total_regional_chunk_count=0.0;
  num_regions=0;
  if (run_num <= (run_thres-1))
    {
          /* BACKWARDS ALGORITHM: (Rabiner 1989, p.263) */
          /* INITIALIZATION: */
      Betasum = 0.0;
      if (run_num == (run_thres-1)) 
	{
	  for (i=0; i < ndonorpops; i++)
	    exp_copy_pop[i]=0.0;
	}
      for(i=0; i < *p_Nhaps; i++)
	{
	  if (newh[(*p_Nloci-1)]!=9) ObsStateProb = (1-MutProb_vec[i]) * (newh[(*p_Nloci-1)] == existing_h[i][(*p_Nloci-1)]) + MutProb_vec[i] * (newh[(*p_Nloci-1)] != existing_h[i][(*p_Nloci-1)]);
	  if (newh[(*p_Nloci-1)]==9) ObsStateProb = 1.0;
	  /*
	  if (ObsStateProb<rounding_val)
	    {
	      printf("Mutation (emission) rate is too low; will likely cause rounding errors. Exiting...\n");
	      exit(1);
	    }
	  */
	  BetavecPREV[i] = 0.0;
	  Betasum = Betasum + TransProb[(*p_Nloci-2)]*copy_prob[i]*ObsStateProb*exp(BetavecPREV[i]);
	  if (run_num == (run_thres-1)) exp_copy_pop[pop_vec[i]]=exp_copy_pop[pop_vec[i]]+exp(BetavecPREV[i]+Alphamat[i][(*p_Nloci-1)]-Alphasum);
	                    // for estimating new mutation rates:
	  expected_differences[i]=expected_differences[i]+exp(Alphamat[i][(*p_Nloci-1)]-Alphasum)*(newh[(*p_Nloci-1)] != existing_h[i][(*p_Nloci-1)]);
	}
      if ((run_num == (run_thres-1)) && (print_file9_ind==1)) 
 	{
	  gzprintf(fout9,"%.0lf",pos[(*p_Nloci-1)]);
	  if (all_versus_all_ind==0)
	    {
	      for (i=0; i < ndonorpops; i++)
		gzprintf(fout9," %lf",exp_copy_pop[i]);
	    }
	  if (all_versus_all_ind==1)
	    {
	      for (i=0; i < ndonorpops; i++)
		{
		  if (i == ind_val) gzprintf(fout9," 0.0");
		  gzprintf(fout9," %lf",exp_copy_pop[i]);
		}
	      if (ind_val == ndonorpops) gzprintf(fout9," 0.0");
	    }
	  gzprintf(fout9,"\n");
	} 

           /* INDUCTION: */
      Betasum = log(Betasum);
          /* CALCULATE EXPECTED NUMBER OF TIMES OF COPYING TO EACH DONOR POP (Rabiner 1989, p.263,265 or Scheet/Stephens 2006 Appendix C): */      
      for (locus = (*p_Nloci-2); locus >= 0; locus--)
	{
	  Betasumnew = 0.0;      
	  large_num = -1.0*Betasum;
	  total_prob=0.0;

	  constant_exclude_i = 0.5;
	  constant_from_i_to_i = 1.0;
	  constant_exclude_i_both_sides = 0.0;
	  expected_chunk_length_sum=0.0;
	  sum_prob=0.0;
	  if (run_num == (run_thres-1)) 
	    {
	      for (i=0; i < ndonorpops; i++)
		exp_copy_pop[i]=0.0;
	    }
 	  for (i = 0; i < *p_Nhaps; i++)
	    {
	      if (newh[locus]!=9) ObsStateProb = (1-MutProb_vec[i]) * (newh[locus] == existing_h[i][locus]) + MutProb_vec[i] * (newh[locus] != existing_h[i][locus]);
	      if (newh[locus]==9) ObsStateProb = 1.0;
	      if (newh[(locus+1)]!=9) ObsStateProbPREV = (1-MutProb_vec[i]) * (newh[(locus+1)] == existing_h[i][(locus+1)]) + MutProb_vec[i] * (newh[(locus+1)] != existing_h[i][(locus+1)]);
	      if (newh[(locus+1)]==9) ObsStateProbPREV = 1.0;
	      /*
	      if ((ObsStateProb<rounding_val) || (ObsStateProbPREV<rounding_val))
		{
		  printf("Mutation (emission) rate is too low; will likely cause rounding errors. Exiting...\n");
		  exit(1);
		}
	      */
	      BetavecCURRENT[i] = log(exp(Betasum+large_num) + (1-TransProb[locus]) * ObsStateProbPREV*exp(BetavecPREV[i] + large_num)) - large_num;
	      if (locus > 0) Betasumnew = Betasumnew + TransProb[(locus-1)]*copy_prob[i]*ObsStateProb*exp(BetavecCURRENT[i] + large_num);
	      if (locus == 0) copy_prob_newSTART[i] = exp(Alphamat[i][0] + BetavecCURRENT[i] - Alphasum);
	      total_prob = total_prob + exp(Alphamat[i][(locus+1)]+BetavecPREV[i]-Alphasum)-exp(Alphamat[i][locus]+BetavecPREV[i]- Alphasum)*ObsStateProbPREV*(1-TransProb[locus]);

	      copy_prob_new[i] = copy_prob_new[i] + exp(Alphamat[i][(locus+1)]+BetavecPREV[i]-Alphasum)-exp(Alphamat[i][locus]+BetavecPREV[i]- Alphasum)*ObsStateProbPREV*(1-TransProb[locus]);

	      total_prob_from_i_to_i = exp(Alphamat[i][locus]+BetavecPREV[i]-Alphasum)*ObsStateProbPREV*(1-TransProb[locus]+TransProb[locus]*copy_prob[i]);
	      total_prob_to_i_exclude_i = exp(Alphamat[i][(locus+1)]+BetavecPREV[i]-Alphasum)-exp(Alphamat[i][locus]+BetavecPREV[i]- Alphasum)*ObsStateProbPREV*(1-TransProb[locus]+TransProb[locus]*copy_prob[i]);
	      total_prob_from_i_exclude_i = exp(Alphamat[i][locus]+BetavecCURRENT[i]-Alphasum) - exp(Alphamat[i][locus]+BetavecPREV[i]-Alphasum)*ObsStateProbPREV*(1-TransProb[locus]+TransProb[locus]*copy_prob[i]);
	      total_prob_from_any_to_any_exclude_i = 1.0-exp(Alphamat[i][locus]+BetavecCURRENT[i]-Alphasum)-exp(Alphamat[i][(locus+1)]+BetavecPREV[i]-Alphasum)+exp(Alphamat[i][locus]+BetavecPREV[i]-Alphasum)*ObsStateProbPREV*(1-TransProb[locus]+TransProb[locus]*copy_prob[i]);

	      regional_chunk_count[i]=regional_chunk_count[i]+(exp(Alphamat[i][(locus+1)]+BetavecPREV[i]-Alphasum)-exp(Alphamat[i][locus]+BetavecPREV[i]- Alphasum)*ObsStateProbPREV*(1-TransProb[locus]));
	      total_regional_chunk_count=total_regional_chunk_count+(exp(Alphamat[i][(locus+1)]+BetavecPREV[i]-Alphasum)-exp(Alphamat[i][locus]+BetavecPREV[i]- Alphasum)*ObsStateProbPREV*(1-TransProb[locus]));
	      ind_snp_sum_vec[pop_vec[i]]=ind_snp_sum_vec[pop_vec[i]]+(exp(Alphamat[i][(locus+1)]+BetavecPREV[i]-Alphasum)-exp(Alphamat[i][locus]+BetavecPREV[i]- Alphasum)*ObsStateProbPREV*(1-TransProb[locus]));

	      //corrected_chunk_count[i]=corrected_chunk_count[i]+(exp(Alphamat[i][(locus+1)]+BetavecPREV[i]-Alphasum)-exp(Alphamat[i][locus]+BetavecPREV[i]- Alphasum)*ObsStateProbPREV*(1-TransProb[locus]))*(1.0+(1.0/(*p_Nhaps))*((p_rhobar*(pos[(locus+1)]-pos[locus])*delta*lambda[locus]/(*p_Nhaps))/(1.0-exp(-1.0*p_rhobar*(pos[(locus+1)]-pos[locus])*delta*lambda[locus]/(*p_Nhaps)))-1.0));
	      corrected_chunk_count[i]=corrected_chunk_count[i]+(exp(Alphamat[i][(locus+1)]+BetavecPREV[i]-Alphasum)-exp(Alphamat[i][locus]+BetavecPREV[i]- Alphasum)*ObsStateProbPREV*(1-TransProb[locus]));
	      if (unlinked_ind==0 && lambda[locus]>=0) expected_chunk_length[i]=expected_chunk_length[i]+100*(pos[locus+1]-pos[locus])*delta*lambda[locus]*(constant_from_i_to_i*total_prob_from_i_to_i+constant_exclude_i*(total_prob_to_i_exclude_i+total_prob_from_i_exclude_i)+constant_exclude_i_both_sides*total_prob_from_any_to_any_exclude_i);  // multiply by 100 to get cM
	      expected_chunk_length_sum=expected_chunk_length_sum+constant_from_i_to_i*total_prob_from_i_to_i+constant_exclude_i*(total_prob_to_i_exclude_i+total_prob_from_i_exclude_i)+constant_exclude_i_both_sides*total_prob_from_any_to_any_exclude_i;
	                    // for estimating new mutation rates:
	      expected_differences[i]=expected_differences[i]+exp(Alphamat[i][locus]+BetavecCURRENT[i]-Alphasum)*(newh[locus] != existing_h[i][locus]);
	      BetavecPREV[i] = BetavecCURRENT[i];

	      if (run_num == (run_thres-1)) exp_copy_pop[pop_vec[i]]=exp_copy_pop[pop_vec[i]]+exp(BetavecCURRENT[i]+Alphamat[i][locus]-Alphasum);

	      sum_prob=sum_prob+total_prob_from_i_to_i+total_prob_to_i_exclude_i+total_prob_from_i_exclude_i;
	    }
	  
	  if ((run_num == (run_thres-1)) && (print_file9_ind==1)) 
	    {
	      gzprintf(fout9,"%.0lf",pos[locus]);
	      if (all_versus_all_ind==0)
		{
		  for (i=0; i < ndonorpops; i++)
		    gzprintf(fout9," %lf",exp_copy_pop[i]);
		}
	      if (all_versus_all_ind==1)
		{
		  for (i=0; i < ndonorpops; i++)
		    {
		      if (i == ind_val) gzprintf(fout9," 0.0");
		      gzprintf(fout9," %lf",exp_copy_pop[i]);
		    }
		  if (ind_val == ndonorpops) gzprintf(fout9," 0.0");
		}
	      gzprintf(fout9,"\n");
	    }
	  

	  expected_transition_prob[locus]=total_prob;
	  if (locus > 0) Betasum = log(Betasumnew) - large_num;

	  if ((total_regional_chunk_count+rounding_val) >= region_size)
	    {
	      for (i = 0; i < *p_Nhaps; i++)
		{
		  regional_chunk_count_sum[pop_vec[i]]=regional_chunk_count_sum[pop_vec[i]]+regional_chunk_count[i];
		  regional_chunk_count[i]=0.0;
		}
	      for (i = 0; i < ndonorpops; i++)
		{
		  regional_chunk_count_sum_final[i]=regional_chunk_count_sum_final[i]+regional_chunk_count_sum[i];
		  regional_chunk_count_sum_squared_final[i]=regional_chunk_count_sum_squared_final[i]+pow(regional_chunk_count_sum[i],2.0);
		  regional_chunk_count_sum[i]=0.0;
		}
	      total_regional_chunk_count=0.0;
	      num_regions=num_regions+1;
	    }
	  total_ind_sum=0.0;
	  for (i = 0; i < ndonorpops; i++)
	    total_ind_sum=total_ind_sum+ind_snp_sum_vec[i];
	  for (i = 0; i < ndonorpops; i++)
	    {
	      snp_info_measure[i]=snp_info_measure[i]+pow((ind_snp_sum_vec[i]/total_ind_sum),2.0);
	      ind_snp_sum_vec[i]=0.0;
	    }
	}     
  for (i=0; i < ndonorpops; i++)
    snp_info_measure[i]=snp_info_measure[i]/(*p_Nloci);

          /* CALCULATE EXPECTED NUMBER OF TOTAL TRANSITIONS, IN ORDER TO ESTIMATE N_e (Scheet/Stephens 2006 Appendix C (C3)): */      
      total_prob=0.0;
      total_gen_dist=0.0;
      for (locus = 0; locus < (*p_Nloci-1); locus++)
	{
	  if (unlinked_ind==0 && lambda[locus] >= 0) total_gen_dist=total_gen_dist+(pos[(locus+1)]-pos[locus])*delta*lambda[locus];
	  if (unlinked_ind==0 && lambda[locus] >= 0) total_prob=total_prob+((p_rhobar*(pos[(locus+1)]-pos[locus])*delta*lambda[locus])/(1.0-exp(-1.0*p_rhobar*(pos[(locus+1)]-pos[locus])*delta*lambda[locus])))*expected_transition_prob[locus];
	}
      if (unlinked_ind==0) N_e_new = total_prob/total_gen_dist;
      if (unlinked_ind==1) N_e_new = 0.0;

              /* CALCULATE SOMETHING ANALAGOUS TO EXPECTED NUMBER OF TIMES EACH HAP i IS VISITED, CONDITIONAL ON THE OBSERVED DATA (I.E  (27) AND PARAGRAPH UNDER (38) IN RABINER 1989, Proceedings of the IEEE 77(2):257-286), BUT -- AS WE'RE ONLY COUNTING CHUNKS -- SUBTRACT OUT TIMES YOU DO NOT SWITCH */
      for (i=0; i < *p_Nhaps; i++) corrected_chunk_count[i]=corrected_chunk_count[i]+copy_prob_newSTART[i];
    }

       /* print-out samples if we've done enough iterations: */
   if (run_num == (run_thres-1))
     {

       N_e_new = p_rhobar;

       for (i=0; i < *p_Nhaps; i++)
	 {
	   copy_prob_new[i] = copy_prob[i];
	   copy_prob_newSTART[i] = copy_probSTART[i];
	 }

           /* calculate Alphasums (for efficient sampling): */
       for (locus=0; locus < *p_Nloci; locus++)
	 {
	   Alphasumvec[locus] = 0.0;
	   large_num = Alphamat[0][locus];
	   for (i = 1; i < *p_Nhaps; i++)
	     {
	       if (Alphamat[i][locus] > large_num)
		 large_num = Alphamat[i][locus];
	     }
	   large_num = -1.0*large_num;
	   for (i = 0; i < *p_Nhaps; i++)
	     Alphasumvec[locus] = Alphasumvec[locus] + exp(Alphamat[i][locus]+large_num);
	   Alphasumvec[locus] = log(Alphasumvec[locus]) - large_num;
	 }

              /* SAMPLING ALGORITHM: (from Falush, Stephens, & Pritchard (2003) Genetics 164:1567-1587) */
       for (j = 0; j < nsampTOT; j++)
	 {
	      /* sample last position: */
	   total_prob = 0.0;
	   large_num = Alphamat[0][(*p_Nloci-1)];
	   for (i = 1; i < *p_Nhaps; i++)
	     {
	       if (Alphamat[i][(*p_Nloci-1)] > large_num)
		 large_num = Alphamat[i][(*p_Nloci-1)];
	     }
	   large_num = -1.0*large_num;
	   random_unif = (double) rand()/RAND_MAX;
	   total_prob = Alphasumvec[(*p_Nloci-1)];
	   prob = 0.0;
	   for (i = 0; i < *p_Nhaps; i++)
	     {
	       prob = prob + exp(Alphamat[i][(*p_Nloci-1)]+large_num);
	       if (random_unif <= exp(log(prob)-large_num-total_prob))
		 {
		   sample_state[(*p_Nloci-1)] = i;
		   break;
		 }
	     }

              /* sample remaining positions: */
	   for (locus = (*p_Nloci-2); locus >= 0; locus--)
	     {
	        // first sample prob you switch and see if you need to
                   // if you do need to switch, you need to go through the below loop to figure out where to switch to
               large_num = -1.0 * Alphasumvec[locus];
	       total_prob = log(exp(Alphasumvec[locus]+large_num)*TransProb[locus]*copy_prob[sample_state[(locus+1)]] + exp(Alphamat[sample_state[(locus+1)]][locus]+large_num)*(1.0-TransProb[locus]))-large_num;
	       no_switch_prob = exp(log(exp(Alphamat[sample_state[(locus+1)]][locus]+large_num)*(1.0-TransProb[locus])) - large_num - total_prob);
               random_unifSWITCH = (double) rand()/RAND_MAX;
	       if (random_unifSWITCH <= no_switch_prob) sample_state[locus] = sample_state[(locus+1)];

               if (random_unifSWITCH > no_switch_prob) 
		 {
		   total_prob = 0.0;
		   large_num = Alphamat[0][locus];
		   for (i = 1; i < *p_Nhaps; i++)
		     {
		       if (Alphamat[i][locus] > large_num)
			 large_num = Alphamat[i][locus];
		     }
		   large_num = -1.0*large_num;

		   random_unif = (double) rand()/RAND_MAX;
		   total_prob = log(exp(Alphasumvec[locus]+large_num)*TransProb[locus]*copy_prob[sample_state[(locus+1)]]) - large_num;
		   prob = 0.0;
		   for (i = 0; i < *p_Nhaps; i++)
		     {
		       prob = prob + exp(Alphamat[i][locus]+large_num)*TransProb[locus]*copy_prob[sample_state[(locus+1)]];
		       if (random_unif <= exp(log(prob)-large_num-total_prob))
			 {
			   sample_state[locus] = i;
			   break;
			 }
		     }
		 }
	     }
	   
	   fprintf(fout,"%d",j+1);
	   for (i = 0; i < *p_Nloci; i++)
	     fprintf(fout," %d",hap_label_vec[sample_state[i]]);
	   fprintf(fout,"\n");
	 }
     }

  for (i = 0; i < *p_Nhaps; i++)
     copy_prob_new_mat[0][i] = copy_prob_new[i];
   for (i = 0; i < *p_Nhaps; i++)
     copy_prob_new_mat[1][i] = copy_prob_newSTART[i];
   for (i = 0; i < *p_Nhaps; i++)
     copy_prob_new_mat[2][i] = corrected_chunk_count[i];
   for (i = 0; i < *p_Nhaps; i++)
     copy_prob_new_mat[3][i] = expected_chunk_length[i];
   for (i = 0; i < *p_Nhaps; i++)
     copy_prob_new_mat[4][i] = expected_differences[i];
   for (i = 0; i < ndonorpops; i++)
     copy_prob_new_mat[5][i] = regional_chunk_count_sum_final[i];
   for (i = 0; i < ndonorpops; i++)
     copy_prob_new_mat[6][i] = regional_chunk_count_sum_squared_final[i];
   for (i = 0; i < ndonorpops; i++)
     copy_prob_new_mat[7][i] = snp_info_measure[i];
   copy_prob_new_mat[0][(*p_Nhaps)] = N_e_new;
   copy_prob_new_mat[1][(*p_Nhaps)] = num_regions;

   for (i=0; i < *p_Nhaps; i++)
     {
       free(Alphamat[i]);
     }
   free(Alphamat);
   free(BetavecPREV);
   free(BetavecCURRENT);
   free(TransProb);
   free(sample_state);
   free(expected_transition_prob);
   free(Alphasumvec);
   free(copy_prob_new);
   free(copy_prob_newSTART);
   free(corrected_chunk_count);
   free(regional_chunk_count);
   free(expected_chunk_length);
   free(expected_differences);
   free(regional_chunk_count_sum);
   free(regional_chunk_count_sum_final);
   free(regional_chunk_count_sum_squared_final);
   free(ind_snp_sum_vec);
   free(snp_info_measure);
   free(exp_copy_pop);

   return(copy_prob_new_mat);
}

int loglik(int nind_tot, int nhaps_startpop, int *p_nloci, int p_nhaps, double N_e_start, double * recom_map, double * MutProb_vec_tot, int nsampTOT, int ndonorpops, int * ndonorhaps_tot, int * include_ind_vec, char ** ind_label_vec, char ** pop_label_vec, char ** donor_pop_vec, int nrecpops, char **rec_pop_vec, double * copy_prob_tot, double * copy_probSTART_tot, int * pop_vec_tot, double region_size, int EMruns, int estimate_copyprob_ind, int estimate_recom_ind, int ne_find, int estimate_mutation_ind, int estimate_mutationALL_ind, int all_versus_all_ind, int prior_donor_probs_ind, int num_rec_ind, int * recipient_ind_vec, char *filename, int donorlist_ind, int haploid_ind, int unlinked_ind, int print_file9_ind, int drift_calc_ind, int num_rec_drift, int rec_ind_topaint_count, FILE *fout, FILE *fout2, FILE *fout3, FILE *fout4, FILE *fout5, FILE *fout6, FILE *fout7, FILE *fout8, FILE *fout9)
{
  double small_copy_val=0.000000000000001; // (!!!) copy props per hap not allowed to go below this value, even if E-M wants to make them lower (!!!)
  int nhaps_condpop, nind_condpop;

  FILE *fd;
  struct data_t *Data;
  int i, j, m, n, r, h;
  int nhaps, num_regions_tot;
  int included_count, included_count_bigloop;
  double sum_total_diff;
  double * total_back_prob = malloc(ndonorpops * sizeof(double));
  double * total_back_probSTART = malloc(ndonorpops * sizeof(double));
  double * total_counts = malloc(ndonorpops * sizeof(double));
  double * total_lengths = malloc(ndonorpops * sizeof(double));
  double * total_differences = malloc(ndonorpops * sizeof(double));
  double * total_region_counts = malloc(ndonorpops * sizeof(double));
  double * total_squared_region_counts = malloc(ndonorpops * sizeof(double));
  double * snp_info_measure_final = malloc(ndonorpops * sizeof(double));
  double N_e_new, N_e;
  double total_prob, total_probSTART;
  double ** copy_prob_pop = malloc(2 * sizeof(double *));
  int * newhap = malloc((*p_nloci) * sizeof(int));
  int drift_malloc_size, drift_malloc_size2, drift_malloc_size3;
  drift_malloc_size=num_rec_drift;
  drift_malloc_size2=num_rec_ind;
  drift_malloc_size3=*p_nloci;
  if (drift_malloc_size==0)
    { 
      drift_malloc_size=1;
      drift_malloc_size2=1;
      drift_malloc_size3=1;
    }
  double * correlated_drift_vec = malloc(drift_malloc_size * sizeof(double));
  nhaps_condpop = p_nhaps-nhaps_startpop;
  nind_condpop=nhaps_condpop/(2-haploid_ind);
  double ** correlated_drift_calc = malloc(drift_malloc_size2 * sizeof(double *));
  double *** RecChromProbArray = malloc(drift_malloc_size * sizeof(double **)); 
  for (i=0; i<drift_malloc_size2; i++)
    correlated_drift_calc[i]=malloc(nind_condpop*sizeof(double));
  for (i=0; i < drift_malloc_size; i++)
    {
      RecChromProbArray[i]=malloc(ndonorpops*sizeof(double *));
      for (j=0; j < ndonorpops; j++)
	RecChromProbArray[i][j]=malloc(drift_malloc_size3 * sizeof(double));
    }
  for (i=0; i < 2; i++)
    copy_prob_pop[i] = malloc(ndonorpops * sizeof(double));
  for (i=0; i < ndonorpops; i++)
    {
      total_back_prob[i] = 0.0;
      total_back_probSTART[i] = 0.0;
    }

  total_prob = 0.0;
  total_probSTART = 0.0;
  int * allelic_type_count_vec = malloc(*p_nloci*sizeof(int));
  int * found_vec = malloc(6*sizeof(int));
  included_count_bigloop=0;
  for (m = 0; m < nind_tot; m ++)
    {
      if (recipient_ind_vec[m]==1)
	{
	  printf(" .....Painting recipient individual %d of %d......\n",included_count_bigloop+1,rec_ind_topaint_count);
	  fprintf(fout3,"%s\n",ind_label_vec[m]);

	  fd = fopen(filename,"r");
	  if (fd == NULL) { printf("error opening %s\n",filename); exit(1);}
	  Data=ReadData(fd,m,2-haploid_ind,include_ind_vec,pop_label_vec,ndonorpops,donor_pop_vec,pop_vec_tot,copy_prob_tot,copy_probSTART_tot,MutProb_vec_tot,ndonorhaps_tot,all_versus_all_ind);

	  for (i=0; i < Data->ndonorhaps; i++)
	    {
	      Data->copy_prob_new[i] = Data->copy_prob[i];
	      Data->copy_prob_newSTART[i] = Data->copy_probSTART[i];
	      if (prior_donor_probs_ind==0)
		{
		  Data->copy_prob_new[i] = Data->copy_prob[i]/Data->ndonorhaps;
		  Data->copy_prob_newSTART[i] = Data->copy_probSTART[i]/Data->ndonorhaps;
		}
	      Data->MutProb_vec_new[i] = Data->MutProb_vec[i];
	    }

                  // find number of alleles per snp (this is NOT every used, but perhaps should be to get default mutation rate correct):
	  for (j=0; j < Data->nsnps; j++)
	    {
	      allelic_type_count_vec[j] = 0;
	      for (i=0; i < 6; i++)
		found_vec[i]=0;
	      for (i=0; i < Data->ndonorhaps; i++)
		{
		  if ((Data->cond_chromosomes[i][j] == 0) && (found_vec[0] == 0))
		    {
		      allelic_type_count_vec[j] = allelic_type_count_vec[j] + 1;
		      found_vec[0] = 1;
		    }
		  if ((Data->cond_chromosomes[i][j] == 1) && (found_vec[1] == 0))
		    {
		      allelic_type_count_vec[j] = allelic_type_count_vec[j] + 1;
		      found_vec[1] = 1;
		    }
		  if ((Data->cond_chromosomes[i][j] == 2) && (found_vec[2] == 0))
		    {
		      allelic_type_count_vec[j] = allelic_type_count_vec[j] + 1;
		      found_vec[2] = 1;
		    }
		  if ((Data->cond_chromosomes[i][j] == 3) && (found_vec[3] == 0))
		    {
		      allelic_type_count_vec[j] = allelic_type_count_vec[j] + 1;
		      found_vec[3] = 1;
		    }
		  if ((Data->cond_chromosomes[i][j] == 4) && (found_vec[4] == 0))
		    {
		      allelic_type_count_vec[j] = allelic_type_count_vec[j] + 1;
		      found_vec[4] = 1;
		    }
		  if ((Data->cond_chromosomes[i][j] == 5) && (found_vec[5] == 0))
		    {
		      allelic_type_count_vec[j] = allelic_type_count_vec[j] + 1;
		      found_vec[5] = 1;
		    }
		}
	    }

	  for (r=0; r < EMruns; r++)
	    {
	      fprintf(fout3,"%d",r);

	      total_prob = 0.0;
	      total_probSTART = 0.0;
	      for (i=0; i < ndonorpops; i++)
		{
		  total_back_prob[i] = 0.0;
		  total_back_probSTART[i] = 0.0;
		  total_counts[i]=0.0;
		  total_lengths[i]=0.0;
		  total_differences[i]=0.0;
		  total_region_counts[i]=0.0;
		  total_squared_region_counts[i]=0.0;
		  snp_info_measure_final[i]=0.0;
		}
	      
	      N_e_new=0.0;
	      num_regions_tot=0;
	      for(h=0; h < (2-haploid_ind); h++)
		{
		  if (r==0)
		    {
		      if (ne_find==0) N_e=N_e_start/Data->ndonorhaps;
		      if (ne_find==1) N_e=N_e_start;
		    }
		  
		  for (n = 0; n < Data->nsnps; n ++)
		    newhap[n] = Data->ind_chromosomes[h][n];

		  nhaps = Data->ndonorhaps;

		  if (r == (EMruns-1))
		    {
		      fprintf(fout,"HAP %d %s\n",h+1,ind_label_vec[m]);
		      if (print_file9_ind==1) gzprintf(fout9,"HAP %d %s\n",h+1,ind_label_vec[m]);
		    }

                         /* SAMPLE FROM PAC CONDITIONAL ON COPY-PROBS: */
		  Data->back_prob = sampler(newhap, Data->cond_chromosomes, p_nloci, &nhaps, &p_nhaps, N_e, Data->MutProb_vec_new, allelic_type_count_vec, recom_map, Data->positions, nsampTOT, Data->copy_prob_new, Data->copy_prob_newSTART, Data->pop_vec, ndonorpops, region_size, r, EMruns, all_versus_all_ind, haploid_ind, unlinked_ind, included_count_bigloop, Data->hap_label_vec, print_file9_ind, fout, fout3, fout9);

		  N_e_new=N_e_new+Data->back_prob[0][Data->ndonorhaps];
		  num_regions_tot=num_regions_tot+Data->back_prob[1][Data->ndonorhaps];

                      /* GET NEW COPY-PROBS BASED ON PAC SAMPLES: */
		  for (i = 0; i < Data->ndonorhaps; i++)
		    {
		      total_back_prob[Data->pop_vec[i]] = total_back_prob[Data->pop_vec[i]] + Data->back_prob[0][i];
		      total_prob = total_prob + Data->back_prob[0][i];
		      total_back_probSTART[Data->pop_vec[i]] = total_back_probSTART[Data->pop_vec[i]] + Data->back_prob[1][i];
		      total_probSTART = total_probSTART + Data->back_prob[1][i];
		      total_counts[Data->pop_vec[i]] = total_counts[Data->pop_vec[i]] + Data->back_prob[2][i];
		      total_lengths[Data->pop_vec[i]] = total_lengths[Data->pop_vec[i]] + Data->back_prob[3][i];
		      total_differences[Data->pop_vec[i]] = total_differences[Data->pop_vec[i]] + Data->back_prob[4][i];
		    }
		  for (i = 0; i < ndonorpops; i++)
		    {
		      total_region_counts[i] = total_region_counts[i] + Data->back_prob[5][i];
		      total_squared_region_counts[i] = total_squared_region_counts[i] + Data->back_prob[6][i];
		      snp_info_measure_final[i] = snp_info_measure_final[i] + Data->back_prob[7][i];
		    }
		}
	      fprintf(fout3," %.10lf %.10lf\n",N_e,Data->MutProb_vec_new[0]);
	      if (estimate_recom_ind==1) N_e=N_e_new/(2.0-haploid_ind);
	      for (i = 0; i < ndonorpops; i++)
		{
		  copy_prob_pop[0][i] = total_back_prob[i]/total_prob;
		  copy_prob_pop[1][i] = total_back_probSTART[i]/total_probSTART;
		}
                       /* RESET COPY-PROBS and MUTATION-PROBS: */
	            // (first check for probabilities of 0:)
	      for (i=0; i < ndonorpops; i++)
		{
		  if (copy_prob_pop[0][i] <= 0)
		    copy_prob_pop[0][i] = small_copy_val*Data->ndonorhaps_vec[i];
		  if (copy_prob_pop[1][i] <= 0)
		    copy_prob_pop[1][i] = small_copy_val*Data->ndonorhaps_vec[i];
		}
	      total_prob = 0.0;
	      total_probSTART = 0.0;
	      for (j=0; j < ndonorpops; j++)
		{
		  total_prob = total_prob + copy_prob_pop[0][j];
		  total_probSTART = total_probSTART + copy_prob_pop[1][j];
		}
	      for (j=0; j < ndonorpops; j++)
		{
		  copy_prob_pop[0][j] = copy_prob_pop[0][j]/total_prob;
		  copy_prob_pop[1][j] = copy_prob_pop[1][j]/total_probSTART;
		}
	      if (estimate_copyprob_ind==1)
		{
		  for (i = 0; i < Data->ndonorhaps; i++)
		    {
		      for (j=0; j < ndonorpops; j++)
			{
			  if (Data->pop_vec[i]==j)
			    {
			      if (Data->ndonorhaps_vec[j]>0)
				{
				  Data->copy_prob_new[i]=copy_prob_pop[0][j]/Data->ndonorhaps_vec[j];
				  Data->copy_prob_newSTART[i]=copy_prob_pop[1][j]/Data->ndonorhaps_vec[j];
				}
			      if (Data->ndonorhaps_vec[j]==0)
				{
				  Data->copy_prob_new[i]=0.0;
				  Data->copy_prob_newSTART[i]=0.0;
				}
			      break;
			    }
			}
		    }
		}
	      if (estimate_mutation_ind==1)
		{
		  for (i = 0; i < Data->ndonorhaps; i++)
		    {
		      for (j=0; j < ndonorpops; j++)
			{
			  if (Data->pop_vec[i]==j)
			    {
			      Data->MutProb_vec_new[i]=total_differences[j]/(*p_nloci*(2-haploid_ind));
			      break;
			    }
			}
		    }
		}
	      
	      if (estimate_mutationALL_ind==1)
		{
		  sum_total_diff=0.0;
		  for (i=0; i < ndonorpops; i++)
		    sum_total_diff=sum_total_diff+total_differences[i]/(*p_nloci*(2-haploid_ind));
		  for (i = 0; i < Data->ndonorhaps; i++)
		    Data->MutProb_vec_new[i] = sum_total_diff;
		}
	      
	      if (r == (EMruns-1))
		{
		     /* print props, lengths, counts, and differences: */
		  fprintf(fout2,"%s",ind_label_vec[m]);
		  fprintf(fout4,"%s",ind_label_vec[m]);
		  fprintf(fout5,"%s",ind_label_vec[m]);
		  fprintf(fout6,"%s",ind_label_vec[m]);
		  fprintf(fout7,"%s",ind_label_vec[m]);
		  fprintf(fout8,"%s",ind_label_vec[m]);
		  fprintf(fout7," %d",num_regions_tot);
		  fprintf(fout8," %d",num_regions_tot);
		  if (all_versus_all_ind==0)
		    {
		      for (j=0; j < ndonorpops; j++)
			{
			  fprintf(fout2," %lf",copy_prob_pop[0][j]);
			  fprintf(fout4," %lf",total_counts[j]);
			  fprintf(fout5," %lf",total_lengths[j]);
			  fprintf(fout6," %lf",total_differences[j]);
			  fprintf(fout7," %lf",total_region_counts[j]);
			  fprintf(fout8," %lf",total_squared_region_counts[j]);
			}
		    }
		  if (all_versus_all_ind==1)
		    {
		      included_count=0;
		      for (j=0; j < nind_tot; j++)
			{
			  if (j==m)
			    {
			      fprintf(fout2," 0.00");
			      fprintf(fout4," 0.00");
			      fprintf(fout5," 0.00");
			      fprintf(fout6," 0.00");
			      fprintf(fout7," 0.00");
			      fprintf(fout8," 0.00");
			    }
			  if (include_ind_vec[j]==1 && j!=m)
			    {
			      fprintf(fout2," %lf",copy_prob_pop[0][included_count]);
			      fprintf(fout4," %lf",total_counts[included_count]);
			      fprintf(fout5," %lf",total_lengths[included_count]);
			      fprintf(fout6," %lf",total_differences[included_count]);
			      fprintf(fout7," %lf",total_region_counts[included_count]);
			      fprintf(fout8," %lf",total_squared_region_counts[included_count]);
			      included_count=included_count+1;
			    }
			}
		    }
		  fprintf(fout2,"\n");
		  fprintf(fout4,"\n");
		  fprintf(fout5,"\n");
		  fprintf(fout6,"\n");
		  fprintf(fout7,"\n");
		  fprintf(fout8,"\n");
		}
	    }
	  DestroyData(Data,2-haploid_ind);
	  fclose(fd);
	}
      if (include_ind_vec[m]==1 && recipient_ind_vec[m]==1)
	included_count_bigloop=included_count_bigloop+1;
    }

  free(newhap);
  for (i=0; i < drift_malloc_size2; i++)
    free(correlated_drift_calc[i]);
  for (i=0; i < drift_malloc_size; i++)
    {
      for (j=0; j < ndonorpops; j++)
	free(RecChromProbArray[i][j]);
      free(RecChromProbArray[i]);
    }
  free(RecChromProbArray);
  free(total_back_prob);
  free(total_back_probSTART);
  free(total_counts);
  free(total_lengths);
  free(total_differences);
  free(total_region_counts);
  free(total_squared_region_counts);
  free(snp_info_measure_final);
  free(allelic_type_count_vec);
  free(found_vec);
  free(correlated_drift_vec);
  free(correlated_drift_calc);

  return(1);
}

void usage()
{
  printf("to run: use './ChromoPainterv2' with following options:\n");
  printf ("       -g <geno.filein>  (REQUIRED; no default)\n");
  printf ("       -r <recommap.filein>  (REQUIRED; no default -- unless using -u switch)\n");
  printf ("       -t <labels.filein> file listing id and population labels for each individual (REQUIRED; no default -- unless using -a switch)\n");
  printf ("       -f <poplist.filein> <f_1> <f_2>  file listing donor and recipient populations (REQUIRED; no default -- unless using -a switch); paint recipient individuals f_1 through f_2 using all donor population haplotypes (use '-f <donorlist.filein> 0 0' to paint all recipient inds)\n");
  printf ("       -i <int>  number of EM iterations for estimating parameters (default=0)\n");
  printf ("       -in  maximize over average switch rate parameter using E-M\n");
  printf ("       -ip  maximize over copying proportions using E-M\n");
  printf ("       -im  maximize over donor population mutation (emission) probabilities using E-M\n");
  printf ("       -iM  maximize over global mutation (emission) probability using E-M\n");
  printf ("       -s <int>  number of samples per recipient haplotype (default=10)\n");
  printf ("       -n <double>  average switch rate parameter start-value (default=400000 divided by total number of donor haplotypes included in analysis)\n");
  printf ("       -p  specify to use prior copying probabilities in population list (-f) file\n");
  printf ("       -m  specify to use donor population mutation (emission) probabilities in population list (-f) file\n");
  printf ("       -M <double>  global mutation (emission) probability (default=Li & Stephen's (2003) fixed estimate)\n");
  printf ("       -k <double>  specify number of expected chunks to define a 'region' (default=100)\n");
  printf ("       -j  specify that individuals are haploid\n");
  printf ("       -u  specify that data are unlinked\n");
  printf ("       -a <a_1> <a_2>  paint individuals a_1 through a_2 using every other individual (use '-a 0 0' to paint all inds)\n");
  printf ("       -b  print-out zipped file with suffix '.copyprobsperlocus.out' containing prob each recipient copies each donor at every SNP (note: file can be quite large)\n");
  printf ("       -o <outfile-prefix>  (default = 'geno.filein')\n");
  printf ("       --help  print this menu\n");
}

int main(int argc, char *argv[])
{
  struct data_t *Data;
  int i,j,k;
  double bpval;
  double totaldonorprobs;
  int log_lik_check;
  int ndonors, ndonorpops, nrecpops;
  int nind, nsites, cond_nhaps, cond_nind, num_rec_drift;
  int donor_count, rec_count, nind_tot, nind_totGENFILE, rec_ind_count, rec_ind_topaint_count, ind_label_find, ind_pop_find;
  int geno_find, recom_find, donorlist_find, EMiter_find, numsamples_find, outfile_find, ne_find, mut_find, num_found, copy_prop_em_find, recom_em_find, mutation_em_find, mutationALL_em_find, all_versus_all_ind, haploid_ind, region_size_find, unlinked_ind, prior_donor_probs_ind, mutation_rate_ind, print_file9_ind, drift_calc_ind, idfile_find;
  char *step;
  char line[2047];
  //char * bigline=malloc(1000000*sizeof(char));
  char waste[2047];
  char waste2[2047];
  char templab[15];
  char donorname[2047];
  char indname[2047];
  FILE *fd, *fd2, *fd3, *fd4, *fout, *fout2, *fout3, *fout4, *fout5, *fout6, *fout7, *fout8, *fout9;
  char * filename = malloc(1000 * sizeof(char));
  char * filenameGEN = malloc(1000 * sizeof(char));
  char * filenameDONORLIST = malloc(1000 * sizeof(char));
  char * filenameID = malloc(1000 * sizeof(char));
  char * filenameOUT = malloc(1000 * sizeof(char));
  srand((unsigned)time(NULL));

  /***********************************************************/
  // DEFAULT VALUES:

  int EMruns = 0;        // number of EMruns 
  int samplesTOT = 10;   // number of final hidden-state samples desired after E-M is finished
  double N_e = 400000;   // scaling constant for recombination rate
  double GlobalMutRate = -9.0;    // global mutation rate per site
  double region_size = 100;    // number of chunks per region -- used to look at variability in copied chunks across regions in order to estimate "c" in Dan Lawson's fineSTRUCTURE

      /* 'NUISSANCE' PARAMETER DETAILS: */
  //double theta = 0.0001;
  //double theta = -9.0;
  double small_recom_val=0.000000000000001;    // lower limit for small genmap rates

  int start_val=0;
  int end_val=0;

  /************************************************************/

  geno_find=0;
  recom_find=0;
  donorlist_find=0;
  idfile_find=0;
  outfile_find=0;
  EMiter_find=0;
  numsamples_find=0;
  ne_find=0;
  mut_find=0;
  region_size_find=0;
  copy_prop_em_find=0;
  recom_em_find=0;
  mutation_em_find=0;
  mutationALL_em_find=0;
  all_versus_all_ind=0;
  haploid_ind=0;
  unlinked_ind=0;
  prior_donor_probs_ind=0;
  mutation_rate_ind=0;  
  print_file9_ind=0;
  drift_calc_ind=0;
  num_found=0;
  for (i=1; i < argc; i++)
    {
      if ((strcmp(argv[i],"-help")==0) || (strcmp(argv[i],"--help")==0) || (strcmp(argv[i],"-h")==0))
	{
	  usage();
	  exit(1);
	}
      if (strcmp(argv[i],"-g")==0)
	{
	  geno_find=1;
	  num_found=num_found+1;
	}
      if (strcmp(argv[i],"-r")==0)
	{
	  recom_find=1;
	  num_found=num_found+1;
	}
      if (strcmp(argv[i],"-f")==0)
	donorlist_find=1;
      if (strcmp(argv[i],"-t")==0)
	{
	  idfile_find=1;
	  num_found=num_found+1;
	}
      if (strcmp(argv[i],"-i")==0)
	{
	  EMiter_find=1;
	  num_found=num_found+1;
	}
      if (strcmp(argv[i],"-s")==0)
	{
	  numsamples_find=1;
	  num_found=num_found+1;
	}
       if (strcmp(argv[i],"-n")==0)
	{
	  ne_find=1;
	  num_found=num_found+1;
	}
       if (strcmp(argv[i],"-M")==0)
	{
	  mut_find=1;
	  num_found=num_found+1;
	}
       if (strcmp(argv[i],"-k")==0)
	{
	  region_size_find=1;
	  num_found=num_found+1;
	}
       if (strcmp(argv[i],"-o")==0)
	{
	  outfile_find=1;
	  num_found=num_found+1;
	}
       if (strcmp(argv[i],"-ip")==0)
	 copy_prop_em_find=1;
       if (strcmp(argv[i],"-in")==0)
	 recom_em_find=1;
       if (strcmp(argv[i],"-im")==0)
	 mutation_em_find=1;
       if (strcmp(argv[i],"-iM")==0)
	 mutationALL_em_find=1;
       //       if (strcmp(argv[i],"-c")==0)
       //	 condition_recipient_inds_find=1;
       if (strcmp(argv[i],"-a")==0)
	 all_versus_all_ind=1;
       if (strcmp(argv[i],"-j")==0)
	 haploid_ind=1;
       if (strcmp(argv[i],"-u")==0)
	 unlinked_ind=1;
       if (strcmp(argv[i],"-p")==0)
	 prior_donor_probs_ind=1;
       if (strcmp(argv[i],"-b")==0)
	 print_file9_ind=1;
       //if (strcmp(argv[i],"-y")==0)
       //indcount_suppress_ind=1;
       if (strcmp(argv[i],"-m")==0)
	 {
	   mutation_rate_ind=1;
	   //num_found=num_found+1;
	 }
    }
  if (argc != (num_found*2+copy_prop_em_find+recom_em_find+mutation_em_find+mutationALL_em_find+mutation_rate_ind+4*donorlist_find+3*all_versus_all_ind+haploid_ind+unlinked_ind+prior_donor_probs_ind+print_file9_ind+1))
    {
      printf("Something wrong with input command line (missing arguments?). Exiting....\n");
      usage();
      exit(1);
    }
  if (donorlist_find==0)
    strcpy(filenameDONORLIST,"NULL");
  if (recom_find==0)
    strcpy(filename,"NULL");
  if (idfile_find==0)
    strcpy(filenameID,"NULL");
  if (((geno_find==0) || (recom_find==0)) && (unlinked_ind==0)) { printf("Error with command line (Each of -g and -r MUST be specified if data are linked). Exiting...\n"); exit(1);}
  if ((geno_find==0) && (unlinked_ind==1)) { printf("Error with command line (-g MUST be specified). Exiting...\n"); exit(1);}
  if ((recom_find==1) && (unlinked_ind==1)) { printf("Data specified as containing unlinked sites (-u). Ignoring supplied recombination rate file....\n");}
  if (all_versus_all_ind==0 && (idfile_find==0 || donorlist_find==0)){ printf("Unless performing all-versus-all ('-a' switch), you MUST specify both the id file ('-t' switch) and file listing donor and recipient populations ('-f' switch). Exiting....\n"); exit(1);}
  if ((mutation_em_find==1) && (mutationALL_em_find==1))
    {
      printf("You have specified to estimate a global mutation (emission) rate and population-specific mutation (emission) rates. Please choose only one of the '-im' and '-iM' switches. Exiting...\n");
      exit(1);
    }
  if ((mutation_rate_ind==1) && (mut_find==1))
    {
      printf("You have provided values for both a global mutation (emission) rate ('-M') and population-specific mutation (emission) rates ('-m'). Please choose only one of the '-m' and '-M' switches. Exiting...\n");
      exit(1);
    }
  if ((mutation_rate_ind==1) && (mutationALL_em_find==1))
    {
      printf("You have specified to estimate a global mutation (emission) rate; will ignore population-specific mutation (emission) rates in %s. If you wish to use donor-specific mutation rates, use the '-im' switch.\n",filenameDONORLIST);
    }
 for (i=1; i < argc; i++)
    {
      if (strcmp(argv[i],"-g")==0)
	{
	  if (argv[(i+1)][0] == '-')
	    {
	      printf("Something wrong with input command line (missing arguments?). Exiting....\n");
	      usage();
	      exit(1);
	    }
 	  strcpy(filenameGEN,argv[(i+1)]);
	  if (outfile_find==0)
	    {
	      strcpy(filenameOUT,argv[(i+1)]);
	      fout = fopen(strcat(filenameOUT,".samples.out"), "w");
	      strcpy(filenameOUT,argv[(i+1)]);
	      fout2 = fopen(strcat(filenameOUT,".priorprobs.out"), "w");
	      strcpy(filenameOUT,argv[(i+1)]);
	      fout3 = fopen(strcat(filenameOUT,".EMprobs.out"), "w");
	      strcpy(filenameOUT,argv[(i+1)]);
	      fout4 = fopen(strcat(filenameOUT,".chunkcounts.out"), "w");
	      strcpy(filenameOUT,argv[(i+1)]);
	      fout5 = fopen(strcat(filenameOUT,".chunklengths.out"), "w");
	      strcpy(filenameOUT,argv[(i+1)]);
	      fout6 = fopen(strcat(filenameOUT,".mutationprobs.out"), "w");
	      strcpy(filenameOUT,argv[(i+1)]);
	      fout7 = fopen(strcat(filenameOUT,".regionchunkcounts.out"), "w");
	      strcpy(filenameOUT,argv[(i+1)]);
	      fout8 = fopen(strcat(filenameOUT,".regionsquaredchunkcounts.out"), "w");
	      strcpy(filenameOUT,argv[(i+1)]);
	      fout9 = gzopen(strcat(filenameOUT,".copyprobsperlocus.out.gz"), "w");
	      strcpy(filenameOUT,argv[(i+1)]);
	    }
	}
      if (strcmp(argv[i],"-r")==0)
	{
	  if (argv[(i+1)][0] == '-')
	    {
	      printf("Something wrong with input command line (missing arguments?). Exiting....\n");
	      usage();
	      exit(1);
	    }
	  strcpy(filename,argv[(i+1)]);
	}
      if (strcmp(argv[i],"-f")==0)
	{
	  if ((argv[(i+1)][0] == '-') || (argv[(i+2)][0] == '-') || (argv[(i+3)][0] == '-'))
	    {
	      printf("Something wrong with input command line (missing arguments?). Exiting....\n");
	      usage();
	      exit(1);
	    }
	  strcpy(filenameDONORLIST,argv[(i+1)]);
	  start_val = atoi(argv[(i+2)]);
	  end_val = atoi(argv[(i+3)]);
	  if ((end_val < start_val) || (start_val < 0) || (end_val < 0))
	    {
	      printf("Invalid start_ind/stop_ind vals ('-f' switch). If you want to paint each recipient individual using every donor individual, use '-f <donorlist.filein> 0 0'. Exiting...\n");
	      exit(1);
	    }
	  if (start_val > 0) start_val=start_val-1;
	}
      if (strcmp(argv[i],"-t")==0)
	{
	  if (argv[(i+1)][0] == '-')
	    {
	      printf("Something wrong with input command line (missing arguments?). Exiting....\n");
	      usage();
	      exit(1);
	    }
	  strcpy(filenameID,argv[(i+1)]);
	}
     if (strcmp(argv[i],"-i")==0)
	{
	  if (argv[(i+1)][0] == '-')
	    {
	      printf("Something wrong with input command line (missing arguments?). Exiting....\n");
	      usage();
	      exit(1);
	    }
	  EMruns = atoi(argv[(i+1)]);
	  if (EMruns < 0)
	    {
	      printf("Number of EM runs must be at least 0. Exiting...\n");
	      exit(1);
	    }
	  if ((EMruns>0) && (copy_prop_em_find==0) && (recom_em_find==0) && (mutation_em_find==0) && (mutationALL_em_find==0))
	    {
	      printf("You have specified to perform E-M iterations, but have not specified which parameter(s) to maximize. If using '-i' switch, please specify at least one of '-in', '-ip', '-im', and/or '-iM'. Exiting...\n");
	      exit(1);
	    }	  
	}
      if (strcmp(argv[i],"-s")==0)
	{
	  if (argv[(i+1)][0] == '-')
	    {
	      printf("Something wrong with input command line (missing arguments?). Exiting....\n");
	      usage();
	      exit(1);
	    }
	  samplesTOT = atoi(argv[(i+1)]);
	  if (samplesTOT < 0)
	    {
	      printf("Number of samples must be >= 0. Exiting...\n");
	      exit(1);
	    }
	}
       if (strcmp(argv[i],"-n")==0)
	{
	  if (argv[(i+1)][0] == '-')
	    {
	      printf("Something wrong with input command line (missing arguments?). Exiting....\n");
	      usage();
	      exit(1);
	    }
	  N_e = atof(argv[(i+1)]);
	  if (N_e <= 0)
	    {
	      printf("Recombination scaling parameter N_e must be > 0. Exiting...\n");
	      exit(1);
	    }
	}
       if (strcmp(argv[i],"-M")==0)
	{
	  if (argv[(i+1)][0] == '-')
	    {
	      printf("Something wrong with input command line (missing arguments?). Exiting....\n");
	      usage();
	      exit(1);
	    }
	  GlobalMutRate = atof(argv[(i+1)]);
	  if (GlobalMutRate==0) GlobalMutRate=-9;
	  if (GlobalMutRate < 0)
	    printf("Mutation (emission) parameter must be > 0. Using Li & Stephens (2003) version of Watterson's estimate instead of user-supplied value...\n");
	}
       if (strcmp(argv[i],"-k")==0)
	{
	  if (argv[(i+1)][0] == '-')
	    {
	      printf("Something wrong with input command line (missing arguments?). Exiting....\n");
	      usage();
	      exit(1);
	    }
	  region_size = atof(argv[(i+1)]);
	  if (region_size < 1)
	    {
	      printf("Region_size must be >= 1. Exiting...\n");
	      exit(1);
	    }
	}
         if (strcmp(argv[i],"-a")==0)
	   {
	     if ((argv[(i+1)][0] == '-') || (argv[(i+2)][0] == '-'))
	       {
		 printf("Something wrong with input command line (missing arguments?). Exiting....\n");
		 usage();
		 exit(1);
	       }
	  start_val = atoi(argv[(i+1)]);
	  end_val = atoi(argv[(i+2)]);
	  if ((end_val < start_val) || (start_val < 0) || (end_val < 0))
	    {
	      printf("Invalid start_ind/stop_ind vals ('-a' switch). If you want to condition each individual on every other individual, use '-a 0 0'. Exiting...\n");
	      exit(1);
	    }
	  if (start_val > 0) start_val=start_val-1;
	}
       if (strcmp(argv[i],"-o")==0)
	 {
	  if (argv[(i+1)][0] == '-')
	    {
	      printf("Something wrong with input command line (missing arguments?). Exiting....\n");
	      usage();
	      exit(1);
	    }
	   strcpy(filenameOUT,argv[(i+1)]);
	   fout = fopen(strcat(filenameOUT,".samples.out"), "w");
  	   strcpy(filenameOUT,argv[(i+1)]);
	   fout2 = fopen(strcat(filenameOUT,".priorprobs.out"), "w");
  	   strcpy(filenameOUT,argv[(i+1)]);
	   fout3 = fopen(strcat(filenameOUT,".EMprobs.out"), "w");
	   strcpy(filenameOUT,argv[(i+1)]);
	   fout4 = fopen(strcat(filenameOUT,".chunkcounts.out"), "w");
	   strcpy(filenameOUT,argv[(i+1)]);
	   fout5 = fopen(strcat(filenameOUT,".chunklengths.out"), "w");
	   strcpy(filenameOUT,argv[(i+1)]);
	   fout6 = fopen(strcat(filenameOUT,".mutationprobs.out"), "w");
	   strcpy(filenameOUT,argv[(i+1)]);
	   fout7 = fopen(strcat(filenameOUT,".regionchunkcounts.out"), "w");
	   strcpy(filenameOUT,argv[(i+1)]);
	   fout8 = fopen(strcat(filenameOUT,".regionsquaredchunkcounts.out"), "w");
	   strcpy(filenameOUT,argv[(i+1)]);
	   fout9 = gzopen(strcat(filenameOUT,".copyprobsperlocus.out.gz"), "w");
	   strcpy(filenameOUT,argv[(i+1)]);
	 }
    }
  if (fout == NULL) {printf("error opening closing file1\n"); exit(1);}
  if (fout2 == NULL) {printf("error opening closing file2\n"); exit(1);}
  if (fout3 == NULL) {printf("error opening closing file3\n"); exit(1);}
  if (fout4 == NULL) {printf("error opening closing file4\n"); exit(1);}
  if (fout5 == NULL) {printf("error opening closing file5\n"); exit(1);}
  if (fout6 == NULL) {printf("error opening closing file6\n"); exit(1);}
  if (fout7 == NULL) {printf("error opening closing file7\n"); exit(1);}
  if (fout8 == NULL) {printf("error opening closing file8\n"); exit(1);}
  if (fout9 == NULL) {printf("error opening closing file9\n"); exit(1);}

  fd = fopen(filenameGEN,"r");
  if (fd == NULL) { printf("error opening %s\n",filenameGEN); exit(1);}
  fgets(line,2047,fd);
  sscanf(line,"%d",&nind_totGENFILE);
  nind_totGENFILE=nind_totGENFILE/(2-haploid_ind);
  fclose(fd);
  if (idfile_find==0)
    {
      // get total number of inds by reading in phase file:
      // fd = fopen(filenameGEN,"r");
      //if (fd == NULL) { printf("error opening %s\n",filenameGEN); exit(1);}
      //fgets(line,2047,fd);
      //sscanf(line,"%d",&nind_tot);
      //nind_tot=nind_tot/(2-haploid_ind);
      /* ??VERY UNCLEAR WHY THIS DOESN'T WORK??
      fgets(line,2047,fd);
      fgets(line,2047,fd);
      nind_tot=0;
      while(!feof(fd))
	{
	  if (fgets(bigline,2047,fd)!=NULL)
	    nind_tot=nind_tot+1;
	}
      fclose(fd);
      printf("%d\n",nind_tot);
      nind_tot=nind_tot/(2-haploid_ind);
      printf("%d\n",nind_tot);
      */
    }
  if (idfile_find==0) nind_tot=nind_totGENFILE;
  if (idfile_find==1)
    {
      // open id file (to get total number of inds)
      fd4 = fopen(filenameID,"r");
      if (fd4 == NULL) { printf("error opening %s\n",filenameID); exit(1);}
      nind_tot=0;
      while(!feof(fd4))
	{
	  if (fgets(line,2047,fd4)!=NULL)
	    nind_tot=nind_tot+1;
	}
      fclose(fd4);
      if (nind_tot != nind_totGENFILE) {printf("number of inds in %s (%d) does not match number of inds in %s (%d)\n. Exiting....\n",filenameID,nind_tot,filenameGEN,nind_totGENFILE); exit(1);}
    }
      // find inds to include (and get population labels for each ind):
  int * pop_vec_tot=malloc(((2-haploid_ind)*nind_tot) * sizeof(int));
  char ** ind_label_vec=malloc(nind_tot * sizeof(char *));
  char ** pop_label_vec=malloc(nind_tot * sizeof(char *));
  for (i=0; i < nind_tot; i++)
    {
      ind_label_vec[i]=malloc(1000*sizeof(char));
      pop_label_vec[i]=malloc(1000*sizeof(char));
    }
  int * include_ind_vec=malloc(nind_tot*sizeof(int));
  if (idfile_find==0)
    {
      for (i=0; i < nind_tot; i++)
	{
	  include_ind_vec[i]=1;
	  sprintf(templab,"%d",i+1);
	  strcpy(indname,"IND");
	  strcat(indname,templab);
	  strcpy(ind_label_vec[i],indname);
	}
      nind=nind_tot;
    }
  if (idfile_find==1)
    {
      fd4 = fopen(filenameID,"r");
      if (fd4 == NULL) { printf("error opening %s\n",filenameID); exit(1);}
      nind=nind_tot;
      for (i=0; i < nind_tot; i++)
	{
	  fgets(line,2047,fd4);
	  step=line;
	  reading(&step,"%s",waste);
	  strcpy(ind_label_vec[i],waste);
	  reading(&step,"%s",waste);
	  strcpy(pop_label_vec[i],waste);
	  reading(&step,"%s",waste);
	  include_ind_vec[i]=1;
	  if (strcmp(waste,"0")==0)
	    {
	      include_ind_vec[i]=0;
	      nind=nind-1;
	    }
	}
      fclose(fd4);
    }
  if (all_versus_all_ind==1)
    {
      if (donorlist_find==1 && idfile_find==1)
	{
	  for (i=0; i < nind_tot; i++)
	    {
	      if (include_ind_vec[i]==1)
		{
		  // open donor-list file (to get information on number of donor/rec pops)
		  fd3 = fopen(filenameDONORLIST,"r");
		  if (fd3 == NULL) { printf("error opening %s\n",filenameDONORLIST); exit(1);}
		  ind_pop_find=0;
		  while(!feof(fd3))
		    {
		      if (fgets(line,2047,fd3)!=NULL)
			{
			  step=line;
			  reading(&step,"%s",waste);
			  reading(&step,"%s",waste2);
			  if (strcmp(waste2,"D")!=0 && strcmp(waste2,"R")!=0){ printf("Second column of %s must contain either a 'D' (for 'donor') or 'R' (for 'recipient'). Exiting....\n",filenameDONORLIST); exit(1);}
			  if (strcmp(pop_label_vec[i],waste)==0)
			    {
			      ind_pop_find=1;
			      break;
			    }
			}
		    }
		  fclose(fd3);
		  if (ind_pop_find==0)
		    { 
		      printf("Population label for ind %s (%s) not found in %s, so excluding this individual....\n",ind_label_vec[i],pop_label_vec[i],filenameDONORLIST);
		      include_ind_vec[i]=0;
		      nind=nind-1;
		    }
		}
	    }
	}

      ndonors=(2-haploid_ind)*nind-2+haploid_ind;
      ndonorpops=ndonors/(2-haploid_ind);
      nrecpops=1;
      donor_count=0;
      for (i=0; i < nind_tot; i++)
	{
	  for (j=0; j < (2-haploid_ind); j++)
	    pop_vec_tot[(2-haploid_ind)*i+j]=-9;
	  if (include_ind_vec[i]==1)
	    {
	      for (j=0; j < (2-haploid_ind); j++)
		pop_vec_tot[(2-haploid_ind)*i+j]=donor_count;
	      donor_count=donor_count+1;
	    }
	}
    }
  if (all_versus_all_ind==0)
    {
      // open donor-list file (to get information on number of donor/rec pops)
      fd3 = fopen(filenameDONORLIST,"r");
      if (fd3 == NULL) { printf("error opening %s\n",filenameDONORLIST); exit(1);}
      ndonorpops=0;
      nrecpops=0;
      while(!feof(fd3))
	{
	  if (fgets(line,2047,fd3)!=NULL)
	    {
	      step=line;
	      reading(&step,"%s",waste);
	      reading(&step,"%s",waste2);
	      if (strcmp(waste2,"D")!=0 && strcmp(waste2,"R")!=0){ printf("Second column of %s must contain either a 'D' (for 'donor') or 'R' (for 'recipient'). Exiting....\n",filenameDONORLIST); exit(1);}
	      if (strcmp(waste2,"D")==0)
		ndonorpops=ndonorpops+1;
	      if (strcmp(waste2,"R")==0)
		nrecpops=nrecpops+1;
	    }
	}
      fclose(fd3);
      if (ndonorpops==0) { printf("No donor populations found in %s! Exiting....\n",filenameDONORLIST); exit(1);}
      if (nrecpops==0) { printf("No recipient populations found in %s! Exiting....\n",filenameDONORLIST); exit(1);}
    }
  int * ndonorhaps_tot=malloc(ndonorpops*sizeof(int));
  int * nrechaps=malloc(nrecpops*sizeof(int));
  double * ndonorprobs = malloc(ndonorpops * sizeof(double));
  double * ndonormutrates = malloc(ndonorpops * sizeof(double));
  char ** donor_pop_vec=malloc(ndonorpops*sizeof(char *));
  for(i=0; i < ndonorpops; i++)
    donor_pop_vec[i]=malloc(1000*sizeof(char));
  char ** rec_pop_vec=malloc(nrecpops*sizeof(char *));
  for(i=0; i < nrecpops; i++)
    rec_pop_vec[i]=malloc(1000*sizeof(char));
  if (all_versus_all_ind==1)
    {
      nrechaps[0]=2-haploid_ind;
      strcpy(rec_pop_vec[0],"RECPOP");
      for (i=0; i < ndonorpops; i++)
	{
	  sprintf(templab,"%d",i+1);
	  strcpy(donorname,"DONOR");
	  strcat(donorname,templab);
	  strcpy(donor_pop_vec[i],donorname);
	  ndonorhaps_tot[i]=2-haploid_ind;
	}
      cond_nhaps = (2-haploid_ind)*nind;
    }
  if (all_versus_all_ind==0)
    {
      for (i=0; i < ndonorpops; i++)
	ndonorhaps_tot[i]=0;
      for (i=0; i < nrecpops; i++)
	nrechaps[i]=0;
      // open donor-list file (to get information on donor/rec pop labels)
      fd3 = fopen(filenameDONORLIST,"r");
      if (fd3 == NULL) { printf("error opening %s\n",filenameDONORLIST); exit(1);}
      donor_count=0;
      rec_count=0;
      totaldonorprobs=0.0;
      while(!feof(fd3))
	{
	  if (fgets(line,2047,fd3)!=NULL)
	    {
	      step=line;
	      reading(&step,"%s",waste);
	      reading(&step,"%s",waste2);
	      if (strcmp(waste2,"R")==0)
		{
		  strcpy(rec_pop_vec[rec_count],waste);
		  rec_count=rec_count+1;
		}
	      if (strcmp(waste2,"D")==0)
		{
		  strcpy(donor_pop_vec[donor_count],waste);
		  if (prior_donor_probs_ind==1)
		    {
		      reading(&step,"%lf",&ndonorprobs[donor_count]);
		      if (ndonorprobs[donor_count]<=0.0000000000000001)
			{
			  printf("Donor copying probabilities in %s must be > 0. Exiting...\n",filenameDONORLIST);
			  exit(1);
			}
		      totaldonorprobs=totaldonorprobs+ndonorprobs[donor_count];
		    }
		  if (mutation_rate_ind==1)
		    {
		      if (prior_donor_probs_ind==0) reading(&step,"%s",waste);
		      reading(&step,"%lf",&ndonormutrates[donor_count]);
		      if (ndonormutrates[donor_count]>1)
			{
			  printf("Donor mutation (emission) probabilities must be <= 1 in %s (use a negative number to specify default). Exiting...\n",filenameDONORLIST);
			  exit(1);
			}
		    }
		  donor_count=donor_count+1;
		}
	    }
	}
      fclose(fd3);
      if (prior_donor_probs_ind==1)
	{
	  if ((totaldonorprobs > 1.00001) || (totaldonorprobs < 0.99999))
	    {
	      printf("Probabilities across all donors in %s does not sum to 1.0 (instead sums to %lf)! Rescaling to sum to 1.0....\n",filenameDONORLIST,totaldonorprobs);
	      for (k=0; k < ndonorpops; k++)
		ndonorprobs[k]=ndonorprobs[k]/totaldonorprobs;
	    }
	}

      // open id file (to get information on donor pop hap numbers and do some checks if necessary)
      fd4 = fopen(filenameID,"r");
      if (fd4 == NULL) { printf("error opening %s\n",filenameID); exit(1);}
      nind=nind_tot;
      for (i=0; i < nind_tot; i++)
	{
	  for (j=0; j < (2-haploid_ind); j++)
	    pop_vec_tot[(2-haploid_ind)*i+j]=-9;
	  fgets(line,2047,fd4);
	  step=line;
	  reading(&step,"%s",waste);
	  reading(&step,"%s",waste);
	  reading(&step,"%s",waste);
	  if (include_ind_vec[i]==0) nind=nind-1;
	  if (include_ind_vec[i]==1)
	    {
	      ind_label_find=0;
	      for (k=0; k < ndonorpops; k++)
		{
		  if (strcmp(pop_label_vec[i],donor_pop_vec[k])==0)
		    {
		      ind_label_find=1;
		      for (j=0; j < (2-haploid_ind); j++)
			pop_vec_tot[(2-haploid_ind)*i+j]=k;
		      ndonorhaps_tot[k]=ndonorhaps_tot[k]+2-haploid_ind;
		      break;
		    }
		}
	      for (k=0; k < nrecpops; k++)
		{
		  if (strcmp(pop_label_vec[i],rec_pop_vec[k])==0)
		    {
		      ind_label_find=1;
		      nrechaps[k]=nrechaps[k]+2-haploid_ind;
		      break;
		    }
		}
	      if (ind_label_find==0)
		{
		  include_ind_vec[i]=0;
		  nind=nind-1;
		  printf("Individual %d (%s) in %s has a population label (%s) not found in %s, so is being excluded.....\n",i+1,ind_label_vec[i],filenameID,pop_label_vec[i],filenameDONORLIST);
		  //printf("Individual %d in %s has a population label (%s) not found in %s, which is not allowed (unless '-a' switch used)! Exiting....\n",ind_count+1;filenameID,pop_label_vec[indcount],filenameDONORLIST);
		  //exit(1);
		}
	    }
	}
      fclose(fd4);
      ndonors=0;
      for (k=0; k < ndonorpops; k++)
	{
	  if (ndonorhaps_tot[k]==0)
	    {
	      printf("No individuals found in %s with population label %s found in %s, which is not allowed (unless '-a' switch used)! Exiting....\n",filenameID,donor_pop_vec[k],filenameDONORLIST);
	      exit(1);
	    }
	  for (j=0; j < nrecpops; j++)
	    {
	      if (strcmp(donor_pop_vec[k],rec_pop_vec[j])==0)
		{
		  printf("WARNING: Population %s is specified as a donor and as a recipient. Will allow 'self-copying' of own population label in recipient individuals....\n",rec_pop_vec[j]);
		  //ndonorhaps_tot[k]=ndonorhaps_tot[k]-2+haploid_ind;
		  break;
		}
	    }
	  ndonors=ndonors+ndonorhaps_tot[k];
	}
      cond_nhaps=0;
      for (k=0; k < nrecpops; k++)
	{
	  if (nrechaps[k]==0)
	    {
	      printf("No individuals found in %s with population label %s found in %s, which is not allowed (unless '-a' switch used)! Exiting....\n",filenameID,rec_pop_vec[k],filenameDONORLIST);
	      exit(1);
	    }
	  cond_nhaps=cond_nhaps+nrechaps[k];
	}
    }
  cond_nind=cond_nhaps/(2-haploid_ind);
              // (0) INITIALIZE copy_prob_tot, copy_probSTART_tot, MutProb_vec_tot:
  double * MutProb_vec_tot = malloc(((2-haploid_ind)*nind_tot) * sizeof(double));
  double * copy_prob_tot = malloc(((2-haploid_ind)*nind_tot) * sizeof(double));
  double * copy_probSTART_tot = malloc(((2-haploid_ind)*nind_tot) * sizeof(double));
  for (i=0; i < nind_tot; i++)
    {
      for (j=0; j < (2-haploid_ind); j++)
	{
	  MutProb_vec_tot[(2-haploid_ind)*i+j]=0;
	  copy_prob_tot[(2-haploid_ind)*i+j]=-9;
	  copy_probSTART_tot[(2-haploid_ind)*i+j]=-9;
	}
      if (include_ind_vec[i]==1 && prior_donor_probs_ind==0)
	{
	  for (j=0; j < (2-haploid_ind); j++)
	    {
	      //copy_prob_tot[(2-haploid_ind)*i+j] = 1.0/ndonors;
	      //copy_probSTART_tot[(2-haploid_ind)*i+j] = 1.0/ndonors;
	      copy_prob_tot[(2-haploid_ind)*i+j] = 1.0;
	      copy_probSTART_tot[(2-haploid_ind)*i+j] = 1.0;
	    }
	}
      if (include_ind_vec[i]==1 && prior_donor_probs_ind==1)
	{
	  for (k=0; k < ndonorpops; k++)
	    {
	      if (pop_vec_tot[(2-haploid_ind)*i]==k)
		{
		  for (j=0; j < (2-haploid_ind); j++)
		    {
		      copy_prob_tot[(2-haploid_ind)*i+j] = ndonorprobs[k]/ndonorhaps_tot[k];
		      copy_probSTART_tot[(2-haploid_ind)*i+j] = ndonorprobs[k]/ndonorhaps_tot[k];
		    }
		  break;
		}
	    }
	}
      if (include_ind_vec[i]==1 && (mutation_rate_ind==0 || mutationALL_em_find==1))
	{
	  for (j=0; j < (2-haploid_ind); j++)
	    MutProb_vec_tot[(2-haploid_ind)*i+j]=GlobalMutRate;
	}
      if (include_ind_vec[i]==1 && mutation_rate_ind==1 && mutationALL_em_find==0)
	{
	  for (k=0; k < ndonorpops; k++)
	    {
	      if (pop_vec_tot[(2-haploid_ind)*i]==k)
		{
		  for (j=0; j < (2-haploid_ind); j++)
		    MutProb_vec_tot[(2-haploid_ind)*i+j] = ndonormutrates[k];
		  break;
		}
	    }
	}
    }

      // open haplotype file (to get information on haplotype numbers)
  fd = fopen(filenameGEN,"r");
  if (fd == NULL) { printf("error opening %s\n",filenameGEN); exit(1);}
  Data = ReadData(fd,0,2-haploid_ind,include_ind_vec,pop_label_vec,ndonorpops,donor_pop_vec,pop_vec_tot,copy_prob_tot,copy_probSTART_tot,MutProb_vec_tot,ndonorhaps_tot,all_versus_all_ind);

  nsites = Data->nsnps;
  if (idfile_find==1 && (nind_tot != (((int) Data->nhaps)/(2-haploid_ind)))) { printf("Total number of individuals in %s does not match that in %s. Exiting....\n",filenameGEN,filenameID); exit(1);}

  if (all_versus_all_ind==1) printf("Will condition each individual on every other individual...\n");
  if (all_versus_all_ind==1 && idfile_find==1) printf("Will use %s for population labels in output (even though performing all-versus-all; i.e. '-a' switch)\n",filenameID);
  //if (all_versus_all_ind==1 && donorlist_find==1) printf("As you have specified '-a' switch for all-versus-all, will ignore '-f' input, including file %s listing donor/recipient populations.\n",filenameDONORLIST);
  if (all_versus_all_ind==1 && donorlist_find==1 && idfile_find==0) printf("As you have specified '-a' switch for all-versus-all but have not specified label file with '-t', will ignore '-f' input, i.e. file %s listing populations to include in all-versus-all analysis.\n",filenameDONORLIST);
  if (all_versus_all_ind==1 && donorlist_find==1 && idfile_find==1) printf("As you have specified '-a' switch for all-versus-all and supplied '-t' and '-f' files, will only consider (non-excluded) individuals in file %s that are from all populations listed in file %s.\n",filenameID,filenameDONORLIST);
  if (prior_donor_probs_ind==1) 
    printf("Using specified prior donor probs from input file....\n");
  if ((mutation_rate_ind==1) && (mutationALL_em_find==0))
    printf("Using specified donor population mutation rates from input file....\n");
  if (copy_prop_em_find==1)
    printf("Running E-M to estimate copying proportions....\n");
  if (recom_em_find==1)
    printf("Running E-M to estimate N_e....\n");
  if (mutation_em_find==1)
    printf("Running E-M to estimate mutation (emission) probabilities....\n");
  if (mutationALL_em_find==1)
    printf("Running E-M to estimate global mutation (emission) probability....\n");
  printf(" Total number of individuals found: %d\n",nind_tot);
  printf(" Total number of individuals included for analysis: %d\n",nind);
  if (all_versus_all_ind==0)
    {
      printf(" Number of donor pops: %d\n",ndonorpops);
      printf(" Number of recipient pops: %d\n",nrecpops);
    }
  if (all_versus_all_ind==1)
    {
      printf(" Number of donor pops: %d (all-versus-all)\n",ndonorpops);
      printf(" Number of recipient pops: %d (all-versus-all)\n",nrecpops);
    }
  if (print_file9_ind==1)
    {
      gzprintf(fout9,"pos");
      if (all_versus_all_ind==0)
	{
	  for (k=0; k < ndonorpops; k++)
	    gzprintf(fout9," %s",donor_pop_vec[k]);
	}
      if (all_versus_all_ind==1)
	{
	  for (j=0; j < nind_tot; j++)
	    {
	      if (include_ind_vec[j]==1)
		gzprintf(fout9," %s",ind_label_vec[j]);
	    }
	}
      gzprintf(fout9,"\n");
    }
  //if (ne_find==0) N_e=N_e/ndonors;

  if (start_val >= cond_nind) 
    {
      printf("Fewer recipient individuals than expected -- either your '-a' switch is mis-specified or there is something wrong with your input file %s? Exiting....\n",filenameGEN);
      exit(1);
    }
  if (drift_calc_ind==0) num_rec_drift=0;
  if (drift_calc_ind==1 && num_rec_drift>cond_nind)
    {
      printf("You have specified to store %d individuals in memory using the '-d' switch, but you have only %d recipient individuals. Will set this value to %d.\n",num_rec_drift,cond_nind,cond_nind);
      num_rec_drift=cond_nind;
    }
  if (drift_calc_ind==1 && num_rec_drift==0) num_rec_drift=cond_nind;
  if (end_val==0 || end_val > cond_nind) end_val=cond_nind;
  int * recipient_ind_vec=malloc(nind_tot*sizeof(int));
  rec_ind_count=0;
  rec_ind_topaint_count=0;
  for (i=0; i < nind_tot; i++)
    {
      recipient_ind_vec[i]=0;
      if (idfile_find==0)
	{
	  if (i>=start_val && i<end_val)
	    {
	      recipient_ind_vec[i]=1;		
	      rec_ind_topaint_count=rec_ind_topaint_count+1;
	    }
	}
      if (idfile_find==1)
	{
	  if (all_versus_all_ind==0 && include_ind_vec[i]==1)
	    {
	      ind_label_find=0;
	      for (k=0; k < nrecpops; k++)
		{
		  if (strcmp(pop_label_vec[i],rec_pop_vec[k])==0)
		    {
		      ind_label_find=1;
		      break;
		    }
		}
	      if (ind_label_find==1)
		{
		  if (rec_ind_count>=start_val && rec_ind_count < end_val)
		    {
		      recipient_ind_vec[i]=1;
		      rec_ind_topaint_count=rec_ind_topaint_count+1;
		    }
		  rec_ind_count=rec_ind_count+1;
		}
	    }
	  if (all_versus_all_ind==1 && include_ind_vec[i]==1)
	    {
	      if (rec_ind_count>=start_val && rec_ind_count<end_val)
		{
		  recipient_ind_vec[i]=1;		
		  rec_ind_topaint_count=rec_ind_topaint_count+1;
		}
	      rec_ind_count=rec_ind_count+1;
	    }
	}
    }

  if (haploid_ind==1) printf("Assuming all inds are haploid....\n");
  if (unlinked_ind==1) 
    {
      printf("Assuming sites are unlinked....\n");
      //EMruns=0;
    }
  //printf(" Number of donor haplotypes = %d\n Number of recipient haplotypes = %d\n",ndonors,(2-haploid_ind)*nind-ndonors);
  printf(" Number of donor haplotypes = %d\n Number of recipient haplotypes = %d\n",ndonors,cond_nhaps);
  if (ne_find==1) printf(" Number of EM-runs = %d\n Number of samples = %d\n N_e value = %lf\n Region size = %lf\n",EMruns,samplesTOT,N_e,region_size);
  if (ne_find==0) printf(" Number of EM-runs = %d\n Number of samples = %d\n N_e value = %lf (divided by number of donor haplotypes)\n Region size = %lf\n",EMruns,samplesTOT,N_e,region_size);
  printf(" Global mutation value = %lf\n",GlobalMutRate);
  printf(" Painting %d of %d recipient individuals....\n",rec_ind_topaint_count,cond_nind);

  if (ne_find==1) fprintf(fout, "EM_iter = %d (N_e = %d / copy_prop = %d / mutation = %d / mutationGLOBAL = %d), nsamples = %d, N_e_start = %lf, region_size = %lf, haplotype dataset = %s, genmap dataset = %s, pop-list dataset = %s, id-label dataset = %s\n", EMruns, recom_em_find, copy_prop_em_find, mutation_em_find, mutationALL_em_find, samplesTOT, N_e, region_size, filenameGEN,filename,filenameDONORLIST,filenameID);
   if (ne_find==0) fprintf(fout, "EM_iter = %d (N_e = %d / copy_prop = %d / mutation = %d / mutationGLOBAL = %d), nsamples = %d, N_e_start = %lf (divided by number of donor haplotypes), region_size = %lf, haplotype dataset = %s, genmap dataset = %s, pop-list dataset = %s, id-label dataset = %s\n", EMruns, recom_em_find, copy_prop_em_find, mutation_em_find, mutationALL_em_find, samplesTOT, N_e, region_size, filenameGEN,filename,filenameDONORLIST,filenameID);

              // (i) GET SAMPLES, RUN E-M:
                  // open recomb map file:
  double * recom_map = malloc((Data->nsnps - 1) * sizeof(double));
  if ((unlinked_ind==1) && (recom_find==0))
    {
      for (j=0; j < (Data->nsnps-1); j++) recom_map[j]=-9.0;
    }
  if (recom_find==1)
    {
      fd2 = fopen(filename,"r");
      if (fd2 == NULL) { printf("error opening recom map input file: %s\n",filename); exit(1);}
      fgets(line,2047,fd2);   // header
      for (j=0; j < (Data->nsnps-1); j++)
	{
	  fgets(line,2047,fd2);
	  step=line;
	  reading(&step,"%lf",&bpval);    // basepair position
	  if (bpval != Data->positions[j])
	    {
	      printf("basepair positions do not match between %s and %s at basepair %d. Exiting....\n",filename,filenameGEN,j+1);
	      exit(1);
	    }
	  reading(&step,"%lf",&recom_map[j]);
	  if (recom_map[j] >= 0 && recom_map[j] <= small_recom_val)
	    {
	      printf("recom rate very low at basepair %lf (%lf). Assuming recomb rate between this snp and next one is %lf....\n",Data->positions[j],recom_map[j],small_recom_val);
	      recom_map[j]=small_recom_val;
	    }
	  if (recom_map[j]<0)
	    {
	      printf("recom rate < 0 at basepair %lf. Assuming transition probability of 1 between this snp and next one....\n",Data->positions[j]);
	      //printf("recom rate must be > 0 (basepair %lf)!! Exiting....\n",Data->positions[j]);
	      //exit(1);
	    }
	}
      fgets(line,2047,fd2);
      step=line;
      reading(&step,"%lf",&bpval);    // basepair position
      if (bpval != Data->positions[(Data->nsnps-1)])
	{
	  printf("basepair positions do not match between %s and %s at basepair %d. Exiting....\n",filename,filenameGEN,Data->nsnps);
	  exit(1);
	}
      fclose(fd2);
    }
           // check ordering of SNPs (only allowed to be less than previous position if recom_map<0 at position -- i.e. suggesting new chromosome):
  for (i=0; i < Data->nsnps; i++)
    {
      if (i > 0)
	{
	  if (Data->positions[i] <= Data->positions[(i-1)] && (recom_map[(i-1)]>=0))
	    {
	      printf("positions in %s are not listed in increasing order at basepairs %lf and %lf. Exiting....\n",filenameGEN,Data->positions[(i-1)],Data->positions[i]);
	      exit(1);
	    }
	}
    }
  DestroyData(Data,2-haploid_ind);
  fclose(fd);

  EMruns = EMruns + 1;

                 /* print-out headers for copy-props, chunk counts, lengths, and differences: */
  fprintf(fout2,"Recipient");
  fprintf(fout4,"Recipient");
  fprintf(fout5,"Recipient");
  fprintf(fout6,"Recipient");
  fprintf(fout7,"Recipient num.regions");
  fprintf(fout8,"Recipient num.regions");
  if (all_versus_all_ind==0)
    {
      for (k=0; k < ndonorpops; k++)
	{
	  fprintf(fout2," %s",donor_pop_vec[k]);
	  fprintf(fout4," %s",donor_pop_vec[k]);
	  fprintf(fout5," %s",donor_pop_vec[k]);
	  fprintf(fout6," %s",donor_pop_vec[k]);
	  fprintf(fout7," %s",donor_pop_vec[k]);
	  fprintf(fout8," %s",donor_pop_vec[k]);
	}
    }
  if (all_versus_all_ind==1)
    {
      for (j=0; j < nind_tot; j++)
	{
	  if (include_ind_vec[j]==1)
	    {
	      fprintf(fout2," %s",ind_label_vec[j]);
	      fprintf(fout4," %s",ind_label_vec[j]);
	      fprintf(fout5," %s",ind_label_vec[j]);
	      fprintf(fout6," %s",ind_label_vec[j]);
	      fprintf(fout7," %s",ind_label_vec[j]);
	      fprintf(fout8," %s",ind_label_vec[j]);
	    }
	}
    }
  fprintf(fout2,"\n");
  fprintf(fout4,"\n");
  fprintf(fout5,"\n");
  fprintf(fout6,"\n");
  fprintf(fout7,"\n");
  fprintf(fout8,"\n");

  log_lik_check = loglik(nind_tot, ndonors, &nsites, (2-haploid_ind)*nind, N_e, recom_map, MutProb_vec_tot, samplesTOT, ndonorpops, ndonorhaps_tot, include_ind_vec, ind_label_vec, pop_label_vec, donor_pop_vec, nrecpops, rec_pop_vec, copy_prob_tot, copy_probSTART_tot, pop_vec_tot, region_size, EMruns, copy_prop_em_find, recom_em_find, ne_find, mutation_em_find, mutationALL_em_find, all_versus_all_ind, prior_donor_probs_ind, rec_ind_count, recipient_ind_vec, filenameGEN, donorlist_find, haploid_ind, unlinked_ind, print_file9_ind, drift_calc_ind, num_rec_drift, rec_ind_topaint_count, fout, fout2, fout3, fout4, fout5, fout6, fout7, fout8, fout9);
  if (log_lik_check != 1)
    {
      printf("Algorithm failed. Check input files and parameters. Exiting....\n");
      exit(1);
    }

  for (i=0; i < nind_tot; i++)
    {
      free(ind_label_vec[i]);
      free(pop_label_vec[i]);
    }
  for (i=0; i < ndonorpops; i++)
    free(donor_pop_vec[i]);
  for (i=0; i < nrecpops; i++)
    free(rec_pop_vec[i]);
  free(ind_label_vec);
  free(pop_label_vec);
  free(donor_pop_vec);
  free(rec_pop_vec);
  free(include_ind_vec);
  free(recipient_ind_vec);
  free(recom_map);
  free(filename);
  free(filenameGEN);
  free(filenameDONORLIST);
  free(filenameID);
  free(copy_prob_tot);
  free(copy_probSTART_tot);
  free(pop_vec_tot);
  free(ndonorhaps_tot);
  free(nrechaps);
  free(ndonorprobs);
  free(ndonormutrates);
  free(MutProb_vec_tot);
  //free(bigline);

  fclose(fout);
  fclose(fout2);
  fclose(fout3);
  fclose(fout4);
  fclose(fout5);
  fclose(fout6);
  fclose(fout7);
  fclose(fout8);
  gzclose(fout9);

  printf("Finished!\n");

  return 0;
}
