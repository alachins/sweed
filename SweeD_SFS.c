/*  
 *  SweeD
 *
 *  Copyright September 2012 by Nikolaos Alachiotis and Pavlos Pavlidis
 *
 *  This program is free software; you may redistribute it and/or modify its
 *  under the terms of the GNU General Public License as published by the Free
 *  Software Foundation; either version 2 of the License, or (at your option)
 *  any later version.
 *
 *  This program is distributed in the hope that it will be useful, but
 *  WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 *  or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
 *  for more details.
 *
 *  For any other enquiries send an email to
 *  Pavlos Pavlidis (pavlidisp@gmail.com) or
 *  Nikolaos Alachiotis (n.alachiotis@gmail.com)
 *  
 */

#include "SweeD.h"


int initSFSfromInput (alignment_struct * alignment)
{
	int i, j, ind, params=0;
	t_sfs sum = 0.0;
	t_sfs accumSFS = 0.0;

	alignment->startSFS = 1;
	alignment->endSFS = alignment->sequences;

	for(i=0;i<alignment->segsites;i++)
	{
		if(alignment->positionsInd[i]!=-1)
		{			
			ind = alignment->x[i];
			assert(ind < alignment->maxx);
		
			if(alignment->n[i] == alignment->sequences)
				alignment->SFS[ind] = alignment->SFS[ind] + 1.0;
			else
				for(j=alignment->sequences - alignment->n[i]; j >= 0 ; --j)
					alignment->SFS[ ind + j ] = alignment->SFS[ind + j] + 1.0;

			if(alignment->folded[i] == 1)
			{
				ind = alignment->n[i] - alignment->x[i];
				assert(ind < alignment->maxx);

				if(alignment->n[i] == alignment->sequences)
					alignment->SFS[ind] = alignment->SFS[ind] + 1.0;
				else
					for(j=alignment->sequences - alignment->n[i]; j >= 0 ; --j)
						alignment->SFS[ ind + j ] = alignment->SFS[ind + j] + 1.0;
			}
		}
	}

	if(alignment->SFS[0] > 0)
	{
		params = 1;
		alignment->startSFS = 0;
		sum += alignment->SFS[0];
	}

	for(i = 1; i < alignment->sequences; i++)
	{
		++params;
		if(alignment->SFS[i] == 0.0)
			alignment->SFS[i] = 1.0;
		sum += alignment->SFS[i];
	}

	if(alignment->SFS[alignment->sequences] > 0)
	{
		++params;
		alignment->endSFS = alignment->sequences+1;
		sum += alignment->SFS[alignment->sequences];
	}

	for(i = 0; i < alignment->maxx; ++i)
	{
		alignment->SFS[i] /= sum;

		accumSFS += alignment->SFS[i];
	}

	assert(accumSFS<=1.0f + 1e-12);

	return params - 1;	
}

void getRelativeLogSFS (t_sfs * logSFS)
{	
	int i;

	assert(alignment->SFS[alignment->startSFS]>0);

	t_sfs norm = log(alignment->SFS[alignment->startSFS]);	

	for(i=alignment->startSFS+1; i<alignment->endSFS;i++)
	{
		assert(alignment->SFS[i] > 0.0);
		logSFS[i-alignment->startSFS-1] = log(alignment->SFS[i])-norm;
	}
}

void printT_SFSarray(t_sfs *a, int size)
{
	int i=0;
	for(i=0; i<size; ++i)
		printf("%e,", a[i]);
	printf("\n");
}

void getSFSfromLogSFS(const t_sfs * logSFS, t_sfs * SFS)
{
	int i;
	t_sfs sum = 0.;

	for(i=0; i<alignment->maxx; ++i)
		SFS[i] = 0.0;

	for(i=alignment->startSFS+1; i<alignment->endSFS; ++i)
		SFS[i] = exp(logSFS[i-alignment->startSFS-1]);

	for(i=alignment->startSFS+1; i<alignment->endSFS; ++i)
		sum += SFS[i];// = exp(logSFS[i-alignment->startSFS-1]));
	
	SFS[alignment->startSFS] = 1.;
	
	sum += 1.;

	sum = 1.0 / sum;

	for(i=0; i<alignment->maxx; ++i)
		SFS[i] = SFS[i] * sum;
}

t_sfs likelihoodSFS_SNP(t_sfs * tmpSFS, int x, int n, unsigned char f, int c)
{
	int minx, maxx, xi, nmax = alignment->sequences;	
	t_sfs probability = 0.;		
	
	minx = x;
	maxx = minx+nmax-n;

	for(xi=minx; xi<=maxx; ++xi)
		probability += tmpSFS[xi]*exp(XchooseY_ln(xi, x)+XchooseY_ln(nmax-xi, n-x)-XchooseY_ln(nmax,n));

	if(f==1)
	{
		x = n-x;
		if(x!=n-x)
		{
			minx = x;
			maxx = minx+nmax-n;

			for(xi=minx; xi<=maxx; ++xi)
				probability += tmpSFS[xi]*exp(XchooseY_ln(xi, x)+XchooseY_ln(nmax-xi, n-x)-XchooseY_ln(nmax,n));
		}		
	}

	return log(probability)*c;
}

void setPatternData (struct pattern * cur, int x, int n, unsigned char f, int c, struct pattern * nxt)
{
	cur->x=x;
	cur->n=n;
	cur->f=f;
	cur->c=c;
	cur->nxt=nxt;
}

void getPatternData (struct pattern * cur, int * x, int * n, unsigned char * f, int * c)
{
	*x=cur->x;
	*n=cur->n;
	*f=cur->f;
	*c=cur->c;
}

#ifdef _USE_PTHREADS
void likelihoodSFS_thread(int tid, int threads)
{
	int i;
	int size = (int)(alignment->patterns / (float)threads);
	int start = tid * size;
	int stop = start + size;

	if(tid==threads-1 && stop!=alignment->patterns)
		stop = alignment->patterns;

	for(i=start;i<stop;i++)
	{
		alignment->threadScoresLIKSFS[tid] +=  likelihoodSFS_SNP(alignment->tmpSFS, alignment->x_pat[i], alignment->n_pat[i], alignment->f_pat[i], alignment->c_pat[i]);
	}
}

t_sfs parallelLikelihoodSFS()
{
	t_sfs likelihood = 0.0;
	int i, threads=threadDataL[0].threadTOTAL;

	for(i=0;i<threads;i++)
		alignment->threadScoresLIKSFS[i]=0.0;

	for(i=1;i<threads;i++)
		threadDataL[i].threadOPERATION = LIKSFS;

	likelihoodSFS_thread(0,threads);

	threadDataL[0].threadBARRIER=1;

	syncThreadsBARRIER();

	for(i=0;i<threads;i++)
		likelihood += alignment->threadScoresLIKSFS[i];

	return likelihood;	
}

#endif
       
t_sfs getLikelihoodSFS(t_sfs *tmpSFS)
{
	t_sfs likelihood=0.0;

#ifdef _USE_PTHREADS
	likelihood = parallelLikelihoodSFS();
#else
	int i;
	for(i=0;i<alignment->patterns;i++)
	{
		likelihood +=  likelihoodSFS_SNP(tmpSFS, alignment->x_pat[i], alignment->n_pat[i], alignment->f_pat[i], alignment->c_pat[i]);
	}
#endif
	return likelihood;
}

t_sfs p_likelihood_freq(const t_sfs * logSFS)
{
	t_sfs likelihood;
  
	getSFSfromLogSFS (logSFS, alignment->tmpSFS);

	likelihood = getLikelihoodSFS(alignment->tmpSFS);

 	return -likelihood;
}

void d_likelihood_freq(const t_sfs * logSFS, t_sfs * outvec)
{
	getgradient(alignment->params, logSFS, alignment->ng, outvec, p_likelihood_freq, alignment->lb, alignment->ub);
}

int stopPatternSearch (int xpat, int npat, unsigned char fpat, int x, int n, unsigned char f)
{
	if(xpat<x)
		return 1;

	if(npat<n)
		return 1;
	
	if(fpat<f)
		return 1;

	return 0;
}

void createPatterns ()
{
	int i=0, x, n, c, match;
	unsigned char f;		

	struct pattern * searchPattern, *lastPattern, *nxtPattern;

	alignment->patternList = malloc(sizeof(struct pattern));
	assert(alignment->patternList!=NULL);

	setPatternData(alignment->patternList, alignment->x[0], alignment->n[0], alignment->folded[0], 1, NULL);

	alignment->patterns=1;

	for(i=1;i<alignment->segsites;i++)
	{
		lastPattern = alignment->patternList;
		searchPattern = alignment->patternList;
		match = 0;

		while(searchPattern!=NULL)
		{
			getPatternData(searchPattern, &x, &n, &f, &c);

			if(alignment->x[i]==x && alignment->n[i]==n && alignment->folded[i]==f)
			{
					searchPattern->c++;
					match = 1;
					break;
			}

			if(stopPatternSearch(x,n,f,alignment->x[i],alignment->n[i],alignment->folded[i]))
				break;				

			lastPattern = searchPattern;
			searchPattern = searchPattern->nxt;
		}

		if(match==0)
		{
			struct pattern * newPattern = malloc(sizeof(struct pattern));
			
			if(searchPattern==alignment->patternList)
			{				
				nxtPattern = alignment->patternList;
				alignment->patternList = newPattern;
			}
			else
			{
				nxtPattern = lastPattern->nxt;
				lastPattern->nxt = newPattern;
			}

			setPatternData(newPattern, alignment->x[i], alignment->n[i], alignment->folded[i], 1, nxtPattern);

			alignment->patterns++;
		}
	}

	alignment->n_pat = malloc(sizeof(al_t)*alignment->patterns);
	assert(alignment->n_pat!=NULL);
	alignment->x_pat = malloc(sizeof(al_t)*alignment->patterns);
	assert(alignment->x_pat!=NULL);
	alignment->f_pat = malloc(sizeof(unsigned char)*alignment->patterns);
	assert(alignment->f_pat!=NULL);
	alignment->c_pat = malloc(sizeof(int)*alignment->patterns);
	assert(alignment->c_pat!=NULL);

	searchPattern = alignment->patternList;
	i=0;

	while(searchPattern!=NULL)
	{
		getPatternData(searchPattern, &x, &n, &f, &c);

		alignment->n_pat[i]=n;
		alignment->x_pat[i]=x;
		alignment->f_pat[i]=f;
		alignment->c_pat[i]=c;
		
		i++;
		searchPattern = searchPattern->nxt;
	}
}

void getFinalSFS (t_sfs * logSFS)
{
	int i;
	double sum=0.0;
		
	for (i=0; i<alignment->maxx; i++)
		alignment->SFS[i]=0.0;

	for (i=0; i<alignment->params; i++)
		sum += (alignment->SFS[i+alignment->startSFS+1] = exp(logSFS[i]));

	sum++;

	alignment->SFS[alignment->startSFS]=1.0;

	for (i=0; i<alignment->maxx; i++)
		alignment->SFS[i] /= sum;

}

int checkSFSind(char * SFSindS, int j)
{
	int i, length = strlen(SFSindS);

	for(i=0;i<length;i++)
	{
		if(SFSindS[i]<48 || SFSindS[i]>57)
		{
			if(SFSindS[i]=='/')
			{	
				if(SFSindS[i+1]=='/')	
					return -1;
				else
				{
					fprintf(stderr, "\n ERROR: Expected // in the input SFS file at the end of the SFS.\n\n");
					exit(0);
				}
			}
	
			fprintf(stderr, "\n ERROR: Invalid character (%c - ASCII: %d) in the input SFS file. %s is not a positive integer. Expected %d.\n\n",SFSindS[i],SFSindS[i],SFSindS,j);
			exit(0);
		}
	}

	return atoi(SFSindS);
}

void computeFactLUT()
{
	int i,scale;
	double tmp;
	factLUT = malloc(sizeof(factorial_LUT));
	assert(factLUT!=NULL);

	factLUT-> val = malloc(sizeof(double)*(alignment->maxn+1));
	assert(factLUT->val!=NULL);

	factLUT-> sc = calloc((alignment->maxn+1),sizeof(int));
	assert(factLUT->sc!=NULL);

	factLUT-> val[0]=0;
	factLUT-> val[1]=1;

	for(i=2;i<alignment->maxn+1;i++)
	{
		scale = factLUT->sc[i-1];
		tmp = factLUT->val[i-1] * i;

		if(tmp>1e6)
		{
			tmp = tmp * 1e-6;
			++scale;				
		}
		
		factLUT->val[i] = tmp;
		factLUT->sc[i] = scale;  
	}

	for(i=1;i<alignment->maxn+1;i++)
		factLUT->val[i] = (factLUT->sc[i]*log(1e6) + log(factLUT->val[i]));

	free(factLUT->sc);
}



void createSFS (FILE * fpSFS, FILE * fpSFSo, int alignmentIndex, int * SFSsize)
{
	int i;
	static int allocSFS = 1;

	if(allocSFS==1)
	{
		alignment->SFS = calloc(alignment->sequences+1, sizeof(t_sfs));
		assert(alignment->SFS!=NULL);
		allocSFS=0;
	}

	if(alignment->sequences+1> (*SFSsize) )
	{
		alignment->SFS = realloc(alignment->SFS, (alignment->sequences+1)*sizeof(t_sfs));
		assert(alignment->SFS!=NULL);
		(*SFSsize) = alignment->sequences+1;
	}

	if(alignment->userSFS==0)
	{

		alignment->params = initSFSfromInput (alignment);

		createPatterns ();

		alignment->tmpSFS = malloc((alignment->maxx)*sizeof(t_sfs));
		assert(alignment->tmpSFS!=NULL);

		alignment->lb=malloc(alignment->params*sizeof(t_sfs));
		assert(alignment->lb!=NULL);

	  	alignment->ub=malloc(alignment->params*sizeof(t_sfs));
		assert(alignment->ub!=NULL);

		alignment->ng = malloc(alignment->params*sizeof(int));
		assert(alignment->ng!=NULL);

		for (i=0; i<alignment->params; i++)
		{
			alignment->lb[i]=-15.0;
			alignment->ub[i]=10.0;
			alignment->ng[i] = 1;
		}

		assert(alignment->params == alignment->endSFS - alignment->startSFS-1);	

		t_sfs * logSFS = malloc(alignment->params*sizeof(t_sfs));
		assert(logSFS!=NULL);

		getRelativeLogSFS (logSFS);

		for (i=0; i<alignment->params; i++)
		{
			if (logSFS[i] < alignment->lb[i]) 
				logSFS[i]=alignment->lb[i];

			if (logSFS[i] > alignment->ub[i])
				logSFS[i]=alignment->ub[i];
		}

		findmax_bfgs(alignment->params, logSFS, p_likelihood_freq, d_likelihood_freq,
		     alignment->lb, alignment->ub, NULL, -1);

		getFinalSFS (logSFS);

		free(logSFS);
		free(alignment->tmpSFS);
		free(alignment->lb);
		free(alignment->ub);
		free(alignment->ng);
		logSFS = NULL;
		alignment->tmpSFS = NULL;
		alignment->lb = NULL;
		alignment->ub = NULL;
		alignment->ng = NULL;

		char SFSvalS[50];
		
		if(alignmentIndex>1 && fpSFSo!=NULL)
			fprintf(fpSFSo,"%d\n",alignmentIndex);

		for(i=0; i<alignment->sequences + 1; ++i)
		{
			sprintf(SFSvalS, "%e", alignment->SFS[i]);
			alignment->SFS[i] = atof(SFSvalS);

			if(fpSFSo!=NULL)
				fprintf(fpSFSo,"%d\t%e\n",i,alignment->SFS[i]);
		}

		if(fpSFSo!=NULL)
			fprintf(fpSFSo,"//");
	}
	else
	{	

	  int *indexes = calloc( alignment->sequences+1, sizeof(int));

	  double sumSFS = 0.;

		alignment->startSFS = 1;
		alignment->endSFS = alignment->sequences;
	
		int index, counter=0;
		char SFSvalS[50]; char SFSindS[50], t;
		
		if(alignmentIndex>1)
		{
			while((t=fgetc(fpSFS))!=EOF)
			{
				if(t=='/')
					counter++;

				if(counter==2)
					break;
			}

			while((t=fgetc(fpSFS))!=EOF)
			{
				if(t=='\n')
					break;
			}			
		}
			

		for(i=0;i<alignment->sequences+1;i++)
		{
			if(fscanf(fpSFS,"%s",SFSindS)==-1)
				break;

		
			index = checkSFSind(SFSindS, i);

			indexes[i] = index;

			if(index==-1)
				break;
			
			if(index!=i)
			{
				fprintf(stderr, "\n ERROR: Wrong index (%d - expected: %d) in the input SFS file.\n\n",index, i);
				exit(0);
			}

			if(fscanf(fpSFS,"%s",SFSvalS)==-1)
			{
				fprintf(stderr, "\n ERROR: Missing SFS value for index %d.\n\n",index);
				exit(0);
			}

			alignment->SFS[index] = atof(SFSvalS);

			sumSFS += alignment->SFS[index];

			if(alignment->SFS[index]<0.0 || alignment->SFS[index]>1.0 )
			{
				fprintf(stderr, "\n ERROR: Invalid SFS value (%e) for index %d in the input SFS file.\n\n",alignment->SFS[index],index);
				exit(0);
			}
		}

		for(i = 0; i < alignment->sequences + 1; ++i)
		  {
		    alignment->SFS[ indexes[i] ] /= sumSFS;
		  }


		if(i<alignment->endSFS)
		{
			if(alignmentIndex==1 || (alignmentIndex>1 && i!=0))
			{
				fprintf(stderr, "\n ERROR: Wrong SFS size (%d - expected: %d).\n\n",i,alignment->endSFS);
				exit(0);
			}
		}		

		if(alignment->SFS[0] > 0)
		{
		  /* TODO why does the next assert exist? */
		  /* if(alignment->startSFS!=0) */
		  /* 	assert(0); */

			alignment->startSFS = 0;
		}

		if(alignment->SFS[alignment->sequences] > 0)
		{
			if(alignment->endSFS!=alignment->sequences+1)
				assert(0);

			alignment->endSFS = alignment->sequences+1;
		}

		if(fpSFSo!=NULL)
		{
			if(alignmentIndex>1)
				fprintf(fpSFSo,"%d\n",alignmentIndex);

			for(i=0; i<alignment->sequences + 1; ++i)
				fprintf(fpSFSo,"%d\t%e\n",i,alignment->SFS[i]);
	
			fprintf(fpSFSo,"//");
		}
		
		free(indexes);
	}
}




#ifdef _ANALYTICAL_SFS

void init_factorialsLUT(mpfr_t * fact, unsigned long int n)
{
	unsigned long int i = 0;

	mpfr_init2(fact[0], ACC);

	mpfr_set_ui(fact[0], 1, GMP_RNDU);


	for(i = 1; i<n+1; ++i)
	{
		mpfr_init2(fact[i], ACC);
		mpfr_mul_ui(fact[i], fact[i-1], i, GMP_RNDU);
	}
}

int is_even(unsigned long int n)
{
	if( (n & 1) == 0)
		return 1;

	return 0;
}


void pochhammer(mpfr_t poch, mpfr_t* factLUT, unsigned long int n, unsigned long int k)
{
	assert(n>=0);
	assert(k>=0);

	//mpfr_init2(poch, ACC);
	mpfr_set(poch, factLUT[n+k-1], GMP_RNDU);
	mpfr_div(poch, poch, factLUT[n-1], GMP_RNDU);
		
}

void binomial(mpfr_t bin, mpfr_t* factLUT, unsigned long int n, unsigned long int k)
{

	if(k > n)
		mpfr_set_ui(bin, 0, GMP_RNDU);
	else if(k == n)
		mpfr_set_ui(bin, 1, GMP_RNDU);
	else if(k < n)
	{
		mpfr_set(bin, factLUT[n], GMP_RNDU);
		mpfr_div(bin, bin, factLUT[k], GMP_RNDU);
		mpfr_div(bin, bin, factLUT[n-k], GMP_RNDU);
	}
}

void get_c (mpfr_t c, mpfr_t* factLUT, unsigned long int n, unsigned long int k)
{
	assert(n >= k);
	mpfr_t bin_nk, poch_nk, poch_kk;

	mpfr_init2(bin_nk, ACC);
	mpfr_init2(poch_nk,ACC);
	mpfr_init2(poch_kk,ACC);
	binomial(bin_nk, factLUT, n, k);
	pochhammer(poch_kk, factLUT, k, k);
	pochhammer(poch_nk, factLUT, n, k);	
	
	//mpfr_init2(c, ACC);
	
	mpfr_set(c, bin_nk, GMP_RNDU);
	mpfr_mul(c, c, poch_kk, GMP_RNDU);
	mpfr_div(c, c, poch_nk, GMP_RNDU);	

	mpfr_clear(bin_nk);
	mpfr_clear(poch_nk);
	mpfr_clear(poch_kk);
}

void get_r (mpfr_t r, mpfr_t* factLUT, unsigned long int k, unsigned long int j)
{
	assert(j <= k);
	long int w = 1;
 	 
  	if(is_even(k-j) == 0)
    		w = -1;
 	
	mpfr_t bin_kj, poch_jkminus1, poch_kkminus1;

	mpfr_init2(bin_kj, ACC);
	mpfr_init2(poch_jkminus1, ACC);
	mpfr_init2(poch_kkminus1, ACC);

	binomial(bin_kj, factLUT, k, j);
	pochhammer(poch_jkminus1, factLUT, j, k-1);
	pochhammer(poch_kkminus1, factLUT, k, k-1);

	//mpfr_init2(r, ACC);
	
	mpfr_set(r, bin_kj, GMP_RNDU);
	mpfr_mul(r, r, poch_jkminus1, GMP_RNDU);
	mpfr_div(r, r, poch_kkminus1, GMP_RNDU);
	
	mpfr_mul_si(r, r, w, GMP_RNDU);

	mpfr_clear(bin_kj);
	mpfr_clear(poch_jkminus1);
	mpfr_clear(poch_kkminus1);
}

void setEvent (event * eventList, int eventIndex, double t, double d)
{
	eventList[eventIndex].t=t;
	eventList[eventIndex].d=d;	
}


void get_ETExpo_summand(mpfr_t s, mpfr_t bin_k2, mpfr_t pieceVal)
{
  
  mpfr_mul(s, bin_k2, pieceVal, GMP_RNDU);
  mpfr_mul_si(s, s, -1, GMP_RNDU);
  mpfr_exp(s, s, GMP_RNDU);
  
}

void get_ETExpo_fraction(mpfr_t s, mpfr_t bin_k2, mpfr_t c, mpfr_t r, unsigned long int n, unsigned long int k, unsigned long int j, double R, int eventsTotal, event *eventList, int storeLUT)
{

  assert(eventsTotal > 1);
  
  if(mpfr_cmp_d(analyticalSFS_LUT[k][j], 0.0) != 0)
    {
      mpfr_set(s, analyticalSFS_LUT[k][j], GMP_RNDU);
      return;
    }
  
  int i;
  
  mpfr_t t, pieceVal, pieceVal_n, tmp, sExpEvent, xmp, EImp;
  mpfr_inits2(ACC, t, pieceVal, pieceVal_n, tmp, sExpEvent, xmp, EImp, (mpfr_ptr)0);
  
  
  mpfr_set_d(tmp, R, GMP_RNDU);
  mpfr_mul_d(tmp, tmp, eventList[1].d, GMP_RNDU);
  
  mpfr_set_d(pieceVal, 1, GMP_RNDU);
  
  mpfr_sub_d(pieceVal, pieceVal, eventList[1].d, GMP_RNDU);
  
  mpfr_div(pieceVal, pieceVal, tmp, GMP_RNDU);
  
  get_ETExpo_summand(t, bin_k2, pieceVal);
  
  mpfr_mul_d(t, t, eventList[1].d, GMP_RNDU);

  mpfr_set(s, t, GMP_RNDU);

  /* printf("first: "); */
  /* mpfr_out_str(NULL, 10, 10, t, GMP_RNDU); */
  /* putchar('\t'); */

    
  for(i=2; i<eventsTotal; ++i)
    {
      
      mpfr_set_d(pieceVal_n, eventList[i].t, GMP_RNDU);
      mpfr_sub_d(pieceVal_n, pieceVal_n, eventList[i-1].t, GMP_RNDU);
      //printf("%.20f %d\n", eventList[i].t, eventsTotal);
      mpfr_div_d(pieceVal_n, pieceVal_n, eventList[i-1].d, GMP_RNDU);
      mpfr_add(pieceVal, pieceVal, pieceVal_n, GMP_RNDU);
      
      get_ETExpo_summand(t, bin_k2, pieceVal);
      
      mpfr_set_d(tmp, eventList[i].d, GMP_RNDU);
      mpfr_sub_d(tmp, tmp, eventList[i-1].d, GMP_RNDU);
      
      mpfr_mul(t, t, tmp, GMP_RNDU);

      /* printf("%ld %ld %ld %d\t", n, k, j, i); */
      /* mpfr_out_str(NULL, 10, 10, t, GMP_RNDU); */
      /* putchar('\t'); */

      mpfr_add(s, s, t, GMP_RNDU);

    }
  
    
  mpfr_div(s, s, bin_k2, GMP_RNDU);

  
  /* mpfr_out_str(NULL, 10, 10, s, GMP_RNDU); */
  /* putchar('\t'); */
  
  //mpfr_eint(sExpEvent, tmp, GMP_RNDU);
  
  /******* SET THE EXPONENTIAL INTEGRAL *********/
  /* long double bin_k2_ld = mpfr_get_ld(bin_k2, GMP_RNDU); */
  /* long double x = (-bin_k2_ld)/R; */
  /* printf("x: %Le\t", x); */
  /* long double EI = -xExponential_Integral_Ei( x ); */
  /* printf("\nEI: %Le\t", EI); */

  /* set the x for the Exponential Integral */
  mpfr_set(xmp, bin_k2, GMP_RNDU);
  mpfr_div_d(xmp, xmp, R, GMP_RNDU);
  mpfr_mul_si(xmp, xmp, -1, GMP_RNDU);
  
  mpfr_Exponential_Integral_Ei(EImp, xmp);
  mpfr_mul_d(EImp, EImp, -1.0, GMP_RNDU);
  

  
  /* printf("k: %lu, d: %f, R: %f, mpfr1: ", k, eventList[1].d, R);  */
  /* mpfr_out_str(NULL, 10, 10, EImp, GMP_RNDU); */
  /* putchar('\n'); */

  mpfr_div_d(xmp, xmp, eventList[1].d, GMP_RNDU);
  mpfr_Exponential_Integral_Ei(sExpEvent, xmp);

  /* printf("k: %lu, d: %f, R: %f, mpfr2: ", k, eventList[1].d, R);  */
  /* mpfr_out_str(NULL, 10, 10, sExpEvent, GMP_RNDU); */
  /* putchar('\n'); */


  mpfr_add(sExpEvent, sExpEvent, EImp, GMP_RNDU);

  
  /* printf("k: %lu, d: %f, R: %f, mpfr: ", k, eventList[1].d, R);  */
  /* mpfr_out_str(NULL, 10, 10, sExpEvent, GMP_RNDU); */
  /* putchar('\n'); */
  /* putchar('\n'); */
  
  
  /* EI += xExponential_Integral_Ei( x/eventList[1].d ); */
  /* printf("EI2: %Le %Le\n", xExponential_Integral_Ei( x/eventList[1].d ), EI); */
  /* mpfr_set_ld(sExpEvent, EI, GMP_RNDU); */
  /********************************************/

  /********** set the exponent *************/
  mpfr_set(tmp, bin_k2, GMP_RNDU);
  mpfr_div_d(tmp, tmp, R, GMP_RNDU);
  mpfr_exp(t, tmp, GMP_RNDU);
  
  /****************************************/
  /* printf("tmp2: "); */
  /* mpfr_out_str(NULL, 10, 10, tmp, GMP_RNDU); */
  /* putchar('\n'); */
  
  mpfr_mul(sExpEvent, sExpEvent, t, GMP_RNDU);
  mpfr_div_d(sExpEvent, sExpEvent, R, GMP_RNDU);

  
  /* printf("first: "); */
  /* mpfr_out_str(NULL, 10, 10, sExpEvent, GMP_RNDU); */
  /* putchar('\n'); */
  
  mpfr_add(s, s, sExpEvent, GMP_RNDU);

  /* if(k == 5 && j==3) */
  /*   { */
  /*     printf("%lu %lu %lu\t", n, k, j); */
  /*     mpfr_out_str(NULL, 10, 20, s, GMP_RNDU); */
  /*   } */
  mpfr_mul(s, s, c, GMP_RNDU); 
  mpfr_mul(s, s, r, GMP_RNDU);
  
  /* mpfr_out_str(NULL, 10, 10, s, GMP_RNDU); */
  /* putchar('\n'); */
  assert(k<=n);
  assert(j<=n);
  
  if(storeLUT==1)
    {
      mpfr_set(analyticalSFS_LUT[k][j], s, GMP_RNDU);
    }
  
  mpfr_clears(t, pieceVal, pieceVal_n, tmp, sExpEvent, xmp, EImp, (mpfr_ptr)0);
  
}


void get_ETpiecewise_summand(mpfr_t s, mpfr_t bin_k2, unsigned int eventIndex, event * eventList)
{

	if(eventIndex==0)
	{
		mpfr_set_ui(s, 1, GMP_RNDU);
		return;
	}

	int i;
	mpfr_t t, tdtot;
	mpfr_init2(t, ACC);
	mpfr_init2(tdtot, ACC);
	
	mpfr_set_d(tdtot, 0., GMP_RNDU);
	for(i=1;i<eventIndex+1;i++)
	{
		mpfr_set_d(t, eventList[i].t, GMP_RNDU);
		mpfr_sub_d(t, t, eventList[i-1].t, GMP_RNDU);
		mpfr_div_d(t, t, eventList[i-1].d, GMP_RNDU);
		
		mpfr_add(tdtot, tdtot, t, GMP_RNDU);
			
	}

	mpfr_mul(s, tdtot, bin_k2, GMP_RNDU);
	mpfr_mul_si(s, s, -1, GMP_RNDU);
	mpfr_exp(s, s, GMP_RNDU);

	mpfr_mul_d(s, s, (eventList[eventIndex].d-eventList[eventIndex-1].d), GMP_RNDU);

	mpfr_clear(t);
	mpfr_clear(tdtot);

}

void get_ETpiecewise_fraction (mpfr_t f, mpfr_t* factLUT, unsigned long int n, unsigned long int k, unsigned long int j, event * eventList, int eventsTotal, int storeLUT)
{
	
   
  if(mpfr_cmp_d(analyticalSFS_LUT[k][j], 0.0) != 0)
    {
      mpfr_set(f, analyticalSFS_LUT[k][j], GMP_RNDU);
      return;
    }

	int i;

	mpfr_t bin_k2, s, c, r;
	mpfr_init2(bin_k2, ACC);
	mpfr_init2(s, ACC);
	mpfr_init2(c, ACC);
	mpfr_init2(r, ACC);


	mpfr_set_d(f, 0., GMP_RNDU);
	binomial(bin_k2, factLUT, k, 2);

	for(i=0;i<eventsTotal;i++)
	{
		get_ETpiecewise_summand(s, bin_k2, i, eventList);
		mpfr_add(f, f, s, GMP_RNDU);
	}

	mpfr_div(f, f, bin_k2, GMP_RNDU);
	get_c(c, factLUT, n, k);
	get_r(r, factLUT, k, j);
	mpfr_mul(f, f, c, GMP_RNDU);
	mpfr_mul(f, f, r, GMP_RNDU);

	assert(k<=n);
	assert(j<=n);
	
	if(storeLUT==1)
	  {
	    mpfr_set(analyticalSFS_LUT[k][j], f, GMP_RNDU);
	  }


	mpfr_clear(bin_k2);
	mpfr_clear(s);
	mpfr_clear(c);
	mpfr_clear(r);
}


void ETpiecewise (mpfr_t ETp, mpfr_t* factLUT, unsigned long int n, unsigned long int j, event* eventList, int eventsTotal, int storeLUT)
{
	unsigned long int k;
	mpfr_t ETfrac;	
	mpfr_init2(ETfrac, ACC);
	
	
	//mpfr_init2(ETp, ACC);
	mpfr_set_d(ETp, 0., GMP_RNDU);

	for(k=j;k<n+1; ++k)
	{
		get_ETpiecewise_fraction (ETfrac, factLUT, n, k, j, eventList, eventsTotal, storeLUT);
		
		//an approximation
		if( mpfr_cmp_d(ETfrac, APPROX) < 0 && mpfr_cmp_d(ETfrac, -APPROX) > 0)
		  {
		    break;
		  }
		
		mpfr_add(ETp, ETp, ETfrac, GMP_RNDU);
	}

	/* mpfr_out_str( NULL, 10, 100, ETp, GMP_RNDU); */
	/* printf("\n"); */
	mpfr_clear(ETfrac);
}


void ETExpo (mpfr_t ETe, mpfr_t* factLUT, unsigned long int n, unsigned long int j, double R, event* eventList, int eventsTotal, int storeLUT)
{

  
	unsigned long int k;
	mpfr_t ETfrac, c, r, bin_k2;	
	mpfr_inits2(ACC, ETfrac, c, r, bin_k2, (mpfr_ptr) 0 );
		
	
	/* mpfr_init2(ETe, ACC); */
	/* mpfr_set_d(ETe, 0., GMP_RNDU); */

	for(k=j;k<n+1; ++k)
	  {
	    get_c(c, factLUT, n, k);
	    get_r(r, factLUT, k, j);

	    

	    binomial(bin_k2, factLUT, k, 2);
	    get_ETExpo_fraction (ETfrac, bin_k2, c, r, n, k, j, R, eventsTotal, eventList, storeLUT);
	    
	    //an approximation
	    if( mpfr_cmp_d(ETfrac, APPROX) < 0 && mpfr_cmp_d(ETfrac, -APPROX) > 0)
	      {
		break;
	      }
	    
	    mpfr_add(ETe, ETe, ETfrac, GMP_RNDU);
	  }
	
	mpfr_clears(ETfrac, c, r, bin_k2, (mpfr_ptr) 0 );
}



void ripiecewise_summand(mpfr_t s, mpfr_t* factLUT, unsigned long int n, unsigned long int k, unsigned long int i,  event* eventList, int eventsTotal, int storeLUT)
{
	mpfr_t t;
	mpfr_init2(t, ACC);

	binomial(s, factLUT, n-i-1, k-2);
	binomial(t, factLUT, n-1, k-1);

	mpfr_div(s, s, t, GMP_RNDU);

	ETpiecewise (t, factLUT, n, k, eventList, eventsTotal, storeLUT);

	mpfr_mul(s, s, t, GMP_RNDU);
	mpfr_mul_ui(s, s, k, GMP_RNDU);
	

	mpfr_clear(t);
}



void riexpo_summand(mpfr_t s, mpfr_t* factLUT, unsigned long int n, unsigned long int k, unsigned long int i,  double R, event* eventList, int eventsTotal, int storeLUT)
{
	mpfr_t t;
	mpfr_init2(t, ACC);
	
	mpfr_set_d(t, 0., GMP_RNDU);
	
	binomial(s, factLUT, n-i-1, k-2);
	binomial(t, factLUT, n-1, k-1);

	mpfr_div(s, s, t, GMP_RNDU);
	
	mpfr_set_d(t, 0., GMP_RNDU);
	ETExpo(t, factLUT, n, k, R, eventList, eventsTotal, storeLUT);
	mpfr_mul(s, s, t, GMP_RNDU);
	mpfr_mul_ui(s, s, k, GMP_RNDU);
	
	mpfr_clear(t);
}


void ripiecewise(mpfr_t rip, mpfr_t* factLUT, unsigned long int n, unsigned long int i,  event* eventList, int eventsTotal, int storeLUT)
{
	unsigned long int k;

	mpfr_t num, den, t;
	mpfr_inits2(ACC, num, den,t, (mpfr_ptr)0);
	
	mpfr_set_ui(num, 0, GMP_RNDU);
	mpfr_set_ui(den, 0, GMP_RNDU);

	for(k=2;k<n-i+1+1;k++)
	{
		ripiecewise_summand(t, factLUT, n, k, i, eventList, eventsTotal, storeLUT);
		mpfr_add(num, num, t, GMP_RNDU);
	}

	for(k=2;k<n+1;k++)
	{
	  /* todo probably here is the reason that we need to initialize t inside the function */
		ETpiecewise (t, factLUT, n, k, eventList, eventsTotal, storeLUT);
		mpfr_mul_ui(t, t, k, GMP_RNDU);
		mpfr_add(den, den, t, GMP_RNDU);
	}

	mpfr_div(rip, num, den, GMP_RNDU);
	
	mpfr_clears( num, den,t, (mpfr_ptr)0);
}


void riexpo(mpfr_t rip, mpfr_t* factLUT, unsigned long int n, unsigned long int i,  double R, event* eventList, int eventsTotal, int storeLUT)
{
	unsigned long int k;

	mpfr_t num, den, t;
	
	mpfr_init2(num, ACC);
	mpfr_init2(den, ACC);
	mpfr_init2(t, ACC);
	mpfr_set_ui(num, 0, GMP_RNDU);
	mpfr_set_ui(den, 0, GMP_RNDU);
	mpfr_set_d(t, 0., GMP_RNDU);
	
	for(k=2;k<n-i+1+1;k++)
	{
	  //double time000 = gettime();
	  //printf("k: %ld\n", k);
	  riexpo_summand(t, factLUT, n, k, i, R, eventList, eventsTotal, storeLUT);
	  mpfr_add(num, num, t, GMP_RNDU);
	  //printf("loop time %.20f\n", gettime() - time000);
	}

	for(k=2;k<n+1;k++)
	{
	  
	  mpfr_set_d(t, 0., GMP_RNDU);
	  ETExpo (t, factLUT, n, k, R, eventList, eventsTotal, storeLUT);
	  mpfr_mul_ui(t, t, k, GMP_RNDU);
	  mpfr_add(den, den, t, GMP_RNDU);
	}
	mpfr_clear(t);
	
	mpfr_div(rip, num, den, GMP_RNDU);
	
	mpfr_clear(num);
	mpfr_clear(den);
	
	
}



void computeAnalyticalSFS(int sequences, FILE * fpSFSo)
{
	unsigned long int n = (unsigned long int) sequences;
	int i,j;
	mpfr_t * factorialsLUT = (mpfr_t*)malloc((2*n+1)*sizeof(mpfr_t));

	init_factorialsLUT(factorialsLUT, 2*n);
	
	analyticalSFS_LUT = malloc( (n+1) * sizeof(mpfr_t*) );
	
	for(i=0; i<n+1; ++i)
	  {
	    //printf("%ld\n", MPFR_PREC_MAX);
	    analyticalSFS_LUT[i] = malloc((n+1) * sizeof(mpfr_t));
	    for(j = 0; j < n+1; ++j)
	      {
		mpfr_init2(analyticalSFS_LUT[i][j], ACC);
		mpfr_set_d(analyticalSFS_LUT[i][j], 0., GMP_RNDU);
	      }
	  }

	alignment->SFS = calloc(sequences+1, sizeof(t_sfs));
	assert(alignment->SFS!=NULL);

	mpfr_t tmp;
	mpfr_init2(tmp, ACC);

	ripiecewise(tmp, factorialsLUT, n, 1, eventList, eventsTotal, 1);
	alignment->SFS[1] = mpfr_get_d(tmp,GMP_RNDU);

	
	for(i=2;i<n;i++)
	{
		ripiecewise(tmp, factorialsLUT, n, i, eventList, eventsTotal, 0);
		alignment->SFS[i] = mpfr_get_d(tmp,GMP_RNDU);
	}

	if(fpSFSo!=NULL)
	{
		for(i=0; i<sequences + 1; ++i)
			fprintf(fpSFSo,"%d\t%e\n",i,alignment->SFS[i]);
	}

	mpfr_clear(tmp);

	
	for(i=0; i<(2*n+1); ++i)
	  mpfr_clear( factorialsLUT[i] );

	free(factorialsLUT);

	for(i=0; i<(n+1); ++i)
	  {
	    for(j=0; j<n+1; ++j)
	      {
		mpfr_clear( analyticalSFS_LUT [i][j] );
	      }
	    free(analyticalSFS_LUT[i]);
	  }

	free(analyticalSFS_LUT);

	mpfr_free_cache();
	
}



void computeAnalyticalSFSExpo(int sequences, double R, FILE * fpSFSo)
{
	unsigned long int n = (unsigned long int) sequences;
	int i,j;
	mpfr_t * factorialsLUT = (mpfr_t*)malloc( (2*n+1) * sizeof(mpfr_t));

	init_factorialsLUT(factorialsLUT, 2*n);

	
	
	analyticalSFS_LUT = malloc( (n+1) * sizeof(mpfr_t*) );
	
	for(i=0; i<n+1; ++i)
	  {
	
	    analyticalSFS_LUT[i] = malloc((n+1) * sizeof(mpfr_t));
	    for(j = 0; j < n+1; ++j)
	      {
		mpfr_init2(analyticalSFS_LUT[i][j], ACC);
		mpfr_set_d(analyticalSFS_LUT[i][j], 0., GMP_RNDU);
	      }
	  }

	
	alignment->SFS = calloc(sequences+1, sizeof(t_sfs));
	assert(alignment->SFS!=NULL);

	mpfr_t tmp;
	mpfr_init2(tmp, ACC);

	riexpo(tmp, factorialsLUT, n, 1, R, eventList, eventsTotal, 1);
	alignment->SFS[1] = mpfr_get_d(tmp,GMP_RNDU);

	
	for(i=2;i<n;i++)
	{
	  //printf("################################## %d\n", i);
	  riexpo(tmp, factorialsLUT, n, i, R, eventList, eventsTotal, 0);
	  alignment->SFS[i] = mpfr_get_d(tmp,GMP_RNDU);
	}

	if(fpSFSo!=NULL)
	{
		for(i=0; i<sequences + 1; ++i)
			fprintf(fpSFSo,"%d\t%e\n",i,alignment->SFS[i]);
	}

	mpfr_clear(tmp);
	
	for(i=0; i<(2*n+1); ++i)
	  mpfr_clear( factorialsLUT[i] );

	free(factorialsLUT);

	for(i=0; i<(n+1); ++i)
	  {
	    for(j=0; j<n+1; ++j)
	      {
		mpfr_clear( analyticalSFS_LUT [i][j] );
	      }
	    free(analyticalSFS_LUT[i]);
	  }

	free(analyticalSFS_LUT);
	
	mpfr_free_cache();
}



#endif
