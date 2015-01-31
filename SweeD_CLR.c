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

void createCLR (int grid) // composite likelihood ratio
{
	int i, startsnp, endsnp;
	float step;

	assert(grid>1);

	clr = (clr_struct *)malloc(sizeof(clr_struct)*grid);

	for(i=0;i<alignment->segsites;i++)
		if(alignment->positionsInd[i]!=-1)
			break;

	assert(i!=alignment->segsites);

	startsnp = alignment->positionsInd[i];

	for(i=alignment->segsites-1;i>-1;i--)
		if(alignment->positionsInd[i]!=-1)
			break;

	assert(i!=-1);

	endsnp = alignment->positionsInd[i];

	step = (double)(endsnp - startsnp)/(grid-1.0);
	
	alignment->sweepStep = step;

	if(step < 1)
	{
		fprintf(stderr,"\n\n WARNING: The requested grid size is too large (%d) for the region between the first and last sites (region size: %d).\n\n", grid, endsnp-startsnp+1);
						
	}

	clr[0].sfRealPos = startsnp;
	clr[0].sfPos = startsnp;
	clr[grid-1].sfRealPos = clr[grid-1].sfPos = endsnp;

	for(i=1;i<grid-1;i++)
	  {
	    clr[i].sfRealPos = (int)(clr[0].sfRealPos + i * step); 
	    clr[i].sfPos = clr[i-1].sfPos + step;
	  }
}

double XchooseY_ln(int x, int y)
{
	assert(x>=0);
	assert(y>=0);

	if(x==y || y==0)
		return 0.0;

		
	assert(y<=x && x!=0);

	return factLUT->val[x]-factLUT->val[x-y]-factLUT->val[y];
}

double XchooseY_ln_2(int x, int y, int * setZero)
{
	assert(x>=0);
	assert(y>=0);

	if(x==y || y==0)
		return 0.0;

	*setZero = 0;
	if(x < y)
	{
		*setZero = 1;
		return 0.;
	}
		
	assert(y<=x && x!=0);

	return factLUT->val[x]-factLUT->val[x-y]-factLUT->val[y];
}


double XchooseY(int x, int y)
{	
	assert(y>=0);

	if (x<y) return 0.0;

	if(x==y || y==0)
		return 1.0;

	return exp(factLUT->val[x]-factLUT->val[x-y]-factLUT->val[y]);
}

double getProb_chooseXfromYgivSFS (int x, int y, double divVal)
{
	int i,setZero=0;
	double prob=0.0, prob1, prob2;
	
	if (x < 0 || x > y)
		return 0.0;

	for (i=x; i<alignment->maxn+1; i++)
		if (y-x>=0 && i<alignment->maxx+1)
		{
			prob1 = XchooseY_ln_2(i,x,&setZero);
			if(setZero == 1)
				continue;
			
			prob2 = XchooseY_ln_2(alignment->sequences-i,y-x,&setZero);
			if(setZero == 1)
				continue;

			//printf("prob1: %e, prob2: %e, divVal: %e\n", prob1, prob2, divVal);

			//prob2tmp =  alignment->SFS[i]*XchooseY(i,x)*XchooseY(alignment->sequences-i,y-x)*divVal;						
			//double tmpprob = alignment->SFS[i]*exp(prob1+prob2+divVal);

			//double tmpprob = alignment->SFS[i]*exp(prob1+prob2)*divVal;

			//printf("tmpprob: %e\n", tmpprob);

			prob += alignment->SFS[i]*exp(prob1+prob2+divVal);
			
			//prob += alignment->SFS[i]*exp(prob1+prob2)*divVal;
		}



	assert(prob <= 1. && prob >= 0);
	

	return prob;
}

//double getProb_Xescapes_ext (int x, double exp_minad)
//{
//	return pow((1.0-exp_minad),x)*pow(exp_minad,(alignment->sequences-x));
//}

double getProb_Xescapes_ext (int x, double exp_minad)
{
	return exp(x*log(1.0-exp_minad) + (alignment->sequences-x)*log(exp_minad));//*pow(exp_minad,(alignment->sequences-x));
}


double getProb_Xescapes (int x, double exp_minad)
{
	return x*log(1.0-exp_minad)+ (alignment->sequences-x)*log(exp_minad);
}
/*
double getProb_Xescapes (double factor1, int x, double exp_minad)
{
	return factor1*pow((1.0-exp_minad),x)*pow(exp_minad,(alignment->sequences-x));
}*/

double getProb_Total (int i, int k, int minx_offset, int maxx_offset)
{
	int j;
	double pr=0.0, sum = 0.0;

	for (j=1; j<alignment->sequences; j++)
		sum += rvLUT[k][j];

	if (maxx_offset == 1)
		sum += rvLUT[k][alignment->sequences];

	if(minx_offset == 0)
		sum += rvLUT[k][0];

	assert(sum!=0.0);

	pr = rvLUT[k][i]/sum;

	assert(pr >= 0.);
	assert(pr <= 1.);	

	return pr;
}

double getY (double * inputVec, int index)
{
	return inputVec[index*2];
}

double getY2 (double * inputVec, int index)
{
	return inputVec[index*2+1];
}

void setY2 (double * inputVec, int index, double value)
{
	inputVec[index*2+1] = value;
}

void spline(double *x, double *y ,int n, double yp1, double ypn, double *y2, double *YY2)
{
	int i,k;
	double p,qn,sig,un,*u, tmp;

	u=malloc(n*sizeof(double));

	if (yp1 > 0.99e30)
	{
		u[0]=0.0;

		setY2(YY2, 0, 0.0);
	}
	else
	{	
		setY2(YY2, 0, -0.5);

		assert(x[1]!=x[0]);
		
		tmp = 1/(x[1]-x[0]);

		u[0]=(3.0*tmp)*((getY(YY2,1)-getY(YY2,0))*tmp-yp1);
	}

	for (i=1;i<=n-2;i++)
	{
		assert(x[i+1]!=x[i-1]);

		sig=(x[i]-x[i-1])/(x[i+1]-x[i-1]);

		p=sig*getY2(YY2,i-1)+2.0;

		assert(p!=0.0);

		setY2(YY2,i, (sig-1.0)/p);

		assert(x[i+1]!=x[i]);
		assert(x[i]!=x[i-1]);
		
		u[i]=(getY(YY2,i+1)-getY(YY2,i))/(x[i+1]-x[i]) - (getY(YY2,i)-getY(YY2,i-1))/(x[i]-x[i-1]);

		u[i]=(6.0*u[i]/(x[i+1]-x[i-1])-sig*u[i-1])/p;		
	}

	//if (ypn > 0.99e30)
		qn=un=0.0;
	//else
	//	exit(-1);

	assert((qn*getY2(YY2,n-2)+1.0)!=0.0);
		
	setY2(YY2,n-1, (un-qn*u[n-2])/(qn*getY2(YY2,n-2)+1.0));

	for (k=n-2;k>=0;k--)
		setY2(YY2,k, getY2(YY2,k)*getY2(YY2,k+1)+u[k]);

	free(u);
}

#ifdef _USE_PTHREADS
void parallelComputationRVLUT_thread(int tid, int threads)
{
	int size = (int)((alignment->sequences+1) / (float)threads);

	assert(size!=0);

	int start = tid * size;
	int stop = start + size;

	if(tid==threads-1 && stop!=alignment->sequences+1)
		stop = alignment->sequences+1;


	int i,k,b, setZero=0;
	double tmpVAL, factor2, factor3, divVal,tmpSFSval;

	for (b=start; b<stop; b++)
	{	
		tmpSFSval = log(alignment->SFS[b]);

		if (b < alignment->maxx)
			for(i=0;i<GRIDSIZE;i++)
				rvLUT[i][b] = exp(tmpSFSval + getProb_Xescapes(alignment->sequences, alignment->gridADs[i]));				

		factor2 = 0.0;

		for (k=0; k<alignment->sequences; k++)
		{			
			factor3 = XchooseY_ln_2(alignment->sequences, k+1, &setZero); 
		  
			divVal = -factor3; 

			tmpVAL = (   getProb_chooseXfromYgivSFS(b+1-alignment->sequences+k, k+1, divVal)*(b+1-alignment->sequences+k)  +   getProb_chooseXfromYgivSFS(b, k+1, divVal)*(k+1-b)   )/(k+1);

			if(tmpVAL==0.0)
			{
				factor2 = factor3;	
				continue;
			}

			tmpVAL = log(tmpVAL)+factor2;

			for(i=0;i<GRIDSIZE;i++)
				rvLUT[i][b] += exp(tmpVAL+getProb_Xescapes(k, alignment->gridADs[i]));

			factor2 = factor3;
		}
	}	
}

void parallelComputationRVLUT()
{

	int i, threads=threadDataL[0].threadTOTAL;


	for(i=1;i<threads;i++)
		threadDataL[i].threadOPERATION = COMPRVLUT;

	parallelComputationRVLUT_thread(0,threads);

	threadDataL[0].threadBARRIER=1;

	syncThreadsBARRIER();
}
#endif

void createPROBS (int probGrid)
{
	int i, j, k, l;
	double prob;	

	alignment->gridProbs = malloc(sizeof(t_sfs**)*(alignment->maxn-alignment->minn+1));
	assert(alignment->gridProbs!=NULL);
	
	int minx_offset = alignment->SFS[0]>0.0?0:1;
	int maxx_offset = alignment->SFS[alignment->sequences]>0.0?1:0;
	
	for(i=0;i<alignment->maxn-alignment->minn+1;i++)
	{
		alignment->gridProbs[i] = malloc(sizeof(t_sfs*)*(alignment->sequences+maxx_offset-minx_offset+1));//(alignment->maxx-alignment->minx+1));
		assert(alignment->gridProbs[i]!=NULL);

		for(j=0;j<alignment->sequences+maxx_offset-minx_offset+1;++j)
		{
			alignment->gridProbs[i][j] = calloc(probGrid, sizeof(t_sfs));
			assert(alignment->gridProbs[i][j]!=NULL);
		}
	}	

	alignment->gridADs = malloc(sizeof(t_sfs)*probGrid);
	assert(alignment->gridADs!=NULL);	

	t_sfs minval=1.0e-8;
	t_sfs maxval=12;
	t_sfs logminval=log(minval);
	t_sfs interval=log(maxval/minval)/(probGrid-1);
	alignment->interval = 1.0/interval;

	rvLUT = malloc(probGrid * sizeof(double*));
	assert(rvLUT!=NULL);

	for(i=0;i<probGrid;i++)
	{
		rvLUT[i]= calloc(alignment->sequences+1,sizeof(double));
		assert(rvLUT[i]!=NULL);
	}

	for(i=0;i<probGrid;i++)
	{	
		alignment->gridADs[i]=exp(logminval+i*interval);
		alignment->gridADs[i]=exp(-alignment->gridADs[i]);
	}


#ifdef _USE_PTHREADS
	parallelComputationRVLUT();
#else
	int b, setZero = 0;
	double tmpVAL, divVal, factor2, factor3, tmpSFSval;

	for (b=0; b<alignment->sequences+1; b++)
	{	
		tmpSFSval = log(alignment->SFS[b]);

		if (b < alignment->maxx)
			for(i=0;i<probGrid;i++)
				rvLUT[i][b] = exp(tmpSFSval + getProb_Xescapes(alignment->sequences, alignment->gridADs[i]));				

		factor2 = 0.0;

		for (k=0; k<alignment->sequences; k++)
		{			
			factor3 = XchooseY_ln_2(alignment->sequences, k+1, &setZero); 
		  
			divVal = -factor3; 

			tmpVAL = (   getProb_chooseXfromYgivSFS(b+1-alignment->sequences+k, k+1, divVal)*(b+1-alignment->sequences+k)  +   getProb_chooseXfromYgivSFS(b, k+1, divVal)*(k+1-b)   )/(k+1);

			if(tmpVAL==0.0)
			{
				factor2 = factor3;	
				continue;
			}

			tmpVAL = log(tmpVAL)+factor2;

			for(i=0;i<probGrid;i++)
				rvLUT[i][b] += exp(tmpVAL+getProb_Xescapes(k, alignment->gridADs[i]));

			factor2 = factor3;
		}
	}
#endif

	for(i=0;i<probGrid;i++)
	{	
		alignment->gridADs[i]=exp(logminval+i*interval);
		
		for(j=alignment->minn;j<alignment->maxn+1;++j)
			for(k=minx_offset;k<j+maxx_offset;k++)
				for(l=k;l<k+alignment->maxn+1 - j;l++)
				{
					prob = exp(XchooseY_ln(l, k)+XchooseY_ln(alignment->maxn-l, j-k)-XchooseY_ln(alignment->maxn, j));
					
					assert(prob>=0.0 && prob<=1.0+1e-6);

					if(prob > 1.0)
						prob = 1.0;
		
					alignment->gridProbs[j-alignment->minn][k-minx_offset][i] += prob * getProb_Total (l, i, minx_offset, maxx_offset);
				}	
	}

//	for(j=alignment->minn;j<alignment->maxn+1;++j)
//		for(k=minx_offset;k<j+maxx_offset;k++)
//			for(i=0;i<probGrid;i++)
//				fprintf(stderr,"1\t%d\t%e\t%e\n",k, alignment->gridADs[i],alignment->gridProbs[j-alignment->minn][k-minx_offset][i]);	
	
	

	alignment->logAD0 = log(alignment->gridADs[0]);
}
