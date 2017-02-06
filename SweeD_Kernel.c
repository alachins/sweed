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


void getMinMaxAlpha (int sweepPosition, t_sfs * minAlpha, t_sfs * maxAlpha)
{
	int i, counter=0;

	double minDist = 0., distDif = 0., totDist = 0.0;

	*maxAlpha = 1300.0;
	
	for(i=0;i<alignment->segsites;i++)
		if(alignment->positionsInd[i]!=-1)
		{
			minDist = abs(alignment->positionsInd[i]-sweepPosition);
			break;
		}

	for(i=0;i<alignment->segsites;i++)
		if(alignment->positionsInd[i]!=-1)
		{
			counter++;

			distDif = (double) abs(alignment->positionsInd[i]-sweepPosition);

			totDist += distDif;

			if(distDif < minDist)
				minDist = distDif;
		}
	

	if (minDist!=0) 
		*maxAlpha = 13.0 / (t_sfs) minDist;

	*minAlpha = 12.0 * counter /totDist;	
}

inline int getGridIndexL (double x)
{
	return (int)((log(x)-log(alignment->gridADs[0]))/alignment->interval);
}

t_sfs splint(int x, int n, int gridSz, double ad)
{  
	int indexn = n-alignment->minn;
	int indexx = x-alignment->startSFS;

	
  
	if(ad<=alignment->gridADs[0])
		return alignment->gridProbs[indexn][indexx][0];
	
	int lInd = (int)((log(ad)-alignment->logAD0)*alignment->interval);
	
	if(lInd >= gridSz -1 )
		return alignment->gridProbs[indexn][indexx][(gridSz-1)*2];

	if (lInd<0)
		lInd = 0;

	int rInd = lInd + 1;

	if(rInd>gridSz-1)
		rInd = gridSz-1;

	assert(rInd - lInd == 1);

	t_sfs gridADDif = alignment->gridADs[rInd]-alignment->gridADs[lInd];

	assert (gridADDif!=0.0);

	//t_sfs gridADDifDiv = 1.0 / gridADDif;	

	//t_sfs a=(alignment->gridADs[rInd]-ad)*gridADDifDiv;
  	//t_sfs b=(ad-alignment->gridADs[lInd])*gridADDifDiv;
	
	assert(indexn >= 0);
	if( indexx < 0 )
	  {
	    fprintf(stderr, "\nx: %d, alignment->startSFS: %d\n", x, alignment->startSFS);
	    assert(indexx >= 0);

	  }
	assert(alignment->minn >= 0);
	assert(alignment->startSFS >= 0);

	//t_sfs y=a*alignment->gridProbs[indexn][indexx][lind2]+
	//	b*alignment->gridProbs[indexn][indexx][rind2]+
	//	((a*a*a-a)*alignment->gridProbs[indexn][indexx][lind2+1]+(b*b*b-b)*alignment->gridProbs[indexn][indexx][rind2+1])*(gridADDif*gridADDif)*div6;

	//if (y==0.0 || y==1.0)
    	//	y=(alignment->gridProbs[indexn][indexx][lind2]*a+b*alignment->gridProbs[indexn][indexx][rind2]);

	//if (y <0.0 || y>1.0)
    	//	y = (alignment->gridProbs[indexn][indexx][lind2]+alignment->gridProbs[indexn][indexx][rind2])*0.5;	

	t_sfs y = alignment->gridProbs[indexn][indexx][lInd] + ((ad - alignment->gridADs[lInd])*(alignment->gridProbs[indexn][indexx][rInd]-alignment->gridProbs[indexn][indexx][lInd]))/gridADDif;

	return y;
}

t_sfs getLikelihood(int sweepPosition, t_sfs alpha, int * sweepWidth, int startPos)
{
	int i, sweepDist, n, x, calcsDone = 0;

	t_sfs ad=0.0, likelihood=0.0, pr = 0.0;

	(*sweepWidth)=0;

	int badFlag = 0;
	

	for(i=startPos;i>-1;i--)
	{

	  badFlag = 0;

		if(alignment->positionsInd[i]!=-1)
		{
			
			sweepDist = abs(alignment->positionsInd[i]-sweepPosition);
			

			if(sweepDist==0)
				ad = alpha * 0.01;
			else
				ad = alpha * sweepDist;
	
		
			if(ad<12.0)
			{

				(*sweepWidth)++;
				
				n = alignment->n[i];
				x = alignment->x[i];	

				if(ad < alignment->gridADs[0])
					pr = alignment->gridProbs[n-alignment->minn][x-alignment->startSFS][0];
				else
					pr = splint(x, n, GRIDSIZE, ad);

				assert(pr>=0.0);

					
				if(pr > 1.0+1e-6)
				  {
				    fprintf(stderr, "TODO.... Warning: A probability measurement is %f > 1.0... for the gridpoint %i\n", pr, i);
				    pr = 1.0;

				    badFlag = 1;
				  }
					
				//assert(pr<=1.0+1e-6);

				if(pr > 1.0)
					pr = 1.0;

				if (alignment->folded[i])
				{					
					pr *= alignment->SFS[x];
						
					x = n-x;
																

					if (x!=n-x)
					{
						if(ad < alignment->gridADs[0])
						  pr += alignment->gridProbs[n-alignment->minn][x-alignment->startSFS][0] * alignment->SFS[x];
						else
						  pr += splint(x, n, GRIDSIZE, ad) * alignment->SFS[x];
						
						pr = pr/(alignment->SFS[x] + alignment->SFS[n-x]);
						
						assert(pr>=0.0);
						
						
						if(pr > 1.0+1e-6)
						  {
						    fprintf(stderr, "TODO.... Warning: A probability measurement is %f > 1.0... for the gridpoint %i\n", pr, i);
						    pr = 1.0;

						    badFlag = 1;
						  }
						
						//assert(pr<=1.0+1e-6);
						
						
						if(pr > 1.0)
						  pr = 1.0;
					}
				}
    				
				assert(i>=0);
				assert(i<alignment->segsites);

				if(badFlag == 0)
				  likelihood += log(pr) - alignment->baseLikelihood[i];				
				else
				  likelihood += 0;

			}
			else
				break;		
		}
	}


	calcsDone = 0;

	for(i=startPos+1;i<alignment->segsites;i++)
	{

	  badFlag = 0;

		if(alignment->positionsInd[i]!=-1)
		{
			sweepDist = abs(alignment->positionsInd[i]-sweepPosition);			

			if(sweepDist==0)
				ad = alpha * 0.01;
			else
				ad = alpha * sweepDist;
	
		
			if(ad<12.0)
			{
				calcsDone = 1;

				(*sweepWidth)++;

				n = alignment->n[i];
				x = alignment->x[i];
	

				if(ad < alignment->gridADs[0])
					pr = alignment->gridProbs[n-alignment->minn][x-alignment->startSFS][0];
				else
					pr = splint(x, n, GRIDSIZE, ad);
			

				assert(pr>=0.0);
				
				
				if(pr > 1.0+1e-6)
				  {
				    fprintf(stderr, "TODO.... Warning: A probability measurement is %f > 1.0... for the gridpoint %i\n", pr, i);
				    pr = 1.0;

				    badFlag = 1;
				  }
				
				//assert(pr<=1.0+1e-6);

				if(pr > 1.0)
					pr = 1.0;

				if (alignment->folded[i])
				{
					pr *= alignment->SFS[x];

					x = n-x;


					if (x!=n-x)
					{
						if(ad < alignment->gridADs[0])
							pr += alignment->gridProbs[n-alignment->minn][x-alignment->startSFS][0] * alignment->SFS[x];
						else
							pr += splint(x, n, GRIDSIZE, ad) * alignment->SFS[x];

						pr = pr/(alignment->SFS[x] + alignment->SFS[n-x]);

						assert(pr>=0.0);

						if(pr > 1.0+1e-6)
						  {
						    fprintf(stderr, "TODO.... Warning: A probability measurement is %f > 1.0... for the gridpoint %i\n", pr, i);
						    pr = 1.0;

						    badFlag = 1;
						  }

						if(pr > 1.0)
							pr = 1.0;
					}

				}
				
				if(badFlag == 0)
				  likelihood += log(pr) - alignment->baseLikelihood[i];				
				else
				  likelihood += 0; 


			}
			else
				if(calcsDone == 1)
					break;		
		}
	}

	return likelihood;
}

int getClosestSNPIndex (int sweepPosition)
{
	int index, leftInd=0, rightInd=alignment->segsites-1;

	while(rightInd - leftInd > 1)
	{
		index = (rightInd + leftInd)>>1;

		if(alignment->positionsInd[index]>sweepPosition)
			rightInd = index;
		else
			leftInd = index;
	}

	return leftInd;
}

#ifdef _USE_PTHREADS
void parallelComputationAlpha_thread(int tid, int threads)
{
	//int size = (int)(gridP / (float)threads);
	//int start = tid * size;
	//int stop = start + size;

	int i;

	//if(tid==threads-1 && stop!=gridP)
	//	stop = gridP;

	//for(i=start;i<stop;i++)
	for(i=0;i<gridP;i++)
	{	
		if(i%threads==tid)
		{
			clr[i].alpha = getAlpha (clr[i].sfRealPos, &clr[i].likelihood);
#ifdef _DO_CHECKPOINTS
			writeCheckpoint();
#endif
		}
	}
}

void computeAlpha_parallel(int grid)
{
	int i, threads=threadDataL[0].threadTOTAL;

	gridP = grid;

	for(i=1;i<threads;i++)
		threadDataL[i].threadOPERATION = COMPALPHA;

	parallelComputationAlpha_thread(0,threads);

	threadDataL[0].threadBARRIER=1;

	syncThreadsBARRIER();
}
#endif

t_sfs getAlpha (int sweepPosition, t_sfs * likelihood)
{
	t_sfs minAlpha, maxAlpha, interval, tol=1.0e-6;

	int localGridSz = 100, i, maxPos=0, startPos;

	t_sfs val[LOCALGRIDSIZE];
	t_sfs lik[LOCALGRIDSIZE];
	int sweepWidth[LOCALGRIDSIZE];

	getMinMaxAlpha (sweepPosition, &minAlpha, &maxAlpha);

	startPos = getClosestSNPIndex (sweepPosition);


	while((maxAlpha-minAlpha)/((maxAlpha+minAlpha)*0.5) > tol)
	{
	
		interval = log(maxAlpha/minAlpha)/(localGridSz-1);

		maxPos = 0;

		for (i=0; i<localGridSz; i++)
		{
			val[i] = exp(log(minAlpha)+i*interval);

			sweepWidth[i] = 0;

			if (i!=0 && sweepWidth[i-1]==0)
				lik[i] = lik[i-1];
			else 
				lik[i] = getLikelihood(sweepPosition, val[i], &sweepWidth[i], startPos);
		
			if (lik[i] > lik[maxPos])
				maxPos = i;
		}

		if (maxPos==0) 
			minAlpha = exp(log(minAlpha)-localGridSz*interval);
		else 
			minAlpha = val[maxPos-1];

		if (maxPos==localGridSz-1) 
			maxAlpha += (val[maxPos-1]-val[maxPos-2]);
		else 
			maxAlpha = val[maxPos+1];

		localGridSz = 5;
	}

	*likelihood = lik[maxPos];

	return val[maxPos];
}

t_sfs computeBaseLikelihood()
{
	int i;
	t_sfs minLikelihood = 0.0;

	alignment->baseLikelihood = malloc(sizeof(t_sfs)*alignment->segsites);

	for(i=0;i<alignment->segsites;i++)
		if(alignment->positionsInd[i]!=-1)
		{
			alignment->baseLikelihood[i] = likelihoodSFS_SNP(alignment->SFS, alignment->x[i], alignment->n[i], alignment->folded[i], 1);
			minLikelihood += alignment->baseLikelihood[i];
		}

	return minLikelihood;
}
