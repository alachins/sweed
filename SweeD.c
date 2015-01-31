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

//#define _PTHREADS_PINNED

#ifdef _DO_CHECKPOINTS
static void dmtcpWriteCheckpoint()
{  
	assert(dmtcpIsEnabled());
    
	int r;
	double t;

	t = gettime();
	r = dmtcpCheckpoint();

	assert(r>=1);

	if(r==1)
		printf(" required %fs)", gettime() - t);

	if(r==2)
	{
		printf(" required %fs) RESTART ", gettime() - t);
	}	
}

void writeCheckpoint()
{
	double 
	delta,
	t = gettime();

	if((delta = (t - lastCheckpointTime)) > checkPointInterval)
	{
		printf("\n\tCheckpoint (after %fs -", delta);
		lastCheckpointTime = t;
		dmtcpWriteCheckpoint();
	}
}
#endif

double gettime(void)
{
	struct timeval ttime;
	gettimeofday(&ttime , NULL);
	return ttime.tv_sec + ttime.tv_usec * 0.000001;
}

void printHeading (FILE * fp)
{
	fprintf(fp,"\n\n                                      _______________");
	fprintf(fp,"\n\n                                           SweeD");
	fprintf(fp,"\n                                      _______________");
	fprintf(fp,"\n\n\n\n SweeD version 3.3.1 released by Nikolaos Alachiotis and Pavlos Pavlidis in January 2015.\n");
}

void printRunInfo (FILE * fp, int argc, char ** argv)
{
	fprintf(fp,"\n\n Command:\n\n\t");

	int i;
			
	for(i=0; i<argc; ++i)
		fprintf(fp," %s",argv[i]);
	
	fprintf(fp,"\n\n\n");
	
}

void introMsg(int argc, char ** argv, FILE * fpInfo)
{
	printHeading (stdout);

	printRunInfo (stdout, argc, argv);
	
	printHeading (fpInfo);

	printRunInfo (fpInfo, argc, argv);
}

void initializeGlobalPointers()
{
	alignment->positions = NULL;
	alignment->positionsInd = NULL;
	alignment->outgroupSequence = NULL;
	alignment->outgroupName = NULL;
	alignment->x = NULL;
	alignment->n = NULL;
	alignment->n_pat = NULL;
	alignment->x_pat = NULL;
	alignment->f_pat = NULL;
	alignment->c_pat = NULL;
	alignment->patternList = NULL;
	alignment->folded = NULL;
	alignment->baseLikelihood = NULL;
	alignment->gridProbs = NULL;
	alignment->gridADs =NULL;

	rvLUT = NULL;
	clr = NULL;
}

#ifdef _USE_PTHREADS

#ifdef _PTHREADS_PINNED
static void pinToCore(int tid)
{
	cpu_set_t cpuset;
         
	CPU_ZERO(&cpuset);    
	CPU_SET(tid, &cpuset);

	if(pthread_setaffinity_np(pthread_self(), sizeof(cpu_set_t), &cpuset) != 0)
	{
		fprintf(stdout, "\n ERROR: Please specify a number of threads that is smaller or equal");
		fprintf(stdout, "\n        to the number of available physical cores (%d).\n\n",tid);
		exit(0);
	}
}
#endif


void initializeThreadData(threadData_t * cur, int i, int threads)
{
	cur->threadID=i;
	cur->threadTOTAL=threads;
	cur->threadBARRIER=0;
	cur->threadOPERATION=BUSYWAIT;
}

void * thread (void * x)
{
	threadData_t * currentThread = (threadData_t *) x;
	
	int tid = currentThread->threadID, threads=currentThread->threadTOTAL;

#ifdef _PTHREADS_PINNED
	pinToCore(tid);
#endif
	
	while(1)
	{
		sleep(0);

		if(threadDataL[tid].threadOPERATION==EXIT)
		{
			return NULL;
		}

		if(threadDataL[tid].threadOPERATION==LIKSFS)
		{
			likelihoodSFS_thread(tid, threads);

			threadDataL[tid].threadBARRIER=1;
	
			while(threadDataL[tid].threadBARRIER==1) sleep(0);			
		}

		if(threadDataL[tid].threadOPERATION==COMPRVLUT)
		{
			parallelComputationRVLUT_thread(tid, threads);

			threadDataL[tid].threadBARRIER=1;
	
			while(threadDataL[tid].threadBARRIER==1) sleep(0);			
		}

		if(threadDataL[tid].threadOPERATION==COMPALPHA)
		{
			parallelComputationAlpha_thread(tid, threads);

			threadDataL[tid].threadBARRIER=1;
	
			while(threadDataL[tid].threadBARRIER==1) sleep(0);			
		}
		
	}
	
	return NULL;				
}

void terminateWorkerThreads()
{
	int i, threads=threadDataL[0].threadTOTAL;
	
	for(i=0;i<threads;i++)
		threadDataL[i].threadOPERATION = EXIT;

	sleep(0);			

	for(i=1;i<threads;i++)
		pthread_join(workerThreadL[i-1],NULL);

	free(threadDataL);
}

void syncThreadsBARRIER()
{
	int i, threads = threadDataL[0].threadTOTAL, barrierS=0;

	while(barrierS!=threads)
	{
		barrierS=0;
		for(i=0;i<threads;i++)
			barrierS += threadDataL[i].threadBARRIER;
	}

	for(i=0;i<threads;i++)
	{
		threadDataL[i].threadOPERATION=BUSYWAIT;
		threadDataL[i].threadBARRIER=0;
	}
}
#endif	

int main(int argc, char** argv)
{
	mainTime0 = gettime();

#ifdef _DO_CHECKPOINTS
  checkPointInterval = 3600.0;
  lastCheckpointTime = mainTime0;
#endif

	int grid = 0, 
	    alignmentLength=0,
	    nxt_alignment=0,
	    alignmentIndex=1,
	    fileFormat=OTHER_FORMAT,
	    resultType=-1,
	  counterRepls = 0,
	  noSeparator = 0,
	  monomorphic = 0, i, threads = 1, curSFSsize=10, sequences=-1, analyticalSFS=0, strictPolymorphic = 0, generateVCFsamplelist=0,
	  generateVCFchromlist=0,
	  minsnps_threshold_user = MINSNPS_THRESHOLD,
	  reports=0;

	unsigned char userSetFolded = 0, userSFS=0, onlySFS=0, onlySF=0;

	unsigned int seed = 0;

	double time0, time1, totalTimeL=.0, totalTimeG0 = gettime(), totalTimeG1, maxLH=0., maxALPHA=0., maxPOS=0., growthRate = 0.;
	
	char inputFileName[INFILENAMESIZE],
	     sfsFileName[INFILENAMESIZE],
	     sfsoFileName[INFILENAMESIZE],
	     sfoFileName[INFILENAMESIZE],
	     infoFileName[INFILENAMESIZE], 
	     clrReportFileName[INFILENAMESIZE],
	     clrReportsFileName[INFILENAMESIZE], 
	     warnFileName[INFILENAMESIZE],
	     sampleVCFfileName[INFILENAMESIZE],
	     chromVCFfileName[INFILENAMESIZE];

	FILE *fpIn=NULL, *fpInfo=NULL, *fpWarnings=NULL, *fpReport=NULL, *fpSFS=NULL, *fpSFSo=NULL, *fpSFo=NULL, *fpVCFsamples=NULL, *fpVCFchroms=NULL;
	
	char outgroupName[SEQNAMESIZE];

	chromList = malloc(sizeof(char**));
	chromList_SZ = -1;

	int alignment_analysis_code=0;
	int processed_chroms=0;

	outgroupName[0] = '\0';
	sfsoFileName[0] = '\0';
	inputFileName[0] = '\0';
	sampleVCFfileName[0]='\0';
	chromVCFfileName[0]='\0';

	div6 = 1.0 / 6.0;

#ifdef _ANALYTICAL_SFS
	eventsTotal = 1;
	eventList = malloc(sizeof(event)*eventsTotal);
	eventList[0].t=0.0;
	eventList[0].d=1.0;
#endif

	

   	commandLineParser(argc, 
			  argv, 
			  inputFileName,
			  sfsFileName,
			  sfsoFileName,
                          sfoFileName, 
			  &grid, 
			  &alignmentLength,	
			  &seed, 
			  &fileFormat, 
			  &resultType,
			  outgroupName,
			  &userSetFolded,
			  &userSFS,
			  &onlySFS,
			  &onlySF,
			  &monomorphic,
			  &strictPolymorphic,
			  &threads,
			  &sequences,
			  &analyticalSFS,
			  &growthRate,
			  &noSeparator,
			  sampleVCFfileName,
			  &generateVCFsamplelist,
			  chromVCFfileName,
			  &generateVCFchromlist,
			  &minsnps_threshold_user,
			  &reports);


	alignment = (alignment_struct *)malloc(sizeof(alignment_struct));

	if(sfsoFileName[0]!='\0')
		fpSFSo = fopen(sfsoFileName,"w");

	strcpy(infoFileName,"SweeD_Info.");
	
	strncat(infoFileName,runName,INFILENAMESIZE-strlen(infoFileName));

	strcpy(warnFileName,"SweeD_Warnings.");
	
	strncat(warnFileName,runName,INFILENAMESIZE-strlen(warnFileName));

	strcpy(clrReportFileName,"SweeD_Report.");

	if(reports==0) // to generate one report file
		strncat(clrReportFileName,runName,INFILENAMESIZE-strlen(clrReportFileName));


	if(inputFileName[0] != '\0')
		fpIn = fopen(inputFileName,"r");
	
	fpInfo = fopen(infoFileName,"w");

	introMsg(argc, argv, fpInfo);

	fprintf(stdout, "\n Input file format (0:ms, 1:fasta, 2:macs, 3:vcf, 4:sf): %d\n", fileFormat);
	fprintf(fpInfo, "\n Input file format (0:ms, 1:fasta, 2:macs, 3:vcf, 4:sf): %d\n", fileFormat);

	if(sampleVCFfileName[0] != '\0' && fileFormat==3)
	{
		if(generateVCFsamplelist==1)
			fpVCFsamples = fopen(sampleVCFfileName,"w");
		else
			fpVCFsamples = fopen(sampleVCFfileName,"r");			
	}

	if(chromVCFfileName[0] != '\0' && fileFormat==3)
	{
		if(generateVCFchromlist==1)
			fpVCFchroms = fopen(chromVCFfileName,"w");
		else
		{
			fpVCFchroms = fopen(chromVCFfileName,"r");

			extract_chromList(fpVCFchroms, fpInfo);
		}			
	}

#ifdef _ANALYTICAL_SFS
	
	if( (analyticalSFS==1 || analyticalSFS==2) && growthRate == 0.)
	  computeAnalyticalSFS(sequences, fpSFSo);
	else if( (analyticalSFS ==1 || analyticalSFS ==2) && growthRate > 0) 
	   computeAnalyticalSFSExpo(sequences, growthRate, fpSFSo); 
	
#endif
	
	if(onlySF!=0)
	{
		fpSFo = fopen(sfoFileName,"w");
		fprintf(fpSFo,"position\tx\tn\tfolded\n");
	}

	if (fileFormat==MS_FORMAT || fileFormat==MACS_FORMAT)
		fpWarnings = fopen(warnFileName, "w");

	if(onlySFS==0 && reports==0)
		fpReport = fopen(clrReportFileName,"w");

	




	alignment->length = alignmentLength;

	alignment->userSetFolded = userSetFolded;
	alignment->userSFS = userSFS;

	if(alignment->userSFS==1)
		fpSFS = fopen(sfsFileName,"r");

	srand(seed);

	initializeGlobalPointers();

	alignment->outgroupIndex = -1;
	alignment->discSites=0;
	
	if(outgroupName[0] != '\0')
	{
		alignment->outgroupName = malloc(SEQNAMESIZE * sizeof(char));
		strcpy(alignment->outgroupName, outgroupName);    
	}

#ifdef _USE_PTHREADS

#ifdef _PTHREADS_PINNED
	pinToCore(0);
#endif

	workerThreadL = malloc (sizeof(pthread_t)*(threads-1));
	assert(workerThreadL!=NULL);

	threadDataL = malloc (sizeof(threadData_t)*threads);
	assert(threadDataL!=NULL);

	for(i=0;i<threads;i++)
		initializeThreadData(&threadDataL[i],i,threads);

	for(i=1;i<threads;i++)
		pthread_create (&workerThreadL[i-1], NULL, thread, (void *) (&threadDataL[i]));

	alignment->threadScoresLIKSFS = malloc (sizeof(t_sfs)*threads);
	assert(alignment->threadScoresLIKSFS);
#endif

	time0 = gettime();

	nxt_alignment = findFirstAlignment(fpIn,fileFormat, fpVCFsamples, generateVCFsamplelist, sampleVCFfileName, fpInfo);		


	if(generateVCFchromlist)
	{	
		fprintf(stdout, " Generating chromosome list...\n");
		fprintf(fpInfo, " Generating chromosome list...\n");
		generate_VCF_chrom_list(fpIn, fpVCFchroms, chromVCFfileName, fpInfo);
	}
	else
	{

		while(nxt_alignment==1)
		{

		counterRepls++;

		


		  initializeAlignmentVariables();
		
		  maxLH=0.;

		  alignment_analysis_code = readAlignment(fpIn, fileFormat, fpInfo, fpSFo, minsnps_threshold_user, alignmentIndex);

		  if(processed_chroms==chromList_SZ)
			break;


		
		  if( alignment_analysis_code == 1 )
		    {


	 		  
		if(fpReport!=NULL && reports==0)
		{

			if(noSeparator == 0)
				fprintf(fpReport,"\n//%d\n",alignmentIndex);
			else
				fprintf(fpReport, "\n");
			  
			if( noSeparator == 0 || counterRepls == 1)
				fprintf(fpReport,"Position\tLikelihood\tAlpha\n");		
		}

		if(reports==1)
		{
			strcpy(clrReportsFileName,clrReportFileName);
			
	
			char str[STRINGLENGTH];
			sprintf(str, "%d", alignmentIndex);

			strncat(clrReportsFileName,runName,INFILENAMESIZE-strlen(clrReportsFileName));
			strncat(clrReportsFileName,".",INFILENAMESIZE-strlen(clrReportsFileName));
			if(fileFormat==3)
				strncat(clrReportsFileName,VCF_alignment_name_prev_to_print,INFILENAMESIZE-strlen(clrReportsFileName));
			else
				strncat(clrReportsFileName,str,INFILENAMESIZE-strlen(clrReportsFileName));

			if(fpReport!=NULL)
			{
				fclose(fpReport);
				fpReport=NULL;
			}

			fpReport = fopen(clrReportsFileName,"w");

			if( noSeparator == 0 || counterRepls == 1)
				fprintf(fpReport,"Position\tLikelihood\tAlpha\n");
		}	  
			  

			 // fprintf(stdout," Alignment %d\n",alignmentIndex);			
			 // fprintf(fpInfo," Alignment %d\n",alignmentIndex);

			 



	#ifdef _ANALYTICAL_SFS
		      if(analyticalSFS==2)
			assert(sequences>=alignment->sequences);
	#endif
		      removeMonomorphicSites (strictPolymorphic, monomorphic, fpInfo);

		      
		      if(alignment -> segsites < MINSNPS_THRESHOLD || alignment -> segsites < minsnps_threshold_user)
			{

			       
			  fprintf(stdout,"\n\n\t\tProcessing ... NOT POSSIBLE (too few SNPs after removing monomorphic)\n\n");
			  fprintf(fpInfo,"\n\n\t\tProcessing ... NOT POSSIBLE (too few SNPs after removing monomorphic)\n\n");
			  
			  if( (onlySF ==0 || onlySF ==2) && onlySFS == 0 )
			    freeAlignment_noSNPs();
		  
			  nxt_alignment = findNextAlignment(fpIn,fileFormat);
		  
			  alignmentIndex++;
		  
			  if(nxt_alignment==1 && fpSFo!=NULL)
			    fprintf(fpSFo, "\n// %d\n",alignmentIndex);
		 
			  continue; 
			}

		      
		      updateSF(fpSFo);

		      
		      
		      if (fileFormat==MS_FORMAT || fileFormat == MACS_FORMAT)
			checkSNIPPositions(fpWarnings, alignment, alignmentIndex);		
		      
		      if(onlySF==0 || onlySF==2 )
			{		
			  
			  fprintf(stdout,"\n\n\t\tProcessing");
			  fprintf(fpInfo,"\n\n\t\tProcessing");
			  fflush(stdout);
			  
			  processed_chroms++;

			  computeFactLUT();
			  
			  if(analyticalSFS==0)
			    createSFS (fpSFS, fpSFSo, alignmentIndex, &curSFSsize);
			  else
			    {
			      alignment->startSFS=1;
			      alignment->endSFS=alignment->sequences;
			    }
			  
			  
			  
			  
	#ifdef _DO_CHECKPOINTS
			  printf("\n");
			  writeCheckpoint();
	#endif
			  
			  if(onlySFS==0)
			    {
			      createCLR (grid);
			      
			      createPROBS (GRIDSIZE);	
			      
			      computeBaseLikelihood();
			      
			      
			      
	#ifdef _DO_CHECKPOINTS
			      writeCheckpoint();
	#endif
			      
	#ifdef _USE_PTHREADS
			      computeAlpha_parallel(grid);
			      
			      for(i=0;i<grid;i++)
				{	
				  if(clr[i].likelihood>maxLH)
				    {
				      maxLH = clr[i].likelihood;
				      maxALPHA = clr[i].alpha;
				      maxPOS = clr[i].sfRealPos;			        					}
				}
	#else
			      for(i=0;i<grid;i++)
				{	
				  clr[i].alpha = getAlpha (clr[i].sfRealPos, &clr[i].likelihood);
				  
				  if(clr[i].likelihood>maxLH)
				    {
				      maxLH = clr[i].likelihood;
				      maxALPHA = clr[i].alpha;
				      maxPOS = clr[i].sfRealPos;					
				    }
				  
	#ifdef _DO_CHECKPOINTS
				  writeCheckpoint();
	#endif
				}
	#endif
			      if(fpReport!=NULL)
				for(i=0;i<grid;i++)
				  fprintf(fpReport,"%.4f\t%e\t%e\n", clr[i].sfPos/* clr[i].sfRealPos*/, clr[i].likelihood,clr[i].alpha);		
			      
			      
			      time1 = gettime();
			      
			      totalTimeL = time1-time0; 
			      
			      fprintf(stdout, ":\t\t%.2f seconds\n\n",totalTimeL);
			      fprintf(fpInfo, ":\t\t%.2f seconds\n\n",totalTimeL);
			      
			      //if(fpInfo!=NULL)
			      fprintf(stdout,"\t\tPosition:\t\t%e\n\t\tLikelihood:\t\t%e\n\t\tAlpha:\t\t\t%e\n\n\n",maxPOS,maxLH,maxALPHA);
			      fprintf(fpInfo,"\t\tPosition:\t\t%e\n\t\tLikelihood:\t\t%e\n\t\tAlpha:\t\t\t%e\n\n\n",maxPOS,maxLH,maxALPHA);
			      
			      time0 = gettime();
			      
			      
			      
			    }
			}

		  if( (onlySF ==0 || onlySF ==2) && onlySFS == 0 )
		    freeAlignment();


		    }
		  else
		    {

			if( alignment_analysis_code == 0 )
			{
		      
			  fprintf(stdout,"\n\n\t\tProcessing ... NOT POSSIBLE (too few SNPs)\n\n");
			  fprintf(fpInfo,"\n\n\t\tProcessing ... NOT POSSIBLE (too few SNPs)\n\n");

			}

			if( (onlySF ==0 || onlySF ==2) && onlySFS == 0 )
			    freeAlignment_noSNPs();
		    }
		  

		
		  nxt_alignment = findNextAlignment(fpIn,fileFormat);
		
		  alignmentIndex++;
		  
		  if(nxt_alignment==1 && fpSFo!=NULL)
		    fprintf(fpSFo, "\n// %d\n",alignmentIndex);
		  
		  
		}
	}
	
	totalTimeG1 = gettime();

#ifdef _USE_PTHREADS
	terminateWorkerThreads();
#endif

	fprintf(stdout, " Total elapsed time %.2f seconds.\n\n\n",totalTimeG1-totalTimeG0);
	fprintf(fpInfo, " Total elapsed time %.2f seconds.\n\n\n",totalTimeG1-totalTimeG0);

	if (fileFormat==MS_FORMAT || fileFormat == MACS_FORMAT) 
		fclose(fpWarnings);

	if(alignment->SFS!=NULL)
	{
		free(alignment->SFS);
		alignment->SFS = NULL;
	}

	if(alignment->VCFsample_valid!=NULL)
	{
		free(alignment->VCFsample_valid);
		alignment->VCFsample_valid = NULL;
	}

	if(chromList!=NULL)
	{
		if(*chromList!=NULL)
		{
			for(i=0;i<chromList_SZ;i++)
				if((*chromList)[i]!=NULL)
					free((*chromList)[i]);

			free(*chromList);
		}
		free(chromList);
	}

	free(alignment);

	if(fpIn!=NULL)
		fclose(fpIn);

	if(fpInfo!=NULL)
		fclose(fpInfo); 

	if(fpReport!=NULL)
		fclose(fpReport);
	
	if(fpSFSo!=NULL)
		fclose(fpSFSo);

	if(fpSFo!=NULL)
		fclose(fpSFo);

	if(fpVCFsamples!=NULL)
		fclose(fpVCFsamples);

	if(fpVCFchroms!=NULL)
		fclose(fpVCFchroms);

#ifdef _ANALYTICAL_SFS
	free(eventList);
#endif

	return 0;
}
