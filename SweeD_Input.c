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

void printVersion (FILE * fp)
{

  
  	/* fprintf(fp,"\n\n"); */
	/* fprintf(fp,"\tVersion:\t\t3.2.6\n\n"); */
	/* fprintf(fp,"\tReleased:\t\tMay 2014\n\n"); */
	/* fprintf(fp, "\tComment:\t\tIt prints only the integer positions\n\n"); */
	/* fprintf(fp,"\n\n"); */

	
  	/* fprintf(fp,"\n\n"); */
	/* fprintf(fp,"\tVersion:\t\t3.2.7\n\n"); */
	/* fprintf(fp,"\tReleased:\t\tMay 2014\n\n"); */
	/* fprintf(fp, "\tComment:\t\tHandles properly the ^M NEWLINE combinations\n\n"); */
	/* fprintf(fp,"\n\n"); */

	
  	/* fprintf(fp,"\n\n"); */
	/* fprintf(fp,"\tVersion:\t\t3.2.8\n\n"); */
	/* fprintf(fp,"\tReleased:\t\tMay 2014\n\n"); */
	/* fprintf(fp, "\tComment:\t\tFix a bug in recognizing the filetype\n\n"); */
	/* fprintf(fp,"\n\n"); */

	
  	/* fprintf(fp,"\n\n"); */
	/* fprintf(fp,"\tVersion:\t\t3.2.9\n\n"); */
	/* fprintf(fp,"\tReleased:\t\tMay 2014\n\n"); */
	/* fprintf(fp, "\tComment:\t\tAllows the distance between grid positions to be less than 1\n\n"); */
	/* fprintf(fp,"\n\n"); */

	
		
//  	fprintf(fp,"\n\n");
//	fprintf(fp,"\tVersion:\t\t3.2.10\n\n");
//	fprintf(fp,"\tReleased:\t\tMay 2014\n\n");
//	fprintf(fp, "\tComment:\t\tskips empty alignments or alignments with no polymorphic sites\n\n");
//	fprintf(fp,"\n\n");


 	/* fprintf(fp,"\n\n"); */
	/* fprintf(fp,"\tVersion:\t\t3.2.11\n\n"); */
	/* fprintf(fp,"\tReleased:\t\tJuly 2014\n\n"); */
	/* fprintf(fp, "\tComment:\t\tCan process only a subset of the VCF samples in the input file\n\n"); */
	/* fprintf(fp,"\n\n"); */


	
      /*fprintf(fp,"\n\n");
	fprintf(fp,"\tVersion:\t\t3.2.12\n\n");
	fprintf(fp,"\tReleased:\t\tNovember 2014\n\n");
	fprintf(fp, "\tComment:\t\tCorrect a bug: We removed a free statement vcf_validsample\n\n");
	fprintf(fp,"\n\n");*/

/* 	fprintf(fp,"\n\n");
	fprintf(fp,"\tVersion:\t\t3.3.0\n\n");
	fprintf(fp,"\tReleased:\t\tDecember 2014\n\n");
	fprintf(fp, "\tComment:\t\tAdded -minsnps flag\n\n");
	fprintf(fp,"\n\n");
*/

 	fprintf(fp,"\n\n");
	fprintf(fp,"\tVersion:\t\t3.3.1\n\n");
	fprintf(fp,"\tReleased:\t\tJanuary 2015\n\n");
	fprintf(fp, "\tComment:\t\tCan process only a subset of the VCF chromosomes in the input file\n\n");
	fprintf(fp,"\n\n");

  
  	/* fprintf(fp,"\n\n"); */
	/* fprintf(fp,"\tVersion:\t\t3.2.4\n\n"); */
	/* fprintf(fp,"\tReleased:\t\tMarch 2014\n\n"); */
	/* fprintf(fp,"\n\n"); */



  	/* fprintf(fp,"\n\n"); */
	/* fprintf(fp,"\tVersion:\t\t3.1\n\n"); */
	/* fprintf(fp,"\tReleased:\t\tJanuary 2013\n\n"); */
	/* fprintf(fp,"\tComments:\t\t1. Fixed bug in the VCF parser that was associated with the handling of missing data\n"); */
	/* fprintf(fp,"\n\n"); */


	
//	fprintf(fp,"\n\n");
//	fprintf(fp,"\tVersion:\t\t3.0\n\n");
//	fprintf(fp,"\tReleased:\t\tNovember 2012\n\n");
//	fprintf(fp,"\tComments:\t\t1. Anlaytical calculation of SFS\n");
//	fprintf(fp,"\n\n");

//	fprintf(fp,"\n\n");
//	fprintf(fp,"\tVersion:\t\t2.6\n\n");
//	fprintf(fp,"\tReleased:\t\tNovember 2012\n\n");
//	fprintf(fp,"\tComments:\t\t1. Scalable for more than 1000 sequences\n");
//	fprintf(fp,"\t\t\t\t2. Parallel version using POSIX threads\n");
//	fprintf(fp,"\t\t\t\t3. Bug fixes\n");
//	fprintf(fp,"\n\n");

//	fprintf(fp,"\n\n");
//	fprintf(fp,"\tVersion:\t\t1.0\n\n");
//	fprintf(fp,"\tReleased:\t\tSeptember 2012\n\n");
//	fprintf(fp,"\tComments:\t\tFirst release\n\n");
//	fprintf(fp,"\n\n");
}

void printHelp (FILE * fp)
{
	fprintf(fp,"\n\n\n");
	fprintf(fp," SweeD | SweeD-P | SweeD-C | SweeD-P-C\n");
	fprintf(fp,"\t -name runName\n");
	fprintf(fp,"\t -input inputFile\n");
	fprintf(fp,"\t -grid gridNumber\n");
	fprintf(fp,"\t[-folded]\n");
	fprintf(fp,"\t[-monomorphic]\n");
	fprintf(fp, "\t[-strictPolymorphic]\n");
	fprintf(fp,"\t[-isfs inputSFS]\n");
	fprintf(fp,"\t[-osfs outputSFS]\n");
	fprintf(fp,"\t[-osf outputSF]\n");
	fprintf(fp,"\t[-threads threadNumber]\n");
	fprintf(fp,"\t[-checkpoint checkpointInterval]\n");
	fprintf(fp,"\t[-eN timeT sizeX]\n");
	fprintf(fp,"\t[-G rateG]\n");	
	fprintf(fp,"\t[-s sequences]\n");
	fprintf(fp,"\t[-h|-help]\n");
	fprintf(fp,"\t[-v|version]\n");
	fprintf(fp,"\t[-noSeparator]\n");
	fprintf(fp,"\t[-sampleList filename]\n");
	fprintf(fp,"\t[-sampleList_out filename]\n");
	fprintf(fp,"\t[-minsnps minimumNumber]\n");
	fprintf(fp,"\t[-chromList filename]\n");
	fprintf(fp,"\t[-chromList_out filename]\n");
	fprintf(fp,"\t[-reports]\n");
	fprintf(fp,"\n\n");
	
	fprintf(fp,"\t-name <STRING>\t\tSpecifies a name for the run and the output files.\n\n");
	fprintf(fp,"\t-input <STRING>\t\tSpecifies the name of the input alignment file.\n");
	fprintf(fp,"\t      \t\t\tSupported file formats: SF (Sweep Finder) format.\n\n");
	fprintf(fp,"\t-grid <INTEGER>\t\tSpecifies the number of positions in the alignment where the CLR will be computed. \n\n");
	fprintf(fp,"\t-folded\t\t\tConsiders the SFS folded (the ancestral and derived states can not be distinguished).\n\n");
	fprintf(fp,"\t-monomorphic\t\tIncludes the monomorphic sites in the analysis. The default action is to discard them.\n\n");
	fprintf(fp, "\t-strictPolymorphic\tDoes not include *potential monomorphic* sites in the analysis. These are sites where for some sequences the state is missing and the remaining are monomorphic.\n\n");
	fprintf(fp,"\t-isfs <STRING>\t\tSpecifies the name of the input SFS file.\n\n");
	fprintf(fp,"\t-osfs <STRING>\t\tSpecifies the name of the output SFS file.\n\n");
	fprintf(fp,"\t-osf <STRING>\t\tSpecifies the name of the output SF file.\n\n");
	fprintf(fp,"\t-threads <INTEGER>\tSpecifies the number of threads.\n\n");
	fprintf(fp,"\t-checkpoint <INTEGER>\tSpecifies the checkpoint interval in seconds (default: 3600).\n\n");
	fprintf(fp,"\t-eN <FLOAT> <FLOAT>\tSets population size to sizeX*N0 at time timeT, where N0 is the present-day population size.\n\n");
	fprintf(fp,"\t-G <FLOAT>\t\tSets the growth rate of the population size at time 0. The growth rate continues to be exponential until\n\t\t\t\tthe -eN command specifies a constant population size.\n\n");	
	fprintf(fp,"\t-s <INTEGER>\t\tSpecifies the number of sequences when no input file is provided.\n\n");
	fprintf(fp,"\t-h|-help\t\tDisplays this help message.\n\n");
	
	fprintf(fp,"\t-v|-version\t\tDisplays version information.\n\n");
	fprintf(fp,"\t-noSeparator\t\tTo suppress printing the // flag that separates datasets. Useful for meta-processing (particularly with R).\n\n");
	fprintf(fp,"\t-sampleList <STRING>\tTo be used with VCF files in order to specify which samples to be included in the analysis.\n\n");
	fprintf(fp,"\t-sampleList_out <STRING> To generate a list of the samples in the input VCF file.\n\n");
	fprintf(fp,"\t-minsnps <INTEGER>\tSpecifies the minimum number of SNPs required to proceed with the analysis (default: 2).\n\n");
	fprintf(fp,"\t-chromList <STRING>\tTo be used with VCF files in order to specify which chromosomes to analyze.\n\n");
	fprintf(fp,"\t-chromList_out <STRING> To generate a list of the chromosomes in the input VCF file.\n\n");
	fprintf(fp,"\t-reports\t\tTo generate each alignment report in a separate file.\n\n");
	fprintf(fp,"\n\n");
}

int flagMatch(FILE *fp, char flag[], int flaglength, char tmp)
{
	int counter = 0;
	while(counter < flaglength)
	{

	
		if(tmp != flag[counter])
		  {
		    break;
		  }
		tmp = fgetc(fp);

		++counter;
	}
	
	return counter;
}

/* void unGetChars(FILE *fp, char* chars, int n) */
/* { */
/* 	int i=0;  */
/* 	fprintf(stderr, "\nput back....\n"); */
/* 	for(i=0; i<n; ++i) */
/* 	  { */
/* 	    fprintf(stderr, "%c", chars[n-1-i]); */
/* 	    ungetc(chars[n-1-i], fp); */
/* 	  } */
/* 	fprintf(stderr, "\n"); */
/* } */



void ignoreLineSpaces(FILE *fp, char *ent)
{
	while(*ent==' '|| *ent == 9) // horizontal tab
		*ent = fgetc(fp);  
}


int getFileFormat (FILE * fp)
{
	char tmp;

	int cnt = 0;

	tmp=fgetc(fp);

	char macsFLAG[] = "COMMAND";

	int flaglength = 7,j;

	char vcfFLAG[] = "##fileformat=VCF";

	int vcf_flaglength = 16;

	char sfFLAG[]="position";

	int sf_flaglength = 8;	  

	int format = -999;

	while(tmp!=EOF)
	{

	  
		if(tmp=='/')
		{
			tmp = fgetc(fp);

			

			if(tmp=='/' && cnt < 3)
			  {
			    format = MS_FORMAT;
			    break; 
			  }
			else
			  tmp = fgetc(fp);				
		}
		else
		{
			if(tmp=='>')
			  {
			    format = FASTA_FORMAT;
			    break;
			  }
			else
			{
			  j = flagMatch(fp, macsFLAG, flaglength, tmp);
			  
			  if(j == flaglength)
			    {
			      format = MACS_FORMAT;
			      break;
			    }
			  else
			    {
			      // unset
			      
			      fseek(fp, -j-1, SEEK_CUR);

			      tmp = fgetc(fp );
			      
			      j = flagMatch(fp, vcfFLAG, vcf_flaglength, tmp);
			      
			      if(j == vcf_flaglength)
				{
				  
				  
				  format = VCF_FORMAT;
				  break;
				}
			      else
				{
				  // unset
				  
				  fseek(fp, -j-1, SEEK_CUR);
				  
				  tmp = fgetc( fp );
				  
				  j = flagMatch(fp, sfFLAG, sf_flaglength, tmp);
				  
				  if( j == sf_flaglength)
				    {
				      tmp = fgetc(fp);
				      
				      ignoreLineSpaces(fp, &tmp);
				      
				      int j1 = flagMatch(fp, "x", 1, tmp);
				      
				      if( j1 == 1)
					{
					  tmp = fgetc( fp );
					  
					  ignoreLineSpaces(fp, &tmp);

					  int j2 = flagMatch(fp, "n", 1, tmp);

					  if( j2 == 1)
					    {					      
					      format = SF_FORMAT;
					      break;
					    }
					}
				    
				    }
				  else
				    {
				      tmp = fgetc(fp);
				      
				      if(tmp == 10 || tmp == 13)
					cnt = 0;
				      
				      if(tmp != ' ')
					cnt++;
				      
				    }
				}
			    }
			}
		}		
	}
	
	if(format == -999)
	  format = OTHER_FORMAT;

	//fprintf(stderr, "\n\nFormat of the data is: %d (0:ms, 1:fasta, 2:macs, 3:vcf, 4:sf)\n", format);
	

	return format;
}

void commandLineParser(int argc, char** argv, 
		       char * infile,
		       char * sfsfile,
		       char * sfsofile,
		       char * sfofile, 
		       int * grid, 
                       int * length,                   
                       unsigned int * seed,
		       int * fileFormat,	
		       int * results,
		       char * outgroupName,
		       unsigned char *userSetFolded,
		       unsigned char * userSFS,
		       unsigned char * onlySFS,
		       unsigned char * onlySF,
		       int * monomorphic,
		       int *strictPolymorphic,
		       int * threads,
		       int * sequences,
		       int * analyticalSFS, 
		       double *growthRate,
		       int *noSeparator,
		       char * samplefile_i,
		       int * generateVCFsamplelist,
		       char * chromfile_i,
		       int * generateVCFchromlist,
		       int * minsnps_threshold_user,
		       int * reports)
{

#ifdef _USE_PTHREADS
	int threadsSet=0;
#endif

	int nameSet = 0, fileSet=0, gridSet=0, i, lengthSet=0;
	
	FILE *fp, *fp2;
	
	strcpy(runName,"x");

	*threads = -1;

	sfofile[0]=0;

	for(i=1; i<argc; ++i)
	{

	  if(!strcmp(argv[i], "-noSeparator") || !strcmp(argv[i], "-noseparator"))
	    {
	      *noSeparator = 1;
	      continue;
	    }

		if(!strcmp(argv[i], "-name")) 
		{ 
			if (i!=argc-1)
			{
				if(argv[i+1][0]=='-')
				{
					fprintf(stderr, "\n ERROR: Run name (string) is missing after argument -name\n\n");
					exit(0);
				}
	
				strcpy(runName,argv[++i]);

				nameSet = 1;
			}
			else
			{
				fprintf(stderr, "\n ERROR: Run name (string) is missing after argument -name\n\n");
				exit(0);
			}
			continue;
		}

		if(!strcmp(argv[i], "-input"))
		{ 
			if (i!=argc-1)
			{	
				if(argv[i+1][0]=='-')
				{
					fprintf(stderr, "\n ERROR: File name (string) is missing after argument -input\n\n");
					exit(0);
				}
					
				strcpy(infile,argv[++i]);

				fp=fopen(infile,"r");

				if (fp==NULL)
				{
					fprintf(stderr, "\n ERROR: File %s does not exist\n\n",infile);
					exit(0);
				}
				else
				{
					*fileFormat = getFileFormat (fp);

					if(* fileFormat == OTHER_FORMAT)
					{
						fprintf(stderr, "\n ERROR: Please check the input file %s. The format could not be recognized!\n\n", infile);
						exit(0);

					}

					if(*fileFormat==FASTA_FORMAT)
						lengthSet=1;
								
					fclose(fp);

					fileSet=1;
				}
			}
			else
			{
				fprintf(stderr, "\n ERROR: File name (string) is missing after argument -input.\n\n");
				exit(0);
			}
			continue;
		}

		if(!strcmp(argv[i], "-isfs"))
		{ 
			if (i!=argc-1)
			{	
				if(argv[i+1][0]=='-')
				{
					fprintf(stderr, "\n ERROR: File name (string) is missing after argument -isfs\n\n");
					exit(0);
				}
		
				strcpy(sfsfile,argv[++i]);

				fp=fopen(sfsfile,"r");

				if (fp==NULL)
				{
					fprintf(stderr, "\n ERROR: File %s does not exist.\n\n",sfsfile);
					exit(0);
				}
				else
				{
					*userSFS = 1;		
					fclose(fp);
				}
			}
			else
			{
				fprintf(stderr, "\n ERROR: File name (string) is missing after argument -isfs\n\n");
				exit(0);
			}	
			continue;
		}

		if(!strcmp(argv[i], "-osfs"))
		{
			if (i!=argc-1)
			{
				if(argv[i+1][0]=='-')
				{
					fprintf(stderr, "\n ERROR: File name (string) is missing after argument -osfs\n\n");
					exit(0);
				}		
				strcpy(sfsofile,argv[++i]);

				*onlySFS = 1;				
			}
			else
			{
				fprintf(stderr, "\n ERROR: File name (string) is missing after argument -osfs\n\n");
				exit(0);
			}
			continue;
		}

		if(!strcmp(argv[i], "-osf"))
		{
			if (i!=argc-1)
			{
				if(argv[i+1][0]=='-')
				{
					fprintf(stderr, "\n ERROR: File name (string) is missing after argument -osf\n\n");
					exit(0);
				}		
				strcpy(sfofile,argv[++i]);

				*onlySF = 1;				
			}
			else
			{
				fprintf(stderr, "\n ERROR: File name (string) is missing after argument -osf\n\n");
				exit(0);
			}
			continue;
		}

    		if(!strcmp(argv[i], "-grid"))
		{
			if (i!=argc-1)
			{
				if(argv[i+1][0]=='-')
				{
					fprintf(stderr, "\n ERROR: Grid number (integer) is missing after argument -grid\n\n");
					exit(0);
				}
				*grid = atoi(argv[++i]);

				if(*grid!=0)
					gridSet=1;
			} 
			else
			{
				fprintf(stderr, "\n ERROR: Grid number (integer) is missing after argument -grid\n\n");
				exit(0);
			}
			continue;
		}

		if(!strcmp(argv[i], "-length"))
		{
			if (i!=argc-1)
			{
				if(argv[i+1][0]=='-')
				{
					fprintf(stderr, "\n ERROR: Length (integer) is missing after argument -length\n\n");
					exit(0);
				}
				
				*length = atoi(argv[++i]);
			
				if(*length!=0)
					lengthSet=1;
			} 
			else
			{
				fprintf(stderr, "\n ERROR: Length (integer) is missing after argument -length\n\n");
				exit(0);
			}
			continue;
		}	
		
		if(!strcmp(argv[i], "-outgroup"))
		{	

			if (i!=argc-1)
			{
				if(argv[i+1][0]=='-')
				{
					fprintf(stderr, "\n ERROR: Outgroup name (string) is missing after argument -outgroup\n\n");
					exit(0);
				}
	
				strcpy(outgroupName, argv[++i]);

			}
			else
			{
				fprintf(stderr, "\n ERROR: Outgroup name (string) is missing after argument -outgroup\n\n");
				exit(0);
			}
			continue;

		}		    

		
		if(!strcmp(argv[i], "-folded"))
		{
			*userSetFolded = 1;
			continue;
		}

		if(!strcmp(argv[i], "-reports"))
		{
			*reports = 1;
			continue;
		}

		if(!strcmp(argv[i], "-monomorphic"))
		{
		    	*monomorphic = 1;
			continue;
		}

		
		if(!strcmp(argv[i], "-strictPolymorphic"))
		{
		    	*strictPolymorphic = 1;
			continue;
		}


    		if(!strcmp(argv[i], "-help")||!strcmp(argv[i], "-h"))
		{ 
			printHeading (stdout);

			printHelp (stdout);

			exit(0);
		}

		if(!strcmp(argv[i], "-version")||!strcmp(argv[i], "-v"))
		{ 
			printVersion (stdout);

			exit(0);
		}

		/*if(!strcmp(argv[i], "-recmap"))
		{ 
			if (i!=argc-1)
			{
				*recfile = malloc(INFILENAMESIZE*sizeof(char));
      				strcpy(*recfile, argv[++i]);
			}
      			continue; 
    		}*/


		if(!strcmp(argv[i], "-sampleList"))
		{ 
			if (i!=argc-1)
			{	
				if(argv[i+1][0]=='-')
				{
					fprintf(stderr, "\n ERROR: File name (string) is missing after argument -sampleList\n\n");
					exit(0);
				}
		
				strcpy(samplefile_i,argv[++i]);

				fp=fopen(samplefile_i,"r");

				if (fp==NULL)
				{
					fprintf(stderr, "\n ERROR: File %s does not exist.\n\n",samplefile_i);
					exit(0);
				}
				else
				{
					fclose(fp);
				}
			}
			else
			{
				fprintf(stderr, "\n ERROR: File name (string) is missing after argument -sampleList\n\n");
				exit(0);
			}	

			continue;
		}


		if(!strcmp(argv[i], "-sampleList_out"))
		{ 
			if (i!=argc-1)
			{	
				if(argv[i+1][0]=='-')
				{
					fprintf(stderr, "\n ERROR: File name (string) is missing after argument -sampleList_out\n\n");
					exit(0);
			}
		
				strcpy(samplefile_i,argv[++i]);

				*generateVCFsamplelist=1;

				// add check for existing file and ask to overwrite it
			}
			else
			{
				fprintf(stderr, "\n ERROR: File name (string) is missing after argument -sampleList_out\n\n");
				exit(0);
			}	

			continue;
		}

		if(!strcmp(argv[i], "-minsnps"))
		{
			if (i!=argc-1)
			{
				if(argv[i+1][0]=='-')
				{
					fprintf(stderr, "\n ERROR: Minimum number of SNPs (integer) is missing after argument -minsnps\n\n");
					exit(0);
				}

				*minsnps_threshold_user = atoi(argv[++i]);
			} 
			else
			{
				fprintf(stderr, "\n ERROR: Minimum number of SNPs (integer) is missing after argument -minsnps\n\n");
				exit(0);	
			}
			continue;
		}

		if(!strcmp(argv[i], "-chromList"))
		{ 
			if (i!=argc-1)
			{	
				if(argv[i+1][0]=='-')
				{
					fprintf(stderr, "\n ERROR: File name (string) is missing after argument -chromList\n\n");
					exit(0);
				}
		
				strcpy(chromfile_i,argv[++i]);

				fp=fopen(chromfile_i,"r");

				if (fp==NULL)
				{
					fprintf(stderr, "\n ERROR: File %s does not exist.\n\n",chromfile_i);
					exit(0);
				}
				else
				{
					fclose(fp);
				}
			}
			else
			{
				fprintf(stderr, "\n ERROR: File name (string) is missing after argument -chromList\n\n");
				exit(0);
			}	

			continue;
		}


		if(!strcmp(argv[i], "-chromList_out"))
		{ 
			if (i!=argc-1)
			{	
				if(argv[i+1][0]=='-')
				{
					fprintf(stderr, "\n ERROR: File name (string) is missing after argument -chromList_out\n\n");
					exit(0);
			}
		
				strcpy(chromfile_i,argv[++i]);

				*generateVCFchromlist=1;

				// add check for existing file and ask to overwrite it
			}
			else
			{
				fprintf(stderr, "\n ERROR: File name (string) is missing after argument -chromList_out\n\n");
				exit(0);
			}	

			continue;
		}	


#ifdef _USE_PTHREADS
    		if(!strcmp(argv[i], "-threads"))
		{
			if (i!=argc-1)
			{
				if(argv[i+1][0]=='-')
				{
					fprintf(stderr, "\n ERROR: Thread number (integer) is missing after argument -threads\n\n");
					exit(0);
				}

				*threads = atoi(argv[++i]);

				threadsSet = 1;
			} 
			else
			{
				fprintf(stderr, "\n ERROR: Thread number (integer) is missing after argument -threads\n\n");
				exit(0);	
			}
			continue;
		}	
#endif

#ifdef _DO_CHECKPOINTS
		if(!strcmp(argv[i], "-checkpoint"))
		{
			if (i!=argc-1)
			{
				if(argv[i+1][0]=='-')
				{
					fprintf(stderr, "\n ERROR: Checkpoint interval (integer) is missing after argument -checkpoint\n\n");
					exit(0);
				}
		
				checkPointInterval = (double)atoi(argv[++i]);				
			} 
			else
			{
				fprintf(stderr, "\n ERROR: Checkpoint interval (integer) is missing after argument -checkpoint\n\n");
				exit(0);	
			}
			continue;
		}
#endif
#ifdef _ANALYTICAL_SFS
		if(!strcmp(argv[i], "-eN"))
		{
			if (i < argc-2)
			{
				eventsTotal++;
				eventList = realloc( eventList, eventsTotal * sizeof(event) );

				if(argv[i+1][0]=='-')
				{
					fprintf(stderr, "\n ERROR: eN time\n\n");
					exit(0);
				}
				
				eventList[eventsTotal-1].t = (double)atof(argv[++i]);
				if(RESCALE)
				  eventList[eventsTotal-1].t *= 2.;

				if(argv[i+1][0]=='-')
				{
					fprintf(stderr, "\n ERROR: eN size\n\n");
					exit(0);
				}

				eventList[eventsTotal-1].d = (double)atof(argv[++i]);
			} 
			else
			{
				fprintf(stderr, "\n ERROR: -eN\n\n");
				exit(0);	
			}
			continue;
		}

		if(!strcmp(argv[i], "-G"))
		  {
		    if (i != argc-1)
			{
				
				if(argv[i+1][0]=='-')
				{
					fprintf(stderr, "\n ERROR: G growthRate\n\n");
					exit(0);
				}
				
				*growthRate = atof(argv[++i]);
				if(RESCALE)
				  (*growthRate) *= 0.5;
			} 
			else
			{
				fprintf(stderr, "\n ERROR: -eN\n\n");
				exit(0);	
			}
			continue;
		}
		
		if(!strcmp(argv[i], "-s"))
		  {
		    if (i!=argc-1)
		      {
			if(argv[i+1][0]=='-')
			  {
			    fprintf(stderr, "\n ERROR: s \n\n");
			    exit(0);
			  }
			*sequences = atoi(argv[++i]);
			
			*analyticalSFS = 1;
		      } 
		    else
		      {
			fprintf(stderr, "\n ERROR: s\n\n");
			exit(0);
		      }
		    continue;
		  }
#endif	


		fprintf(stderr, "\n ERROR: Invalid argument %s\n\n",argv[i]);
		exit(0);
	}

#ifdef _ANALYTICAL_SFS

	assert(*growthRate >= 0.);

	if(*analyticalSFS == 0 &&  *growthRate > 0)
	{
		fprintf(stdout, "\n ERROR: Please specify the sample size for the analytical SFS calculation with -s\n\n");
		exit(0);
	}

	
	if(*growthRate > 0 && eventsTotal == 1)
	{
		fprintf(stdout, "\n ERROR: Please specify a population size change event with -eN\n\n");
		exit(0);
	}
	

	if(*analyticalSFS == 0 &&  eventsTotal > 1)
	{
		fprintf(stdout, "\n ERROR: Please specify the sample size for the analytical SFS calculation with -s\n\n");
		exit(0);
	}

	if(growthRate > 0 && eventsTotal > 1)
	  {

	    eventList = realloc(eventList, (eventsTotal+1) * sizeof(event));
	    for(i=eventsTotal; i>=2; --i)
	      {
		eventList[i].d = eventList[i-1].d;
		eventList[i].t = eventList[i-1].t;
	      }
	    eventList[1].t = eventList[2].t - TOLTIME; 
	    eventList[1].d = exp(-(*growthRate) * eventList[1].t);
	    
	    eventsTotal++;


	    
	  }

	if(*growthRate > 0 && eventsTotal <= 1)
	  {
	    assert(*growthRate > 0 && eventsTotal > 1);
	  }

	if(*analyticalSFS==0 && gridSet==1 && nameSet==1 && fileSet==1)
	{
		*onlySFS=0;

		if(*onlySF==1)
			*onlySF = 2;
	}

	if(*analyticalSFS==1)
	{
		if(gridSet==1 && nameSet==1 && fileSet==1)
		{
			*analyticalSFS=2;
		}
		else
		{
			if(gridSet==1 && fileSet==0)
			{
				fprintf(stdout, "\n ERROR: Please specify an alignment with -input\n\n");
				exit(0);
			}

			if(gridSet==0 && fileSet==1)
			{
				fprintf(stdout, "\n TODO ERROR: Please specify a grid size with -grid \n\n");
				exit(0);
			}

			
		}		
	}

	if(*analyticalSFS==1)
	{
		if(sfsofile[0] == '\0')
		{
			fprintf(stdout, "\n ERROR: Please specify -osfs\n\n");
			exit(0);
		}
	}

		if(*analyticalSFS==2)
			*onlySFS=0;
	

	if (nameSet==0 && *onlySFS==0 && *onlySF!=1)
	{
		fprintf(stdout, "\n ERROR: Please specify a name for this run with -name\n\n");
		exit(0);
	}

	if (fileSet==0 && *analyticalSFS==0)
	{
		fprintf(stdout, "\n ERROR: Please specify an alignment with -input\n\n");
		exit(0);
	}

	if (gridSet==0 && *onlySFS==0 && *onlySF!=1 && *analyticalSFS==0 && *generateVCFsamplelist==0)
	{
		fprintf(stdout, "\n ERROR: Please specify the number of positions to compute the CLR at with -grid\n\n");
		exit(0);
	}

	if(lengthSet==0 && (*fileFormat==MS_FORMAT || *fileFormat==MACS_FORMAT))
	{
		fprintf(stdout, "\n ERROR: Please specify the alignment length with -length\n\n");
		exit(0);
	}

#else

	if(*analyticalSFS==0 && gridSet==1 && nameSet==1 && fileSet==1)
	{
		*onlySFS=0;

		if(*onlySF==1)
			*onlySF = 2;
	}


	if (nameSet==0 && *onlySFS==0 && *onlySF!=1)
	{
		fprintf(stdout, "\n ERROR: Please specify a name for this run with -name\n\n");
		exit(0);
	}

	if (fileSet==0)
	{
		fprintf(stdout, "\n ERROR: Please specify an alignment with -input\n\n");
		exit(0);
	}

	if (gridSet==0 && *onlySFS==0 && *onlySF!=1 && *generateVCFsamplelist==0 && *generateVCFchromlist==0)
	{
		fprintf(stdout, "\n ERROR: Please specify the number of positions to compute the CLR at with -grid\n\n");
		exit(0);
	}

	if(lengthSet==0 && (*fileFormat==MS_FORMAT || *fileFormat==MACS_FORMAT))
	{
		fprintf(stdout, "\n ERROR: Please specify the alignment length with -length\n\n");
		exit(0);
	}

#endif


#ifdef _USE_PTHREADS
	if (threadsSet==0)
	{
		fprintf(stderr, "\n ERROR: Please specify the number of threads with -threads\n\n");
		exit(0);
	}
	else
	{
		if(*threads<=1)
		{
			fprintf(stderr, "\n ERROR: Please specify a number of threads greater than 1\n\n");
			exit(0);
		}
	}
#endif


	if(*fileFormat==VCF_FORMAT)
	{
		if(*generateVCFsamplelist==1 && *generateVCFchromlist==0)
		{
			fp=fopen(samplefile_i,"r");

			if (fp!=NULL)
			{
				fclose(fp);
				
				fprintf(stdout, " \nFile %s already exists! Do you want to overwrite it? (y:yes or any other key to exit)\n", samplefile_i);
	
				if(getchar()!='y')
				{
					exit(0);
				}					
				
				
			}

		}

		if(*generateVCFchromlist==1 && *generateVCFsamplelist==0)
		{
			fp=fopen(chromfile_i,"r");

			if (fp!=NULL)
			{
				fclose(fp);
				
				fprintf(stdout, " \nFile %s already exists! Do you want to overwrite it? (y:yes or any other key to exit)\n", chromfile_i);
	
				if(getchar()!='y')
				{
					exit(0);
				}					
				
				
			}

		}

		if(*generateVCFchromlist==1 && *generateVCFsamplelist==1)
		{
			fp=fopen(chromfile_i,"r");
			fp2=fopen(samplefile_i,"r");

			if(!strcmp(chromfile_i, samplefile_i))
			{
				fprintf(stderr, " \nERROR: Filenames %s and %s are the same. \n\n", chromfile_i, samplefile_i);

				exit(0);
			}

			if (fp!=NULL && fp2==NULL)
			{
				fclose(fp);
				
				fprintf(stdout, " \nFile %s already exists! Do you want to overwrite it? (y:yes or any other key to exit)\n", chromfile_i);
	
				if(getchar()!='y')
				{
					exit(0);
				}					
				
				
			}

			if (fp2!=NULL && fp==NULL)
			{
				fclose(fp2);
				
				fprintf(stdout, " \nFile %s already exists! Do you want to overwrite it? (y:yes or any other key to exit)\n", samplefile_i);
	
				if(getchar()!='y')
				{
					exit(0);
				}					
				
				
			}

			if (fp2!=NULL && fp!=NULL)
			{
				fclose(fp);
				fclose(fp2);
				
				fprintf(stdout, " \nFiles %s and %s already exist! Do you want to overwrite them? (y:yes or any other key to exit)\n", samplefile_i, chromfile_i);
	
				if(getchar()!='y')
				{
					exit(0);
				}					
				
				
			}



		}


	}
}

void ignoreNewLineChars(FILE *fp, char *ent)
{

  int brk = 0;
  
  while( *ent == 10 || *ent == 13)
    {
      *ent = fgetc(fp);
      
      if( *ent != 10 && *ent != 13)
	{
	  ungetc(*ent, fp);
	  
	  brk = 1;
	  
	  break;
	}
    }
  if(brk == 1)
    *ent = 10;

}

int isEndOfLine2(char ent)
{
	if(ent == 10 || ent == 13)
		return 1;

	return 0;
}

int getNextString(FILE *fp, char ** word, int *readEOL, int *readEOF, int *wordLength)
{
  
  /* int iseol = 0; */
  
  *readEOL = *readEOF = 0;

  char ent = fgetc(fp);

  /* if( ent == 10 || ent == 13) */
  /*   iseol++; */

  int i=0;

  //ignoreNewLineChars(fp, &ent);

  ignoreLineSpaces(fp, &ent);
	
  if(ent == EOF)
    {
      *readEOF = 1;
	    
      (*word)[0] = '\0';
	    
      return 0;
    }
	
  if(isEndOfLine2(ent))
    {
      /* iseol++; */
      *readEOL = 1;
      (*word)[0] = '\0';
      
      /* if(iseol == 2) */
      /* 	return 2; */
      
      return 0;
    }
	
  while(isSpace(ent)==0 && isEndOfLine2(ent)==0 && ent != EOF)
    {

      if( i+1 >= (*wordLength) )
	{
	  (*wordLength) = (*wordLength) << 1; 

	  (*word) = realloc( (*word), (*wordLength) * sizeof(char) );
		    
	  //checkAlloc(*word);
	}
		
      (*word)[i++] = ent;
      
      ent = fgetc(fp);
    }
  
  (*word)[i] = '\0';
  
  if(isEndOfLine2(ent))
    *readEOL = 1;
  	
  return 1;
}

int getNextString_all_lines(FILE *fp, char ** word, int *readEOL, int *readEOF, int *wordLength)
{
  
  *readEOL = *readEOF = 0;

  char ent = fgetc(fp);

  int i=0;

  while(isSpace(ent) || isEndOfLine2(ent))
  {
	ent = fgetc(fp);
  }

	
  if(ent == EOF)
    {
      *readEOF = 1;
	    
      (*word)[0] = '\0';
	    
      return 0;
    }

	
  while(isSpace(ent)==0 && isEndOfLine2(ent)==0 && ent != EOF)
    {

      if( i+1 >= (*wordLength) )
	{
	  (*wordLength) = (*wordLength) << 1; 

	  (*word) = realloc( (*word), (*wordLength) * sizeof(char) );
		    
	}
		
      (*word)[i++] = ent;
      
      ent = fgetc(fp);
    }
  
  (*word)[i] = '\0';
  
  if(isEndOfLine2(ent))
    *readEOL = 1;
  	
  return 1;
}

int findFirstAlignment(FILE *fp, int format, FILE *fpVCFsamples, int generateVCFsamplelist, char * vcf_samples_filename, FILE * fpInfo)
{
	char tmp;
	int temp=-1, tip = 0, i;

	if(format==FASTA_FORMAT)
	{
		tmp=fgetc(fp);

		while(tmp!=EOF)
		{
			if(tmp=='>') // first sequence
				return 1;
			else
				tmp = fgetc(fp);
		}
	}
	else if(format == MS_FORMAT)
	{
	    	tmp=fgetc(fp);

		while(tmp!=EOF)
		{
			if(tmp=='/')
			{
				tmp = fgetc(fp);

				if(tmp=='/')
					return 1;
				else
					tmp = fgetc(fp);				
			}
			else
				tmp = fgetc(fp);
		}
	}
	else if(format == MACS_FORMAT)
	{
		char word[1000];
		int sitenumber;

		while( fscanf(fp, "%s", word) )
		{
			if(!strcmp(word, "SITE:") )
			{
				temp = fscanf(fp, "%d", &sitenumber);
				
				assert(temp==1);
				
				if(sitenumber == 0)
					return 1;
			}
		}
	}
	else if(format == SF_FORMAT)
	{
		tmp=fgetc(fp);	    

		ignoreAll(fp, &tmp);
		    
		goToLineEnd(fp);

		return 1;
	}
	else if(format == VCF_FORMAT)
	{

		char ** string = (char **) malloc (sizeof(char*));
		
		(*string) = (char *) malloc(sizeof(char)*STRINGLENGTH);
		
		int status=-1, eol=0, eof=0, maxLength=STRINGLENGTH;

		VCF_header_lines = 0;

		while(1)
		{
			tip = status = getNextString (fp, string, &eol, &eof, &maxLength);
			

			if( eol == 1)
			  VCF_header_lines++;
			
			if(status==1)
				if(strcmp((*string),"#CHROM")==0)
				{
					VCF_header_lines++;
					break;
				}

			if(eof==1)
				assert(0);		
		}

		
		char * headerFields[VCF_HLENGTH];
		int fieldInd=1;
		int VCFsamples=0;
		
		headerFields[0] = "#CHROM";
		headerFields[1] = "POS";
		headerFields[2] = "ID";
		headerFields[3] = "REF";
		headerFields[4] = "ALT";
		headerFields[5] = "QUAL";
		headerFields[6] = "FILTER";
		headerFields[7] = "INFO";
		headerFields[8] = "FORMAT";


		
		while((tip = getNextString (fp, string, &eol, &eof, &maxLength) ) ==1 )
		{

			if(fieldInd<VCF_HLENGTH)
				if(strcmp(headerFields[fieldInd],(*string))!=0)
				{
					fprintf(stderr, "\n\n ERROR: VCF header field %s is missing.\n\n",headerFields[fieldInd]);
					assert(0);
	
				}
			
			fieldInd++;

			if(fieldInd>=VCF_HLENGTH)	
				break;

			if(eol==1 || eof==1)
				assert(0);			

		}
		
		int sampleList_size=1, sampleList_index=0;
		char ** sampleList = (char **) malloc (sizeof(char *)*sampleList_size);
		
		while(getNextString (fp, string, &eol, &eof, &maxLength)==1)
		{
			if(sampleList_index==sampleList_size)
			{
				sampleList_size++;
				sampleList = realloc(sampleList, sampleList_size * sizeof(char*));
			}

			sampleList[sampleList_index] = malloc(sizeof(char)*STRINGLENGTH);
			strcpy(sampleList[sampleList_index++], *string);

			VCFsamples++;

//			printf("Sample: %s \n",string[0]);	

			if(alignment->outgroupName!=NULL)
				if(strcmp(alignment->outgroupName, *string)==0)
					alignment->outgroupIndex = VCFsamples;
		
			if(eol==1)
				break;

			if(eof==1)
				assert(0);			
		
		}

		if(generateVCFsamplelist==1)
		{
			for(i=0;i<VCFsamples;i++)
			{
				fprintf(fpVCFsamples,"%s\n",sampleList[i]);
			}

			fprintf(stdout, "\n\n A list of VCF samples has been stored in file %s!\n\n",vcf_samples_filename);

			return 0;
		}

		
		fprintf(stdout, " Total number of samples in the VCF:\t%d\n",VCFsamples);
		fprintf(fpInfo, " Total number of samples in the VCF:\t%d\n",VCFsamples);

		alignment->VCFsamples = VCFsamples;
		alignment->VCFsample_valid = malloc(sizeof(int)*VCFsamples);

		for(i=0;i<VCFsamples;i++)
			alignment->VCFsample_valid[i]=0; // set all samples to invalid


		if(fpVCFsamples==NULL)
		{
			for(i=0;i<VCFsamples;i++)
				alignment->VCFsample_valid[i]=1; // set all samples to valid


			fprintf(stdout, " Samples excluded from the analysis:\t0\n\n\n");
			fprintf(fpInfo, " Samples excluded from the analysis:\t0\n\n\n");
		}
		else
		{
			int validVCFsamples_matched=0;

			while(getNextString_all_lines (fpVCFsamples, string, &eol, &eof, &maxLength)==1)
			{
				//printf("%s\n",*string);

				for(i=0;i<VCFsamples;i++)
				{
					if(!strcmp(*string, sampleList[i]))
					{
						//printf("Match at %d %s %s\n",i, *string, sampleList[i]);
						alignment->VCFsample_valid[i]=1;
						validVCFsamples_matched++;
						break;
					}
				}

				if(i==VCFsamples)
				{
					fprintf(stdout, " *** NOT FOUND: Sample \"%s\" does not exist in the VCF file ***\n",*string);
					fprintf(fpInfo, " *** NOT FOUND: Sample \"%s\" does not exist in the VCF file ***\n",*string);
				}
			
			}

				fprintf(stdout, " Samples excluded from the analysis:\t%d\n\n\n",VCFsamples-validVCFsamples_matched);
				fprintf(fpInfo, " Samples excluded from the analysis:\t%d\n\n\n",VCFsamples-validVCFsamples_matched);
		}
	
		assert(VCFsamples!=0);

		tmp = fgetc(fp);
		
		ungetc(tmp, fp);
		
		ignoreNewLineChars(fp, &tmp);
		

		tip = getNextString (fp, string, &eol, &eof, &maxLength);

		strncpy(VCF_alignment_name, *string, MAX_CHROM_NAME_VCF);

		assert(strlen(VCF_alignment_name)!=0);

		free(*string);
		free(string);
			
		return 1;
	}


	return 0;		
}




int findNextAlignment(FILE *fp, int fileFormat)
{
	char stop,tmp;
	int temp=-1;

	if (fileFormat == MS_FORMAT || fileFormat == SF_FORMAT)
		stop = '/';
	else
		stop = '>';
	
	
	if(fileFormat == MS_FORMAT || fileFormat == FASTA_FORMAT || fileFormat == SF_FORMAT)
	  {
	    tmp=fgetc(fp);
	    
	    while(tmp!=EOF)
	      {
		if(tmp==stop)
		  return 1;	
		else
		  {
		    tmp = fgetc(fp);
		    
		  }
	      }	    	    
	  }
	
	if(fileFormat == MACS_FORMAT)
	  {
	    char word[100];
	    int sitenumber;
	    
	    int nextAl=1;
	    
	    while( (nextAl = fscanf(fp, "%s", word)))
	      {
		
		if( nextAl <= 0)
		  return 0;
		
		if(!strcmp(word, "SITE:") )
		  {
		    temp = fscanf(fp, "%d", &sitenumber);

		    assert(temp==1);
		    
		    if(sitenumber == 0)
		      return 1;
		  }
	      }
	  }

	
	if(fileFormat == VCF_FORMAT)
	{
		if(nxtVCFalignment==0)
			return 1;
	}

	
	return 0;
}



int mapNuclToInt(char a)
{
  if(a == AD || a == ade) 	return 0;
  if(a == CY || a == cyt) 	return 1;
  if(a == GU || a == gua) 	return 2;
  if(a == TH || a == thy) 	return 3;
  return -1;
}


int mapCharToInt(char a)
{
	if(a == ZERO) 	return 0;
	if(a == ONE) 	return 1;
	if(a=='X' || a=='K' || 
	   a=='M' || a=='R' || 
	   a=='Y' || a=='S' || 
	   a=='W' || a=='B' || 
	   a=='V' || a=='H' || 
	   a=='D' || a=='.') return 2;
	if(a == GAP || a == UN ) return 2;
	if(a == AD || a == ade) 	return 3;
	if(a == CY || a == cyt) 	return 4;
	if(a == GU || a == gua) 	return 5;
	if(a == TH || a == thy) 	return 6;
	return 7;
}




char getCharBIN (char input, char state0, char state1)
{
	if(input==GAP)
 		return GAP;
	
	if(input == UN)
	  	return UN;

	if(input==state0)
		return ZERO;

	if(input==state1)
		return ONE;

	return 'X';
}



void initializeAlignmentVariables()
{
  alignment->positionsInd = NULL;
  alignment->x = NULL;
  alignment->n = NULL;
  alignment->folded = NULL;
  alignment->outgroupName = NULL;
  alignment->outgroupSequence = NULL;
  alignment->n_pat = NULL;
  alignment->x_pat = NULL;
  alignment->f_pat = NULL;
  alignment->c_pat = NULL;
  alignment->baseLikelihood = NULL;
  if(alignment->patternList!=NULL)
    {
      struct pattern * cur = alignment->patternList;
      
      while(cur!=NULL)
	{
	  alignment->patternList = cur->nxt;
	  
	  free(cur);
	  
	  cur = alignment->patternList;
	}
      
      alignment->patternList = NULL;	
    }
  
  alignment->gridProbs = NULL;
  
  alignment->gridADs = NULL;
  
  rvLUT = NULL;

  clr=NULL;
    
  //  factLUT->val = NULL;
  
  factLUT = NULL;
}


void freeAlignment()
{
	if(alignment->positionsInd != NULL)
	{
		free(alignment->positionsInd);
		alignment->positionsInd = NULL;
	}

	if(alignment->x != NULL)
	{
		free(alignment->x);
		alignment->x = NULL;
	}

	if(alignment->n != NULL)
	{
		free(alignment->n);
		alignment->n = NULL;
	}

	if(alignment->folded != NULL)
	{
		free(alignment->folded);
		alignment->folded = NULL;
	}

	if(alignment->outgroupName != NULL)
	{
		free(alignment->outgroupName);
		alignment->outgroupName = NULL;
	}

	if(alignment->outgroupSequence != NULL)
	{
		free(alignment->outgroupSequence);
		alignment->outgroupSequence = NULL;
	}

	if(alignment->n_pat != NULL)
	{
		free(alignment->n_pat);
		alignment->n_pat = NULL;
	}

	if(alignment->x_pat != NULL)
	{
		free(alignment->x_pat);
		alignment->x_pat = NULL;
	}

	if(alignment->f_pat != NULL)
	{
		free(alignment->f_pat);
		alignment->f_pat = NULL;
	}

	if(alignment->c_pat != NULL)
	{
		free(alignment->c_pat);
		alignment->c_pat = NULL;
	}

	if(alignment->baseLikelihood!=NULL)
	{
		free(alignment->baseLikelihood);
		alignment->baseLikelihood = NULL;
	}

	if(alignment->patternList!=NULL)
	{
		struct pattern * cur = alignment->patternList;
	
		while(cur!=NULL)
		{
			alignment->patternList = cur->nxt;

			free(cur);

			cur = alignment->patternList;
		}

		alignment->patternList = NULL;	
	}
	
	int i,j;
	
	int minx_offset = alignment->SFS[0]>0.0?0:1;
	int maxx_offset = alignment->SFS[alignment->sequences]>0.0?1:0;
	
	if(alignment->gridProbs!=NULL)
	{
		for(i=0;i<alignment->maxn-alignment->minn+1;++i)
		{
			if(alignment->gridProbs[i]!=NULL)
			{
				for(j=0;j<alignment->sequences+maxx_offset-minx_offset+1;++j)
				{
					if(alignment->gridProbs[i][j]!=NULL)
					{
						free(alignment->gridProbs[i][j]);
						alignment->gridProbs[i][j] = NULL;
					}
				}

				free(alignment->gridProbs[i]);
				alignment->gridProbs[i] = NULL;
			}
		}

		free(alignment->gridProbs);
		alignment->gridProbs = NULL;
	}

	if(alignment->gridADs!=NULL)
	{
		free(alignment->gridADs);
		alignment->gridADs = NULL;
	}

	if(rvLUT!=NULL)
	{
		for(i=0;i<GRIDSIZE;i++)
		{
			if(rvLUT[i]!=NULL)
			{
				free(rvLUT[i]);
				rvLUT[i] = NULL;
			}
		}	

		free(rvLUT);
		rvLUT = NULL;
	}

	if(clr!=NULL)
	{
		free(clr);
		clr=NULL;
	}

	if(factLUT != NULL && factLUT->val!=NULL)
	{
		free(factLUT->val);
		factLUT->val = NULL;
	}

	if(factLUT!=NULL)
	{
		free(factLUT);
		factLUT = NULL;
	}
}



void freeAlignment_noSNPs()
{

  int i;
	if(alignment->positionsInd != NULL)
	{
		free(alignment->positionsInd);
		alignment->positionsInd = NULL;
	}

	if(alignment->x != NULL)
	{
		free(alignment->x);
		alignment->x = NULL;
	}

	if(alignment->n != NULL)
	{
		free(alignment->n);
		alignment->n = NULL;
	}

	if(alignment->folded != NULL)
	{
		free(alignment->folded);
		alignment->folded = NULL;
	}

	if(alignment->outgroupName != NULL)
	{
		free(alignment->outgroupName);
		alignment->outgroupName = NULL;
	}

	if(alignment->outgroupSequence != NULL)
	{
		free(alignment->outgroupSequence);
		alignment->outgroupSequence = NULL;
	}

	if(alignment->n_pat != NULL)
	{
		free(alignment->n_pat);
		alignment->n_pat = NULL;
	}

	if(alignment->x_pat != NULL)
	{
		free(alignment->x_pat);
		alignment->x_pat = NULL;
	}

	if(alignment->f_pat != NULL)
	{
		free(alignment->f_pat);
		alignment->f_pat = NULL;
	}

	if(alignment->c_pat != NULL)
	{
		free(alignment->c_pat);
		alignment->c_pat = NULL;
	}

	if(alignment->baseLikelihood!=NULL)
	{
		free(alignment->baseLikelihood);
		alignment->baseLikelihood = NULL;
	}

	if(alignment->patternList!=NULL)
	{
		struct pattern * cur = alignment->patternList;
	
		while(cur!=NULL)
		{
			alignment->patternList = cur->nxt;

			free(cur);

			cur = alignment->patternList;
		}

		alignment->patternList = NULL;	
	}
	

	if(alignment->gridProbs!=NULL)
	{

	  free(alignment->gridProbs);
	  alignment->gridProbs = NULL;
	}
	
	if(alignment->gridADs!=NULL)
	  {
		free(alignment->gridADs);
		alignment->gridADs = NULL;
	  }

	if(rvLUT!=NULL)
	  {
		for(i=0;i<GRIDSIZE;i++)
		{
		  if(rvLUT[i]!=NULL)
			{
			  free(rvLUT[i]);
			  rvLUT[i] = NULL;
			}
		}	
		
		free(rvLUT);
		rvLUT = NULL;
	  }
	
	if(clr!=NULL)
	  {
	    free(clr);
	    clr=NULL;
	}
	
	if(factLUT != NULL && factLUT->val!=NULL)
	  {
	    free(factLUT->val);
	    factLUT->val = NULL;
	}
	
	if(factLUT!=NULL)
	{
		free(factLUT);
		factLUT = NULL;
	}
}


char getCharacterImputeBIN (int *input)
{
  	char states[2] = {'0', '1'};
  
  	int total = input[0] + input[1];
  
  	int rn = rand() % total + 1;

 	int j = 0;
  
  	int current = input[j];
  
  	while(rn > current)
      		current += input[++j];
  
  	return states[j];
}




char getCharacterImputeDNA(int *input)
{
  	char states[4] = {'A', 'C', 'G', 'T'};
  
  	// For DNA, the numbers of A, C, G, T are in positions 3,4,5,6 
  	int total = input[3] + input[4] + input[5] + input[6];
  
  	int rn = rand() % total + 1;

  	int j = 3;

  	int current = input[j];  

  	while(rn > current)
      		current += input[++j];
  
  	return states[j-3];
}



void switchValues (int * input, char * sortedStates, int index0, int index1)
{
	int tmp;
        char tmpC;

	if(input[index1]>input[index0])
	{
		tmp = input[index0];
		tmpC = sortedStates[index0];
		input[index0] = input[index1];
		sortedStates[index0] = sortedStates[index1];
		input[index1] = tmp;
   		sortedStates[index1] = tmpC;
	}
}

void sort(int * freqs, char * sortedStates)
{
	switchValues (freqs, sortedStates, 0, 1);
	switchValues (freqs, sortedStates, 2, 3);
	switchValues (freqs, sortedStates, 0, 2);
	switchValues (freqs, sortedStates, 1, 3);
	switchValues (freqs, sortedStates, 1, 2);	
}

int minorToMajor(int* inp, char* outp)
{
 	char major, alpha[4] = {'A', 'C', 'G', 'T'};
   
 	int total = 0, rv, i;
  
	sort(inp,alpha); 
  
	assert(inp[0] >= inp[1] && inp[1] >= inp[2] && inp[2] >= inp[3]);
	
  	if(inp[2] + inp[3] == 0)
    		return -1;

  	total = inp[0] + inp[1];
  
  	for(i=0; i<2; ++i)
    	{
	      rv= rand() % total + 1;
	      major = (rv <= inp[0] ) ? alpha[0] : alpha[1];
	      outp[i] = alpha[ i+2 ];
	      outp[i+2] = major;
    	}

    	return 1;
}




int isEOFnear(FILE *fp, char *ent)
{
	while(*ent < 33 && *ent != EOF)
		*ent = fgetc(fp);

	if(*ent == EOF)
		return 1;

	return 0;
}

void ignoreSpaces(FILE * fp, char * ent)
{
	while(*ent != EOF && (*ent==10 || *ent==13 || *ent==32 || *ent==9) )
		*ent = fgetc(fp);
}

void print2DIntArray(int **array, int n1, int n2)
{
  int i,j;
  for(i=0; i<n1; ++i)
    {
      for(j=0; j<n2; ++j)
	printf("%d,", array[i][j]);
      printf("\n");
    }
}

void printIntArray(int *array, int n)
{
	int i;
	for( i = 0; i<n; ++i)
	{
	  printf("%d,", array[i]);
	}
}

void printUCharArray(unsigned char *array, int n)
{
	int i;
	for( i = 0; i<n; ++i)
	{
	  printf("%u,", array[i]);
	}
}


void printFloatArray(float *array, int n)
{
	int i;
	for( i = 0; i<n; ++i)
	{
	  printf("%f,", array[i]);
	}
}


void printCharArray(char *array, int n)
{
	int i;
	for( i = 0; i<n; ++i)
	{
		printf("%c,", array[i]);
	}
}



int compare(const void * a, const void * b)
{
  return ( *(int*)a - *(int*)b );
}



int matchIndex(char *s, char *q)
{
	int len1 = strlen(s);
	int len2 = strlen(q);

	int i,j = 0;

	for(i=0; i<len1; ++i)
	{
		for(j=0; j<len2; ++j)
		{
			fprintf(stderr, "%c - %c\n", s[i], q[j]);

			if( s[i] == q[j] && i < len1 - 1)
				i++;
			else
				break;
		}

		//fprintf(stderr, "j is: %d and length is: %d\n", j, len2);

		if(j == len2 )
			return i;
	}


	return -1;
}


int parseInfoFreq(char *word, float *freq)
{
	int len = strlen(word);
	int i = 0, j = 0;

	char * p = "AF=";
	char freqString[100];

	int ind = matchIndex(word, p);

	//fprintf(stderr, "%d\n", ind);

	if(ind > 0)
	{
		i = ind;
		while(i < len && word[i] != ';')
		{
			freqString[j++] = word[i++];
//			fprintf(stderr, "%c.", word[i-1]);
		}
	}
	else
		{
			return 0;
		}

	*freq = atof(freqString);

	return 1;

}

int parseInfoVar(char *word, char *type, int start)
{
	int len = strlen(word);
	int i = start, j = 0;

	char * p = "VT=";
	//char freqString[100];

	int ind = matchIndex(word, p);

	//fprintf(stderr, "%d\n", ind);

	if(ind > 0)
	{
		i = ind;
		while(i < len && word[i] != ';')
		{
			type[j++] = word[i++];
//			fprintf(stderr, "%c.", word[i-1]);
		}
	}
	else
	{
		return 0;
	}

	return 1;

}























//
//void parseInfo(FILE *fp, char * flag, char *word)
//{
//	fscanf(fp, "%99cAF=%[^;];", word);
//}

void goToLineEnd(FILE *fp)
{
	char ent;
	while( (ent = fgetc(fp) ) != EOF )
	  {
	    
	    if(isEndOfLine(ent))
	      break;
	  }
}


int exactMatch(char * w1, char *w2)
{

	int i;

	int len = strlen(w1);
	if(len != strlen(w2))
		return 0;

	for(i=0; i<len; ++i)
	{
		if(w1[i]!=w2[i])
			return 0;
	}
	return 1;

}



int getStates(char *word, char* states, int start, int *end)
{
	int len = strlen(word);
	int i, j=start, sep = 1;

	for(i=0; i<len; ++i)
	{
		if(word[i] == ',')
		{
			sep = 1;
			++j;
		}
		else
			if(sep == 1 && isValidDNACharacter(word[i]) == 1)
			{
				sep = 0;
				states[j] = word[i];
			}
			else
				return 0;
	}

	*end = j+1;

	return 1;
}

int isValidDNACharacter(char input)
{
	if(input==AD) return 1;
	if(input==CY) return 1;
	if(input==GU) return 1;
	if(input==TH) return 1;
	if(input==GAP) return 1;
	if(input==UN) return 1;
	if(input == ade) return 1;
	if(input == cyt) return 1;
	if(input == gua) return 1;
	if(input == thy) return 1;
	

	if(input=='X' || input=='K' || 
	   input=='M' || input=='R' || 
	   input=='Y' || input=='S' || 
	   input=='W' || input=='B' || 
	   input=='V' || input=='H' || 
	   input=='D' || input=='.') return 1;

	return 0;
} 

int isACGT(char input)
{
  	if(input==AD || input == ade) return 1;
	if(input==CY || input == cyt) return 1;
	if(input==GU || input == gua) return 1;
	if(input==TH || input == thy) return 1;

	return 0;
}

int isBinary(char input)
{
  	if(input==ONE) return 1;
	if(input==ZERO) return 1;
	
	return 0;
}



char filterAmb (char input)
{
	if(input=='X' || input=='K' || 
	   input=='M' || input=='R' || 
	   input=='Y' || input=='S' || 
	   input=='W' || input=='B' || 
	   input=='V' || input=='H' || 
	   input=='D' || input=='.') return UN;

	return input;
}



void removeEndGaps (char * line)
{
	int i, len = strlen(line);

	for(i=len-1;i>-1;i--)
	{
		if(isSpace(line[i]) || isEndOfLine(line[i]))
			line[i] = '\0';
      		else
			break;
	}
}



int isOutgroup (char * seqName, char * outgroup)
{
	if(outgroup==NULL)
		return 0;

	if(!strcmp(seqName,outgroup))
		return 1;

	return 0;
}


void updateStateCounters (alignment_struct * alignment, char ent, int index)
{
	if(ent==AD || ent == ade)
	{
		++alignment->stateA[index];
		return;
	}

	if(ent==CY || ent == cyt)
	{
		++alignment->stateC[index];
		return;
	}

	if(ent==GU || ent == gua)
	{
		++alignment->stateG[index];
		return;
	}

	if(ent==TH || ent == thy)
	{
		++alignment->stateT[index];
		return;
	}		
}

int readFASTAsequenceLight(FILE * fp, alignment_struct * alignment, char * seqName, int outgroup)
{	
	int DIM = 10, i=0, y=0;
	char ent;

	if(alignment->segsites==-1)
	{	
		alignment->outgroupSequence = calloc(DIM, sizeof(char));
		assert(alignment->outgroupSequence!=NULL);

		alignment->stateA = calloc (DIM, sizeof(al_t));
		assert(alignment->stateA!=NULL);

		alignment->stateC = calloc (DIM, sizeof(al_t));
		assert(alignment->stateC!=NULL);

		alignment->stateG = calloc (DIM, sizeof(al_t));
		assert(alignment->stateG!=NULL);

		alignment->stateT = calloc (DIM, sizeof(al_t));
		assert(alignment->stateT!=NULL);
	}
	else
	{
		DIM = alignment->segsites;
	}


	ent = toupper(fgetc(fp));

	ignoreSpaces(fp, &ent);

	while(ent != EOF && ent != '>' && ent != '/')
	{
		if(isValidDNACharacter(ent))
		{		
	
			if(outgroup==1)
			{
				alignment->outgroupSequence[i] = ent;
				//updateStateCounters (alignment, ent, i);
			}		
			else
				updateStateCounters (alignment, ent, i); // TODO: check whether we need to call this for the outgroup as well
		
			++i;			
		}
		else
		{		
			if(ent >= 33 && ent <= 126)
			{
				fprintf(stderr,"\n ERROR: Invalid character (%c) at position %d of sequence %s.\n\n", ent, i, seqName);
				exit(0);				
			}
		}
		
		ent = toupper(fgetc(fp));

		if(i==DIM)
		{
			DIM = DIM + 10;

			alignment->outgroupSequence = realloc(alignment->outgroupSequence, DIM * sizeof(char));
			assert(alignment->outgroupSequence != NULL);

			alignment->stateA = realloc(alignment->stateA, DIM * sizeof(al_t));
			assert(alignment->stateA != NULL);

			alignment->stateC = realloc(alignment->stateC, DIM * sizeof(al_t));
			assert(alignment->stateC != NULL);

			alignment->stateG = realloc(alignment->stateG, DIM * sizeof(al_t));
			assert(alignment->stateG != NULL);

			alignment->stateT = realloc(alignment->stateT, DIM * sizeof(al_t));
			assert(alignment->stateT != NULL);

			for(y=DIM-10;y<DIM;y++)
			{
				alignment->outgroupSequence[y] = 0;
				alignment->stateA[y]=0;
				alignment->stateC[y]=0;
				alignment->stateG[y]=0;
				alignment->stateT[y]=0;
			}
		}
	}

	if(alignment->segsites==-1)
		alignment->segsites = i;
	
	if(alignment->segsites != -1 && i != alignment->segsites)
	{
		fprintf(stderr, "\n ERROR: Length of sequence %s is %d. Expected %d.\n\n", seqName, i, alignment->segsites);
		exit(0);	
	}

	alignment->length = alignment->segsites;	

	if(ent == '>')
		return 1;

	return 0;
}

void createXNFvectors (alignment_struct * alignment)
{
	alignment->n = malloc(sizeof(al_t)*alignment->segsites);
	assert(alignment->n!=NULL);

	alignment->x = malloc(sizeof(al_t)*alignment->segsites);
	assert(alignment->x!=NULL);

	alignment->folded = malloc(sizeof(al_t)*alignment->segsites);
	assert(alignment->folded!=NULL);

	int i;

	for(i=0;i<alignment->segsites;i++)
	{
		//printf("A: %d C: %d G: %d T: %d -- %c\n",alignment->stateA[i],alignment->stateC[i],alignment->stateG[i],alignment->stateT[i], alignment->outgroupSequence[i]);
		alignment->n[i] = alignment->stateA[i] + alignment->stateC[i] + alignment->stateG[i] + alignment->stateT[i];

		// no outgroup
		if(alignment->outgroupSequence[i]==0) 
		{
			alignment->folded[i] = 1;

			if(alignment->stateA[i]!=0)
			{
				alignment->x[i] = alignment->stateA[i];
				continue;
			}
			if(alignment->stateC[i]!=0)
			{
				alignment->x[i] = alignment->stateC[i];
				continue;
			}
			if(alignment->stateG[i]!=0)
			{
				alignment->x[i] = alignment->stateG[i];
				continue;
			}
			if(alignment->stateT[i]!=0)
			{
				alignment->x[i] = alignment->stateT[i];
				continue;
			}
			alignment->x[i]=0;
			continue;			
		}

		// outgroup
		if(alignment->outgroupSequence[i]==AD || alignment->outgroupSequence[i]== ade)
		{
			if(alignment->stateA[i]!=0)
			{
				alignment->folded[i] = 0;
				alignment->x[i] = alignment->n[i]-alignment->stateA[i];
				/*printf("-- %d --   segsites,x,n,stateA,position: %d %d %d %d\n", i, alignment->segsites, alignment->x[i], alignment->n[i], alignment->stateA[i]); */
			}
			else
			{
				alignment->folded[i] = 1;

				if(alignment->stateC[i]!=0)
				{
					alignment->x[i] = alignment->stateC[i];
					continue;
				}
				if(alignment->stateG[i]!=0)
				{
					alignment->x[i] = alignment->stateG[i];
					continue;
				}
				if(alignment->stateT[i]!=0)
				{
					alignment->x[i] = alignment->stateT[i];
					continue;
				}
				alignment->x[i]=0;
			}

			continue;
		}
		if(alignment->outgroupSequence[i]==CY || alignment->outgroupSequence[i]==cyt)
		{

			if(alignment->stateC[i]!=0)
			{
				alignment->folded[i] = 0;
				alignment->x[i] = alignment->n[i]-alignment->stateC[i];
			}
			else
			{
				alignment->folded[i] = 1;

				if(alignment->stateA[i]!=0)
				{
					alignment->x[i] = alignment->stateA[i];
					continue;
				}
				if(alignment->stateG[i]!=0)
				{
					alignment->x[i] = alignment->stateG[i];
					continue;
				}
				if(alignment->stateT[i]!=0)
				{
					alignment->x[i] = alignment->stateT[i];
					continue;
				}
				alignment->x[i]=0;
			}

			continue;
		}
		if(alignment->outgroupSequence[i]==GU || alignment->outgroupSequence[i]==gua)
		{
			//--alignment->n[i];

			if(alignment->stateG[i]!=0)
			{
				alignment->folded[i] = 0;
				alignment->x[i] = alignment->n[i]-alignment->stateG[i];
			}
			else
			{
				alignment->folded[i] = 1;

				if(alignment->stateC[i]!=0)
				{
					alignment->x[i] = alignment->stateC[i];
					continue;
				}
				if(alignment->stateA[i]!=0)
				{
					alignment->x[i] = alignment->stateA[i];
					continue;
				}
				if(alignment->stateT[i]!=0)
				{
					alignment->x[i] = alignment->stateT[i];
					continue;
				}
				alignment->x[i]=0;
			}
			continue;
		}

		if(alignment->outgroupSequence[i]==TH || alignment->outgroupSequence[i]== thy)
		{
			//--alignment->n[i];
		
			if(alignment->stateT[i]!=0)
			{
				alignment->folded[i] = 0;
				alignment->x[i] = alignment->n[i]-alignment->stateT[i];
			}
			else
			{
				alignment->folded[i] = 1;

				if(alignment->stateC[i]!=0)
				{
					alignment->x[i] = alignment->stateC[i];
					continue;
				}
				if(alignment->stateG[i]!=0)
				{
					alignment->x[i] = alignment->stateG[i];
					continue;
				}
				if(alignment->stateA[i]!=0)
				{
					alignment->x[i] = alignment->stateA[i];
					continue;
				}
				alignment->x[i]=0;
			}

			continue;
		}

		// Ambiguity character

		alignment->folded[i] = 1;

		if(alignment->stateA[i]!=0)
		{
			alignment->x[i] = alignment->stateA[i];
			continue;
		}
		if(alignment->stateC[i]!=0)
		{
			alignment->x[i] = alignment->stateC[i];
			continue;
		}
		if(alignment->stateG[i]!=0)
		{
			alignment->x[i] = alignment->stateG[i];
			continue;
		}
		if(alignment->stateT[i]!=0)
		{
			alignment->x[i] = alignment->stateT[i];
			continue;
		}
		alignment->x[i]=0;		
	}

	/* if(alignment->stateA!=NULL) */
	/* 	free(alignment->stateA); */
	/* if(alignment->stateC!=NULL) */
	/* 	free(alignment->stateC); */
	/* if(alignment->stateG!=NULL) */
	/* 	free(alignment->stateG); */
	/* if(alignment->stateT!=NULL) */
	/* 	free(alignment->stateT); */


}

void skipToNextAlignmentDelimiter(FILE *fp, char d)
{
  char ent = fgetc( fp );
  
  while( ent != d && ent != EOF)
    {
      ent = fgetc( fp );
    }
  
  
}

/* int skipLine (FILE * fp) */
/* { */
/*   char tmp; */

/*   while( (tmp = fgetc(fp) ) != EOF) */
/*     if(tmp == 10 || tmp == 13) */
/*       break; */
  
/*   if(tmp == EOF) */
/*     return 1; */

/*   return 0; */
/* } */



int readAlignmentFASTA(FILE *fp, alignment_struct *alignment, FILE *fpInfo, FILE *fpSFo, int minsnps_threshold_user, int alignmentIndex)
{
	fprintf(stdout," Alignment %d\n",alignmentIndex);			
	fprintf(fpInfo," Alignment %d\n",alignmentIndex);

	int nxt_seq = 1, outgroup=0, outgroupCounter=0, i;

	char seqName[SEQNAMESIZE];
	
	char * sucseqName = NULL;
	
	alignment->segsites = -1;
	alignment->sequences = 0;

	alignment->maxx = -1;
	alignment->maxn = -1;

	alignment->minx = MAXINT;
	alignment->minn = MAXINT;
	
	alignment->stateA = NULL; 
	alignment->stateC = NULL;
	alignment->stateG = NULL;
	alignment->stateT = NULL;
	
	while(nxt_seq==1)
	{
		sucseqName = fgets(seqName, SEQNAMESIZE, fp);
		assert(sucseqName!=NULL);
		removeEndGaps (seqName);

		alignment->sequences++;
		
		outgroup = 0;
		if(isOutgroup(seqName, alignment->outgroupName)==1)
		{
			outgroup = 1;
			--alignment->sequences;
			outgroupCounter++;
		}		

		nxt_seq = readFASTAsequenceLight(fp, alignment, seqName, outgroup);
	}

	assert(outgroupCounter<=1);
	//printf("segsites: %d length: %d\n",alignment->segsites, alignment->length);
	assert(alignment->segsites<=alignment->length);

	createXNFvectors (alignment);

	alignment->positionsInd = malloc(sizeof(int)*alignment->segsites);

	for(i=0;i<alignment->segsites;i++)
		alignment->positionsInd[i]=i+1;

	if(outgroupCounter==0)
	{
		fprintf(stdout,"\t\tOutgroup:\t\tnone\n");
		fprintf(fpInfo,"\t\tOutgroup:\t\tnone\n");
	}
	else
	{
		fprintf(stdout,"\t\tOutgroup:\t\t%s\n", alignment->outgroupName);
		fprintf(fpInfo,"\t\tOutgroup:\t\t%s\n", alignment->outgroupName);
	}

	//fprintf(stdout,"\t\tSequences:\t\t%d\n\t\tSites:\t\t\t%d", alignment->sequences,alignment->segsites);

//	fprintf(fpInfo,"\t\tSequences:\t\t%d\n\t\tSites:\t\t\t%d", alignment->sequences,alignment->segsites);
//			fprintf(stderr,"\n ERROR: Invalid character (%c) at position %d of sequence %s.\n\n", ent, i, seqName);
//				exit(0);

	for(i=0;i<alignment->segsites;i++)
	{
		if(alignment->n[i]<alignment->minn)
			alignment->minn = alignment->n[i];

		if(alignment->x[i]<alignment->minx)
			alignment->minx = alignment->x[i];

		if(alignment->n[i] > alignment->sequences)
			alignment->sequences = alignment->n[i];

		alignment->maxn = alignment->sequences;		

		alignment->maxx = alignment->sequences;
	}

	for(i=0;i<alignment->segsites;i++)
		if((alignment->x[i]==alignment->n[i]) || (alignment->folded[i]==1 && alignment->x[i]==0))
		{
			alignment->maxx++;
			break;
		}


	//alignment->positionsInd = malloc(alignment->segsites*sizeof(int));

	//for(i=0;i<alignment->segsites;i++)
	//  alignment->positionsInd[i] = (int)(alignment->positions[i] * (float)alignment->length);

	//free(alignment->positions);


	fprintf(stdout,"\n\t\tSequences:\t\t%d",alignment->sequences);	
	fprintf(fpInfo,"\n\t\tSequences:\t\t%d",alignment->sequences);

	fprintf(stdout,"\n\t\tSites:\t\t\t%d",alignment->segsites);	
	fprintf(fpInfo,"\n\t\tSites:\t\t\t%d",alignment->segsites);

	return 1;
}



void checkSNIPPositions (FILE* fp, alignment_struct * alignment, int index)
{
	int i,j=0;

	int t1 = alignment->positionsInd[0],
	    t2;
	
	fprintf(fp, "\n// Alignment %d\n\n",index);
	
	for(i=1;i<alignment->segsites;i++)
	{
		t2 = alignment->positionsInd[i];
	
		if (t2==t1)	
			fprintf(fp, " SNPs %d and %d correspond to the same alignment position: %d\n",j,i,alignment->positionsInd[i]);	
		
		t1 = t2;
                
		j=i;		
	}
	fprintf(fp,"\n");	
}

int isEndOfLine(char ent)
{
	if(ent == 10|| ent == 13 || ent == 32)
		return 1;

	return 0;
}



int isSpace(char ent)
{
	if(ent == 9 || ent == 32)
		return 1;

	return 0;
}



void ignoreAll(FILE * fp, char * ent)
{
	while(*ent < 33 && *ent != EOF)
		*ent = fgetc(fp);
}

int skipLine (FILE * fp)
{
	char tmp;

	while( (tmp = fgetc(fp) ) != EOF)
	  if(tmp == 10 || tmp == 13)
	    break;
	
	if(tmp == EOF)
	  return 1;
	
	return 0;	
}


int isValidVCFBase(char input)
{
	if(input=='A' || input=='C' || 
	   input=='G' || input=='T' || 
	   input=='N' || input=='a' || 
	   input=='c' || input=='g' || 
	   input=='t' || input=='n' || input == '.') return 1;

	return 0;
}

int scanStateVector(char * stateVector, char X)
{
	int i;

	for(i=0;i<MAX_STATES_VCF;i++)
		if(stateVector[i]==X)
			return 1;

	return 0;
}

int getStatesNum (char * stateVector)
{
  int i, states=0, l = strlen(stateVector);
  
  for(i=0;i<l;i++)
	{
		if(stateVector[i]!='X')
			states++;
		else
			return states;
	}

	return states;
}

int getStatesREF (char * string, char * stateVector, int line)
{
	int i, j, index=0, elen=0, slen=strlen(string);

	char CharToStore = 'X';

	for(j=0;j<MAX_STATES_VCF;j++)
		stateVector[j] = CharToStore;
	
	stateVector[MAX_STATES_VCF] = 0;

	for(i=0;i<slen;i++)
	{

		if(string[i]!=',')
		{
		
			if(string[i]=='<')
			{
				for(j=0;j<MAX_STATES_VCF;j++)
					stateVector[j] = 'X';

				return 0;
			}
			
			assert(isValidVCFBase(string[i])==1);

			elen++;

			CharToStore = string[i];
		}

		if(string[i]==',' || i == slen-1)
		{
			if(elen==1)
			{
				if(scanStateVector(stateVector, CharToStore)==0)
				{
					stateVector[index++] = CharToStore;
					elen = 0;	
				}
				else
				{
					fprintf(stderr, "\n\n ERROR: Nucleotide %c (field REF) in line %d appears twice.\n\n", CharToStore, line);
					assert(0); 
				}			
			}
			else
			{
				for(j=0;j<MAX_STATES_VCF;j++)
					stateVector[j] = 'X';

				return 0;
			}
		}
	}
	return 1;
}

int getStatesALT (char * string, char * stateVector, int line)
{
	int i, j, index=0, elen=0, slen=strlen(string);

	char CharToStore = 'X';

	for(j=0;j<MAX_STATES_VCF;j++)
		if(stateVector[j] == 'X')
		{
			index = j;
			break;
		}

	assert(index!=0);

	for(i=0;i<slen;i++)
	{

		if(string[i]!=',')
		{

			if(string[i]=='<')
			{
				for(j=0;j<MAX_STATES_VCF;j++)
					stateVector[j] = 'X';

				return 0;
			}

			if( isValidVCFBase( string[i] ) != 1 ) 
			  {
			    fprintf(stderr, "\n\nBase %c at line %d is not valid\n\n", string[i], line);
			    assert(isValidVCFBase(string[i])==1);
			  }

			elen++;

			CharToStore = string[i];
		}

		if(string[i]==',' || i == slen-1)
		{
			if(elen==1)
			{
				if(scanStateVector(stateVector, CharToStore)==0)
				{
					stateVector[index++] = CharToStore;
					elen = 0;	
				}
				else
				{
					fprintf(stderr, "\n\n ERROR: Nucleotide %c (field ALT) in line %d appears twice (in REF or ALT).\n\n", CharToStore, line);
					assert(0);
				}			
			}
			else
			{
				for(j=0;j<MAX_STATES_VCF;j++)
					stateVector[j] = 'X';

				return 0;
			}
		}
	}
	return 1;
}

float getValueAF (char * string, int line)
{
  int i, j, len = strlen(string);
  char AF_s [AFMAXLENGTH + 1];
  AF_s[0] = 0;
  
  float AF=-1.0;
  
  if(strcmp(".",string)==0 || len<4)
    return -1;
  
  for(i=0;i<len-3;i++)
    {	
      if(string[i]=='A' && string[i+1]=='F' && string[i+2]=='=')
	{
	  if(i==0)
	    break;
	  else
	    if(string[i-1]==';')
	      break;
	}
    }
  
  if(i==len-3)
    return -1.0;
  
  i = i+3;
  
  for(j=0;j<AFMAXLENGTH;j++)
    {
      if(string[i]==';' || i==len)
	break;
      
      if(string[i]=='.' || (string[i]>=48 && string[i]<=57) || string[i] == 'e' || string[i] == '-' || string[i] == '+')
	{
	  AF_s[j] = string[i];
	  i++;
	  AF_s[j+1] = 0;
	}
      else
	{
	  fprintf(stderr, "\n\n ERROR: Invalid character (%c) in AF of field INFO in line %d.\n\n", string[i], line );
	  assert(0);
	}
    }
  
  if(strlen(AF_s)==0)
    {
      fprintf(stderr, "\n\n ERROR: AF in field INFO in line %d has no value.\n\n", line);
      assert(strlen(AF_s)!=0);
    }
  
  AF = atof(AF_s);
  
  if(AF < 0.0 || AF > 1.0)
    {
      fprintf(stderr, "\n\n ERROR: AF (%f) in line %d should be >= 0.0 and <= 1.0.\n\n", AF, line);
      assert(AF>=0.0 && AF<=1.0);
    }
  
  return AF;
}

int checkVTisSNP (char * string)
{
	int i, len = strlen(string);

	if(strcmp(".",string)==0 || len<4)
		return -1;

	for(i=0;i<len-3;i++)
	{
		if(string[i]=='V' && string[i+1]=='T' && string[i+2]=='=')
		{
			if(i==0)
				break;
			else
				if(string[i-1]==';')
					break;
		}	
	}

	if(i==len-3)
		return -1;

	i = i+3;

	assert(i<=len-3);

	if(string[i]=='S' && string[i+1]=='N' && string[i+2]=='P')
		return 1;
	else
		return 0;
}

int getGTpos(char * string, int line)
{
	int i, len = strlen(string), GTposition = 0;

	for(i=0;i<len-1;i++)
	{
		if(string[i]==':')
			GTposition ++;

		if(string[i]=='G' && string[i+1]=='T')
		{
			if(i==0)
			{				
				if(i+2==len)
					return GTposition;

				else
					if(string[i+2]==':')
						return GTposition;
			}
			else
			{
				if(string[i-1]==':')
				{
					if(i+2==len)
						return GTposition;
					else
						if(string[i+2]==':')
							return GTposition;
				}
			}
		}			
	}

	GTposition = -1;

	return GTposition;
}

int getGTfield (char * string, int GTpos)
{
	int i=0, pos=GTpos, len = strlen(string),j=0, counter=0;

	assert(pos!=-1);	

	for(i=0;i<len;i++)
	{
		if(string[i]!=':')
		{	
			if(pos==0)
			{
				string[j++]=string[i];

				if(string[i]=='|' || string[i]=='/')
					counter++;
			}
		}
		else
		{
			pos--;
			
			if(pos==-1)
			{
				string[j]=0;
				break;
			}
		}
	}

	assert(strlen(string) > 0);	

	return counter+1;
}

void dataShuffleKnuth(char * data, int startIndex, int endIndex)
{
	if(startIndex == endIndex)
		return;

	int i, index;
	char tmp;

	for (i = endIndex; i > startIndex; i--)
	{
		index = startIndex + (rand() % (i - startIndex + 1));

		tmp = data[index];
		data[index] = data[i];
		data[i] = tmp;
	}
}

void getGTdata (char * string, char * stateVector, int statesTotal, char * sampleData)
{
	int i, j=0, index=0, start=0, end=0, len = strlen(string);

	for(i=0;i<len;i++)
	{	
		if(string[i]>=48 && string[i]<=57)
		{
			index = string[i]-48;

			assert(index<statesTotal);
			
			sampleData[j++] = stateVector[index];
		}
		else
		{
			if(string[i]=='.')
			{
				sampleData[j++] = 'N';
			}
			
			if(string[i]=='/')
			{
				end++;
			}

			if(string[i]=='|')
			{
				dataShuffleKnuth(sampleData, start, end);
				start = j;
				end = j;
			}			
		}
	}

	dataShuffleKnuth(sampleData, start, end);
}

void processSampleVCF (char * string, int GTpos, char ** sampleData, int * sampleDataMemSize, char * stateVector, int statesTotal)
{
	
	int i, dataSize = getGTfield (string,GTpos);

	if(dataSize+1>*sampleDataMemSize)
	{
		*sampleDataMemSize = dataSize+1;

		*sampleData = realloc((*sampleData), (*sampleDataMemSize) * sizeof(char));		
	}

	for(i=0;i<(*sampleDataMemSize);i++)
		(*sampleData)[i]=0;

	getGTdata (string, stateVector, statesTotal, *sampleData);

	(*sampleData)[dataSize]=0;
}


int readLine_VCF (FILE * fp, char ** string, int lineIndex, alignment_struct * alignment, int * DIM, char ** sampleData, int * sampleDataMemSize, int * SNP_SZ, int * tmpSNPtableLineIndex, int outgroup, int * errors)
{

  int eol=0, eof=0, maxLength=STRINGLENGTH;
  char tmp;
  int elementIndex = -1, lineSkipped=0, isEOF = 1, i,j=0, y, s=0, tip = 0;

	char stateVector[MAX_STATES_VCF+1]; 

	
	
	int position=-1, statesREF=0, statesALT=0, VTisSNP=0, GTpos=-1, sampleIndex;
	//float AF=0.5;


	if(lineIndex==0)
		elementIndex++;	

	

	while(getNextString (fp, string, &eol, &eof, &maxLength)==1)
	{
		
		elementIndex++;

		

		switch(elementIndex)
		{

			

			case 0: // CHROM
				if(strcmp(VCF_alignment_name,*string)!=0)
				{
					strncpy(VCF_alignment_name, *string, MAX_CHROM_NAME_VCF);

					assert(strlen(VCF_alignment_name)!=0);

					return 0;
				}
				break;

			case 1: // POS
				
				position = atoi(*string);

				break;


			case 2: // ID
	
				break;


			case 3: // REF
				
				getStatesREF (*string, stateVector, lineIndex+1+VCF_header_lines);
				
				statesREF = getStatesNum (stateVector);	

				if(statesREF==0)
				{
					isEOF = skipLine(fp);

					if (isEOF == 1)
						return -1;

					lineSkipped=1;

					alignment->discSites += lineSkipped;

					break;
				}

				alignment->outgroupSequence[alignment->segsites] = stateVector[0];
	
				break;


			case 4: // ALT
			  
			  getStatesALT (*string, stateVector, lineIndex+1+VCF_header_lines);			
			  statesALT = getStatesNum (stateVector);
			  
			  
				
			  if(statesALT != 0 && (stateVector[ statesALT - 1] == '.') )
			    for( s = statesALT - 1; s < statesALT-1 + statesREF; ++s)
			      stateVector[s] = stateVector[s-statesREF];
			  
				
				if(statesALT==0)
				{
					isEOF = skipLine(fp);

					if (isEOF == 1)
						return -1;

					lineSkipped=1;

					alignment->discSites += lineSkipped;

					break;
				}

				break;


			case 5: // QUAL
				
				break;


			case 6: // FILTER
			  
			  if((strcmp("PASS",*string)!=0) && (strcmp(".", *string)!=0) )
				{
					isEOF = skipLine(fp);

					if (isEOF == 1)
						return -1;

					lineSkipped=1;

					alignment->discSites += lineSkipped;

					break;
				}

				break;


			case 7: // INFO
				
//				AF = getValueAF(*string, lineIndex+1+VCF_header_lines); // returns -1.0 if AF not found

/*				if(AF==0.0 || AF==1.0)
				{
					isEOF = skipLine(fp);

					if (isEOF == 1)
						return -1;

					lineSkipped=1;

					break;
				}
*/
				VTisSNP = checkVTisSNP(*string); // returns -1.0 if VT not found

				if(VTisSNP==0)
				{
					isEOF = skipLine(fp);

					if (isEOF == 1)
						return -1;

					lineSkipped=1;

					alignment->discSites += lineSkipped;

					break;
				}

				break;


			case 8: // FORMAT

				GTpos = getGTpos(*string, lineIndex+1+VCF_header_lines); // returns -1 if GT not found

				if(GTpos==-1)		
				{
					isEOF = skipLine(fp);

					if (isEOF == 1)
						return -1;

					lineSkipped=1;
	
					alignment->discSites += lineSkipped;

					break;
				}
			
				break;

			default:

				break;

		}

		if(elementIndex == VCF_HLENGTH-1 || eol==1 || eof==1 || lineSkipped==1)		
			break;
	}

	if (eof== 1 || position == -1)
		return -1;

	if(lineSkipped!=1)
	{

		(*tmpSNPtableLineIndex)++;

		alignment->segsites++;
	
		if(alignment->segsites>=2)
			assert(position>=alignment->positionsInd[alignment->segsites-2]);

		if(alignment->segsites>=*DIM)
		{
			*DIM = *DIM + 10;

			alignment->positionsInd = (int *) realloc(alignment->positionsInd,*DIM*sizeof(int));
			assert(alignment->positionsInd != NULL);			

			alignment->outgroupSequence = realloc(alignment->outgroupSequence, *DIM * sizeof(char));
			assert(alignment->outgroupSequence != NULL);

			alignment->stateA = realloc(alignment->stateA, *DIM * sizeof(al_t));
			assert(alignment->stateA != NULL);

			alignment->stateC = realloc(alignment->stateC, *DIM * sizeof(al_t));
			assert(alignment->stateC != NULL);

			alignment->stateG = realloc(alignment->stateG, *DIM * sizeof(al_t));
			assert(alignment->stateG != NULL);

			alignment->stateT = realloc(alignment->stateT, *DIM * sizeof(al_t));
			assert(alignment->stateT != NULL);

			for(y=*DIM-10;y<*DIM;y++)
			{
				alignment->stateA[y]=0;
				alignment->stateC[y]=0;
				alignment->stateG[y]=0;
				alignment->stateT[y]=0;
			}
				
		}

		alignment->positionsInd[alignment->segsites-1] = position;

		sampleIndex = 0;

	//	for(i=0;i<alignment->VCFsamples;i++)
	//		printf("loc %d : %d\n",i, alignment->VCFsample_valid[i]);

//		assert(0);
		while( ( tip = getNextString (fp, string, &eol, &eof, &maxLength) ) == 1 )
		{
		  
			//printf("sample %d: %s\n",sampleIndex,*string);
			assert(sampleIndex<alignment->VCFsamples);

			
			if(alignment->VCFsample_valid[sampleIndex]==1)
			{
			//	printf("procssing sample: %d\n",sampleIndex);
				processSampleVCF(*string, GTpos, sampleData, sampleDataMemSize, stateVector, statesALT);

				if(*tmpSNPtableLineIndex==0 && sampleIndex+1!=outgroup)
					*SNP_SZ += strlen(*sampleData);
		
				for(i=0;i<strlen(*sampleData);i++)
				{
					if(sampleIndex+1==outgroup)			
						alignment->outgroupSequence[alignment->segsites-1] = (*sampleData)[0];		
					else
					{
						if(j>=*SNP_SZ)
						{
						    fprintf(stderr, "\n\n ERROR: There are more than %d nucleotides in line %d. Expected %d (according to the first SNP in line %d).\n\n", j, lineIndex+1+VCF_header_lines, *SNP_SZ, VCF_first_SNP );

						    (*errors)++;
						    assert(j<*SNP_SZ);
						}
						updateStateCounters (alignment, (*sampleData)[i], alignment->segsites-1);
						++j;
					}				
				}

				
			}
			sampleIndex++;
			if(eol==1 || eof==1)
				break;
		}

		//printf("processed samples: %d\n",sampleIndex);
		
		if(eol == 1)
		  {
		    tmp = fgetc(fp);
		    ungetc(tmp, fp);
		    ignoreNewLineChars(fp, &tmp);
		  }

		
		if(j!=*SNP_SZ)
		{
			fprintf(stderr, "\n\n ERROR: There are %d nucleotides in line %d. Expected %d (according to the first SNP in line %d).\n\n", j, lineIndex+1+VCF_header_lines, *SNP_SZ, VCF_first_SNP );
			(*errors)++;
			//return 0;
			assert(j==*SNP_SZ);
		}

		if (eof== 1)
			return -1;
	}

	return 1;
}

int scanListfor(char * chromName)
{
	int i;
	assert(chromList_SZ>=1);

	for(i=0;i<chromList_SZ;i++)
	{
		
		if(!strcmp((*chromList)[i],chromName))
		{
			return 1;
		}
	}
	return 0;
}

void skipAlignment(FILE *fp, char * VCF_alignment_name)
{
	int inAlignment = 1, lineIndex=-1;
	int eol=0, eof=0, maxLength=STRINGLENGTH;
	char ** string = (char **) malloc (sizeof(char*));
	(*string) = (char *) malloc(sizeof(char)*STRINGLENGTH);

	while(inAlignment==1)
	{
		lineIndex++;

		skipLine( fp );
		getNextString (fp, string, &eol, &eof, &maxLength);
		if(strcmp(VCF_alignment_name,*string)!=0)
		{
			strncpy(VCF_alignment_name, *string, MAX_CHROM_NAME_VCF);

			assert(strlen(VCF_alignment_name)!=0);

			inAlignment=0;
		}
	}

	VCF_header_lines = lineIndex + VCF_header_lines;

	free(*string);
	free(string);
	
}
int readAlignmentVCF(FILE *fp, alignment_struct * alignment, FILE *fpInfo, FILE *fpSFo, int minsnps_threshold_user, int alignmentIndex)
{

	if(chromList_SZ!=-1)
	{
		if(scanListfor(VCF_alignment_name)==0)
		{
			skipAlignment(fp, VCF_alignment_name);
			return 2;
		}
	}

	int errors = 0;

	alignment->maxx = -1;
	alignment->maxn = -1;
	alignment->minx = MAXINT;
	alignment->minn = MAXINT;

	char ** string = (char **) malloc (sizeof(char*));
	(*string) = (char *) malloc(sizeof(char)*STRINGLENGTH);

	int DIM=2;
	alignment->positionsInd = (int *) malloc(DIM*sizeof(int));

	alignment->outgroupSequence = malloc(DIM * sizeof(char));
	assert(alignment->outgroupSequence!=NULL);

	alignment->stateA = calloc (DIM, sizeof(al_t));
	assert(alignment->stateA!=NULL);

	alignment->stateC = calloc (DIM, sizeof(al_t));
	assert(alignment->stateC!=NULL);

	alignment->stateG = calloc (DIM, sizeof(al_t));
	assert(alignment->stateG!=NULL);

	alignment->stateT = calloc (DIM, sizeof(al_t));
	assert(alignment->stateT!=NULL);

	alignment->segsites = 0;

	fprintf(stdout," Alignment %d\n",alignmentIndex);			
	fprintf(fpInfo," Alignment %d\n",alignmentIndex);

	fprintf(stdout,"\n\t\tChromosome:\t\t%s", VCF_alignment_name);
	fprintf(fpInfo,"\n\t\tChromosome:\t\t%s", VCF_alignment_name);

	strcpy(VCF_alignment_name_prev_to_print, VCF_alignment_name);	

	int inAlignment = 1, lineIndex=-1, i;

	int sampleDataMemSize = 2;
	char ** sampleData = (char **) malloc(sizeof(char*));
	(*sampleData) = (char *) malloc(sizeof(char)*sampleDataMemSize);

	int SNP_SZ = 0;
	//int SNP_NUM = 10;

	int prevtmpSNPtableLineIndex,tmpSNPtableLineIndex=-1;

	int outgroup = alignment->outgroupIndex;

	while(inAlignment==1)
	{
		lineIndex++;

		prevtmpSNPtableLineIndex = tmpSNPtableLineIndex;

		inAlignment = readLine_VCF (fp, string, lineIndex, alignment, &DIM, sampleData, &sampleDataMemSize, &SNP_SZ, &tmpSNPtableLineIndex, outgroup, &errors);		

		if(inAlignment==1)
		{
			if(prevtmpSNPtableLineIndex!=tmpSNPtableLineIndex)
			{
				if(prevtmpSNPtableLineIndex == -1)
					VCF_first_SNP = lineIndex + VCF_header_lines + 1;
			}
		}
	}

	nxtVCFalignment = inAlignment;

	if(nxtVCFalignment==0)
		VCF_header_lines = lineIndex + VCF_header_lines;

	alignment->sequences = SNP_SZ;

	alignment->length = alignment->positionsInd[alignment->segsites-1];

	assert(alignment->segsites<=alignment->length);

	free(*sampleData);
	free(sampleData);
	free(*string);
	free(string);

	//for(i=0;i<DIM;i++)
	//	alignment->outgroupSequence[i]=0;

	createXNFvectors (alignment);

	for(i=0;i<alignment->segsites;i++)
	{
		if(alignment->n[i]<alignment->minn)
			alignment->minn = alignment->n[i];

		if(alignment->x[i]<alignment->minx)
			alignment->minx = alignment->x[i];

		if(alignment->n[i] > alignment->sequences)
			alignment->sequences = alignment->n[i];

		alignment->maxn = alignment->sequences;		

		alignment->maxx = alignment->sequences;
	}

	for(i=0;i<alignment->segsites;i++)
		if((alignment->x[i]==alignment->n[i]) || (alignment->folded[i]==1 && alignment->x[i]==0))
		{
			alignment->maxx++;
			break;
		}

	if(alignment->outgroupSequence!=NULL)
	{
		free(alignment->outgroupSequence);
		alignment->outgroupSequence=NULL;
	}
	
	fprintf(stdout,"\n\t\tSequences:\t\t%d",alignment->sequences);	
	fprintf(fpInfo,"\n\t\tSequences:\t\t%d",alignment->sequences);

	fprintf(stdout,"\n\t\tSites:\t\t\t%d",alignment->segsites);	
	fprintf(fpInfo,"\n\t\tSites:\t\t\t%d",alignment->segsites);

	if(errors!=0)
		return 0;

	return 1;
}

void updateSF(FILE * fpSFo)
{
	int i;

	if(fpSFo!=NULL)
		for(i=0;i<alignment->segsites;i++)
			if(alignment->positionsInd[i]!=-1)
				fprintf(fpSFo,"%d\t%d\t%d\t%d\n",alignment->positionsInd[i],alignment->x[i],alignment->n[i],alignment->folded[i]);
}


void removeMonomorphicSites (int strictPolymorphic, int monomorphic, FILE * fp)
{
	int i, counter=0;

	for(i=0;i<alignment->segsites;i++)
		if(alignment->n[i]==0)
		{
			alignment->positionsInd[i]=-1;
			counter++;
		}

	if( monomorphic == 1 && strictPolymorphic == 1)
	  {
	    fprintf(stderr, "Cannot be both monomorphic == 1 && strictPolymorphic == 1\n\n");
	    assert( monomorphic == 0 || strictPolymorphic == 0);
	  }

	if (monomorphic==0)
	{
	  
	  for(i=0;i<alignment->segsites;i++)
	    {
	      if(alignment->n[i] == 0)
		continue;
	      
	      if( (alignment->n[i]==alignment->sequences && (alignment->x[i]==0 || alignment->x[i]==alignment->n[i]) ) ||
		  (strictPolymorphic == 1 && alignment->n[i] < alignment->sequences &&  alignment->x[i] == 0  ) ||
		  (strictPolymorphic == 1 && alignment->n[i] < alignment->sequences &&  alignment->x[i] == alignment->n[i]  ))
		{
		  alignment->positionsInd[i]=-1;			
		  
		  counter++;
		}
	    }
	
	  if( alignment->segsites - counter < 0)
	    {
	      fprintf(stderr, "counter: %d, alignment->segsites: %d\n", counter, alignment->segsites);
	      assert(alignment->segsites-counter >= 0);
	    }
	  
		if(counter != 0)
		{
			
			alignment->p_t = malloc(sizeof(int)*(alignment->segsites-counter));
			alignment->f_t = malloc(sizeof(unsigned char)*(alignment->segsites-counter));
			alignment->x_t = malloc(sizeof(al_t)*(alignment->segsites-counter));
			alignment->n_t = malloc(sizeof(al_t)*(alignment->segsites-counter));

			int j=0;

			for(i=0;i<alignment->segsites;i++)
				if(alignment->positionsInd[i]!=-1)
				{
					alignment->p_t[j] = alignment->positionsInd[i];
					alignment->f_t[j] = alignment->folded[i];
					alignment->x_t[j] = alignment->x[i];
					alignment->n_t[j] = alignment->n[i];
					j++;
				}

			alignment->segsites = alignment->segsites-counter;

			free(alignment->positionsInd);
			free(alignment->folded);
			free(alignment->x);
			free(alignment->n);

			alignment->positionsInd = alignment->p_t;
			alignment->folded = alignment->f_t;
			alignment->x = alignment->x_t;
			alignment->n = alignment->n_t;
		}

	


		alignment->maxx = -1;
		alignment->maxn = -1;
		alignment->minx = MAXINT;
		alignment->minn = MAXINT;


		for(i=0;i<alignment->segsites;i++)
		{
			if(alignment->positionsInd[i]!=-1)
			{
				if(alignment->n[i]<alignment->minn)
					alignment->minn = alignment->n[i];

				if(alignment->x[i]<alignment->minx)
					alignment->minx = alignment->x[i];

				if(alignment->n[i] > alignment->sequences)
					alignment->sequences = alignment->n[i];

				alignment->maxn = alignment->sequences;		

				alignment->maxx = alignment->sequences;
			}
		}

		for(i=0;i<alignment->segsites;i++)
			if(alignment->positionsInd[i]!=-1)
			{
				if((alignment->x[i]==alignment->n[i]) || (alignment->folded[i]==1 && alignment->x[i]==0))
				{
					alignment->maxx++;
					break;
				}
			}
	}	

	fprintf(stdout,"\n\t\tDiscarded sites:\t%d",counter + alignment->discSites);	
	fprintf(fp,"\n\t\tDiscarded sites:\t%d",counter + alignment->discSites);
}

int readAlignmentMACS(FILE *fp, alignment_struct *alignment, FILE *fpInfo, FILE *fpSFo, int minsnps_threshold_user, int alignmentIndex)
{

	fprintf(stdout," Alignment %d\n",alignmentIndex);			
	fprintf(fpInfo," Alignment %d\n",alignmentIndex);

	int  i, DIM = 2, prevDIM;
	char ent;
	int nsnp = 0;
	int sitevar;
	char siteflag[100];

	alignment->maxx = -1;
	alignment->maxn = -1;
	
	alignment->n = calloc(DIM, sizeof(al_t));
	assert(alignment->n != NULL);

	alignment->x = calloc(DIM, sizeof(al_t));
	assert(alignment->x != NULL);

	alignment->folded = calloc(DIM, sizeof(unsigned char));
	assert(alignment->folded != NULL);

	alignment->positions = calloc(DIM, sizeof(float));
	assert(alignment->positions != NULL);

	
	// read first line 
	int temp = fscanf(fp, "%f", &alignment->positions[0]);
	assert(temp==1);

	alignment->sequences = 0;


	// read first SNP
	while((ent = fgetc(fp)))
	  {
		ignoreLineSpaces(fp, &ent);

		if(isEndOfLine(ent))
		  break;
				
		if(isBinary(ent))
		  {
		    ++alignment->n[0];
		    alignment->x[0] += (al_t) (ent - '0');
		  }

		alignment->folded[0] = alignment->userSetFolded;
		alignment->sequences++;
	}



	nsnp = 1;

	// now rest of the lines
	while(1)
	  {
	    if(nsnp >= DIM)
	      {
		prevDIM = DIM;
		
		DIM = DIM + 10;
		
		alignment->n = realloc(alignment->n, DIM*sizeof(al_t));
		assert(alignment->n != NULL);
		
		alignment->x = realloc(alignment->x, DIM*sizeof(al_t));
		assert(alignment->x != NULL);
		
		alignment->folded = realloc(alignment->folded, DIM*sizeof(al_t));
		assert(alignment->folded != NULL);
		
		alignment->positions = realloc(alignment->positions, DIM*sizeof(float));
		assert(alignment->positions != NULL);
		
		
		for(i=prevDIM; i<DIM; ++i)
		  {
		    alignment->x[i] = 0;
		    alignment->n[i] = 0;
		    alignment->folded[i] = alignment->userSetFolded;
		  }
	      }
	    
	    temp = fscanf(fp, "%s",siteflag);
		assert(temp==1);
	    
	    if(!strcmp(siteflag,"SITE:"))
	      {
		temp = fscanf(fp, "%d %f",&sitevar, &alignment->positions[nsnp]);
		assert(temp==2);

		ent = fgetc(fp);

		ignoreLineSpaces(fp, &ent);

				
		while(!isEndOfLine(ent) && !isSpace(ent))
		  {
		    if(isBinary(ent))
		      {
			++alignment->n[nsnp];
			alignment->x[nsnp] += (al_t) (ent - '0');
		      }
		    ent = fgetc(fp);		    
		  }
		
			
		alignment->folded[nsnp] = alignment->userSetFolded;
		
		nsnp++;
	      }
	    else
	      break;
	  }
	
	
	alignment->segsites = nsnp;	
	assert(alignment->segsites<=alignment->length);


	for(i=0;i<alignment->segsites;i++)
	{
		if(alignment->n[i]<alignment->minn)
			alignment->minn = alignment->n[i];

		if(alignment->x[i]<alignment->minx)
			alignment->minx = alignment->x[i];

		if(alignment->n[i] > alignment->sequences)
			alignment->sequences = alignment->n[i];

		alignment->maxn = alignment->sequences;		

		alignment->maxx = alignment->sequences;
	}

	for(i=0;i<alignment->segsites;i++)
		if((alignment->x[i]==alignment->n[i]) || (alignment->folded[i]==1 && alignment->x[i]==0))
		{
			alignment->maxx++;
			break;
		}


	alignment->positionsInd = malloc(alignment->segsites*sizeof(int));

	for(i=0;i<alignment->segsites;i++)
	  alignment->positionsInd[i] = (int)(alignment->positions[i] * (float)alignment->length);

	free(alignment->positions);


	fprintf(stdout,"\n\t\tSequences:\t\t%d",alignment->sequences);	
	fprintf(fpInfo,"\n\t\tSequences:\t\t%d",alignment->sequences);

	fprintf(stdout,"\n\t\tSites:\t\t\t%d",alignment->segsites);	
	fprintf(fpInfo,"\n\t\tSites:\t\t\t%d",alignment->segsites);

	return 1;

}

int readAlignmentMS(FILE *fp, alignment_struct *alignment, FILE * fpInfo, FILE *fpSFo, int minsnps_threshold_user, int alignmentIndex)
{

	fprintf(stdout," Alignment %d\n",alignmentIndex);			
	fprintf(fpInfo," Alignment %d\n",alignmentIndex);

  char ent, stringtemp[100];
  
	
  /* get rid of the first line information */
  while( (ent = fgetc(fp) ) != '\n');

  int i, temp = fscanf(fp,"%s %d", stringtemp, &alignment->segsites); 

  if(strcmp(stringtemp, "segsites:") != 0 )
    {
      fprintf(stderr, "ERROR! \"segsites:\" are expected but %s has been found\n", stringtemp);
      assert(0);
    }

  
  if( skipLine( fp ) == 1)
    {
      fprintf(stderr, "Unexpected end of file after segsites....\n\n");
      assert(0);
    }

	  
  alignment->positions = malloc(sizeof(float)*alignment->segsites); 
  
  alignment->positionsInd = malloc(sizeof(int)*alignment->segsites);

  if( alignment->segsites > 0)
    {
      temp = fscanf(fp, "%s", stringtemp);
      if(strcmp( stringtemp, "positions:") != 0)
	{
	  fprintf(stderr, "ERROR: \"positions:\" is expected\n\n");
	  assert(0);
	}
    }

  
  alignment->maxx = -1;
  alignment->maxn = -1;
  alignment->minx = MAXINT;
  alignment->minn = MAXINT;
  
  alignment->sequences = 0;
  
  for(i=0;i<alignment->segsites;i++)
    {
      temp = fscanf(fp,"%f",&alignment->positions[i]);
      assert(temp==1);		
    }

  for(i=0;i<alignment->segsites;i++)
    alignment->positionsInd[i] = (int)(alignment->positions[i] * alignment->length);
  
  free(alignment->positions);

  
  
  alignment->n = calloc(alignment->segsites, sizeof(al_t)); 
  assert(alignment->n != NULL);
  
  alignment->x = calloc(alignment->segsites, sizeof(al_t));
  assert(alignment->x != NULL);
  
  alignment->folded = calloc(alignment->segsites, sizeof(unsigned char)); 
  assert(alignment->folded != NULL);
  
  alignment->sequences=0;
  
  if(alignment -> segsites == 0 )
    skipToNextAlignmentDelimiter(fp, '/');
  


  while(alignment->segsites > 0)
    {
      ent = fgetc(fp);
      
      while(ent=='\n' || ent==' ')
	ent = fgetc(fp);		
      
      if(ent=='/' || ent==EOF)
	break;		
      
      i=0;
      
      while(ent!='\n' && ent!=' ')
	{
	  if(isBinary(ent))
	    {
	      alignment->x[i] += (al_t)(ent - '0');
	      ++alignment->n[i];
	    }
	  
	  alignment->folded[i] = alignment->userSetFolded; 
	  
	  ++i;
	  
	  ent = fgetc(fp);
	}
    }
  
  for(i=0;i<alignment->segsites;i++)
    {
      if(alignment->n[i]<alignment->minn)
	alignment->minn = alignment->n[i];
      
      if(alignment->x[i]<alignment->minx)
	alignment->minx = alignment->x[i];
      
      if(alignment->n[i] > alignment->sequences)
	alignment->sequences = alignment->n[i];
      
      alignment->maxn = alignment->sequences;		
      
      alignment->maxx = alignment->sequences;
    }
  
  for(i=0;i<alignment->segsites;i++)
    if((alignment->x[i]==alignment->n[i]) || (alignment->folded[i]==1 && alignment->x[i]==0))
      {
	alignment->maxx++;
	break;
      }
  
  assert(alignment->segsites<=alignment->length);
  
  fprintf(stdout,"\n\t\tSequences:\t\t%d",alignment->sequences);	
  fprintf(fpInfo,"\n\t\tSequences:\t\t%d",alignment->sequences);
  
  fprintf(stdout,"\n\t\tSites:\t\t\t%d",alignment->segsites);	
  fprintf(fpInfo,"\n\t\tSites:\t\t\t%d",alignment->segsites);

  if(alignment -> segsites < MINSNPS_THRESHOLD || alignment -> segsites < minsnps_threshold_user)
    return 0;

  return 1;
}


int readAlignmentSF(FILE *fp, FILE * fpInfo, FILE *fpSFo, int minsnps_threshold_user, int alignmentIndex)
{

	fprintf(stdout," Alignment %d\n",alignmentIndex);			
	fprintf(fpInfo," Alignment %d\n",alignmentIndex);

  int i = 0, j = 0, DIM=1, x = 0, n = 0, folded = 0, cnt = 0;
  char *fgetsValue = NULL; 
  double position = 0.;
  char tmpString[100];
  char ent;
  char line[LEN];

  alignment->x = calloc(DIM, sizeof(al_t));
  assert(alignment->x != NULL);

  alignment->n = calloc(DIM, sizeof(al_t));
  assert(alignment->n != NULL);

  alignment->folded = calloc(DIM, sizeof(unsigned char));
  assert(alignment->folded != NULL);

  alignment->positionsInd = calloc(DIM, sizeof(int));
  assert(alignment->positionsInd != NULL);

  alignment->segsites = -1;

  alignment->maxx = -1;
  alignment->maxn = -1;
  alignment->minx = MAXINT;
  alignment->minn = MAXINT;

  alignment->sequences = -1;

  while(1)
    {

      cnt++;

      char charTMP = fgetc(fp);
      
      ungetc(charTMP, fp);

      //fprintf(stderr, "CHAR: %c--%d\n", charTMP, charTMP);
      
      fgetsValue = fgets( line, LEN, fp);

      //fprintf(stderr, "%d -- %s -- %d-- %p\n", cnt, line, strlen(line), fgetsValue);
 	    
      if( fgetsValue != NULL) 
	{
	  for(j = 0; j < LEN; ++j)
	    {
	      if(line[0] == 9 || line[0] == 32)
		continue;
	      else
		break;
	    }
	      
	  if(j < LEN && line[j] == '/')
	    {
	      ungetc('/', fp);
	      break;
	    }

	  if(j < LEN && (line[j] == 10 || line[j] == 13 || line[j] == 0))
	    {
	      continue;
	    }
	 
	  if(sscanf(line, "%s\t%d\t%d\t%d", tmpString, &x, &n, &folded)!=4)
	    {
	      
	      ent = fgetc(fp);
	      fprintf(stderr,"\n ERROR: Invalid character (%c - ASCII: %d) in the input file. (Alignment separator: //), last: %c,%d\n\n", ent, ent, line[j], line[j]);
	      
	      assert(0);
	    }
	  
	  position = atof(tmpString);
	  
	  if(i >= DIM)
	    {
	      DIM = DIM + 10;
	      
	      alignment->x = realloc(alignment->x, DIM * sizeof(al_t));
	      assert(alignment->x != NULL);
	      
	      alignment->n = realloc(alignment->n, DIM * sizeof(al_t));
	      assert(alignment->n != NULL);

	      alignment->folded = realloc(alignment->folded, DIM * sizeof(unsigned char));
	      assert(alignment->folded != NULL);

	      alignment->positionsInd = realloc(alignment->positionsInd, DIM * sizeof(int));
	      assert(alignment->positionsInd != NULL);
	    }

	  if(x>n)
	    {
	      fprintf(stderr,"\n ERROR: Invalid input values for x (%d) and n (%d) at position %.0f. x must be less than or equal to n.\n\n", x, n, position);
	      exit(0);
	    }
	  
	  alignment->x[i] = x;
	  alignment->n[i] = n;
	  alignment->positionsInd[i] = (int)position;
	  
	  if(alignment->userSetFolded==1)
	    alignment->folded[i] = 1;
	  else
	    alignment->folded[i] = folded;
		
	  alignment->segsites = i+1;
	  
	  if(i>0)
	    {	
	      if(alignment->positionsInd[i]<alignment->positionsInd[i-1])
		printf("%d %d\n", alignment->positionsInd[i],alignment->positionsInd[i-1]);

	      assert(alignment->positionsInd[i]>=alignment->positionsInd[i-1]);
	    }

	  if(n<alignment->minn)
	    alignment->minn = n;

	  if(x<alignment->minx)
	    alignment->minx = x;

	  if(n > alignment->sequences)
	    alignment->sequences = n;
	  
	  alignment->maxn = alignment->sequences;		

	  alignment->maxx = alignment->sequences;

	  ++i;
	}
      
      else
	{
	  ent = EOF;
	  break;
	}
      
    }
  
  
  for(i=0;i<alignment->segsites;i++)
    if((alignment->x[i]==alignment->n[i]) || (alignment->folded[i]==1 && alignment->x[i]==0))
      {
	alignment->maxx++;
	break;
      }
  
  fprintf(stdout,"\n\t\tSequences:\t\t%d",alignment->sequences);	
  fprintf(fpInfo,"\n\t\tSequences:\t\t%d",alignment->sequences);
  
  fprintf(stdout,"\n\t\tSites:\t\t\t%d",alignment->segsites);	
  fprintf(fpInfo,"\n\t\tSites:\t\t\t%d",alignment->segsites);

  if(alignment -> sequences < 2 || alignment -> segsites < MINSNPS_THRESHOLD || alignment -> segsites < minsnps_threshold_user)
    return 0;
  
  return 1;
  
}

void generate_VCF_chrom_list(FILE * fpin, FILE * fpVCFchroms, char * chromVCFfileName, FILE * fpInfo)
{
	char ** string = (char **) malloc (sizeof(char*));
	(*string) = (char *) malloc(sizeof(char)*STRINGLENGTH);
	
	int eol=0, eof=0, maxLength=STRINGLENGTH;

	char prev_VCF_alignment_name [MAX_CHROM_NAME_VCF];
	char cur_VCF_alignment_name [MAX_CHROM_NAME_VCF];

	strcpy(prev_VCF_alignment_name,VCF_alignment_name);
	strcpy(cur_VCF_alignment_name,VCF_alignment_name);

	fprintf(fpVCFchroms, "%s\n", cur_VCF_alignment_name);

	while(skipLine (fpin)!=1) // while not EOF
	{
		getNextString (fpin, string, &eol, &eof, &maxLength);

		strcpy(cur_VCF_alignment_name, *string);

		if(strcmp(prev_VCF_alignment_name,cur_VCF_alignment_name))
		{
			fprintf(fpVCFchroms, "%s\n", cur_VCF_alignment_name);

			strcpy(prev_VCF_alignment_name, cur_VCF_alignment_name);
		}

		
	}

	fprintf(stdout, "\n A list of VCF chromosomes has been stored in file %s\n\n",chromVCFfileName);
	fprintf(fpInfo, "\n A list of VCF chromosomes has been stored in file %s\n\n",chromVCFfileName);

	free(*string);
	free(string);
}

void extract_chromList(FILE * fpVCFchroms, FILE * fpInfo)
{
	char ** string = (char **) malloc (sizeof(char*));
	(*string) = (char *) malloc(sizeof(char)*STRINGLENGTH);

	int eol=0, eof=0, maxLength=STRINGLENGTH;
	
	int chromIndex=-1;
	
	chromList_SZ = 0;

	while(getNextString_all_lines (fpVCFchroms, string, &eol, &eof, &maxLength)==1)
	{
		chromIndex++;
		chromList_SZ++;

		if(chromIndex==0)
			*chromList = malloc(sizeof(char*)*chromList_SZ);
		else
			*chromList = realloc (*chromList, (chromList_SZ)*sizeof(char*));

		(*chromList)[chromIndex] = malloc(sizeof(char)*MAX_CHROM_NAME_VCF);
		strcpy((*chromList)[chromIndex], *string);
	}

	free(*string);
	free(string);

	fprintf(stdout, " Number of chromosomes to be analyzed:\t%d\n",chromList_SZ);
	fprintf(fpInfo, " Number of chromosomes to be analyzed:\t%d\n",chromList_SZ);
}

int readAlignment(FILE *fp, int format, FILE * fpInfo, FILE *fpSFo, int minsnps_threshold_user, int alignmentIndex)
{
  if(format == SF_FORMAT)
    return readAlignmentSF(fp, fpInfo, fpSFo, minsnps_threshold_user, alignmentIndex);

  if(format == MS_FORMAT)
    return readAlignmentMS(fp, alignment, fpInfo, fpSFo, minsnps_threshold_user, alignmentIndex);

  if(format == MACS_FORMAT)
    return readAlignmentMACS(fp, alignment, fpInfo, fpSFo, minsnps_threshold_user, alignmentIndex);

  if(format==FASTA_FORMAT)
    return readAlignmentFASTA(fp, alignment, fpInfo, fpSFo, minsnps_threshold_user, alignmentIndex);

  if(format==VCF_FORMAT)
    return readAlignmentVCF(fp, alignment, fpInfo, fpSFo, minsnps_threshold_user, alignmentIndex);

  return 0;
} 
