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


#ifndef _SweeD_H
#define _SweeD_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <limits.h>
#include <ctype.h>

#include "SweeD_BFGS.h"

#ifdef _USE_PTHREADS
#include <pthread.h>
#include <unistd.h>

#define EXIT 127
#define BUSYWAIT 0
#define LIKSFS 1
#define COMPRVLUT 2
#define COMPALPHA 3

#endif

#ifdef _DO_CHECKPOINTS
#include "dmtcpaware.h"
#endif

#ifdef _ANALYTICAL_SFS
#include <mpfr.h>
#include <float.h>
#define ACC 10000
#define APPROX 0
#define TOLTIME 1e-7
#define RESCALE 1
#endif

#define ENTRIES_PER_INT 32
#define DENOMINATOR_OFFSET 0.00001
#define REGISTER_WIDTH 32
#define MAXSIZE 10000
#define INFILENAMESIZE 1000
#define SEQNAMESIZE 1000
#define DECREASE 0.99
#define MINSNPPERWINDOW 5
#define STRINGLENGTH 1000
#define VCF_HLENGTH 9 // number of fields in the VCF header line
#define MAX_CHROM_NAME_VCF 100
#define MAX_STATES_VCF 5
#define AFMAXLENGTH 9
#define RESULTS_DEFAULT 0
#define RESULTS_ALL 1

#define RSQUARE 0
#define DOM 1
#define ABSDOM 2
#define JUSTD 3
#define ABSD 4
#define ABSDOM2 5

#define ZERO '0'
#define ONE  '1'
#define GAP '-'
#define AD 'A'
#define CY 'C'
#define GU 'G'
#define TH 'T'
#define UN 'N'
#define ade 'a'
#define cyt 'c'
#define gua 'g'
#define thy 't'

#define BINARY 2
#define BINARY_WITH_GAPS 3
#define DNA 4
#define DNA_WITH_GAPS 5
#define STATESALL 8

#define OTHER_FORMAT -1
#define MS_FORMAT 0
#define FASTA_FORMAT 1
#define MACS_FORMAT 2
#define VCF_FORMAT 3
#define SF_FORMAT 4

#define INTOKS 8

#define N_REF 2
#define N_ALT 3
#define N_INFO 6
#define MAX_STATES 100
#define MAXINT 2147483647
#define LEN 5000
#define MINSNPS_THRESHOLD 3

char bits_in_16bits [0x1u << 16];

int linkage_disequilibrium;

int borderTol;

int VCF_header_lines;
char VCF_alignment_name [MAX_CHROM_NAME_VCF];
char VCF_alignment_name_prev_to_print [MAX_CHROM_NAME_VCF];
int VCF_first_SNP;
char *** chromList;
int chromList_SZ;

int nxtVCFalignment;

char runName[INFILENAMESIZE];

typedef float cor_t;



cor_t ABS (cor_t input);

double mainTime0;

int map[INTOKS];

int nofTokens;

int nofSamples;

typedef int al_t;
typedef double t_sfs;

struct pattern
{
	int x;
	int n;
	unsigned char f;
	int c;
	struct pattern * nxt;		
};

typedef struct 
{
 	al_t	*n; // number of valid bases in every site
	al_t	*n_t; // number of valid bases in every site
	al_t *x; // number of derived bases in every site
	al_t *x_t; // number of derived bases in every site
	float 	       * positions;
	int 	       * positionsInd;
	int 	       * p_t;
	unsigned char *folded;
	unsigned char *f_t;
	unsigned char userSetFolded;
	unsigned char userSFS;
	char* outgroupSequence;
  	char* outgroupName;
	int outgroupIndex;
	int sequences;
  	int foldedFormat; // if outgroup exists we set this to zero
	int 		 segsites;	
	int		 length;

	int minn;
	int minx;	
  	int maxn;
	int maxx;

	al_t * stateA;
	al_t * stateC;
	al_t * stateG;
	al_t * stateT;


	t_sfs * SFS;
	t_sfs * tmpSFS;
	int startSFS;
	int endSFS;

	int * ng;

	t_sfs * lb;
	t_sfs * ub;

	int params;

	al_t * n_pat;
	al_t * x_pat;
	unsigned char * f_pat;
	int *  c_pat;
	int patterns;

	t_sfs *** gridProbs;
	t_sfs * gridADs;
	t_sfs interval;
	t_sfs logAD0;	

	t_sfs * baseLikelihood;

	struct pattern * patternList;
	
	float sweepStep;

	int discSites;

	int VCFsamples;
	int * VCFsample_valid;

#ifdef _USE_PTHREADS
	t_sfs * threadScoresLIKSFS;
#endif

} alignment_struct;

alignment_struct * alignment;

double * rv;
double ** rvLUT;
double ** rvLUT2;


#define GRIDSIZE 300
#define LOCALGRIDSIZE 100

double div6; 


typedef struct
{
	double * val;
	int * sc;

} factorial_LUT;

factorial_LUT * factLUT;


typedef struct
{
  double        sfPos;
  int 		sfRealPos; 
  t_sfs		likelihood; 
  t_sfs		alpha; 
 
} clr_struct;

clr_struct * clr;

#ifdef _ANALYTICAL_SFS
typedef struct
{
	double t;
	double d;
} event;

int eventsTotal;
event * eventList;

mpfr_t ** analyticalSFS_LUT;

void computeAnalyticalSFS(int sequences, FILE * fp);
void computeAnalyticalSFSExpo(int sequences, double R, FILE * fpSFSo);
long double xExponential_Integral_Ei( long double x );
void mpfr_Exponential_Integral_Ei( mpfr_t EImp, mpfr_t xmp );
#endif

void printHeading (FILE * fp);
void computeCorrelationMatrixPairwise(alignment_struct * alignment, clr_struct * omega, int omegaIndex, int firstRowIndex, void * threadData, cor_t ** myCorrelationMatrix, char * lookuptable);
void applyCorrelationMatrixAdditions (clr_struct * omega, int omegaIndex, int firstRowIndex, cor_t ** correlationMatrix);
void overlapCorrelationMatrixAdditions (alignment_struct * alignment, clr_struct * omega, int lvw_i, int cvw_i, int * firstRowToCopy, int * firstRowToCompute, int * firstRowToAdd);
void shiftCorrelationMatrixValues (clr_struct * omega, int lvw_i, int cvw_i, int firstRowToCopy, cor_t ** correlationMatrix);


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
		       int* noSeparator,
		       char * samplefile_i,
		       int * generateVCFsamplelist,
		       char * chromfile_i,
		       int * generateVCFchromlist,
		       int * minsnps_threshold_user,
		       int * reports
);

void removeMonomorphicSites (int strictPolymorphic, int monomorphic, FILE * fp);

t_sfs computeBaseLikelihood();
t_sfs getAlpha (int sweepPosition, t_sfs * likelihood);


int isBinary(char input);

void printIntArray(int *array, int n);

void printUCharArray(unsigned char *array, int n);

void printFloatArray(float *array, int n);

void printCharArray(char *array, int n);

int isEndOfLine(char ent);

void goToLineEnd(FILE *fp);

int findFirstAlignment(FILE *fp, int format, FILE *fpVCFsamples, int generateVCFsamplelist, char * vcf_samples_filename, FILE * fpInfo);

int findNextAlignment(FILE *fp, int fileFormat);

void freeAlignment();

void freeAlignment_noSNPs();

void initializeAlignmentVariables();

int readAlignment(FILE *fp, int format, FILE * fpInfo, FILE * sfO, int minsnps_threshold_user, int alignmentIndex);

void compressAlignment(alignment_struct *alignment);

void knuthShuffle(int orgArray[], int arraySize, int start, int end);

int parseSample(FILE *fp, int targetTok, char **word, int *wordLength, int curTok);

void unGetChars(FILE *fp, char* chars, int n);

int isValidDNACharacter(char input);

void parseIndividual(char *word, char** genotype, char delimiter, int token, char* states, char phased, char unphased, int *ploidy, int *max_ploidy);

void ignoreAll(FILE * fp, char * ent);

int isSpace(char ent);

double XchooseY_ln(int x, int y);

void checkSNIPPositions (FILE* fp, alignment_struct * alignment, int index);

void createSFS (FILE * fpSFS, FILE * fpSFSo, int alignmentIndex, int * SFSsize);
void createCLR (int grid);
void createPROBS (int probGrid);
void createPatterns(void);
void computeFactLUT(void);

double gettime(void);

t_sfs likelihoodSFS_SNP(t_sfs * tmpSFS, int x, int n, unsigned char f, int c);

void updateSF(FILE * fpSFo);

void generate_VCF_chrom_list(FILE * fpin, FILE * fpVCFchroms, char * chromVCFfileName, FILE *fpInfo);
void extract_chromList(FILE * fpVCFchroms, FILE * fpInfo);

#ifdef _USE_PTHREADS

typedef struct
{
	int threadID;
	int threadTOTAL;

	int threadBARRIER;
	int threadOPERATION;

}threadData_t;

pthread_t * workerThreadL;
threadData_t * threadDataL;

int gridP;
void likelihoodSFS_thread(int tid, int threads);
void parallelComputationRVLUT_thread(int tid, int threads);
void computeAlpha_parallel(int grid);
void parallelComputationAlpha_thread(int tid, int threads);
void syncThreadsBARRIER();


#endif

#ifdef _DO_CHECKPOINTS
double lastCheckpointTime;
double checkPointInterval;
extern void writeCheckpoint();
#endif

#endif
