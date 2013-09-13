/*
	GADA v1.0 Genome Alteration Detection Algorithm 
    Copyright (C) 2008  Childrens Hospital of Los Angeles

	GADA is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    any later version.

    GADA is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with GADA.  If not, see <http://www.gnu.org/licenses/>.

	Author: 
		Roger Pique-Regi    piquereg@usc.edu

*/

#include <stdio.h>
#include <malloc.h>
#include "BaseGADA.h"
#include <math.h>


// Global variables with algorithm parameters

double T=5.0; //Backward elimination threshold
//double T2=5.0; //Segment collapse to base Non-Alteration level threshold
double BaseAmp=0.0;  //Base-level
double a=0.2; //SBL hyperprior parameter
double sigma2=-1; //Variance observed, if negative value, it will be estimated by the mean of the differences
			  // I would recommend to be estimated on all the chromosomes and as a trimmed mean. 
int MinLen=0; //Lenght in number of probes for a CNA segment to be called significan. 
int SelectClassifySegments=0; //Classify segment into altered state (1), otherwise 0
int SelectEstimateBaseAmp=1; //Estimate Neutral hybridization amplitude.
char *InputFile;
char *OutputFile;

void help_message(FILE *fd)
{
	//fprintf(fd,"# Welcome to GADA 1.0 \n");
	fprintf(fd,"# Usage:\n");
	fprintf(fd,"# \t GADA [-T 5] [-a 0.8] [-s -1] [-M 3] < input.txt > output.txt\n");
	fprintf(fd,"# \t input.txt is a single column text file with no header\n");
	fprintf(fd,"# Possible options:\n");
	//fprintf(fd,"# \t -i\t Input file. Otherwise standard input assumed\n");
	//fprintf(fd,"# \t -o\t Output file. Otherwise standard output assumed\n");
	fprintf(fd,"# \t -a\t is the SBL hyperprior parameter for a breakpoint. It is the \n");
	fprintf(fd,"# \t\t shape parameter of the Gamma distribution. Default value %g.\n",a);
	fprintf(fd,"# \t\t Higher (lower) value more (less) breakpoints expected a priori\n");
	fprintf(fd,"# \t -T\t is the backward elimination critical value for a breakpoint. \n");
	fprintf(fd,"# \t\t i.e. the minimum difference between the segment means divided\n");
	fprintf(fd,"# \t\t by sigma. The default value for T is %g\n",T);
	fprintf(fd,"# \t -M\t is the minimum size in number of probes for a segment to be \n");
	fprintf(fd,"# \t\t deemed significant. Default value is %d.\n",MinLen);	
	fprintf(fd,"# \t -s\t The variance estimate for the noise. If not provided or \n");
	fprintf(fd,"# \t\t negative it will be estimated from the provided data. \n");
	fprintf(fd,"# \t\t We recomend to estimate this on the entire data and not \n");
    fprintf(fd,"# \t\t separately for each chromosome.\n");
	fprintf(fd,"# \t -c\t Classify segments into altered state (L)oss -- (N)eutral \n");
	fprintf(fd,"# \t\t (G)ain). If c option is not specified, segments are returned \n");
	fprintf(fd,"# \t\t with their mean.\n");
	fprintf(fd,"# \t -b\t Mean amplitude associated to the Neutral state. If not \n");
	fprintf(fd,"# \t\t provided, and c option is used, then it is estimated as the\n");
	fprintf(fd,"# \t\t median value of all probes hybridization values after running \n");
	fprintf(fd,"# \t\t the algorithm. We recomend to estimate this on chromosomes\n");
	fprintf(fd,"# \t\t that are known to have a Neutral state on most areas. In some\n");
	fprintf(fd,"# \t\t cases this value may be known if we have been applied some \n");
	fprintf(fd,"# \t\t normalization, preprocessing or using another sample as ref.\n");
	fprintf(fd,"# \t -h \t Prints this help message.\n");
}

void help_and_exit(FILE *fd,int code)
{	
	fprintf(fd, "Invalid syntax. Use GADA -h for help\n");
	exit(code);
}


void Configure(int ac, char *av[])
{
	int CLcount,i;
	FILE *fd;

	// Parse the command line
	CLcount=1;


	fd=stdout;
	fprintf(fd,"# NumArgs = %d \n",ac);
	fprintf(fd,"# CallStr = ");

	for(i=0;i<ac;i++)
		fprintf(fd,"%s ",av[i]);
	fprintf(fd,"\n# Parsing Arguments: \n");


	while (CLcount < ac)
	{
		if (0 == strncmp (av[CLcount], "-h", 2))
		{
			help_message(stderr);
			exit(0);
		}

		else if (0 == strncmp (av[CLcount], "-c", 2))
		{
			SelectClassifySegments = 1;			
			CLcount+=1;
			fprintf(fd,"# -c option activated to classify segments into altered states\n");
		}
/*		else if (0 == strncmp (av[CLcount], "-i", 2)) //! Input file
		{
			if(0 == strncmp (av[CLcount+1], "-", 1))help_and_exit(stderr,1);
			strcpy(InputFile,av[CLcount+1]);
			CLcount += 2;
			fprintf(fd,"# Input file: %s \n",InputFile);
		}
		else if (0 == strncmp (av[CLcount], "-o", 2)) //! Output File
		{
			if(0 == strncmp (av[CLcount+1], "-", 1))help_and_exit(stderr,1);
			strcpy(OutputFile,av[CLcount+1]);
			CLcount += 2;
			fprintf(fd,"# Output file: %s \n",OutputFile);
		}
*/		else if (0 == strncmp (av[CLcount], "-a", 2)) //! a parameter
		{
			if(0 == strncmp (av[CLcount+1], "-", 1))help_and_exit(stderr,1);
			sscanf (av[CLcount+1], "%lf", &a);
			CLcount += 2;
			fprintf(fd,"# a= %g \n",a);
		}
		else if (0 == strncmp (av[CLcount], "-T", 2)) //! T parameter
		{
			if(0 == strncmp (av[CLcount+1], "-", 1))help_and_exit(stderr,1);
			sscanf (av[CLcount+1], "%lf", &T);
			CLcount += 2;
			fprintf(fd,"# T= %g \n",T);
		}
		else if (0 == strncmp (av[CLcount], "-M", 2)) //! MinLen parameter
		{
			if(0 == strncmp (av[CLcount+1], "-", 1))help_and_exit(stderr,1);
			sscanf (av[CLcount+1], "%d", &MinLen);
			CLcount += 2;
			fprintf(fd,"# MinLen= %d \n",MinLen);
		}
		else if (0 == strncmp (av[CLcount], "-b", 2)) //! BaseAmp parameter
		{
			//if(0 == strncmp (av[CLcount+1], "-", 1))help_and_exit(stderr,1);
			sscanf (av[CLcount+1], "%lf", &BaseAmp);
			CLcount += 2;
			fprintf(fd,"# BaseAmp= %g \n",BaseAmp);
			SelectEstimateBaseAmp=0;
		}
		else if (0 == strncmp (av[CLcount], "-s", 3)) //! sigma2 parameter
		{
			//if(0 == strncmp (av[CLcount+1], "-", 1))help_and_exit(stderr,1);
			sscanf (av[CLcount+1], "%lf", &sigma2);
			CLcount += 2;
			fprintf(fd,"# sigma2= %g \n",sigma2);
		}		
		else
		{			
			help_and_exit(stderr,1);
		}
	}
}

int main(int argc, char *argv[])
{
	int M=1000;
	int i;
	int *Iext;
	int K;
	double *tn;
	int *SegLen;
	double *SegAmp;
	double *SegState;
	double *Wext;

//	FILE *fin,*fout;

	tn=calloc(M,sizeof(double));

	fprintf(stdout,"#GADA v1.0 Genome Alteration Detection Algorithm\n");
	fprintf(stdout,"#Copyright (C) 2008  Childrens Hospital of Los Angeles\n");
	fprintf(stdout,"# author: Roger Pique-Regi piquereg@usc.edu\n");
    
	Configure(argc,argv);

	fprintf(stdout,"#Parameter setting: a=%g,T=%g,MinLen=%d,sigma2=%g,BaseAmp=%g\n",a,T,MinLen,sigma2,BaseAmp);
	
	i=0;
	while(!feof(stdin))
	{
		fscanf(stdin,"%lf",&tn[i++]);
		if(i>=M){
			M=M+1000;		
			tn=realloc(tn,M*sizeof(double));
		}
	}
	M=i-1;
	tn=realloc(tn,M*sizeof(double));

	fprintf(stdout,"# Reading M=%d probes in input file\n",M);
	
	K=SBLandBE(tn,M,&sigma2,a,0,0,&Iext,&Wext);

	//K=SBLandBE(tn,M,&sigma2,a,T,MinLen,&Iext,&w);
    fprintf(stdout,"# Overall mean %g\n",Wext[0]);
	fprintf(stdout,"# Sigma^2=%g\n",sigma2);
	fprintf(stdout,"# Found %d breakpoints after SBL\n",K);

	
	BEwTandMinLen(Wext,Iext,&K,sigma2,T,MinLen); 
    fprintf(stdout,"# Kept %d breakpoints after BE\n",K);

	SegLen=calloc(K+1,sizeof(int));
	SegAmp=calloc(K+1,sizeof(double));
	IextToSegLen(Iext,SegLen,K);
	IextWextToSegAmp(Iext,Wext,SegAmp,K);
    fprintf(stdout,"# Making segments\n");

	//Collapse Segments
	if(SelectClassifySegments==1)
	{
		if(SelectEstimateBaseAmp==1)
		{
			BaseAmp=CompBaseAmpMedianMethod(SegLen,SegAmp,K);
			fprintf(stdout,"# Estimating BaseAmp\n");
		}
		fprintf(stdout,"# BaseAmp=%g \n",BaseAmp);
		fprintf(stdout,"# Classify Segments \n",BaseAmp);
		

		SegState=calloc(K+1,sizeof(double));
		for(i=0;i<=K;i++)SegState[i]=SegAmp[i];
		CollapseAmpTtest(SegState,SegLen,K,BaseAmp,sigma2,T);
	}


	if(SelectClassifySegments==0)
	{
    	fprintf(stdout,"Start\tStop\tLength\tAmpl\n");
		for(i=0;i<K+1;i++)
			fprintf(stdout,"%d\t%d\t%d\t%g\n",Iext[i]+1,Iext[i+1],SegLen[i],SegAmp[i]);
	}
	else if(SelectClassifySegments==1)
	{
	   	fprintf(stdout,"Start\tStop\tLenght\tAmpl\tState\n");
		for(i=0;i<K+1;i++)
		{
			fprintf(stdout,"%d\t%d\t%d\t%g\t",Iext[i]+1,Iext[i+1],SegLen[i],SegAmp[i]);
			if(SegState[i]>BaseAmp)
				fprintf(stdout,"G");
			else if(SegState[i]<BaseAmp)
				fprintf(stdout,"L");
			else 
				fprintf(stdout,"N");
			fprintf(stdout,"\n");
		}
		

	}




}