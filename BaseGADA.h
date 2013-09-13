/*=================================================================
 * BaseGenomeBreaks.h 
 *
 *=================================================================*/
/*
	This File is part of GADA

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

#ifndef _BaseGADA_H_
#define _BaseGADA_H_

//#include "matlabdefines.h"

void reconstruct (double *wr,int M,double *aux_vec);
void BubbleSort (int *I,int L);
void doubleBubbleSort (double *D,int *I,int L);
void TrisolveREG(double *t0,double *tu,double *tl,double *coef,double *sol,int sizeh0);
void DiagOfTriXTri(double *ll,double *l0,double *lu,double *rl,double *r0,double *ru,double *d,int N);
void tridiagofinverse(double *t0,double *tl,double *itl,double *it0,double *itu,int N,double *d,double *e);
void ForwardElimination(double *A,int N);
void BackSubstitution(double *A,int N);
void BackwardElimination(double *A,int N);
void TriSolveINV(double *AA,int M, int N, double *x,double *d,double *e);
void ComputeH(double *h0, double *h1, int M);
void ComputeFdualXb(int M, double *b);
// 20080119 REMOVED void ComputeHs(int *s,double *a,int M,int Ms,double *h0,double *h1);
void ComputeHs(int *s,int M,int Ms,double *h0,double *h1);
void TriSymGaxpy(double *t0, double *t1, double *x, int M, double *y);
void ComputeT(double *h0,double *h1,int M,double *alfa,double sigma,double *t0,double *tl,double *tu);
int findminus(double *alpha,int Ms,double maxalpha,int *sel);
int simpletresholding(double *inputvector,int N,double thres,double *disc);
void computesegmentmeans(double *inputvector,int N,double *disc,int numdisc,double *amp);
void reconstructoutput(double *rec,int N,double *disc,int numdisc,double *amp);
int SBL(
    double *y, //I -- 1D array with the input signal
    int *I, //IO -- 1D array with the initial (final) candidate breakpoints
    double *alpha, //I -- 1D array with the initial (final) hyperparameter inv. varainces.
    double *w, //O -- 1D array containing the breakpoint weigths or posterior mean. 
    double *sigw, //O -- Posterior variances, I would not really need them since they can be computed from alpha and H
    int M, //Initial size of the array in y
    int *K, //Size of the I alpha w
    
    //Algorithm parameters:
    double sigma2, //Noise estimated 
    double a,      //
    double b,       
    double maxalpha,  //Basis reduction parameter 
    int    maxit,     //Max number of iterations
    double tol,       //Tolerance for convergence
    int debug       //verbosity... set equal to 1 to see messages  0 to not see them
    );

int BEthresh( //To eliminate...   
    double *Scores,
    int Nscores,
    double *wr,
    int *indsel,
    int *pointNumRem,
    double *pointTau
    );

int SBLandBE( //Returns breakpoint list lenght.
    double *tn,
    int M,  //length of the noisy signal tn
    double *sigma2, //If sigma2 < 0, compute sigma2 (Input/Output)
    double a,      // SBL parameter
    double T,      // Threshold to prune
    int MinSegLen, //Minimum length of the segment.
    int **pI,   //Returns breakpoint positions
    double **pw //Returns breakpoint weights.
    //int *pK    
    );


void Project(
    double *y,
    int M,
    int *I,
    int L,
    double *xI,
    double *wI
    );
void IextToSegLen(
	int *Iext, // Enters Iext
	int *SegLen,		// Outputs SegLen.. can be the same as Iext?
	int K			// Length Iext - 1
	);
void IextWextToSegAmp(
	int *Iext, 
	double *Wext,
	double *SegAmp,
	int K 
	);
void CompZ(// computes z=F'y for entire possilbe breakpoint positions (normalized PWC)
	double *y,
	double *z,
	int M
	);
void ComputeHsIext(
    //input variables:
    int *Iext,     // Indices of selection,
    int K,     // Length of the indices, 
    double *h0, // Returning diagonal of H,
    double *h1  // Returning upper diagonal of H
    );
void ProjectCoeff ( //IextYobs2Wext
    double *y,
    int M,
    int *Iext,
    int K,
    double *Wext
	);
void 
CollapseAmpTtest(//Uses a T test to decide which segments collapse to neutral
	double *SegAmp, //Segment Amplitudes (input output)
	const int *SegLen, //Segment Lengths
	int K, //Number of segments.
	double BaseAmp, //Reference amplitude to compare
	double sigma2,  //Reference noise 
	double T		//Critical value that decides when to colapse
	);
double // Returns BaseAmp corresponding to the base level.
CompBaseAmpMedianMethod( //Computes the median recontruction level, as baseline level.
	const int *SegLen,    //Lengths corresponding to the amplitudes
	const double *SegAmp, //Amplitudes !!! assumed already ordered...
	int K
	);

void 
ClassifySegments(
	double *SegAmp,  
	int *SegLen,
	double *SegState,
	int K,
	double BaseAmp,
	double ploidy,
    double sigma2,  //Reference noise 
	double T		//Critical value that decides when to colapse
	);

void
ComputeTScores(
	const double *Wext,
	const int *Iext,
	double *Scores,
	int K,
	int start,
	int end
	);

int BEwTscore(
    double *Wext,  //IO Breakpoint weights extended notation...
	int *Iext,     //IO Breakpoint positions in extended notation...
	double *tscore,
    int *pK,       //IO Number breakpoint positions remaining.
    double T      //IP  Threshold to prune
	);

int BEwTandMinLen( //Returns breakpoint list lenght. with T and MinSegLen
    double *Wext,  //IO Breakpoint weights extended notation...
	int *Iext,     //IO Breakpoint positions in extended notation...
    int *pK,       //IO Number breakpoint positions remaining.
	double sigma2, //IP If sigma2 
    double T,      //IP  Threshold to prune,  T=T*sqrt(sigma2);
    int MinSegLen  //IP Minimum length of the segment.
);
int
RemoveBreakpoint(
	double *Wext,
	int *Iext,
	int K,
	int jrem
	);


#endif //_BaseGADA_H_
