/*=================================================================
 * BaseGADA.c 
 * All basic functions for SBL, BE, etc. 
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
		Jordi Monso-Varona

*/
#include "BaseGADA.h"
#include <math.h>

#ifdef _MATLAB_
    #include "matlabdefines.h"
#endif //_MATLAB_
#ifndef _MATLAB_
	#include "consoledefines.h"
#endif //_MATLAB_

#define min(x,y) x<y?x:y

/* General Notation definition
K		= number of breakpoints, segments is K+1 
M		= number of probes in chromosome / unit of analysis (can be a chr arm)
I Ir Ired	Breakpoint set, double array I[0] contains first breakpoint I[K-1] contains Kth breakpoint length 0 if no breakpoints
Ie Iext		Breakpoint set, extended notation I[0] contains 0 I[K] contains Kth breakpoint I[K+1]=M
We		Aligned with Ie We[0] contains the overall mean and We[K+1] nothing
SegAmp  Segment amplitudes lenght K+1
SegLen  Segment longitudes lenght K+1 and should sum up to M?
*/



void IextToSegLen(
	int *Iext, // Enters Iext
	int *SegLen,		// Outputs SegLen.. can be the same as Iext?
	int K			// Length Iext - 1
	)
{
	int k;
	for(k=1;k<=K+1;k++)
		SegLen[k-1]=Iext[k]-Iext[k-1];
}


void IextWextToSegAmp(
	int *Iext, 
	double *Wext,
	double *SegAmp,
	int K 
	)
{
	int k;
	double M,TotalMean,AuxMean;
	M=Iext[K+1];

	TotalMean=Wext[0];
	SegAmp[0]=0;
	for(k=1;k<=K;k++)
		SegAmp[k]=Wext[k]/sqrt((double)(M-Iext[k])*(double)Iext[k]/M)+SegAmp[k-1];
	AuxMean=0;
	for(k=0;k<=K;k++)
		AuxMean=AuxMean+SegAmp[k]*(double)(Iext[k+1]-Iext[k]);
	AuxMean=AuxMean/M;
	for(k=0;k<=K;k++)
		SegAmp[k]=SegAmp[k]-AuxMean+TotalMean;

}



// Implements: Reconstruct.m -- Recontructs x from w with zero mean
void reconstruct (
    double *wr,
    int M,
    double *aux_vec
    )
{
    int i=0;
    double aux_double=0;
    
    for (i=0;i<M;i++){
        aux_vec[i]=0;
    }
    for (i=0;i<M;i++){
        aux_vec[i]=(-1)*sqrt(((double)(M-i)*(double)(i+1))/(double)(M+1));
    }
    for (i=0;i<M;i++){
        wr[i]=(wr[i]/aux_vec[i]);
    }
    for (i=0;i<M-1;i++){
        wr[M-(i+2)]=wr[M-(i+2)]+wr[M-(i+1)];
    }
    aux_double=0;
    for (i=0;i<M;i++)
    {
        aux_double=aux_double+wr[i];
    }
    aux_double=(aux_double/M);
    for (i=0;i<M;i++)
    {
        wr[i]=wr[i]-aux_double; // I have y mean in aux_double
    }
}
    
    


/* Bubble Sort */
void BubbleSort (
    int *I,
    int L
    )
{
    int i=0;
    int aux=0;
    
    for (i=0;i<L-1;i++) {
        if (I[i]>I[i+1]){
            aux=I[i+1];
            I[i+1]=I[i];
            I[i]=aux;
            i=-1;
        }
    }
}
void doubleBubbleSort (
    double *D,
	int *I,
    int L
    )
{
    int i=0;
    double Daux=0;
	int Iaux=0;

    for (i=0;i<L-1;i++) {
        if (D[i]>D[i+1]){
            Daux=D[i+1];
            D[i+1]=D[i];
            D[i]=Daux;
			
			if(I!=NULL)
			{
				Iaux=I[i+1];
				I[i+1]=I[i];
				I[i]=Iaux;
			}

            i=-1;
        }
    }
}


/////////////  TRIDIAGONAL MATRIX OPERATION FUNCTIONS ////////////////////
void 
TrisolveREG(
    //input variables:
    double *t0,
    double *tu,
    double *tl,
    double *coef,
    double *sol,
    int sizeh0
    )
{
    double m=0;
    int i=0;
    
    for (i=0;i<sizeh0-1;i++){
        //forward elimination
        m=tl[i]/t0[i];          
        t0[i+1]=t0[i+1]-m*tu[i];
        coef[i+1]=coef[i+1]-coef[i]*m;
    }
    sol[sizeh0-1]=coef[sizeh0-1]/t0[sizeh0-1];
    for (i=0;i<sizeh0-1;i++){
        sol[sizeh0-2-i]=(coef[sizeh0-2-i]-(tu[sizeh0-2-i]*sol[sizeh0-1-i]))/t0[sizeh0-2-i];
    }
}
            
// DiagOfTriXTri  Diagonal of L*R where L,R tridiagonal (MexDiagOfTriXTri)
void   
DiagOfTriXTri(
    //Input variables:
    double *ll, // ll Lower diagonal
    double *l0, // l0 Central diagonal
    double *lu, // lu Upper diagonal
    double *rl, // rl Lower diagonal
    double *r0, // r0 Central diagonal
    double *ru, // ru Upper diagonal
    double *d,  // d=DiagOfTriXTri(ll,l0,lu,rl,r0,ru)
    int N  //Number of variables or length of l0, 
    )
{
    register int i;
    
    d[0]=l0[0]*r0[0]+lu[0]*rl[0];
    for (i=1;i<N-1;i++){
        d[i]=ll[i-1]*ru[i-1]+l0[i]*r0[i]+lu[i]*rl[i];
    }
    d[N-1]=ll[N-2]*ru[N-2]+l0[N-1]*r0[N-1];
    
}

/* Tridiagonal de la inversa imput parameters: (tl,t0,d,e)
 * 
 */
void   
tridiagofinverse(
    //Input variables:
    double *t0,
    double *tl,
    double *itl,
    double *it0,
    double *itu,
               
    int N,  //Number of variables
    double *d, //Array cointaining the upper diagonal of pseudo-inverse
    double *e //Array cointaining the lower diagonal of pseudo-inverse
    )
{   
    int i;    
    it0[0]=1/(t0[0]*(1-e[0]*d[0]));
    
    for(i=1;i<(N-1);i++)
    {
        it0[i]=1/((t0[i]-(tl[i-1]*d[i-1]))*(1-e[i]*d[i]));
    }
    it0[N-1]=1/(t0[N-1]-tl[N-2]*d[N-2]);
    
    for (i=0;i<(N-1);i++)
    {
        itl[i]=-e[i]*it0[i];
        itu[i]=-d[i]*it0[i+1];
    }    
}


//////////               MEXTRISOLVE            ////////////////////
void   
ForwardElimination(
    //Input variables:
    double *A, //2D Array containing the [tu';t0';tl';b']
               // Retuns [d';t0';tl';z']
    int N  //Number of variables
    )
{
    register int i=0;
    register double *d,*tl,*t0,*z,*b;
    register double m,m2;
    
    
    d=&A[0];
    t0=&A[1];
    tl=&A[2];
    b=&A[3];
    z=b;
    
    //First iteration 
    //%normalization of the first row
    m2=*t0;
    *d/=m2; //d=d/m
    *b/=m2; //b=b/m
    //Pointer relocation
    b=b+4;
    t0=t0+4;
    
    for(i=1;i<(N-1);i++)
    {
        m=*tl;
        tl+=4;
        m2=(*t0)-(*d)*m;
        t0+=4;        
        d+=4;
        *d/=m2; //tu=tu/m2      
        *b=(*b-m*(*z))/m2;
        b+=4;
        z+=4;
    }
    
    //Last step    
    m=*tl; 
    m2=(*t0)-(*d)*m;
    *b=(*b-m*(*z))/m2;
}

void   
BackSubstitution(
    //Input variables:
    double *A, //2D Array containing the [d';!t0';!tl';z'] from Forward Elimination
               // Retuns [d';t0';tl';x']

    int N  //Number of variables
    )
{
    register int i=0;
    register double *d,*z,*x;
       
    // +4 jumps pointer to the element in next column of A
    // Pointer initialization
    z=&A[N*4-1];
    x=z-4;
    d=&A[(N-2)*4];
    
    for(i=0;i<(N-1);i++)
    {
        x[0]-=x[4]*x[-3];
        x-=4;
    }
    
}   
void   
BackwardElimination(
    //Input variables:
    double *A, //2D Array containing the [tu';t0';tl';!b']
               // Retuns [tu';t0';e';b']
    int N  //Number of variables
    )
{
    register int i=0;
    double *e,*tu;
    register double *t0;
    e=&A[4*N-6];   
    t0=&A[4*N-3];    
    tu=&A[4*(N-2)];   
    *e/=*t0;
    t0-=4;
    
    for(i=1;i<(N-1);i++)
    {        
        t0[-3]/=t0[0]-(t0[-1])*(t0[+1]);
        t0-=4;        
    }
    
}


void
TriSolveINV(double *AA,int M /*rows*/, int N/*columns*/, double *x,double *d,double *e)
{
    int i;
    
    double aux;
    
    
    /* Save backup copy of tu d*/
    for(i=0;i<N-1;i++) 
    {
        d[i]=AA[(i<<2)];
//        x[i]=A[(i<<2)+3];
    }     
    
    /*Run algorithm*/
    ForwardElimination(AA,N);
    BackSubstitution(AA,N);
    
    /*Store solution and Recompose A*/
    for(i=0;i<N-1;i++) 
    {
        aux=AA[(i<<2)];
        AA[(i<<2)]=d[i];
        d[i]=aux;
        
        x[i]=AA[(i<<2)+3];
    }     
    x[i]=AA[(i<<2)+3];   //x[N-1]=A[N<<2-1];

    
    BackwardElimination(AA,N);
        
    /*Store solution on the output vectors*/
    for(i=0;i<N-1;i++){
        e[i]=AA[(i<<2)+2];
    }    
}

////////////////////////////////////////////////////////////////////////

void 
ComputeFdualXb(
    //imput variables:
    int M,
    double *b
    )
{
    int i=0;
    double a=0;
    
    /* diff */
    for (i=0;i<M-1;i++){
        b[i]=b[i+1]-b[i];
        
    }
    b[M-1]=0;
    for (i=0;i<M-1;i++){
        a=(double)(M-1-i)*(double)(i+1)/(double)M;
        b[i]=b[i]*sqrt(a);
    }
}

void CompZ(// computes z=F'y for entire possilbe breakpoint positions (normalized PWC)
	double *y,
	double *z,
	int M
	)
{
	double LeftSum;
	double RightSum;
	double *aux;
	int m;

	aux=myCalloc(M-1,sizeof(double));
	
	LeftSum=0;
	RightSum=0;
	for(m=1;m<M;m++)
	{
		LeftSum=LeftSum+y[m-1];
		aux[m-1]=(-1)*sqrt( (double)(M-m) / ((double)M*(double)m) )*LeftSum;	
	}
	for(m=M-1;m>=1;m--)
	{
		RightSum=RightSum+y[m];
		z[m-1]=aux[m-1]+sqrt( (double)m / ((double)(M-m)*(double)(M)) )*RightSum;	
	}
//	myPrintf("\n CHECKING: Rsum=%g,Lsum=%g,M=%d\n",RightSum,LeftSum,M);
//	myPrintf("\n z CHECKING: z[0]=%g,z[1]=%g,z[2]=%g,z[M-3]=%g,z[M-2]=%g\n",z[0],z[1],z[2],z[M-3],z[M-2]);
//	myPrintf("\n z CHECKING: z[0]=%g,z[1]=%g,z[2]=%g,z[M-1]=%g,z[M]=%g\n",aux[0],aux[1],aux[2],aux[M-1],aux[M]);


	myFree(aux);
}


/* Compute H   fuction [h0,h1]=CompH(dim) */
void   
ComputeH(
    //Input variables:
    double *h0,
    double *h1,
    int M  //Number of variables
    )
{
    int i;

    for (i=0;i<M-1;i++){
        h0[i]=((double)(M-1-i)*(double)(i+1))/(double)M;
    }
    for (i=0;i<M-2;i++){
        h1[i]=-sqrt((h0[i+1]*h0[i]));
    }
    for (i=0;i<M-1;i++){
        h0[i]=2*h0[i];
    }      
}


// Computes the H at the vector of indices...
void 
ComputeHs(
    //input variables:
    int *s,     // Indices of selection,
    int M,      // Length of the chormosome, 
    int K,     // Length of the indices, 
    double *h0, // Returning diagonal of H,
    double *h1  // Returning upper diagonal of H
    )
{
    int i;
    double iC,iL,iR;
    //double M;
    //M=(double)MM;
    
    i=0;
    
    if (K==1)
    {   
        iL=0;iC=(double)(s[i]+1);iR=M;    
        h0[i]=((M-iC)*iC/M)*(iR-iL)/((iR-iC)*((iC-iL)));
        //h1[i]=0;   
    }
    else
    {
        iL=0;iC=(double)(s[i]+1);iR=(double)(s[i+1]+1);
        h0[i]=((M-iC)*iC/M)*(iR-iL)/((iR-iC)*((iC-iL)));
        h1[i]=-sqrt((M-iC)*iC*(M-iR)*iR)/(M*(iR-iC));
        for (i=1;i<K-1;i++){
            //        h0[i]=((double)((M-s[i]-1)*(s[i]+1))/(double)M)*(double)(s[i+1]-s[i-1])/(double)((s[i+1]-s[i])*((s[i]-s[i-1])));
            //        h1[i]=sqrt((double)(((M-s[i]-1)*(s[i]+1))*(M-s[i+1]-1)*(s[i+1]+1)))/(double)(M*(s[i+1]-s[i]));
            iL=iC;iC=iR;iR=(double)(s[i+1]+1);
            h0[i]=((M-iC)*iC/M)*(iR-iL)/((iR-iC)*((iC-iL)));
            h1[i]=-sqrt((M-iC)*iC*(M-iR)*iR)/(M*(iR-iC));
        }
        //i=K-1
        iL=iC;iC=iR;iR=M;
        h0[i]=((M-iC)*iC/M)*(iR-iL)/((iR-iC)*((iC-iL)));
        //    h0[i]=((double)((M-s[i]-1)*(s[i]+1))/(double)M)*(double)(M-s[i-1]-1)/(double)((M-s[i]-1)*((s[i]-s[i-1]))); //s[K]=M;
    }
       
}

void 
ComputeHsIext(
    //input variables:
    int *Iext,     // Indices of selection,
    int K,     // Length of the indices, 
    double *h0, // Returning diagonal of H,
    double *h1  // Returning upper diagonal of H
    )
{
    int i,M;
    double iC,iL,iR;
    //double M;
    //M=(double)MM;
    
	M=Iext[K+1];

    //iL=0;iC=(double)Iext[1];iR=(double)s[2];
    for (i=1;i<K;i++){
       //iL=iC;iC=iR;iR=(double)(Iext[i+1]);
	   iL=Iext[i-1];iC=Iext[i];iR=(double)Iext[i+1];
       h0[i-1]=((M-iC)*iC/M)*(iR-iL)/((iR-iC)*((iC-iL)));
       h1[i-1]=-sqrt((M-iC)*iC*(M-iR)*iR)/(M*(iR-iC));
    }
    //i=K
    //iL=iC;iC=iR;iR=M;
    iL=Iext[i-1];iC=Iext[i];iR=(double)Iext[i+1];
    h0[i-1]=((M-iC)*iC/M)*(iR-iL)/((iR-iC)*((iC-iL)));       
}

void 
TriSymGaxpy(
    //input variables:
    double *t0,
    double *t1,
    double *x,
    int M,
    double *y
    
    )
{
    int i;
    
    if (M==1){
        y[0]=t0[0]*x[0];
    }
    else
    {        
        y[0]=t0[0]*x[0]+t1[0]*x[1];        
        for (i=1;i<M-1;i++){
            y[i]=t0[i]*x[i]+t1[i]*x[i+1]+t1[i-1]*x[i-1];
        }
        y[M-1]=t0[M-1]*x[M-1]+t1[M-2]*x[M-2];
    }
}

void        
ComputeT(
    double *h0,
    double *h1,
    int M,
    double *alfa, 
    double sigma, /*pass 1 if scale not available*/
    double *t0,
    double *tl,
    double *tu
    /* in theory there is a parameter called 'scale' but we don't use it */
    )
{
    int i=0;
    
    for (i=0;i<M;i++){
        t0[i]=(h0[i]*alfa[i]*sigma)+1;
        if (i<M-1){
            tu[i]=h1[i]*alfa[i+1]*sigma;
            tl[i]=h1[i]*alfa[i]*sigma;
        }
    }
}        
        
int 
findminus(/*allocates all those lower-than-maxalpha alpha's INDEXES in sel vector  */
    //Input variables:
    double *alpha,
    int K,
    double maxalpha,
    int *sel //index vector
    )
{
    int i,j=0;
    for (i=0;i<K;i++){
        if (alpha[i]<maxalpha){
            sel[j]=i;
            j++;
        }
    }
    return (j);
}
    

/* simpletresholding 
 * Applies a simple tresholding algorithm to find the discontinuities
 */
int   //Returns the number of discontinuities
simpletresholding(
    //Input variables:
    double *inputvector, //1D Array containing the input vector
    int N,  //Vector length
    double thres, //Treshold value
    //Output variables:
    double *disc //1D empty array, with memory already allocated for finding up to N discontinuities positions
    )
{
    int i=0;
    int numdisc=0; //Number of discontinuities
    int state=0; //To keep up the state 0(<tresh),1(>tresh)
    
    //First step: Find discontinuities
    if(inputvector[0]>=thres) state=1;    
    for(i=0;i<N;i++){
        if (state==0)
        {
            if (inputvector[i]>=thres) 
            {
                disc[numdisc]=(double)(i+1); //Matlab uses 1..N indices instead of 0..(N-1)
                numdisc++;
                state=1;
//                myPrintf("i%d numdisc%d State change state%d.\n",i,numdisc,state);
            }
        }
        else
        {
            if (inputvector[i]<thres) 
            {
                disc[numdisc]=(double)(i+1); //Matlab uses 1..N indices instead of 0..(N-1)
                numdisc++;
                state=0;
  //              myPrintf("i%d numdisc%d State change state%d.\n",i,numdisc,state);
            }            
        }
    }    
    
    return numdisc;
}


/* computesegmentmeans 
 * Computed segment means.
 */
void
computesegmentmeans(
    //Input variables:
    double *inputvector, //1D Array containing the input vector
    int N,  //Vector length
    double *disc, //1D empty array, with memory already allocated for finding up to N discontinuities positions
    int numdisc, //Length of the discontinuity vector
    //Output variables:    
    double *amp  //1D empty array, containing the amplitudes of the segments between discotinuities, 
    )
{
    int i,j;
    int prevdisc=0;

    for(i=0;i<numdisc;i++)
    {
        amp[i]=0;
        for(j=prevdisc;j<(int)(disc[i]-1);j++)
            amp[i]+=inputvector[j];
        amp[i]=amp[i]/((double)(disc[i]-1-prevdisc));
        prevdisc=(int)(disc[i]-1);
    }
    amp[i]=0;
    for(j=prevdisc;j<N;j++)
        amp[i]+=inputvector[j];
    amp[i]=amp[i]/(double)(N-prevdisc);   
}

/* reconstructoutput 
 * Reconstruct signal from the discontinuity location, and segment 
 * amplitudes
 */
void
reconstructoutput(
    //Output variables:
    double *rec, //1D Array returning the reconstructed vector
    //Input variable:
    int N,  //Vector length
    double *disc, //1D empty array, with memory already allocated for finding up to N discontinuities positions
    int numdisc, //Length of the discontinuity vector
    //Output variables:    
    double *amp  //1D empty array, containing the amplitudes of the segments between discotinuities, 
    )
{
    int i,j;
    
    j=0;
    for(i=0;i<numdisc;i++)
        for(;j<(int)(disc[i]-1);j++)
            rec[j]=amp[i];
    //If no discontinuities or from the last to the end 
    for(;j<N;j++)
        rec[j]=amp[i];    
}

/* +++++++++++============================================++++++++++++++++++++++++++++++************************************** */

/* SBL function 
 * returns number of EM iterations
 */
int SBL(
    double *y, //I -- 1D array with the input signal
    int *I, //IO -- 1D array with the initial (final) candidate breakpoints
    double *alpha, //I -- 1D array with the initial (final) hyperparameter inv. varainces.
    double *w, //O -- 1D array containing the breakpoint weigths or posterior mean. 
    double *sigw, //O -- Posterior variances, I would not really need them since they can be computed from alpha and H
    int M, //Initial size of the array in y
    int *pK, //Size of the I alpha w
    
    //Algorithm parameters:
    double sigma2, //Noise estimated 
    double a,      //
    double b,       
    double maxalpha,  //Basis reduction parameter 
    int    maxit,     //Max number of iterations
    double tol,       //Tolerance for convergence
    int debug         //verbosity... set equal to 1 to see messages  0 to not see them
    ){
    int n,i,K;
    int M0,sizesel;
    double *w0;
    double mynorm,myaux;
   
    //Extra memory necessary to store variables during the algorithm.... (require memory initialization)
    double *h0;
    double *h1;
    double *xx;
    double *z;
    double *t0;
    double *tl;
    double *tu;
    double *wpred;
    double *d;
    double *e;
    int *sel;
    double *AA;
    double *yy;
    
    K=*pK;
    M0=M-1;
    
    //non discontinuities case
    if (K==0){
        w[0]=-1;
        *pK=K;
        sigw[0]=-1;
        return 0;
    }
    //Memory initialization of the outputs. (Already with mem assigned)
//    w=myCalloc(K,sizeof(double));
//    sigw=myCalloc(K,sizeof(double));
    
    //Memory initialization (internal to be freed)
    yy=myCalloc(M,sizeof(double));
    t0=myCalloc(M,sizeof(double));
    tl=myCalloc(M-1,sizeof(double));
    tu=myCalloc(M-1,sizeof(double));
    AA=myCalloc(4*M,sizeof(double));
    d=myCalloc(M0-1,sizeof(double));
    e=myCalloc(M0-1,sizeof(double));
    h0=myCalloc(M0,sizeof(double));
    h1=myCalloc(M0-1,sizeof(double));
    wpred=myCalloc(K,sizeof(double));
    sel=myCalloc(K,sizeof(int));
    z=myCalloc(M0,sizeof(double));    
    xx=myCalloc(K,sizeof(double));//myDoubleMAlloc(K);

    //Create a copy of the input
    for(i=0;i<M;i++)
        yy[i]=y[i];
    
   //myPrintf("\n\nOPERATIONS BEFORE LOOP:\n");
    
   // myPrintf("\nCOMPUTE H\n");
    
    ComputeH(h0,h1,M);
    //myPrintf("H0\n\nh0[0]:%g\nh0[1]:%g\nh0[2]:%g\nh0[3]:%g\nh0[%d]:%g\n",h0[0],h0[1],h0[2],h0[3],M0-1,h0[M0-1]);
    //myPrintf("H1\n\nh1[0]:%g\nh1[1]:%g\nh1[2]:%g\nh1[3]:%g\nh1[%d]:%g\n",h1[0],h1[1],h1[2],h1[3],M0-2,h1[M0-2]); 
    //myPrintf("\nCOMPUTE F DUAL\n");   
    
    ComputeFdualXb(M,yy);//checked
    
    w0=yy;//w0 now has M-1 dimmension
    //myPrintf("W0\n\nw0[0]:%g\nw0[1]:%g\nw0[2]:%g\nw0[3]:%g\nw0[%d]:%g\n",w0[0],w0[1],w0[2],w0[3],M0-1,w0[M0-1]);
    for (i=0;i<K;i++){
        xx[i]=w0[i];
        t0[i]=h0[i];
    }

    //myPrintf("H0\n\nh0[0]:%g\nh0[1]:%g\nh0[2]:%g\nh0[3]:%g\nh0[%d]:%g\n",t0[0],t0[1],t0[2],t0[3],M0-1,t0[M0-1]);

    //myPrintf("\nTRISOLVE REGULAR\n");    
    TrisolveREG(t0,h1,h1,xx,z,M0); //checked
    //myPrintf("Z\n\nz[0]:%g\nz[1]:%g\nz[2]:%g\nz[3]:%g\nz[%d]:%g\n",z[0],z[1],z[2],z[3],M0-1,z[M0-1]); 
    //myPrintf("Z\n\nz[0]:%g\nz[1]:%g\nz[2]:%g\nz[3]:%g\nz[%d]:%g\n",xx[0],xx[1],xx[2],xx[M0-2],M0-1,xx[M0-1]); 
    
    //myPrintf("\n\nSUBSELECTION STARTS\n");
    /* initialize if subselection */
    if (K<M0){
        //myPrintf("\nCOMPUTE HS SUBSELECTION\n");
        ComputeHs(I,M,K,h0,h1);//checked

        //myPrintf("h0s\n\nh0[0]:%g\nh0[1]:%g\nh0[2]:%g\nh0[3]:%g\nh0[%d]:%g\n",h0[0],h0[1],h0[2],h0[3],K-1,h0[K-1]); 
        //myPrintf("h1s\n\nh1[0]:%g\nh1[1]:%g\nh1[2]:%g\nh1[3]:%g\nh1[%d]:%g\n",h1[0],h1[1],h1[2],h1[3],K-2,h1[K-2]); 
        //z(I)
        //myPrintf("I[0]:%d\n",I[0]);
        
        for (i=0;i<K;i++)
            t0[i]=z[I[i]];        
        
        for (i=0;i<K;i++)
            xx[i]=0;
       
        //myPrintf("\nTRISYMGAXPY SUBSELECTION\n");
        TriSymGaxpy(h0,h1,t0,K,xx);//checked
        w0=xx; // error??? nono I think it's perfecly fine!!!
        //myPrintf("W0 TRISYMGAXPY\n\nw0[0]:%g\nw0[1]:%g\nw0[2]:%g\nw0[3]:%g\nd[%d]:%g\n",w0[0],w0[1],w0[2],w0[3],K-1,w0[K-1]); 
    }
    //myPrintf("\n\n*****************\nITERATION STARTS\n*****************\n\n");
    
/********************************/
/******* EM loop         ********/
/********************************/
    
    
    for (n=0;n<maxit;n++){
      
       for (i=0;i<K;i++){
           wpred[i]=w[i];
       }
       ComputeT(h0,h1,K,alpha,sigma2,t0,tl,tu);//checked  
//        if (n==0){
//        //myPrintf("COMPUTE T FIRST ITERATION:\n");
//        //myPrintf("T0\n\nt0[0]:%g\nt0[1]:%g\nt0[2]:%g\nt0[3]:%g\nt0[%d]:%g\n",t0[0],t0[1],t0[2],t0[3],K-1,t0[K-1]); 
//        //myPrintf("TL\n\ntl[0]:%g\ntl[1]:%g\ntl[2]:%g\ntl[3]:%g\ntl[%d]:%g\n",tl[0],tl[1],tl[2],tl[3],K-2,tl[K-2]); 
//        //myPrintf("TU\n\ntu[0]:%g\ntu[1]:%g\ntu[2]:%g\ntu[3]:%g\ntu[%d]:%g\n",tu[0],tu[1],tu[2],tu[3],K-2,tu[K-2]); 
//        }
       
       if (K==1){
           w[0]=w0[0]/t0[0];
           sigw[0]=1/t0[0]*h0[0]*sigma2;           
       }else{
           for(i=0;i<K;i++){
               if (i!=K-1){
                   AA[(i<<2)]=tu[i];
                   AA[(i<<2)+2]=tl[i];
               }
               AA[(i<<2)+1]=t0[i];
               AA[(i<<2)+3]=w0[i];
           }
           
           TriSolveINV(AA,4,K,w,d,e);//checked
//            if (n==0){
//            myPrintf("TRISOLVE INVERSE FIRST ITERATION:\n");
//            myPrintf("W TRISOLVE\n\nw[0]:%g\nw[1]:%g\nw[2]:%g\nw[3]:%g\nw[%d]:%g\n",w[0],w[1],w[2],w[3],K-1,w[K-1]); 
//            myPrintf("D TRISOLVE\n\nd[0]:%g\nd[1]:%g\nd[2]:%g\nd[3]:%g\nd[%d]:%g\n",d[0],d[1],d[2],d[3],K-2,d[K-2]); 
//            myPrintf("E TRISOLVE\n\ne[0]:%g\ne[1]:%g\ne[2]:%g\ne[3]:%g\ne[%d]:%g\n",e[0],e[1],e[2],e[3],K-2,e[K-2]); 
//            }
           
           
           
           tridiagofinverse(t0,tl,tl,t0,tu,K,d,e);// itl it0 itu on the sames without i to save memory can be recycled.
           // now itl itu it0 cointain the tridiagonal of the inverse....
//            if (n==0){
//            myPrintf("TRISOLVE FIRST ITERATION:\n");
//            myPrintf("IT0\n\nit0[0]:%g\nit0[1]:%g\nit0[2]:%g\nit0[3]:%g\nit0[%d]:%g\n",it0[0],it0[1],it0[2],it0[3],K-1,it0[K-1]); 
//            myPrintf("ITL\n\nitl[0]:%g\nitl[1]:%g\nitl[2]:%g\nitl[3]:%g\nitl[%d]:%g\n",itl[0],itl[1],itl[2],itl[3],K-2,itl[K-2]); 
//            myPrintf("ITU\n\nitu[0]:%g\nitu[1]:%g\nitu[2]:%g\nitu[3]:%g\nitu[%d]:%g\n",itu[0],itu[1],itu[2],itu[3],K-2,itu[K-2]); 
//            }
           
           
           DiagOfTriXTri(tl,t0,tu,h1,h0,h1,sigw,K);//checked
//            if (n==0){
//            myPrintf("TRIDIAGOFINVERSE FIRST ITERATION:\n");
//            myPrintf("DIAG\n\ndiag[0]:%g\ndiag[1]:%g\ndiag[2]:%g\ndiag[3]:%g\ndiag[%d]:%g\n",diag[0],diag[1],diag[2],diag[3],K-1,diag[K-1]); 
//            myPrintf("DIAGOFTRIXTRI ITERATION:%d\n",n);
//            }
           
           for(i=0;i<K;i++){
               sigw[i]=sigma2*sigw[i];
           }
           
       }
       for (i=0;i<K;i++){
           alpha[i]=(1+2*a)/(w[i]*w[i]+sigw[i]+2*b);
       }
       
//        if (n==0){
//            myPrintf("ALPHA & DIAGSIGMA FIRST ITERATION:\n");
//            myPrintf("\n\nMS:%d\n\n",K);
//            myPrintf("ALPHA\n\nalpha[0]:%g\nalpha[1]:%g\nalpha[2]:%g\nalpha[3]:%g\nalpha[%d]:%g\n",alpha[0],alpha[1],alpha[2],alpha[3],K-1,alpha[K-1]); 
//            myPrintf("DIAGSIGMA\n\nd[0]:%g\ndiagsigma[1]:%g\ndiagdigma[2]:%g\ndiagsigma[3]:%g\ndiagsigma[%d]:%g\n",sigw[0],sigw[1],sigw[2],sigw[3],K-1,sigw[K-1]); 
//        }

       //euclidean norm of wpred-w  CRITERIUM of CONVERGENCE
       /*
       mynorm=0;
       for (i=0;i<K;i++){
           mynorm=mynorm+(wpred[i]-w[i])*(wpred[i]-w[i]);
       }       
       mynorm=sqrt(mynorm);
       */

       //max diff of wpred-w  CRITERIUM of CONVERGENCE
       mynorm=0;
       for (i=0;i<K;i++){
           //myaux=abs(wpred[i]-w[i]);
           myaux=(wpred[i]-w[i]);
           if (myaux<0)myaux=-myaux;
           if (myaux>mynorm)mynorm=myaux;          
       }             
       
//        if (n==0){
//            myPrintf("EUCLIDEAN NORM FIRST ITERATION:\n");
//            myPrintf("\n\neuclidean_norm:%g\n\n",mynorm);
//        }
       if (mynorm<tol){
           if (debug>0){ 
			   myPrintf("#      SBL: Converged after %d iterations within tolerance %g, M=%d \n",n,tol,K);
           }
           break;
       }
       
       sizesel=findminus(alpha,K,maxalpha,sel);
//        if (n==0){
//            myPrintf("\n\nSIZESEL FIRST ITERATION:\n");
//            myPrintf("SEL\n\nsel[0]:%d\nsel[1]:%d\nsel[2]:%d\nsel[3]:%d\nsel[%d]:%d\n",sel[0],sel[1],sel[2],sel[3],K-1,sel[K-1]); 
//        }
       if (sizesel==0){
           K=0;
           w[0]=-1;
           I[0]=-1;
           sigw[0]=-1;
           alpha[0]=-1;
		   myPrintf("#      SBL: After %d iterations, No disconinuities found M=%d \n",n,K);        
           break;
       }
       if (sizesel<K){
           K=sizesel;           
           
           //I(sel)
           for (i=0;i<sizesel;i++)
               I[i]=I[sel[i]];
               //myPrintf("I[%d]=%d\n",i,I[i]);
                      
           ComputeHs(I,M,K,h0,h1);

//            if (k==0){
//            myPrintf("\n\nCOMPUTE Hs FIRST TIME REDUCTION:\n");
//            myPrintf("h0s sizesel<K\n\nh0[0]:%g\nh0[1]:%g\nh0[2]:%g\nh0[3]:%g\nh0[%d]:%g\n",h0[0],h0[1],h0[2],h0[3],K-1,h0[K-1]); 
//            myPrintf("h1s sizesel<K\n\nh1[0]:%g\nh1[1]:%g\nh1[2]:%g\nh1[3]:%g\nh1[%d]:%g\n",h1[0],h1[1],h1[2],h1[3],K-2,h1[K-2]); 
//            }
           
           //alpha(sel)
           for (i=0;i<K;i++)
               alpha[i]=alpha[sel[i]];
                            
          //myPrintf("I[0]:%d\n",I[0]);
           for (i=0;i<K;i++)
                t0[i]=z[I[i]];           
          
           for (i=0;i<K;i++)
               xx[i]=0;
           
           TriSymGaxpy(h0,h1,t0,K,xx);//checked until here
           w0=xx;
//            if (k==0){
//                myPrintf("\n\nTRISYMGAXPY FIRST TIME REDUCTION:\n");
//                myPrintf("W0 trisymgaxpy\n\nw0[0]:%g\nw0[1]:%g\nw0[2]:%g\nw0[3]:%g\nw0[%d]:%g\n",w0[0],w0[1],w0[2],w0[3],K-1,w0[K-1]); 
//            }
           
//            if (k==0){
//                myPrintf("W FIRSt TIME REDUCTION\n\nw[0]:%g\nw[1]:%g\nw[2]:%g\nw[3]:%g\nw[%d]:%g\n",w[0],w[1],w[2],w[3],K-1,w[K-1]); 
//                myPrintf("SEL FIRSt TIME REDUCTION\n\nsel[0]:%d\nsel[1]:%d\nsel[2]:%d\nsel[3]:%d\nsel[%d]:%d\n",sel[0],sel[1],sel[2],sel[3],K-1,sel[K-1]); 
//            }
           for (i=0;i<K;i++){
               w[i]=w[sel[i]];
               wpred[i]=wpred[sel[i]];               
           }
//             if (k==0){
//                myPrintf("W FIRST TIME REDUCTION AFTER W[SEL]\n\nw[0]:%g\nw[1]:%g\nw[2]:%g\nw[3]:%g\nw[%d]:%g\n",w[0],w[1],w[2],w[3],K-1,w[K-1]); 
//                k=1;
//             }
       
       }//End if of column elimination
    }//End loop MAXIT
    /*********************/
    
       if (debug>0){
           if (i>=maxit){              
			   myPrintf("#      SBL: Converged??? Stopped after %d iterations with change %g, M=%d\n",n,mynorm,K);
           }
       }
       
       #ifdef _DEBUG_SBL_
            myPrintf("_SBL_ Reffiting K=%d\n",K);
       #endif
       
       // Projection onto Breakpoint set... (REFITTING)
       if (K!=0){
           
           //myPrintf("I[0]:%d\n",I[0]);
          for (i=0;i<K;i++)
              t0[i]=z[I[i]]; 
          
          for (i=0;i<K;i++)
              sigw[i]=h0[i];

          for (i=0;i<K;i++)
              w[i]=0;          
          
          TriSymGaxpy(h0,h1,t0,K,w);
       }else{
           h0[0]=-1;
       }
    
    //Memory freeing....
    myFree(yy);
    myFree(h0);
    myFree(h1);
    myFree(xx);
    myFree(z);
    myFree(t0);
    myFree(tl);
    myFree(tu);
    myFree(wpred);
    myFree(d);
    myFree(e);
    myFree(sel);
    myFree(AA);
    
    
    *pK=K;
    return n;
}


/**************************************************************************************************************************/

    /*deppending on the Tau value this function will behave in two different 
    ways. First, if Tau is lower than zero (-1 will be our agreement) the 
    algorithm will remove the less-important break-point left. On the other
    hand, if it is greater than zero, it will work as a regular threshold 
    algorithm i.e.,it will remove all those break-point who don't reach the 
    selected value TAU.*/
int BEthresh(    
    double *Scores,
    int Nscores,
    double *wr,
    int *indsel,
    int *pointNumRem,
    double *pointTau
    )
{
    int M=Nscores+1; //genome length
    double vmin=-1e100; //Forces to enter at least once for tau<0
    int imin=-1;
    int inrem=-1;
    int i=0;
    int aux=0;
    int iC=0;
    int iR=0;
    int iL=0;
    double h0R=0;
    double h0L=0;
    double Tau=*pointTau;
    int NumRem=*pointNumRem;
    
    #ifdef _DEBUG_
        //checked-> inputs are OK
        int n=0;
        myPrintf("INPUT CHECKING:\n\nNscores=%d\n\nSCORES[indsel[0]]=%g\nSCORES[indsel[1]]=%g\nSCORES[indsel[2]]=%g\nSCORES[indsel[%d]]=%g\n",Nscores,Scores[indsel[0]],Scores[indsel[1]],Scores[indsel[2]],NumRem-1,Scores[indsel[NumRem-1]]);
        myPrintf("\nwr[indsel[0]]=%g\nwr[indsel[1]]=%g\nwr[indsel[2]]=%g\nwr[indsel[%d]]=%g\n\nindsel[0]=%d\nindsel[1]=%d\nindsel[2]=%d\nindsel[%d]=%d\n",wr[indsel[0]],wr[indsel[1]],wr[indsel[2]],NumRem,wr[indsel[NumRem-1]],indsel[0],indsel[1],indsel[2],NumRem-1,indsel[NumRem-1]);
        myPrintf("\nNumRem=%d\n\nTau=%g\n\nNscores=%d\n\n",NumRem,Tau,Nscores);
        
    #endif
    
    
    
    while (NumRem>0)
    {
        #ifdef _DEBUG_
            n++;
            myPrintf("ITERATION NUMBER:\n\nn=%d\n\n",n);
        #endif
        
        //search for the minimum score>0
        imin=0;
        vmin=Scores[indsel[imin]];
        for (i=1;i<NumRem;i++) //search of the lowest Score (less important)
        {
            aux=indsel[i];
            if (Scores[aux]<vmin)
            {
                imin=i;
                vmin=Scores[aux];
            }
        }
        inrem=indsel[imin]; // Returns the absolute position of the disc to remove.
        
        #ifdef _DEBUG_
            myPrintf("NEXT COMPONENT OUT:\n\nimin=%d\nvmin=%g\n\n",imin,vmin);
        #endif
        
        if ((vmin>Tau)&&(Tau>=0)) //if the lowest score is greater than TAU, we don't do anything
        {
            
            #ifdef _DEBUG_
                myPrintf("Threshold is not enough --> THERE IS NO MORE TO EXTRACT\n\n");
            #endif
        
            break;
        }
        if (imin==-1) 
        {
            
            #ifdef _DEBUG_
                myPrintf("No discontinuities left --> THERE IS NO MORE TO EXTRACT\n\n");
            #endif
        
            //there is nothing left
            NumRem=0;
            vmin=0;
            break;
        }
        /*If I'm here then imin is what I want to remove... because there is still
         *something to remove (non zero scores) and they are lower than my threshold*/
        
        iC=indsel[imin];
        
        #ifdef _DEBUG_
            myPrintf("IC Component out:\n\niC=%d\n\n",iC);
        #endif
        
        //only one left
        if (NumRem==1)
        {
            aux=indsel[imin];
            Scores[aux]=0;
            wr[aux]=0;
            NumRem=0;
            
            #ifdef _DEBUG_
                myPrintf("There is one left to extract:\n\nindsel[imin]=%d\n\n",aux);
            #endif
            
            break;
        }
        /*is it efficient to look for them, or will it be better to mantain an
        indsel like in the matlab version??? Mantain indsel
        Update left weight if it exist, and right weight if it exist.*/     
        
        if (imin==0) //if I should remove the leftmost w
        {
            #ifdef _DEBUG_
                myPrintf("LEFTMOST REMOVING:\n\n",imin,vmin);
            #endif
            
            iR=indsel[imin+1];//indsel[1]
            iL=0;
            //recompute weights
            wr[iR]=wr[iR]+(double)sqrt((double)(M-iR-1)/(M-iC-1)*(double)(iR+1)/(iC+1))*(double)(iC+1-iL)/(iR+1-iL)*(double)wr[iC];
            wr[iC]=0;
            
            #ifdef _DEBUG_
                myPrintf("wr[iR]=%g:\niR=%d\n",wr[iR],iR);
            #endif
            
            //I don't modify wr[iL] coz there is no left breakpoint
            //recompute Scores
            Scores[iC]=0;
            iC=indsel[imin+1];
            if (imin+1==NumRem-1)////imin+1 or imin+2???????
            {
                /*if there are just two discontinuitues so now we will 
                 work with the rightmost w (or Score)*/
                iR=M-1;
            }else
            {
                iR=indsel[imin+2];
            }
            h0R=(double)(M-iC-1)*(iC+1)/M*(iR+1-iL)/(iR-iC)/(iC+1-iL);
            if (wr[iC]>=0)
            {
                Scores[iC]=(double)(wr[iC]/sqrt(h0R));
            }else
            {
                Scores[iC]=(double)(-wr[iC]/sqrt(h0R));
            }
            
            #ifdef _DEBUG_
                myPrintf("h0R=%g\niC=%d\nScores[iC]=%g\n\n",h0R,iC,Scores[iC]);
            #endif
            
        }
        else if(imin==NumRem-1)
        {
            #ifdef _DEBUG_
                myPrintf("RIGHTMOST REMOVING:\n\n",imin,vmin);
            #endif
            
            iL=indsel[imin-1];
            iR=M;
            //recompute weights
            wr[iL]=wr[iL]+sqrt((double)(M-iL-1)/(M-iC-1)*(double)(iL+1)/(iC+1))*(double)(iR-iC-1)/(iR-iL-1)*wr[iC]; 
            wr[iC]=0;
            #ifdef _DEBUG_
                myPrintf("wr[iL]=%g:\niL=%d\n",wr[iL],iL);
            #endif

            //Scores
            Scores[iC]=0;
            /*new iC is gonna be placed in imin-1 so look if is the most
             left breakpoint*/
            if (imin-1==0)
            {
                iL=-1;
            }else
            {
                iL=indsel[imin-2];
            }
            iC=indsel[imin-1];
            h0L=(double)(M-iC-1)*(iC+1)/M*(iR-iL-1)/(iR-iC-1)/(iC-iL);
            if (wr[iC]>=0)
            {
                Scores[iC]=(double)(wr[iC]/sqrt(h0L));
            }else
            {
                Scores[iC]=(double)(-wr[iC]/sqrt(h0L));
            }
            #ifdef _DEBUG_
                myPrintf("h0L=%g\niC=%d\nScores[iC]=%g\n\n",h0L,iC,Scores[iC]);
            #endif
        }else
        {
            #ifdef _DEBUG_
                myPrintf("LEFT&RIGHT REMOVING:\n\n",imin,vmin);
            #endif
            //removing left and right
            iR=indsel[imin+1];
            iL=indsel[imin-1];
            wr[iR]=wr[iR]+sqrt((double)(M-iR-1)/(M-iC-1)*(double)(iR+1)/(iC+1))*(double)(iC-iL)/(iR-iL)*wr[iC];
            wr[iL]=wr[iL]+sqrt((double)(M-iL-1)/(M-iC-1)*(double)(iL+1)/(iC+1))*(double)(iR-iC)/(iR-iL)*wr[iC];
            wr[iC]=0;
            #ifdef _DEBUG_
                myPrintf("M-iR-1=%d:\nM-iC-1=%d\niR+1/iC+1=%g:\nsqrt(all)=%g\n",M-iR-1,M-iC-1,(iR+1)/(iC+1),sqrt(((M-iR-1)/(M-iC-1))*((iR+1)/(iC+1))));
                myPrintf("wr[iL]=%g:\niL=%d\nwr[iR]=%g:\niR=%d\n",wr[iL],iL,wr[iR],iR);
            #endif

            //Scores
            Scores[iC]=0;
            iL=indsel[imin-1];
            iC=indsel[imin+1];
            if (imin+1==NumRem-1)//iC is now the most right point
            {
                iR=M-1;
            }else
            {
                iR=indsel[imin+2];
            }
            h0R=(double)(M-iC-1)*(iC+1)/M*(iR-iL)/(iR-iC)/(iC-iL);
            if (wr[iC]>=0)
            {
                Scores[iC]=(double)(wr[iC]/sqrt(h0R));
            }else
            {
                Scores[iC]=(double)(-wr[iC]/sqrt(h0R));
            }
            #ifdef _DEBUG_
                myPrintf("h0R=%g\niCfirst=%d\nScores[iC]=%g\n\n",h0R,iC,Scores[iC]);
            #endif
            //scores
            if (imin-1==0)//iC is now the most left point
            {
                iL=-1;
            }
            else
            {
                iL=indsel[imin-2];
            }
            iC=indsel[imin-1];
            iR=indsel[imin+1];
            h0L=(double)(M-iC-1)*(iC+1)/M*(iR-iL)/(iR-iC)/(iC-iL);
            if (wr[iC]>=0)
            {
                Scores[iC]=(double)(wr[iC]/sqrt(h0L));
            }else
            {
                Scores[iC]=(double)(-wr[iC]/sqrt(h0L));
            }
            #ifdef _DEBUG_
                myPrintf("h0L=%g\niCsecond=%d\nScores[iC]=%g\n\n",h0L,iC,Scores[iC]);
            #endif
        }
        NumRem=NumRem-1;
        for (i=imin;i<NumRem;i++)
        {
            indsel[i]=indsel[i+1];
            
        }
        #ifdef _DEBUG_
                myPrintf("INDSEL RECOMPUTE:\n\nindsel[0]=%d\nindsel[1]=%d\nindsel[2]=%d\nindsel[%d]=%d\n\n",indsel[0],indsel[1],indsel[2],NumRem-1,indsel[NumRem-1]);
        #endif
        if (Tau<0)
        {   
            #ifdef _DEBUG_
                myPrintf("EXIT BECAUSE TAU IS LOWER THAN ZERO\n\n");
            #endif
            break;
        }
    }
    
    #ifdef _DEBUG_
        myPrintf("OUTPUT Parameters:\n\nvmin=%g\nNumRem=%d\nimin=%d",vmin,NumRem,imin);
    #endif
    
    *pointTau=vmin;
    *pointNumRem=NumRem;
    return(inrem);
}
        
/*************************************************************************/
/*
% CollapseAmp -- Collapses the amplitudes of the non significantly altered segments into a base level. 
% [OutAmp,State]=CollapseAmp(SegAmp,SegLen,BaseAmp,sigma2,T);
*/

void 
CollapseAmpTtest(//Uses a T test to decide which segments collapse to neutral
	double *SegAmp, //Segment Amplitudes (input output)
	const int *SegLen, //Segment Lengths
	int L, //Number of segments.
	double BaseAmp, //Reference amplitude to compare
	double sigma2,  //Reference noise 
	double T		//Critical value that decides when to colapse
	)
{
	int k;

//	if(L==1)
//	{
//		if( fabs(SegAmp[0]-BaseAmp)/sqrt(sigma2/(double)SegLen[0]) < T)
//			SegAmp[0]=BaseAmp;
//	}
	for(k=0;k<L;k++)
		if( fabs(SegAmp[k]-BaseAmp)/sqrt(sigma2/(double)SegLen[k]) < T)
		{
			//but do it only if one of the neigboring ones have been collapsed
			if((k>0)&&(fabs(SegAmp[k-1]-BaseAmp)/sqrt(sigma2/(double)SegLen[k-1]) < T))
				SegAmp[k]=BaseAmp;
			if((k<L-1)&&(fabs(SegAmp[k+1]-BaseAmp)/sqrt(sigma2/(double)SegLen[k+1]) < T))
				SegAmp[k]=BaseAmp;
			//or it is the initial segment or the final segment of the unit, since we assume that the unseen neighbors where in collapsed state
			if((k==0)||(k==L-1))
				SegAmp[k]=BaseAmp;
            //or it has larger size than the neighboring segments
			if((k>0)&&(SegLen[k]>SegLen[k-1]))
				SegAmp[k]=BaseAmp;
			if((k<L-1)&&(SegLen[k]>SegLen[k+1]))
				SegAmp[k]=BaseAmp;

		}

/*	
	//Classify amplitudes...
	for(k=0;k++;k<K)
	{
		if(SegAmp[k]>BaseAmp)
			SegAmp=+1; //Gain
		else if(SegAmp[k]<BaseAmp)
			SegAmp=-1;
		else
			SegAmp=0;		
	}
*/
}
/*************************************************************************/
/*
% ClassifySegments -- Classifies/Collapses reconstructed segments into Altered/Gain/Loss
% Gain uses positive numbers log2(3)-log2(2)
*/

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
	)
{
	int k;
	double c;
	double aux;
	
	for(k=0;k<=K;k++)
	{
		if( fabs(SegAmp[k]-BaseAmp)/sqrt(sigma2/(double)SegLen[k]) < T)
		{
			SegState[k]=ploidy;
		}
		else if((SegAmp[k]-BaseAmp)>0)
		{
			for(c=ploidy;c<100;c++)
			{
				aux=(SegAmp[k]-BaseAmp-log2(c)+log2(ploidy))/sqrt(sigma2/(double)SegLen[k]);
				if((SegAmp[k]-BaseAmp-log2(c)+log2(ploidy))/sqrt(sigma2/(double)SegLen[k]) < T )
				{
//					c=c-1;
					break;
				}
			}
			SegState[k]=c;
		}
		else if((SegAmp[k]-BaseAmp)<0)
		{
			for(c=ploidy;c>0;c--)
			{
				if((SegAmp[k]-BaseAmp-log2(c)+log2(ploidy))/sqrt(sigma2/(double)SegLen[k]) > -T )
				{
					break;
				}
			}
			SegState[k]=c;
		}
	}
}

/*	
	//Classify amplitudes...
	for(k=0;k++;k<K)
	{
		if(SegAmp[k]>BaseAmp)
			SegAmp=+1; //Gain
		else if(SegAmp[k]<BaseAmp)
			SegAmp=-1;
		else
			SegAmp=0;		
	}
*/


double // Returns BaseAmp corresponding to the base level.
CompBaseAmpMedianMethod( //Computes the median recontruction level, as baseline level.
	const int *SegLen,    //Lengths corresponding to the amplitudes
	const double *SegAmp, //Amplitudes !!! assumed already ordered...
	int K
	)
{	
	int M,k,RunLen;
	double BaseAmp=0;
	
	//If they need to be sorted use the following...
	double *D;
	int *I;
	D=myCalloc(K+1,sizeof(double));
	I=myCalloc(K+1,sizeof(double));
	for(k=0;k<K+1;k++)
	{
		D[k]=SegAmp[k];
		I[k]=SegLen[k];
	}
	doubleBubbleSort(D,I,K+1); //I need indexes of the sort
	SegAmp=D;
	SegLen=I;

	M=0;
	for(k=0;k<=K;k++)
		M=M+SegLen[k];
#ifdef _DebugCompBaseAmpMedianMethod_
	myPrintf("_DebugCompBaseAmpMedianMethod_: M%d K%d M/2%d\n",M,K,M/2);
#endif
	
	RunLen=0;
	k=0;
	while(RunLen<M/2)RunLen=RunLen+SegLen[k++];
#ifdef _DebugCompBaseAmpMedianMethod_
	myPrintf("_DebugCompBaseAmpMedianMethod_: k%d RunLen=%d SegAmp[k-1]%g SegAmp[k]%g\n",k,RunLen,SegAmp[k-1],SegAmp[k]);
#endif

   BaseAmp=SegAmp[k-1];

#ifdef _DebugCompBaseAmpMedianMethod_
   printf("_DebugCompBaseAmpMedianMethod_: BaseAmp=%g\n",BaseAmp);
#endif

   return BaseAmp;
}






/**************************************************************************************************************************/
int SBLandBE( //Returns breakpoint list lenght.
    double *y, 
    int M,  //length of the noisy signal tn
    double *Psigma2, //If sigma2 < 0, compute sigma2
    double a,      // SBL parameter
    double T,      // Threshold to prune
    int MinSegLen, //Minimum length of the segment.
    int **pIext,   //Returns breakpoint positions
    double **pWext //Returns breakpoint weights.
    //int *pK    
)
{
    int i,K;
    double myaux;
	double sigma2;
    double ymean;
    double tol,maxalpha,b;
    double *alpha,*Wext,*tn,*aux;
    int *Iext;
    int maxit;
    int NumEMsteps;

	sigma2=*Psigma2;
		
	tn=myCalloc(M,sizeof(double));
	for(i=0;i<M;i++)
		tn[i]=y[i];
    
    //If sigma2 < 0, compute sigma2
    if(sigma2<0){
        sigma2=0;
        for(i=1;i<M;i++){
            myaux=tn[i]-tn[i-1];
            sigma2+=(0.5*myaux*myaux);
        }
        sigma2=sigma2/(M-1);
    }      
    
    //Mean removal
    ymean=0;
    for(i=0;i<M;++i)
        ymean+=tn[i];
    ymean=ymean/M;
    for(i=0;i<M;++i)
        tn[i]=tn[i]-ymean;
    
    // SBL optimization parameters
    tol=1E-10;     //1E-10 or 1E-8 seems to work well for this parameter. -- => ++ conv time
    maxalpha=1E8;  //1E8 better than 1E10 seems to work well for this parameter. -- => -- conv time
    maxit=10000;   //Maximum number of iterations to reach convergence...
    b=1E-20;       //
    
    //Call to SBL 
    #ifdef _DEBUG_SBLBE_
        myPrintf("_SBLBE_ Memory initialization\n");
    #endif

    K=M-1;
    Iext=myCalloc(M,sizeof(int)); //R
    Wext=myCalloc(M,sizeof(double)); //R
    alpha=myCalloc(M,sizeof(double));   
    aux=myCalloc(M,sizeof(double));

    
    #ifdef _DEBUG_SBLBE_
        myPrintf("_SBLBE_ Breakpoint Initialization\n");
    #endif
    //Initialize breakpoints
    for (i=0;i<M;i++)
        alpha[i]=0.0;
    for (i=0;i<M;i++)
        Wext[i]=0.0;
    for (i=0;i<M;i++)
        Iext[i]=i;
    
    #ifdef _DEBUG_SBLBE_
        myPrintf("_SBLBE_ SBL starts\n");
    #endif
    NumEMsteps=SBL(tn,Iext,alpha,Wext+1,aux,M,&K,sigma2,a,b,maxalpha,maxit,tol,1);

    //Freeing memory
    myFree(alpha);
	myFree(tn);
	myFree(aux);


	//Convert Iext and Wext to the extended notation.
	Iext=realloc(Iext,(K+2)*sizeof(int));
	for(i=(K+1);i>0;i--)
		Iext[i]=Iext[i-1]+1;
	Iext[0]=0;
	Iext[K+1]=M;

	Wext=realloc(Wext,(K+1)*sizeof(double));
	//for(i=K;i>0;i++)
	//	Wext[i]=Wext[i-1];
	Wext[0]=ymean;
    
    #ifdef _DEBUG_SBLBE_
	myPrintf("_SBLBE_ Backward Elimination T=%g MinLen=%d\n",T,MinSegLen);
    #endif
	BEwTandMinLen(Wext,Iext,&K,sigma2,T,MinSegLen); 

	Iext=realloc(Iext,(K+2)*sizeof(int));
	Wext=realloc(Wext,(K+1)*sizeof(double));

	#ifdef _DEBUG_SBLBE_
        myPrintf("_SBLBE_ After BE K=%d \n",K);
    #endif
    
    
    *pIext=Iext;
    *pWext=Wext;
    *Psigma2=sigma2;

    return K;    
}


/**************************************************************************************************************************/
int BEwTandMinLen( //Returns breakpoint list lenght. with T and MinSegLen
    double *Wext,  //IO Breakpoint weights extended notation...
	int *Iext,     //IO Breakpoint positions in extended notation...
    int *pK,       //IO Number breakpoint positions remaining.
	double sigma2, //IP If sigma2 
    double T,      //IP  Threshold to prune,  T=T*sqrt(sigma2);
    int MinSegLen  //IP Minimum length of the segment.
)
{
    int i,K,imin,M; 
//    double myaux;
    double vmin;
    double *tscore; //Statistical Scores,  
    int *L; //Vector with the smallest of the two neighboring segments of each breakpoint.
    int smallinside; //Variable that indicates that there are still small segments to eliminate.

	K=*pK;       //Number of breakpoints
	M=Iext[K+1]; //Total length		
    T=T*sqrt(sigma2); //Adjusting T to the noise power
    
    if(MinSegLen>0)
        smallinside=1; 
    else
        smallinside=0;
        
    L=myCalloc(K+1,sizeof(int)); //Bring from outside?
    tscore=myCalloc(K+1,sizeof(double)); //
    
	//Computing scores
    #ifdef _DEBUG_BEwTandMinLen_
        myPrintf("_BEwTandMinLen_ Computing scores\n");
    #endif    
    for (i=0;i<K+1;i++)
        tscore[i]=0;    

    ComputeTScores(Wext,Iext,tscore,K,1,K);

	#ifdef _DEBUG_BEwTandMinLen_
       	myPrintf("\n_BEwTandMinLen_:K%d w[0]=%g,w[1]=%g,w[2]=%g,w[K-1]=%g,w[K]=%g\n",K,Wext[0],Wext[1],Wext[2],Wext[K-1],Wext[K]); 
       	myPrintf("\n_BEwTandMinLen_:K%d t[0]=%g,t[1]=%g,t[2]=%g,t[K-1]=%g,t[K]=%g\n",K,tscore[0],tscore[1],tscore[2],tscore[K-1],tscore[K]); 
    #endif
    #ifdef _DEBUG_BEwTandMinLen_
        myPrintf("_BEwTandMinLen_ Backward Elimination sigma2=%g T=%g MinSegLen %d\n",sigma2,T,MinSegLen);
    #endif

	BEwTscore(Wext,Iext,tscore,&K,T);

	#ifdef _DEBUG_BEwTandMinLen_
       	myPrintf("\n_BEwTandMinLen_:K%d w[0]=%g,w[1]=%g,w[2]=%g,w[K-1]=%g,w[K]=%g\n",K,Wext[0],Wext[1],Wext[2],Wext[K-1],Wext[K]); 
       	myPrintf("\n_BEwTandMinLen_:K%d t[0]=%g,t[1]=%g,t[2]=%g,t[K-1]=%g,t[K]=%g\n",K,tscore[0],tscore[1],tscore[2],tscore[K-1],tscore[K]); 
    #endif


    #ifdef _DEBUG_BEwTandMinLen_
        myPrintf("_BEwTandMinLen_ Small segment prunning ML=%d SI=%d \n",MinSegLen,smallinside);
    #endif    
    while(smallinside)
    {
        //Compute segment lengths of flanking breakpoints
        for(i=1;i<K+1;i++)
            L[i]=min(Iext[i]-Iext[i-1],Iext[i+1]-Iext[i]);
        
        //Find smallest segment with lowest score to remove
        imin=-1;
        vmin=1E100;
        for(i=1;i<K+1;i++)
            if((tscore[i]<vmin)&&(L[i]<MinSegLen)){
                vmin=tscore[i];
                imin=i;
            }

        if(imin<0)
            smallinside=0;
        else
        {
			#ifdef _DEBUG_BEwTandMinLen_
				myPrintf("_BEwTandMinLen_ Removing %d %d %g Numrem(old)=%d",Iext[imin],imin,vmin,K);
			#endif
			tscore[imin]=-1;
			BEwTscore(Wext,Iext,tscore,&K,T);

			#ifdef _DEBUG_BEwTandMinLen_
				myPrintf("_BEwTandMinLen_ K = %d \n",K);
			#endif
        }
        if(K<1)smallinside=0;
    }    
    #ifdef _DEBUG_BEwTandMinLen_
        myPrintf("T=%g K=%d\n",T,K);
        myPrintf("_BEwTandMinLen_ Small segment prunning ends imin=%d vmin=%g\n",imin,vmin);
    #endif    
         
  // for (i=1;i<K;i++)
  //     Wext[i]=0.0;
  // for (i=1;i<K;i++)
  //     Wext[i]=w2[Iext[i]];

    #ifdef _DEBUG_BEwTandMinLen_
        myPrintf("_BEwTandMinLen_ After BE K=%d \n",K);
    #endif
    
    //Freeing memory
    myFree(tscore);
    myFree(L);

	*pK=K;
    return K;    
}

/******************************************************/
    //BEwTscore(Iext,Wext,h0,h1,tscore,&K,T);  //Need to update BEthres to operate on the Iext Wext notation...
int BEwTscore(
    double *Wext,  //IO Breakpoint weights extended notation...
	int *Iext,     //IO Breakpoint positions in extended notation...
	double *tscore,
    int *pK,       //IO Number breakpoint positions remaining.
    double T      //IP  Threshold to prune
	)
{
	int i,K,jmin,M; 
    double vmin;

	K=*pK;       //Number of breakpoints
	M=Iext[K+1]; //Total length		

#ifdef _DebugBEwTscore_
	myPrintf("_DebugBEwTscore_ BE starts K=%d M=%d T=%g\n",K,M,T);
#endif    


	vmin=0;
	while((vmin<T)&&(K>0))
	{
		//Find breakpoint with lowest score to remove
        jmin=-1;
        vmin=1E100;
        for(i=1;i<K+1;i++)
            if(tscore[i]<vmin){
                vmin=tscore[i];
                jmin=i;
            }
	    #ifdef _DebugBEwTscore_
			myPrintf("_DebugBEwTscore_ Smallest breakpoint imin=%d vmin=%g T=%g\n",jmin,vmin,T);
		#endif    
		if(vmin<T) //Remove breakpoint at imin
		{
#ifdef _DebugBEwTscore_
			myPrintf("_DebugBEwTscore_ Removing imin=%d vmin=%g T=%g\n",jmin,vmin,T);
#endif    

#ifdef _DebugBEwTscore_
			myPrintf("_DebugBEwTscore_ Removing imin=%d vmin=%g T=%g\n",jmin,vmin,T);
#endif    

#ifdef _DebugBEwTscore_
			myPrintf("_DebugBEwTscore_ %d, W[jmin-2]=%g W[jmin-1]=%g W[jmin]=%g W[jmin+1]=%g W[jmin+2]=%g T=%g\n",jmin,Wext[jmin-2],Wext[jmin-1],Wext[jmin],Wext[jmin+1],Wext[jmin+2],T);
			myPrintf("_DebugBEwTscore_ %d, I[jmin-2]=%d I[jmin-1]=%d I[jmin]=%d I[jmin+1]=%d I[jmin+2]=%d T=%g\n",jmin,Iext[jmin-2],Iext[jmin-1],Iext[jmin],Iext[jmin+1],Iext[jmin+2],T);
			myPrintf("_DebugBEwTscore_ %d, t[jmin-2]=%g t[jmin-1]=%g t[jmin]=%g I[jmin+1]=%g t[jmin+2]=%g T=%g\n",jmin,tscore[jmin-2],tscore[jmin-1],tscore[jmin],tscore[jmin+1],tscore[jmin+2],T);
#endif    
			for(i=jmin;i<K;i++)
				tscore[i]=tscore[i+1];
	   		K=RemoveBreakpoint(Wext,Iext,K,jmin);
	        ComputeTScores(Wext,Iext,tscore,K,jmin-1,jmin);

#ifdef _DebugBEwTscore_
			myPrintf("_DebugBEwTscore_ %d, W[jmin-2]=%g W[jmin-1]=%g W[jmin]=%g W[jmin+1]=%g W[jmin+2]=%g T=%g\n",jmin,Wext[jmin-2],Wext[jmin-1],Wext[jmin],Wext[jmin+1],Wext[jmin+2],T);
			myPrintf("_DebugBEwTscore_ %d, I[jmin-2]=%d I[jmin-1]=%d I[jmin]=%d I[jmin+1]=%d I[jmin+2]=%d T=%g\n",jmin,Iext[jmin-2],Iext[jmin-1],Iext[jmin],Iext[jmin+1],Iext[jmin+2],T);
			myPrintf("_DebugBEwTscore_ %d, t[jmin-2]=%g t[jmin-1]=%g t[jmin]=%g I[jmin+1]=%g t[jmin+2]=%g T=%g\n",jmin,tscore[jmin-2],tscore[jmin-1],tscore[jmin],tscore[jmin+1],tscore[jmin+2],T);
#endif    

		}
	}
	#ifdef _DebugBEwTscore_
			myPrintf("_DebugBEwTscore_ BE ends imin=%d vmin=%g T=%g K=%d\n",jmin,vmin,T,K);
	#endif  

	*pK=K;
	return K;
}

int
RemoveBreakpoint(
	double *Wext,
	int *Iext,
	int K,
	int jrem
	)
{
	int j;
	double iC,iL,iR,M;

	M=(double)Iext[K+1];
	iL=(double)Iext[jrem-1];
	iC=(double)Iext[jrem];
	iR=(double)Iext[jrem+1];

	//Change coefficients
	if(jrem>1)
		Wext[jrem-1] = Wext[jrem-1] + sqrt((M-iL)/(M-iC)*iL/iC)*(iR-iC)/(iR-iL) * Wext[jrem]; 
	if(jrem<K)
		Wext[jrem+1] = Wext[jrem+1] + sqrt((M-iR)/(M-iC)*iR/iC)*(iC-iL)/(iR-iL) * Wext[jrem]; 
	Wext[jrem]=0; //

	//Shorten list
	for(j=jrem;j<K;j++)
		Wext[j]=Wext[j+1];
	for(j=jrem;j<K+1;j++)
		Iext[j]=Iext[j+1];
	return K-1;
}

void
ComputeTScores(
	const double *Wext,
	const int *Iext,
	double *Scores,
	int K,
	int start,
	int end
	)
{
	int j;
	double h0,M;
	
	M=(double)Iext[K+1];

	for(j=start;j<=end;j++)
	{
		h0=(double)(M-Iext[j]) * (double)Iext[j] / M * (double)(Iext[j+1]-Iext[j-1]) / (double)(Iext[j+1]-Iext[j]) / (double)(Iext[j]-Iext[j-1]);
        Scores[j]=fabs(Wext[j])/sqrt(h0);
	}
}


void Project (
    double *y,
    int M,
    int *I,
    int L,
    double *xI,
    double *wI
    )
{
    // Intern variables declaration
    double aux_double=0;
    double ymean=0;
    int i=0;
    double *h0;
    double *h1;
    double *w0;
    double *z;
    double *temp;
    double *wr;
    double *aux_vec;
    
    // Variables inizialitation
    h0=myCalloc(M-1,sizeof(double));
    h1=myCalloc(M-2,sizeof(double));
    z=myCalloc(M-1,sizeof(double));
    temp=myCalloc(L,sizeof(double));
    wr=myCalloc(M,sizeof(double));
    aux_vec=myCalloc(M,sizeof(double));
    
//     h0=mxGetPr(mxCreateDoubleMatrix(1,M-1,mxREAL)); 
//     h1=mxGetPr(mxCreateDoubleMatrix(1,M-2,mxREAL));
//     z=mxGetPr(mxCreateDoubleMatrix(1,M-1,mxREAL));
//     temp=mxGetPr(mxCreateDoubleMatrix(1,L,mxREAL));
//     wr=mxGetPr(mxCreateDoubleMatrix(1,M,mxREAL));    
//     aux_vec=mxGetPr(mxCreateDoubleMatrix(1,M,mxREAL));
    
    // Remove y's mean
    aux_double=0;
    
    for (i=0;i<M;i++)
    {
        aux_double=aux_double+y[i];
    }
    ymean=(double)(aux_double/M);
    for (i=0;i<M;i++)
    {
        y[i]=y[i]-ymean; // I have y mean in aux_double
    }
    #ifdef _DEBUG_
        //checked-> inputs are OK
        myPrintf("MEAN CHECKING:\n\nymean=%g\n",ymean);
        myPrintf("\ny[0]=%g\ny[1]=%g\ny[2]=%g\ny[3]=%g\ny[%d]=%g\nM=%d\n",y[0],y[1],y[2],y[3],M-1,y[M-1],M);
        
    #endif
    
    //start reconstruction
    if (L>0)
    {
        //I sort -> maybe desorded
        BubbleSort(I,L);
        #ifdef _DEBUG_
        //checked-> inputs are OK
        myPrintf("BubbleSort CHECKING:\n");
        myPrintf("\nI[0]=%d\nI[1]=%d\nI[2]=%d\nI[3]=%d\nI[%d]=%d\nL=%d\n",I[0],I[1],I[2],I[3],L-1,I[L-1],L);
        
        #endif
        ComputeH(h0,h1,M);
        #ifdef _DEBUG_
        //checked-> inputs are OK
        myPrintf("COMPUTE H CHECKING:\n\nh0[0]=%g\nh0[1]=%g\nh0[2]=%g\nh0[3]=%g\nh0[%d]=%g\nsizeh0=%d\n",h0[0],h0[1],h0[2],h0[3],M-2,h0[M-2],M-1);
        myPrintf("\nh1[0]=%g\nh1[1]=%g\nh1[2]=%g\nh1[3]=%g\nh1[%d]=%g\nsizeh0-1=%d\n",h1[0],h1[1],h1[2],h1[3],M-3,h1[M-3],M-2);
        
        #endif
        ComputeFdualXb(M,y);
        w0=y;  //careful I just erased y's value
        #ifdef _DEBUG_
        //checked-> inputs are OK
        myPrintf("Compute F dual CHECKING:\n\nw0[0]=%g\nw0[1]=%g\nw0[2]=%g\nw0[3]=%g\nw0[%d]=%g\nM=%d\n",w0[0],w0[1],w0[2],w0[3],M-1,w0[M-2],M-1);
        
        #endif
        TrisolveREG(h0,h1,h1,w0,z,M-1);
        #ifdef _DEBUG_
        //checked-> inputs are OK
        myPrintf("Trisolve Checking:\n\nz[0]=%g\nz[1]=%g\nz[2]=%g\nz[3]=%g\nz[%d]=%g\nM=%d\n",z[0],z[1],z[2],z[3],M-1,z[M-2],M-1);
        
        #endif
        ComputeHs(I,M,L,h0,h1);
        #ifdef _DEBUG_
        //checked-> inputs are OK
        myPrintf("Compute Hs2 CHECKING:\n\nh0[0]=%g\nh0[1]=%g\nh0[2]=%g\nh0[3]=%g\nh0[%d]=%g\nL=%d\n",h0[0],h0[1],h0[2],h0[3],L-1,h0[L-1],L);
        myPrintf("\nh1[0]=%g\nh1[1]=%g\nh1[2]=%g\nh1[3]=%g\nh1[%d]=%g\nL-1=%d\n",h1[0],h1[1],h1[2],h1[3],L-2,h1[L-2],L-1);
        
        #endif
        for (i=0;i<L;i++){
            temp[i]=z[I[i]];
        }
        
        for (i=0;i<L;i++){
        y[i]=0;
        }
        #ifdef _DEBUG_
        //checked-> inputs are OK
        myPrintf("Z(I) CHECKING:\n\ntemp[0]=%g\ntemp[1]=%g\ntemp[2]=%g\ntemp[3]=%g\ntemp[%d]=%g\nL=%d\n",temp[0],temp[1],temp[2],temp[3],L-1,temp[L-1],L);
        
        #endif
        TriSymGaxpy(h0,h1,temp,L,wI);
        #ifdef _DEBUG_
        //checked-> inputs are OK
        myPrintf("TriSymGaxpy CHECKING:\n\nwI[0]=%g\nwI[1]=%g\nwI[2]=%g\nwI[3]=%g\nwI[%d]=%g\nL=%d\n",wI[0],wI[1],wI[2],wI[3],L-1,wI[L-1],L);
        
        #endif
        for (i=0;i<L;i++){
            wr[I[i]]=wI[i];
            #ifdef _DEBUG_
            //checked-> inputs are OK
            myPrintf("wr CHECKING:\n\nwr[I[%d]]=%g\nI[i]=%d\n",i,wr[I[i]],I[i]);
        
            #endif
        }
        #ifdef _DEBUG_
        //checked-> inputs are OK
        myPrintf("wr CHECKING:\n\nwr[I[0]]=%g\nwr[I[1]]=%g\nwr[%d]=%g\nM=%d\n",wr[I[0]],wr[I[1]],I[L-1],wr[I[L-1]],M);        
        #endif
        
        reconstruct(wr,M,aux_vec);
        #ifdef _DEBUG_
        //checked-> inputs are OK
        myPrintf("RECONSTRUCT CHECKING:\n\nwr[0]=%g\nwr[1]=%g\nwr[2]=%g\nwr[3]=%g\nwr[%d]=%g\nM=%d\naux_double=%g\nymean=%g\n",wr[0],wr[1],wr[2],wr[3],M-1,wr[M-1],M,aux_double,ymean);
        
        #endif
        for (i=0;i<M;i++)
        {
            xI[i]=ymean+wr[i];
        }
    }
    else
    {
        L=0;
        for (i=0;i<M;i++)
        {
            xI[i]=ymean;
        }
    }
    
    myFree(h0);
    myFree(h1);
    myFree(z);
    myFree(temp);
    myFree(wr);
    myFree(aux_vec);
    
}

void ProjectCoeff ( //IextYobs2Wext
    double *y,
    int M,
    int *Iext,
    int K,
    double *Wext
	)
{
    // Intern variables declaration
    double ymean=0;
    int i=0;
    double *h0;
    double *h1;
    double *z;
    
    // Variables inizialitation
    h0=myCalloc(K,sizeof(double));
    h1=myCalloc(K-1,sizeof(double));
    z=myCalloc(M-1,sizeof(double));
    
    
	ymean=0;    
    for (i=0;i<M;i++)
        ymean=ymean+y[i];
    ymean=ymean/M;

	// Remove y's mean, Not necessary
//    for (i=0;i<M;i++)
//        y[i]=y[i]-ymean;	

    if (K>0)
    {
        //I sort -> assumed that I is already ordered
        // BubbleSort(Iext,K+2);
	
		CompZ(y,z,M);
//		myPrintf("\n CHECKING: ymean=%g,M=%d,K=%d\n",ymean,M,K);//
//		myPrintf("\n z CHECKING: z[0]=%g,z[1]=%g,z[2]=%g,z[M-3]=%g,z[M-2]=%g\n",z[0],z[1],z[2],z[M-3],z[M-2]);

		for (i=1;i<=K;i++)
            z[i-1]=z[Iext[i]-1];
//    	myPrintf("\n z CHECKING: z[0]=%g,z[1]=%g,z[2]=%g,z[K-2]=%g,z[K-1]=%g\n",z[0],z[1],z[2],z[K-2],z[K-1]);   

		ComputeHsIext(Iext,K,h0,h1);
        //myPrintf("\n h0 CHECKING: h0[0]=%g,h0[1]=%g,h0[2]=%g,h0[K-2]=%g,h0[K-1]=%g,h0[K]=%g\n",h0[0],h0[1],h0[2],h0[K-2],h0[K-1]);   
        //myPrintf("\n h1 CHECKING: h1[0]=%g,h1[1]=%g,h1[2]=%g,h1[K-2]=%g,h1[K-1]=%g\n",h1[0],h1[1],h1[2],h1[K-2]);   
		
		for (i=0;i<K+1;i++)
			Wext[i]=0.0;

        TriSymGaxpy(h0,h1,z,K,Wext+1);
 //   	myPrintf("\n w CHECKING: w[0]=%g,w[1]=%g,w[2]=%g,w[K-1]=%g,w[K]=%g\n",Wext[0],Wext[1],Wext[2],Wext[K-1],Wext[K]);   
    }
    Wext[0]=ymean;
    
    myFree(h0);
    myFree(h1);
    myFree(z);
   
}


