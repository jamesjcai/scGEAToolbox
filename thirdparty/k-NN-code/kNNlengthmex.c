/**************************************************************************/
/*  kNNlengthmex.c: C-mex code, written for Matlab to build k-NN graph    */
/*  Yields cumulative graph length of power weighted edges                */
/*  Code Copyright (c): Huzefa Neemuchwala, hneemuch@umich.edu            */
/*  (All Rights Reserved)                                                 */
/*  while a PhD student, with advising from  Alfred Hero, EECS            */
/*  University of Michigan, Ann Arbor                                     */
/*  Departments of Biomedical Engineering, EECS and Radiology             */
/*  Part of this code has been obtained from the web,                     */
/*  Thank you, Anonymous                                                  */
/*  This is the final version.                                            */
/*  Usage:                                                                */
/*   kNNlength = kNNlengthmex(rrw, N, dim, kneightbors, gamma)            */
/*   Variables and input formats are best explained with accompanying     */
/*  Matlab Code                                                           */
/*  Contact author (hneemuch@umich.edu) for questions,                    */
/*  comments, bugs or concerns                                            */
/**************************************************************************/

/**************************************************************************/
/* Program to perform orthogonal range searches and nearest neighbor      */
/* querys in a more sophisticated k-d tree.  In this implementation the,  */
/* nodes on any given level of the tree do not have the same              */
/* discriminating dimension as the discrimiator is chosen based on the    */
/* dimension with   the "maxspead."                                       */
/*                                                                        */
/* References:  J.H. Friedman, J.L. Bentley, R.A. Finkel  "An Algorithm   */
/* for Finding Best Matches in Logarithmic Expected Time."                */
/* ACM Transactions on Mathematical Software, Vol 3 No. 3 Sept. 1977      */
/* pp. 209-226.                                                           */
/**************************************************************************/

#include <stdio.h>
#include <math.h>
#include "mex.h"
#include "optkd.h"
#include "nn.h"


/* Used to create a new tree in the k-d tree */
#define TESTTREE(PP)  ((PP) = (optkdNode *)mxMalloc(sizeof(optkdNode)))
#define NEWTREE(PP)  if (TESTTREE(PP)==NULL) \
                         {printf("memory error\n");return;}
#define min(X, Y)  ((X) < (Y) ? (X) : (Y))
static void principal (double, double, double *, double, double * , /*double *,*/ double);
void Selection(double **, int, int, int, int);
int findmaxspread(int, int, int, double **);
optkdNode *BuildkdTree(double**, int, int, int);
optkdNode *BuildOptTree(double **, int, int);
void rnnEuclidean(optkdNode *, double *, double **, int, int, int*);
int *kdOptNNQuery(double **, int, double *, int, int, optkdNode *, int);
void KillOptTree(optkdNode *);
void PQupheap(double *, int *, int);
void PQInsert(double, int, double *, int *);
void PQdownheap(double *, int *, int, int);
void PQreplace(double, double *, int *, int);
/*double EuclidDist2(double **, int, double *, int);*/

int *perm;  /* permutation array */
/*int *optfound; /*This one will be killed by KillOptTree()*/
static int Count=0;
double *nndist;

/*double (*Distance)();
extern double fabs();*/



void mexFunction(int nlhs, mxArray *plhs[], int nrhs,
                 const mxArray *prhs[])
{

  double numpoints, dimension, kneighbors;
  double* prLxo;
  double *rrw;
  /*double *Graph;*/
  double gamma;
  

  /* Check for proper number of arguments. */
  if (nrhs != 5) {
    mexErrMsgTxt("Five inputs required : kNNgraphmex (rrw, N, dimension, k_as_in_kNN, gamma)");
  } else if (nlhs > 1) {
    mexErrMsgTxt("Only one output argument");
  }

  /* Assign pointers to each input */
  /* Use mex syntax */
  /* matrices for each time point */
  rrw = (double *)mxGetPr(prhs[0]);

  /* parameters for analysis */
  numpoints = *mxGetPr(prhs[1]);
  dimension = *mxGetPr(prhs[2]);
  kneighbors= *mxGetPr(prhs[3]);
  kneighbors = kneighbors + 1;  /* Because first neighbor is the point itself */
  gamma = *mxGetPr(prhs[4]);
  if (kneighbors > numpoints) 
    mexErrMsgTxt("Input Error 1: Number of NN cant be larger than total number of points ");
  /* Create matrix for the return argument. */
  /* this is basically useless */
  plhs[0] = mxCreateDoubleMatrix(1,1,mxREAL);
  /*plhs[1] = mxCreateDoubleMatrix( numpoints*kneighbors, 1, mxREAL);*/
  
  /* Assign pointer to the output */
  prLxo = (double *)mxGetPr(plhs[0]); 
  /*Graph = (double *)mxGetPr(plhs[1]);*/
  /*printf("%g %g %g\n",numpoints, dimension,kneighbors);*/
  principal (numpoints, dimension, rrw, kneighbors, prLxo , /*Graph,*/ gamma);
  /*free(nndist);
  free(perm);*/
}


/************Need***************************************************************/

void Selection(a, l,N, k,discrim)

/***************************************************************************/
/* Makes the perm partition the array Values along the element k.          */
/* Adapted from Sedgewick's Algorithms in C (p. 128)                       */
/***************************************************************************/

double **a;
int l,N,k,discrim;

{
double v;
int t,i,j,r;

r=N;

while(r>l) {
  v=a[perm[r]][discrim]; i=l-1; j=r;
  for (;;) {
    while (a[perm[++i]][discrim] < v);
    while (a[perm[--j]][discrim] > v && j>l); 
    if (i >= j) break;
    t=perm[i]; perm[i] = perm[j]; perm[j]=t;
  }
  t=perm[i]; perm[i] = perm[r]; perm[r]=t;
  if (i>=k) r=i-1;
  if (i<=k) l=i+1;
}
}

/*************Need***************************************************************/

int findmaxspread(l,u,dimension,points)

/****************************************************************************/

int l,u,dimension;
double **points;

{
int i,j,maxdim;
double max       =-999999999.0,
       min       = 999999999.0,
       maxspread =-999999999.0;
for (i=0; i < dimension; i++) {
  max =-999999999.0;
  min = 999999999.0;
  for (j=l; j <= u; j++) {
    if (max < points[perm[j]][i]) { 
      max = points[perm[j]][i];
    }
    if (min > points[perm[j]][i]) { 
      min = points[perm[j]][i];
    }
    if (maxspread < fabs(max-min)) {
      maxspread = fabs(max-min);
      maxdim = i;
    }
  }
}
return(maxdim);
}

/********Need***********************************************************************/

optkdNode *BuildkdTree(points,l,u,dimension)

/*******************************************************************************/

int l,u;
double **points;

{
optkdNode *p;
int m;

  NEWTREE(p);
  if (u-l+1 <= BUCKETSIZE) {
    p->bucket = 1;
    p->lopt = l;
    p->hipt = u;
    p->loson = NULL;
    p->hison = NULL;
  } else {
    p->bucket =0;
    p->discrim = findmaxspread(l,u,dimension,points);
    m=(l+u)/2;
    Selection(points,l,u,m,p->discrim);
    p->cutval = points[perm[m]][p->discrim];
    p->loson = BuildkdTree(points,l,m,dimension);
    p->hison = BuildkdTree(points,m+1,u,dimension);
  }
  return(p);
}

/***Need****************************************************************************/

optkdNode *BuildOptTree(points,numPoints,dimension)

/*******************************************************************************/
int dimension,numPoints;
double **points;
{

int j;

  /* initialize perm array */
  perm = (int *) mxMalloc(numPoints*sizeof(int));
  for (j=0; j < numPoints; j++) {
    perm[j]=j;
  }
  return(BuildkdTree(points,0,numPoints-1,dimension));
}


/***********Need********************************************************************/

void rnnEuclidean(p,querpoint,points,dimension,numpoints,optfound)

/*******************************************************************************/

/* special searching algorithm to take advantage of the fact that square roots
   do not need to be evaulated */

optkdNode *p;
double *querpoint;
double **points;
int dimension,numpoints;
int * optfound;

{
  int i,j,k;
  double d,thisdist,val,thisx;

  if (p->bucket) {
    for (i=p->lopt; i <= p->hipt; i++) {
      thisdist=0.0;
      for (j=0; j<dimension; j++) {
        d=(querpoint[j]-points[perm[i]][j]);
        thisdist=thisdist+d*d;
      }        

      if (optfound[0] < numpoints) {
        PQInsert(thisdist,perm[i],nndist,optfound);
      } else {
        PQreplace(thisdist,nndist,optfound,perm[i]);
      }
    }
  } else {
    val = querpoint[p->discrim] - p->cutval;
    if (val < 0) {
      rnnEuclidean(p->loson,querpoint,points,dimension,numpoints,optfound);
      if (nndist[1] >= val*val) {
        rnnEuclidean(p->hison,querpoint,points,dimension,numpoints,optfound);
      }
    } else {
      rnnEuclidean(p->hison,querpoint,points,dimension,numpoints,optfound);
      if (nndist[1] >= val*val) {
        rnnEuclidean(p->loson,querpoint,points,dimension,numpoints,optfound);
      }
    }
  }
}


/*****************Need**************************************************************/

int *kdOptNNQuery(points,dimension, querpoint,numNN,Metric,root,MinkP)

/*******************************************************************************/

optkdNode *root;
double *querpoint, **points;
int dimension,numNN,MinkP;

{
  int j;
  int *optfound;
  /* set up found array */
  optfound = (int *) mxMalloc((numNN+1)*sizeof(int));
  optfound[0]=1;  /* for now */

  /* nndist is a priority queue of the distances of the nearest neighbors found */
  nndist = (double *)mxMalloc((numNN+1)*sizeof(double));
  for (j=0; j < numNN+1; j++) {
    nndist[j] = 99999999999.0;
  }

  switch(Metric) {
    case EUCLIDEAN :  rnnEuclidean(root,querpoint,points,dimension,numNN,optfound);
                      break;
    /*case MANHATTAN : Distance=ManhattDist;
                      rnnGeneral(root,querpoint,points,dimension,numNN,MinkP);
                      break;
    case L_INFINITY: Distance=LInfinityDist;
                      rnnGeneral(root,querpoint,points,dimension,numNN,MinkP);
                      break;
    case L_P       : Distance=LGeneralDist;
                      rnnGeneral(root,querpoint,points,dimension,numNN,MinkP);
                      break;*/
  }
  /*for (j=0;j<numNN;j++)
  printf("%g\n",nndist[j]);
  free(nndist);*/
  return(optfound);
  /*free(optfound);*/
}


/******Need*********************************************************************/

void KillOptTree(P)

/***************************************************************************/

/*  Kills a kd-tree to avoid memory holes.   */


optkdNode *P;

{
  /*if (perm != NULL) {
    free(perm);
  }  /* free permutation array */

  if (P==NULL) {
    return;
  } /* just to be sure */
  if (P->loson != NULL) {
    KillOptTree(P->loson);
  }

  if (P->hison != NULL) {
    KillOptTree(P->hison);
  }

  free(P);

}



/*************************************************************************************/
/*  Code to implement the abstract data type priority queue for use in j nearest     */
/*  neighbor searching.  Actual implementation is done using heaps.                  */
/*                                                                                   */
/*  Adapted from Sedgewick's: Algorithms in C p. 148-160.                            */
/*************************************************************************************/

/* 
   The heap data structure consists of two priority queues.  One for the j-smallest
   distances encountered, one to keep the indexes into the points array of the  
   points corresponding to the j-smallest distances. 
*/



/*********Need****************************************************************************/

void PQupheap(DistArr,FoundArr,k)

/*************************************************************************************/

double *DistArr;  /* j-smallest distances encountered */

int *FoundArr,k;

{
double v;
int j;

v=DistArr[k]; DistArr[0] = 999999999999999.0;
j=FoundArr[k];

while(DistArr[k/2] <= v) {
  DistArr[k] = DistArr[k/2];
  FoundArr[k] = FoundArr[k/2];
  k=k/2;
}
DistArr[k] = v;
FoundArr[k] = j;
}

/***********Need**************************************************************************/

void PQInsert(distance,index,DistArr,FoundArr)

/*************************************************************************************/

double distance,*DistArr;
int   index, *FoundArr;

{
  FoundArr[0]=FoundArr[0]+1;
  DistArr[FoundArr[0]] = distance;
  FoundArr[FoundArr[0]] = index;
  PQupheap(DistArr,FoundArr,FoundArr[0]);
}



/************Need*************************************************************************/

void PQdownheap(DistArr,FoundArr,k,index)

/*************************************************************************************/

double *DistArr;  /* j-smallest distances encountered */

int *FoundArr,k,index;

{

int j,l,N;
double v;

v=DistArr[k];

N = FoundArr[0];  /* tricky patch to maintain the data structure */
FoundArr[0]=index;

while (k <= N/2) {
  j=k+k;
  if (j < N && DistArr[j] <DistArr[j+1]) j++;
  if (v>=DistArr[j]) break;
  DistArr[k]=DistArr[j]; 
  FoundArr[k]=FoundArr[j];
  k=j;
}

DistArr[k] = v;
FoundArr[k]= index;
FoundArr[0]=N;  /* restore data struct */


}

/*************Need************************************************************************/

void PQreplace(distance,DistArr,FoundArr,index)

/*************************************************************************************/

double *DistArr,distance;
int *FoundArr;

{
  DistArr[0]=distance;
  PQdownheap(DistArr,FoundArr,0,index);
}


/******************************************************************************************************/

void principal (double numpoints, double dimension, double *rrw, double kneighbors, double *prLxo /*, double *Graph*/, double gamma) 

/******************************************************************************************************/
{
    double **Points;
    double * querpoint;
    optkdNode *OptkdRoot;
    int *Found, MinkP;
    int Metric;
    int i, j, k;
    Metric = EUCLIDEAN;
    MinkP = 0;
	numpoints = (int) numpoints;
	dimension = (int) dimension;
	kneighbors = (int) kneighbors;
    querpoint=(double *) mxMalloc(sizeof(double)*(dimension));
    
    Points = (double **)mxMalloc((numpoints)*sizeof(double *));

    for (k = 0; k < numpoints; k++) {
        Points[k] = (double *)mxMalloc((dimension)*sizeof(double));
    }
    for (k = 0, i = 0; k < numpoints; k++) {
        for (j = 0; j < dimension; j++, i++) {
            Points[k][j] = rrw[i];
        }
    }
    *prLxo = 0.0;
    /*Here we are going points first, then dimension
    for (j = 0, i = 0; j < dimension; j++) {
	for (k = 0; k < numpoints; k++, i++) {
		Points[k][j] = rrw[i];
        }
    }*/

    OptkdRoot = BuildOptTree(Points, numpoints, dimension);
    
    for (i = 0, k = 0; i < numpoints; i++) {        
            
	    for (j = 0; j < dimension; j++) 
               querpoint[j]=Points[i][j];
            
            Found = kdOptNNQuery (Points, dimension, querpoint, kneighbors, Metric, OptkdRoot, MinkP); 
            
          for (j = kneighbors; j >= 1; j--, k++) /*Found[0] is a dummy. */
	    {
		*prLxo += pow(sqrt( nndist[j] ),gamma);
/*		Graph[k] = Found[j] + 1;*/
	    }
    }
   /*KillOptTree(OptkdRoot);
    OptkdRoot=NULL;*/
    /*free(Found);*/
    /*for (k=0;k<numpoints;k++)
    free(Points[k]);
    free(querpoint);*/
}

