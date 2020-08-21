#include "mex.h"
#include <iostream>
#include <stdlib.h>
#include <fstream>
#include <math.h>
#include <string>
#include <getopt.h>
#include "nr.h"

using namespace std;

#define TRUE 1
#define FALSE 0


typedef struct 
{
  double *t;
  double *y;
  double *sig;
  double maxLambda;
  long nKnots;
  int gcv;
} glob;


// Input
#define DATA   prhs[0]
#define ALPHA  prhs[1]
#define MAXA   prhs[2]
#define NOGCV  prhs[3]

// Output

#define SPLINE plhs[0]
#define ALPHA_OUT plhs[1]

void mexFunction( int nlhs,
		  mxArray *plhs[],
		  int nrhs,
		  const mxArray *prhs[]
		  )
{
  
  glob globs;
  long i,j,n,col;
  double *g,*gam,alpha;
  //mex input
  double *data_,*alpha_,*maxa_,*nogcv_;
  //mex output
  double *spline_,*alpha_out_;
   
  //initialise mexstuff
  data_=mxGetPr(DATA);
  alpha_=mxGetPr(ALPHA);
  maxa_=mxGetPr(MAXA);
  nogcv_=mxGetPr(NOGCV);

  globs.nKnots=mxGetM(DATA);
  col=mxGetN(DATA);  

  n=globs.nKnots;
  g=dvector(1,n);
  gam=dvector(1,n);

  globs.t=dvector(1,n);
  globs.y=dvector(1,n);
  globs.sig=dvector(1,n);

  if(col==2)
    {
      for(i=1;i<=n;i++)
	{
	  globs.t[i]=data_[i-1];
	  globs.y[i]=data_[n+i-1];
	  globs.sig[i]=1.;
	}
    }
  else if(col==3)
    {
      for(i=1;i<=n;i++)
	{
	  globs.t[i]=data_[i-1];
	  globs.y[i]=data_[n+i-1];
	  globs.sig[i]=data_[2*n+i-1];
	}
    }
  else
    {
      mexWarnMsgTxt("Insufficient data structure. Exit.\n");
      return;
    }
  
  alpha=fabs(alpha_[0]);
  globs.maxLambda=fabs(maxa_[0]);
  globs.gcv=(int)nogcv_[0];


  if(alpha > globs.maxLambda)
    alpha=globs.maxLambda;

  if(globs.gcv==TRUE)
    splines_gcv(globs.t,globs.y,globs.sig,n,&alpha,g,gam,globs.maxLambda);
  else
    splines(globs.t,globs.y,globs.sig,n,alpha,g,gam);

  //writing output

  SPLINE=mxCreateDoubleMatrix(n,3,mxREAL);
  spline_=mxGetPr(SPLINE);
  ALPHA_OUT=mxCreateDoubleMatrix(1,1,mxREAL);
  alpha_out_=mxGetPr(ALPHA_OUT);
  alpha_out_[0] = alpha;

  spline_[0]=globs.t[1];
  spline_[n]=g[1];
  spline_[2*n]=0.;
  for(i=2;i<=n-1;i++)
    {
      spline_[i-1]=globs.t[i];
      spline_[n+i-1]=g[i];
      spline_[2*n+i-1]=gam[i-1];
    }
  spline_[n-1]=globs.t[n];
  spline_[2*n-1]=g[n];
  spline_[3*n-1]=0.; 

  //free memory
  free_dvector(globs.y,1,n);
  free_dvector(globs.t,1,n);
  free_dvector(globs.sig,1,n);
  free_dvector(g,1,n);
  free_dvector(gam,1,n);
}

