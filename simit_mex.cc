#include "mex.h"
#include<iostream>
#include<fstream>
#include<stdlib.h>
#include<math.h>

#include "def.h"
#include "model.h"
#include "nr.h"

using namespace std;

//Begin definition of module prototypes

GlobExp *parseopts(int argc, char *argv[],Glob *globs,char *outstr);

void initialise(GlobExp ex[],Glob *globs,int simit);

void tabulateValues (Glob *globs,GlobExp *ex,double t0, double t1, double *t, long n, double *state,
		     double **y, double ***dmds, double ***dmdp, double **dyds,double **dydp);

void outSimit(GlobExp ex[],Glob *globs,double *t,long n,double **y);
void freeMem(GlobExp *ex,Glob *globs,int simit);

//End definition of module prototypes


// DEBUG stream
ofstream *dbg;
// Input
#define IN   prhs[0]

// Output

#define N    plhs[0]
#define YOUT plhs[1]
#define DMDS plhs[2]
#define DMDP plhs[3]


void mexFunction( int nlhs,
		  mxArray *plhs[],
		  int nrhs,
		  const mxArray *prhs[]
		  )
{
  long k,i,j,n,ldum;
  char outstr[100];
  Glob globs;
  GlobExp *ex;
  double time,**y,*state;
  double t0; //init. time
  double *t,*y0,*p;
  double **yout,*yout_,**ydum;
  double ***dmds,***dmdp,*dmds_,*dmdp_;
  double **dyds,**dydp;
  double *n_,*spline_;
  double *tdum=dvector(1,1);
  mxArray *splinedat;
  mxArray *tmp;
  mxArray *SPLINE;
  double *tmp_;
  int hidden;

  // input -> initial time
  tmp=mxGetField(IN,0,"t");
  tmp_=mxGetPr(tmp);
  n=mxGetNumberOfElements(tmp);
  t=dvector(1,n);
  for(i=1;i<=n;i++)
    t[i]=tmp_[i-1];
 
  //allocation of init. state and parameters
  y0=dvector(1,NEQNS);
  p=dvector(1,NPARAMS);

  // input -> initial time
  tmp=mxGetField(IN,0,"t0");
  tmp_=mxGetPr(tmp);
  t0=tmp_[0];

  // input -> y0
  tmp=mxGetField(IN,0,"y0");
  tmp_=mxGetPr(tmp);
  if(mxGetNumberOfElements(tmp)!=NEQNS)
    {
      if(mxGetNumberOfElements(tmp)!=0)
	mexWarnMsgTxt("Insufficient number of initial values.\nUse predefined.\n");
      for(i=1;i<=NEQNS;i++)
	y0[i]=DefYValues[i-1];
    }
  else
    {
      for(i=1;i<=NEQNS;i++)
	y0[i]=tmp_[i-1];
    }
  
  // input -> p
  tmp=mxGetField(IN,0,"p");
  tmp_=mxGetPr(tmp);
  if(mxGetNumberOfElements(tmp)!=NPARAMS)
    {
      if(mxGetNumberOfElements(tmp)!=0)
	mexWarnMsgTxt("Insufficient number of parameters.\nUse predefined.\n");
      for(i=1;i<=NPARAMS;i++)
	p[i]=DefParameters[i-1];
    }  
  else
    {
      for(i=1;i<=NPARAMS;i++)
	p[i]=tmp_[i-1];
    }

  // input -> eps
  tmp=mxGetField(IN,0,"eps");
  tmp_=mxGetPr(tmp);
  globs.eps=fabs(tmp_[0]);
  
  // input -> integrator
  tmp=mxGetField(IN,0,"int");
  tmp_=mxGetPr(tmp);
  globs.integrator=abs((int)tmp_[0]);
  globs.stiff=TRUE;
  // input -> spline
  SPLINE=mxGetField(IN,0,"spline");
  
  // input -> hidden
  tmp=mxGetField(IN,0,"hidden");
  tmp_=mxGetPr(tmp);
  hidden=(int)tmp_[0];

  //allocation of trajectory vect. for output
  if(hidden==TRUE)
    YOUT=mxCreateDoubleMatrix(n,NOBS+NEQNS,mxREAL);
  else
    YOUT=mxCreateDoubleMatrix(n,NOBS,mxREAL);
  yout_=mxGetPr(YOUT);  
  yout=dmatrix(1,n,1,NOBS+NEQNS);
  ydum=dmatrix(1,2,1,NOBS+NEQNS);

  //initialise some global parameters in globs
  globs.noGnu=TRUE;
  globs.npar=NPARAMS;
  globs.noMeasurements=FALSE;
  globs.doP=ivector(1,globs.npar);
  for(k=1;k<=globs.npar;k++)
    globs.doP[k]=TRUE;
  globs.maxit=1000;
  globs.gnuFp=NULL;
  globs.wait=FALSE;
  globs.usesig=FALSE;
  globs.stiff=TRUE;
  globs.maxstp=5000;
  globs.minimp=0.05;
  globs.nowait=FALSE;
  globs.elastic=1.;
  globs.reg=FALSE;
  globs.epsilon=1e-10;
  globs.lambda=1e6;
  globs.dt=0.1;
  globs.sig=0.;
  globs.nrExp=1;  

  ex=new GlobExp[globs.nrExp+1];   
  //initialise some global parameters

  ex[1].nobs=NOBS;
  ex[1].nvar=NEQNS;
  ex[1].par=dvector(1,NPARAMS);
  ex[1].y0=dvector(1,NEQNS);


  for(k=1;k<=NPARAMS;k++)
    ex[1].par[k]=p[k];

 
  ex[1].fitstart=t[1];
  ex[1].fitend=t[n];
  //initial values

  state=dvector(1,ex[1].nvar);
  for(k=1;k<=NEQNS;k++)
    state[k]=y0[k];


  N=mxCreateDoubleMatrix(3,1,mxREAL);
  n_=mxGetPr(N);
  
  n_[0]=NEQNS;
  n_[1]=NOBS;
  n_[2]=NPARAMS;
  
  // spline section
  if(NSPLINES==0)
    {
      globs.initSpline=TRUE;
      initialise(ex,&globs,TRUE);
    }
  else
    {
      globs.initSpline=FALSE;      
      initialise(ex,&globs,TRUE);  
      ex[1].splineNodes=(double **) malloc((size_t)(NSPLINES+1)*sizeof(double*));
      ex[1].splineY=(double **) malloc((size_t)(NSPLINES+1)*sizeof(double*));
      ex[1].splineGam=(double **) malloc((size_t)(NSPLINES+1)*sizeof(double*));
      ex[1].nNodes=(long unsigned*) malloc((size_t)(NSPLINES+1)*sizeof(long unsigned*));

      for(i=1;i<=NSPLINES;i++)
	{
	  if(mxGetNumberOfElements(SPLINE) < NSPLINES)
	    {
	      cerr << "Please specify " << NSPLINES << " spline(s).\n";
	      return;
	    }

	  splinedat=mxGetCell(SPLINE,i-1);
	  spline_=mxGetPr(splinedat);
	  ex[1].nNodes[i]=mxGetM(splinedat);
	  if(mxGetN(splinedat)!=3)
	    {
	      cerr << "Invalid data for spline " << i << ".\n";
	      return;
	    }
	  ex[1].splineNodes[i]=(double *) malloc((size_t)(ex[1].nNodes[i]+1)*sizeof(double));
	  ex[1].splineY[i]=(double *) malloc((size_t)(ex[1].nNodes[i]+1)*sizeof(double));
	  ex[1].splineGam[i]=(double *) malloc((size_t)(ex[1].nNodes[i]+1)*sizeof(double));	      

	  //copy spline data
	  for(j=1;j<=ex[1].nNodes[i];j++)
	     {
	       ex[1].splineNodes[i][j]=spline_[j-1];
	       ex[1].splineY[i][j]=spline_[ex[1].nNodes[i]+j-1];
	       ex[1].splineGam[i][j]=spline_[2*ex[1].nNodes[i]+j-1];	    
	     }
	}
    }

  //... and integrate
  try
    {
      if(nlhs!=2)
	{
	  dmds=d3tensor(1,n,1,NEQNS,1,NEQNS);
	  dmdp=d3tensor(1,n,1,NEQNS,1,NPARAMS);
	  dyds=dmatrix(1,NEQNS,1,NEQNS);
	  dydp=dmatrix(1,NEQNS,1,NPARAMS);
	  
	  tabulateValues(&globs,&ex[1],t0,t[n],t,n,state,yout,dmds,dmdp,dyds,dydp);
	  
	  DMDS=mxCreateDoubleMatrix(NOBS*NEQNS,n,mxREAL);
	  dmds_=mxGetPr(DMDS);
	  DMDP=mxCreateDoubleMatrix(NOBS*NPARAMS,n,mxREAL);
	  dmdp_=mxGetPr(DMDP);
	  

	  //copying 
	  
	  for(k=1;k<=NEQNS;k++) 
	    for(i=1;i<=NOBS;i++)
	      for(j=1;j<=n;j++)
		{
		  dmds_[(j-1)*NOBS*NEQNS+(k-1)*NOBS+i-1]=dmds[j][i][k];
		}
	  
	  
	  for(k=1;k<=NPARAMS;k++)
	    for(i=1;i<=NOBS;i++)
	      for(j=1;j<=n;j++)
		{
		  dmdp_[(j-1)*NOBS*NPARAMS+(k-1)*NOBS+i-1]=dmdp[j][i][k];
		}
	  
	  free_d3tensor(dmds,1,n,1,NEQNS,1,NEQNS);
	  free_d3tensor(dmdp,1,n,1,NEQNS,1,NPARAMS);
	  free_dmatrix(dyds,1,NEQNS,1,NEQNS);
	  free_dmatrix(dydp,1,NEQNS,1,NPARAMS);
	}
      else
	{
	  tabulateValues(&globs,&ex[1],t0,t[n],t,n,state,yout,NULL,NULL,NULL,NULL);
	  
	}

      if(hidden==TRUE)
	{
	  ldum=NOBS+NEQNS;
	  for(i=1;i<=NEQNS;i++)
	    state[i]=y0[i];
	  tdum[1]=t[1];
	  if(t0!=t[1])
	    tabulateValues(&globs,&ex[1],t0,t[1],tdum,1,state,ydum,NULL,NULL,NULL,NULL);
	  for(j=1;j<=NEQNS;j++)
	    yout[1][NOBS+j]=state[j];

	  for(i=2;i<=n;i++)
	    {
	      tdum[1]=t[i]; 
	      tabulateValues(&globs,&ex[1],t[i-1],t[i],tdum,1,state,ydum,NULL,NULL,NULL,NULL);
	      for(j=1;j<=NEQNS;j++)
		yout[i][NOBS+j]=state[j];
	     
	    }
	}
      else
	ldum=NOBS;

      for(i=1;i<=ldum;i++)
	for(j=1;j<=n;j++)
	  yout_[(i-1)*n+j-1]=yout[j][i];
    }
  catch (int i)
    {
      return;
    }
  //free mem
  freeMem(ex,&globs,TRUE);
  free_dvector(state,1,ex[1].nvar);
  free_dmatrix(yout,1,n,1,NEQNS+NOBS);
  free_dvector(t,1,n);
  free_dvector(y0,1,NEQNS);
  free_dvector(p,1,NPARAMS);
  free_dvector(ex[1].par,1,NPARAMS);
  free_dvector(ex[1].y0,1,NEQNS);
  free_dvector(tdum,1,1);
  free_dmatrix(ydum,1,2,1,NOBS+NEQNS);

  delete ex;
}
