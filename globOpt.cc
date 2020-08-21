#include<iostream>
#include<fstream>
#include<stdlib.h>
#include<stdio.h>

#include "libSRES/sharefunc.h"
#include "libSRES/ESSRSort.h"
#include "libSRES/ESES.h"

#include "def.h"
#include "model.h"
#include "nr.h"

using namespace std;

//Begin definition of module prototypes

void intODE(GlobExp *ex,Glob *globs,int doDerivs,int doPlot,long expNr);

double computeRight(Glob *globs,GlobExp *ex);

void solvLin(Glob *globs,GlobExp *ex,int computeCovar);

double dampIt(Glob *globs,GlobExp *ex,double **P0, double ***S0,
              double **dP, double ***dS,double *uS);


//Pointer to the structures ex and globs
Glob *globsPtr;
GlobExp *exPtr;


//fitness function for libSRES

void fitness(double *x,double *f,double *g)
{
  long nExp,j;
  long pos=0;
  double sum=0;
  
  for (nExp=1;nExp<=globsPtr->nrExp;++nExp)
    {
      for(j=1;j<=NEQNS;j++)
	{
	  if(globsPtr->y0fix[j]==TRUE)
	    {
	      exPtr[nExp].yTry[1][j]=x[pos];
	      pos++;
	    }
	} 
      for(j=1;j<=globsPtr->npar;j++)
	{
	  if(globsPtr->doP[j]=='L')
	    {
	      exPtr[nExp].par[j]=x[pos];
	      pos++;
	    }
	} 
    }
  for(j=1;j<=globsPtr->npar;j++)
    {
      if(globsPtr->doP[j]==TRUE)
	{
	  for (nExp=1;nExp<=globsPtr->nrExp;++nExp)
	    exPtr[nExp].par[j]=x[pos];
	  pos++;
	}
    } 
  
  
  for (nExp=1;nExp<=globsPtr->nrExp;++nExp)
        {
	  try
	    {
	      intODE(&exPtr[nExp],globsPtr,FALSE,FALSE,nExp);
	      // set up internal data (residues, etc.)
	      sum+=computeRight(globsPtr,&exPtr[nExp]);
	    }
	  catch(int i)
	    {
	      sum+=1e10;
	    }
	}
  globsPtr->chisq=sum;
  (*f)=sum;
}


void globOpt(GlobExp ex[],Glob *globs)
{
  globsPtr=globs;
  exPtr=ex;


  globsPtr->nIter=0;
  long i,j;
  //variables for libSRES
  ESParameter *param;
  ESPopulation *population;
  ESStatistics *stats;
  ESfcnTrsfm *trsfm;
  unsigned int seed;
  int es;
  int dim=0;
  double *ub, *lb;
  int miu, lambda, gen;
  double gamma, alpha, varphi;
  int retry;
  double pf; 

  //keeping old values
  long unsigned *nms=lvector(1,globsPtr->nrExp);
  for(i=1;i<=globsPtr->nrExp;i++)
    {
      nms[i]=exPtr[i].nPoints;
      exPtr[i].nPoints=2;
    }

  seed = shareDefSeed;
  gamma = esDefGamma;
  alpha = esDefAlpha;
  varphi = esDefVarphi;
  retry = esDefRetry;
  pf = essrDefPf;
  es = esDefESSlash;
  miu = 30;
  lambda = 200;
  gen = 1750;
  trsfm = NULL;

  //determine # of decision variables
  for(i=1;i<=globsPtr->nrExp;i++)
    {
      for(j=1;j<=NEQNS;j++)
	{
	  if(globsPtr->y0fix[j]==TRUE)
	    dim++;
	}
      for(j=1;j<=globsPtr->npar;j++)
	{
	  if(globsPtr->doP[j]=='L')
	    dim++;
	}
    }      
  for(j=1;j<=globsPtr->npar;j++)
    {
      if(globsPtr->doP[j]==TRUE)
	dim++;
    }


  //bounds ...for test only

  lb=dvector(0,dim);
  ub=dvector(0,dim);

  for(i=0;i<dim;i++)
    {
      lb[i]=0;
      ub[i]=2;
    }

  ESInitial(seed, &param, trsfm, fitness, es,0, dim, ub, lb, \
	    miu, lambda, gen, gamma, alpha, varphi,  \
            retry, &population, &stats);
  
  for(int i=1;i<=globsPtr->maxit;i++)
    { 
      globsPtr->nIter++;
      ESStep(population, param, stats, pf);
      if(stats->curgen >= param->gen)
	break;
    }
  
  //copying old stuff back
  for(i=1;i<=globsPtr->nrExp;i++)
    {
      exPtr[i].nPoints=nms[i];
    }

  cout << endl;

  free_dvector(lb,0,dim);
  free_dvector(ub,0,dim);
  free_lvector(nms,1,globsPtr->nrExp);
}
