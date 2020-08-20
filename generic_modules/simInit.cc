#include<iostream>
#include<fstream>
#include<math.h>
#include<string.h>
#include<stdio.h>

#include "../def.h"
#include "../model.h"
#include "../nr.h"

using namespace std;

//modules

//Begin definition of module prototypes

void tabulateValues (Glob *globs,GlobExp *ex,double t0, double t1, double *t, 
		     long n, double *state, double **y, double ***dmds, 
		     double ***dmdp, double **dyds,double **dydp);

void simInit(GlobExp *ex,Glob *globs)
{
  long i,j;
  long nvar=ex->nvar;
  long nPoints=ex->nPoints;
  long idum=-42;
  double *t=dvector(1,2);
  double **y=dmatrix(1,2,1,nvar);
  //RCL 08/02/19
  double xi, d, w1, w2;
  long nobs=ex->nobs, nMeasure=ex->nMeasure, k;


  if (ex->y0) 
    {
      //read from command line
      for (i=1; i<=nvar; ++i) 
	ex->yTry[1][i]=ex->y0[i];
    } 
  else
    {
      for (i=1; i<=nvar; ++i) 
	{
	  //read from model.cc
	  ex->yTry[1][i]=DefYValues[i-1];
	}
    }

  //RCL
  k=1;

  //integrating 

  for(i=2;i<=nPoints;i++)
    {
      for(j=1;j<=nvar;j++)
	ex->yTry[i][j]=ex->yTry[i-1][j];
      t[1]=ex->mesh[i-1];
      t[2]=ex->mesh[i];

      //defined in intODE.cc 
      tabulateValues(globs,ex,t[1],t[2],t,2,ex->yTry[i],y,NULL,NULL,NULL,NULL);
      for(j=1;j<=nvar;j++)
	//add Gaussian jitter
	ex->yTry[i][j]+=globs->pert*ex->yTry[i][j]*gasdev(&idum);

      //RCL 07/29/2019 Testing a hypothesis that the observable values at MS boundaries are not being used

      //RCL 08/02/19 Hints taken from setInitialValues.cc
      xi=ex->mesh[i];
      if (ex->xMeasure[k]>=ex->mesh[i])
        {
          for (j=1; j<=nvar; ++j)
            if (j<=nobs)
              //this assumes that the first nobs are the first variables
              ex->yTry[i][j]=ex->yMeasure[k][j];
        }
      else if (k==nMeasure) // we have reached the last point, as per increment of k below
        {
          // after last measurement: take last measured value
          for (j=1; j<=nvar; ++j)
            if (j<=nobs)
              ex->yTry[i][j]=ex->yMeasure[k][j];
            //else
            //    ex->yTry[i][j]=1;
        }
      else
        {
        // interpolate linearly
          while (k<nMeasure && ex->xMeasure[k]<xi) //terminates either at the  last index
            //or at the first one that is at least equal to xi
            ++k;
          d=ex->xMeasure[k]-ex->xMeasure[k-1];
          w1=(ex->xMeasure[k]-xi)/d;
          w2=(xi-ex->xMeasure[k-1])/d;
          for (j=1; j<=nvar; ++j)
            if (j<=nobs)
              ex->yTry[i][j]=w1*ex->yMeasure[k-1][j] + w2*ex->yMeasure[k][j];
        }

    }
   


  free_dvector(t,1,2);
  free_dmatrix(y,1,2,1,nvar);
}
