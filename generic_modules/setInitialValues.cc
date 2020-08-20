#include<iostream>
#include<fstream>
#include <math.h> //122118: added due to hack below

#include "../def.h"
#include "../model.h" //RMM added for model specific unobserved var initliazation (using steady-state eqs)

using namespace std;

void setInitialValues(GlobExp *ex, int nP, double *parameters, double *yValues)
{
  long i,j,k;
  double xi, d, w1, w2;
  long nobs=ex->nobs, nvar=ex->nvar, nMeasure=ex->nMeasure;
  long nPoints=ex->nPoints;
  
  /* set initial values */
  //Note: it appears that yValues will be NULL because the call from initialize.cc does not yet have y0 defined.
  if (yValues) 
    {
      //for the first shooting interval, initialize as set in the initial conditions
      for (j=1; j<=nvar; ++j) //if observed vars are not the first several vaars, can just pull j's from list
	ex->yTry[1][j]=yValues[j];
      i=2; // index for multiple shooting interval
    } 
  else
    i=1;

  k=1;
  
  while (i<=nPoints) //loop over the multiple shooting intervals boundary points
    {
      xi=ex->mesh[i]; //current mesh time point
      
      // make educated guess for yTry(xi)
      // could be improved (polynomial interpolation/extrapolation etc.)
      
      if (ex->xMeasure[k]>=xi) 
	{
	  // filling in observed measurements directly, this loop could be problematic if first vars are not measured ones (RMM)
	  for (j=1; j<=nobs; ++j)
	      //this assumes that the first nobs are the first variables
	      ex->yTry[i][j]=ex->yMeasure[k][j];

	  //model specific initialization of none observed variables
	  unobsinit(ex->yTry[i], parameters);
	  
	} else if (k==nMeasure) // we have reached the last point, as per increment of k below 
	  {
	    // after last measurement: take last measured value
	    for (j=1; j<=nvar; ++j)
	      if (j<=nobs)
		ex->yTry[i][j]=ex->yMeasure[k][j];
	      else
		ex->yTry[i][j]=1;
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
	  for (j=1; j<=nobs; ++j) 
	    ex->yTry[i][j]=w1*ex->yMeasure[k-1][j] + w2*ex->yMeasure[k][j];
	  unobsinit(ex->yTry[i], parameters);
	}
      for (j=1; j<=nvar; ++j)
	      *dbg << "yTry Me # " << ex->yTry[i][j] << '\n';
      ++i;
    }  
}
