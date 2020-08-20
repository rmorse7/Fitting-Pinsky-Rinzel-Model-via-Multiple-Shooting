#include<iostream>
#include<math.h>
#include<stdlib.h>

#include "../def.h"
#include "../model.h"
#include "../nr.h"

using namespace std;

// objective Function
double objectiveF(Glob *globs,GlobExp *ex)
{
  int i,j;
  double df,f;
  
  f=0; //this will be the value of the objective function
  for (i=1; i<=ex->nMeasure; ++i) //loop over time measurements
    {
      for (j=1; j<=ex->nobs; ++j)  //loop over observations
	{
	  //the residues are calculated in the function computeRight; square them and normalize by variance (sigma^2)
	  df = (ex->residues[i][j]*ex->residues[i][j])/
	    (ex->sigma[i][j]!=0 ? (ex->sigma[i][j]*ex->sigma[i][j]) : 1.0); //use 1 if sigma is zero
#ifdef DEBUGF
	  if (df>0.1) 
	    *dbg << "f increment " << df << " for i = " << i << ", j = " << j << "\n";
#endif
	  f += df;
	}
    }
  ex->objF=f;
#ifdef DEBUGF
  *dbg << "nMeasure = " << ex->nMeasure << ", nobs = " << ex->nobs << ", f = " << f << "\n";
#endif
  return(f);
}

double computeRight(Glob *globs,GlobExp *ex)
{
    long i,j,k;

    // calculate discrepancies (deviation from continuous trajectory), residues and constraints
    long nPoints=ex->nPoints, nvar=ex->nvar;
    long nMeasure=ex->nMeasure, nobs=ex->nobs;
    long nP=globs->npar;

    // no continuity constraints for first and last mesh point; set discrepancy h to zero
    for (j=1;j<=nvar;++j) {
	ex->h[1][j]=0;
	ex->h[nPoints][j]=0;
    }
    for (i=2;i<nPoints;++i) //compute discrepancy
      {
	for (j=1;j<=nvar;++j)
	  ex->h[i][j]=ex->yComp[i][j]-ex->yTry[i][j];
      }
    for (i=1;i<=nMeasure;++i) //compute residues; loop over time measurements
      {
	for (j=1;j<=nobs;++j) //loop over observations
	  ex->residues[i][j]=ex->yPred[i][j]-ex->yMeasure[i][j];
      }

    //extra equality and inequality constraints
    R2(globs,ex, TRUE);
    R3(globs,ex ,TRUE);

#ifdef PRINTDERIVS
    for (i=1; i<=nPoints; ++i) {
	*dbg << "dyds[" << i << "] = (";
	for (j=1; j<=nvar; ++j) {
	    for (k=1; k<=nvar; ++k)
		*dbg << " " << ex->dyds[i][j][k];
	    *dbg << ", ";
	}
	*dbg << ")\ndydp[" << i << "] = (";
	for (j=1; j<=nvar; ++j) {
	    for (k=1; k<=nP; ++k)
		*dbg << " " << ex->dydp[i][j][k];
	    *dbg << ", ";
	}
	*dbg << ")\n";
    }
    for (i=1; i<=nMeasure; ++i) {
	*dbg << "dmds[" << i << "] = (";
	for (j=1; j<=nobs; ++j) {
	    for (k=1; k<=nvar; ++k)
		*dbg << " " << ex->dmds[i][j][k];
	    *dbg << ", ";
	}
	*dbg << ")\ndmdp[" << i << "] = (";
	for (j=1; j<=nobs; ++j) {
	    for (k=1; k<=nP; ++k)
		*dbg << " " << ex->dmdp[i][j][k];
	    *dbg << ", ";
	}
	*dbg << ")\n";
    }
    for (i=1; i<=ex->me; ++i) {
	*dbg << "dR2ds[" << i << "] = (";
	for (j=1; j<=nPoints; ++j) {
	    for (k=1; k<=nvar; ++k)
		*dbg << " " << ex->dR2ds[i][j][k];
	    *dbg << ", ";
	}
	*dbg << ")\n";
	*dbg << "dR2dp[" << i << "] = (";
	for (j=1; j<=nP; ++j)
	    *dbg << " " << ex->dR2dp[i][j];
	*dbg << ")\n";
    }
    for (i=1; i<=ex->mg; ++i) {
	*dbg << "dR3ds[" << i << "] = (";
	for (j=1; j<=nPoints; ++j) {
	    for (k=1; k<=nvar; ++k)
		*dbg << " " << ex->dR3ds[i][j][k];
	    *dbg << ", ";
	}
	*dbg << ")\n";
	*dbg << "dR3dp[" << i << "] = (";
	for (j=1; j<=nP; ++j)
	    *dbg << " " << ex->dR3dp[i][j];
	*dbg << ")\n";
    }
    dbg->flush();
#endif

    return(objectiveF(globs,ex));
}
