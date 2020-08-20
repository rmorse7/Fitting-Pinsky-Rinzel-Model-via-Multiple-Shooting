#include<math.h>
#include<stdlib.h>
#include<iostream>

#include "../def.h"
#include "../model.h"
#include "../nr.h"

using namespace std;

void freeMem(GlobExp *ex,Glob *globs,int simit)
{
  long nExp,i;
  long nrExp=globs->nrExp;


  delete globs->gnuFp;
  
  if(simit==FALSE)
    {
      for (nExp=1;nExp<=nrExp;++nExp) 
	{ 
	  long nPoints=ex[nExp].nPoints;
	  long nP=globs->npar;
	  long nvar=ex[nExp].nvar, nobs=ex[nExp].nobs;
	  long nMeas=ex[nExp].nMeasure;
	  long  me=ex[nExp].me, mg=ex[nExp].mg;
	  free_dvector(ex[nExp].errP,1,nP);
	  free_dvector(ex[nExp].errY0,1,nvar);
	  free_dmatrix(ex[nExp].yTry,1,nPoints,1,nvar);
	  free_dmatrix(ex[nExp].yTrySave,1,nPoints,1,nvar);
	  free_dmatrix(ex[nExp].yComp,1,nPoints,1,nvar);
	  free_dmatrix(ex[nExp].yPred,1,nMeas,1,nobs);
	  free_dmatrix(ex[nExp].h,1,nPoints,1,nvar);   
	  free_dmatrix(ex[nExp].residues,1,nMeas,1,nobs);
	  free_d3tensor(ex[nExp].dyds,1,nPoints,1,nvar,1,nvar);
	  free_d3tensor(ex[nExp].dydp,1,nPoints,1,nvar,1,nP);
	  free_d3tensor(ex[nExp].dmds,1,nMeas,1,nobs,1,nvar);
	  free_d3tensor(ex[nExp].dmdp,1,nMeas,1,nobs,1,nP);
	  if (me>0) 
	    {
	      free_dvector(ex[nExp].r2,1,me);
	      free_d3tensor(ex[nExp].dR2ds,1,me,1,nPoints,1,nvar);
	      free_dmatrix(ex[nExp].dR2dp,1,me,1,nP);
	    }
	  if (mg>0) 
	    {
	      free_dvector(ex[nExp].r3,1,mg);
	      free_d3tensor(ex[nExp].dR3ds,1,mg,1,nPoints,1,nvar);
	      free_dmatrix(ex[nExp].dR3dp,1,mg,1,nP);
	    }
	  
	  free_dvector(ex[nExp].ua,1,nMeas*nobs);
	  free_dmatrix(ex[nExp].Ea,1,nMeas*nobs,1,nvar);
	  free_dmatrix(ex[nExp].Pa,1,nMeas*nobs,1,nP);
	  free_dvector(ex[nExp].ue,1,me);
	  free_dmatrix(ex[nExp].Ee,1,me,1,nvar);
	  free_dmatrix(ex[nExp].Pe,1,me,1,nP);
	  free_dvector(ex[nExp].ug,1,mg);
	  free_dmatrix(ex[nExp].Eg,1,mg,1,nvar);
	  free_dmatrix(ex[nExp].Pg,1,mg,1,nP);
	  free_dmatrix(ex[nExp].dS,1,nPoints,1,nvar);
	  free_dvector(ex[nExp].dP,1,nP);
	  
	}
    }
  for (nExp=1;nExp<=nrExp;++nExp) 
    {
      free(ex[nExp].splineNodes);
      free(ex[nExp].splineY);
      free(ex[nExp].splineGam);
      free(ex[nExp].nNodes);
      for(i=1;i<=NSPLINES;i++)
	{
	  free(ex[nExp].splineNodes[i]);
	  free(ex[nExp].splineY[i]);
	  free(ex[nExp].splineGam[i]);
	}
    }
  if(simit==FALSE && globs->covar!=NULL)
    free_dmatrix(globs->covar,1,globs->fitdim,1,globs->fitdim);
}
