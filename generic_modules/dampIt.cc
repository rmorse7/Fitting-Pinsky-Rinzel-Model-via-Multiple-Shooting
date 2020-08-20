#include<iostream>
#include<fstream>
#include<math.h>
#include<string.h>
#include<stdio.h>

#include "../def.h"
#include "../model.h"
#include "../nr.h"

using namespace std;



//Begin definition of module prototypes

//Module modules/setInitialValues.cc

void intODE(GlobExp *ex,Glob *globs,int doDerivs,int doPlot,long expNr);

double computeRight(Glob *globs,GlobExp *ex);

void solvLin(Glob *globs,GlobExp *ex,int computeCovar);

//End definition of module prototypes

//norm of dx

double Norm(Glob *globs, GlobExp *ex,double ***dS, double **dP,double uS) 
{
  double norm=0;
  if (!dS) 
    {
    cerr <<"dS==NULL in Norm\n";
    throw 1;
    }
  long i,ivar,iexp,ip;
  for (iexp=1;iexp<=globs->nrExp;++iexp)
    {
      for (i=1; i< ex[iexp].nPoints; ++i)
	for (ivar=1; ivar<=ex[iexp].nvar; ++ivar)
	  norm+=pow(uS*dS[iexp][i][ivar],2);

      for (ip=1;ip<=globs->npar; ++ip)
	{
	  if(globs->doP[ip]=='L')
	    norm+=pow(dP[iexp][ip],2);
	}
    } // end of: for (iexp=1;iexp<=globs->nrExp;++iexp)
  for (ip=1;ip<=globs->npar; ++ip)
    {
      if(globs->doP[ip]==TRUE)
	norm+=pow(dP[1][ip],2);
    }
  return sqrt(norm);
  
} // Norm


// perform one damping subiteration

double subIt(Glob *globs, GlobExp *ex, double t, double normdx,
	     double **p0, double ***s0, double **dp, double ***ds,double uS) 
{
  long iexp,i,j;
  long maxPoints=0,maxVar=0;
  long nP=globs->npar;

  for (i=1;i<=globs->nrExp;++i) 
    {
      if (ex[i].nPoints>maxPoints) 
	maxPoints=ex[i].nPoints;
      if (ex[i].nvar>maxVar)       
	maxVar=ex[i].nvar;
    }
  
  double **newdP=dmatrix(1,globs->nrExp,1,nP);                      // parameter updates
  double ***newdS=d3tensor(1,globs->nrExp,1,maxPoints,1,maxVar); // yTry updates
  
  
  for(iexp=1;iexp<=globs->nrExp;iexp++) 
    {
      GlobExp *expi=&ex[iexp];
      for(i=1; i< expi->nPoints; i++)
	{
	  for(j=1; j<=expi->nvar; j++)
	    {
	      expi->yTry[i][j] = s0[iexp][i][j]+t*uS*ds[iexp][i][j];
	      // positiveness of initial values --> test
	      if(expi->yTry[i][j] < 0.0)
		expi->yTry[i][j] = 0.0;
	      expi->h[i][j]*=uS;
	    }
	}
      for(i=1; i<=nP; i++)
	{
	  expi->par[i] = p0[iexp][i]+t*dp[iexp][i];
	}
    }

  //integrate and compute value of objective function
  globs->chisq=0.;
  for (iexp=1;iexp<=globs->nrExp;iexp++) 
    {
      intODE(&ex[iexp],globs,FALSE,FALSE,iexp);
      globs->chisq+=computeRight(globs,&ex[iexp]);
      
    }
  *dbg  <<"subtotal chisq=" << globs->chisq << endl;
  solvLin(globs,ex,FALSE);

  for (iexp=1;iexp<=globs->nrExp;iexp++) 
    {      
      GlobExp *expi=&ex[iexp];
      for (i=1; i< expi->nPoints; i++)
	{
	  for (j=1; j<=expi->nvar; j++)
	    {
	      newdS[iexp][i][j]=ex[iexp].dS[i][j];
	    }
	}
      for(i=1; i<=nP; i++)
	{
	newdP[iexp][i]=ex[iexp].dP[i];
	}
    }
  //new increments in newdP and newdS

  //w-estimator
  double sum=0;
  for (iexp=1;iexp<=globs->nrExp;++iexp) 
    {
      GlobExp *expi=&ex[iexp];
      for (i=1; i< expi->nPoints; ++i)
	for (j=1; j<=expi->nvar; ++j)
	  sum+=pow(newdS[iexp][i][j]-(1-t)*uS*ds[iexp][i][j],2);
      for (i=1; i<=nP; ++i)
	{
	  if(globs->doP[i]=='L')
	    sum+=pow(newdP[iexp][i]-(1-t)*dp[iexp][i],2);
	}
    }
  for (i=1; i<=nP; ++i)
    {
      if(globs->doP[i]==TRUE)
	sum+=pow(newdP[1][i]-(1-t)*dp[1][i],2);
    }
	
  free_d3tensor(newdS,1,globs->nrExp,1,maxPoints,1,maxVar);
  free_dmatrix(newdP,1,globs->nrExp,1,nP);
  return 2*sqrt(sum)/(pow(t*normdx,2));
} // subIt


double dampIt(Glob *globs,GlobExp *ex,double **P0, double ***S0, 
	      double **dP, double ***dS,double *uS)
{
  int etaOk=FALSE,i=1;
  double tau=0.95;
  double t;
  double tau_min=0.1;
  double eta0=1.;
  double eta2=1.8; //should be parsed

  double normdx=Norm(globs,ex,dS,dP,*uS);
  // damping algorithm

  if(globs->wquer<0)
    {
      globs->wquer=1.*eta0/normdx;
    }
  
  else
    t=eta0/(globs->wquer*normdx);
  
  if (t > tau) 
    {
      t=1;
    }
  else if(t <= tau_min)
    {
      t=tau_min;
    }
 
  globs->wquer=subIt(globs, ex, t,normdx,P0, S0, dP, dS,*uS);
  if(globs->silent!=TRUE)
    {
      cout << "Start damping:\n";
      cout << "#" << i << " (Predicted)  t = " << t << "\t";
      cout << "chisq=" << globs->chisq;
      cout << endl;
    }

  etaOk=((globs->wquer*t*normdx)<eta2);
  i++;

  while (!etaOk) 
    {
      t=eta0/(normdx*globs->wquer);
        
      if (t > tau) 
	{
	  t=1;
	  etaOk=TRUE;
	  continue;
	}
      else if(t < tau_min)
	{
	      t=tau_min;
	      etaOk=TRUE;
	      if(globs->silent!=TRUE)
		{
		  cout << "#" << i <<  " (Corrected)  t = " << t << "\t";
		  cout << "\nMinimal damping factor tau_min= " << tau_min << " reached.\n";
		}
	      continue;
	}
      globs->wquer=subIt(globs, ex, t, normdx,P0, S0, dP, dS,*uS);
      if(globs->silent!=TRUE)
	{
	  cout << "#" << i <<  " (Corrected)  t = " << t << "\t";
	  cout << "chisq=" << globs->chisq;
	  cout << endl;
	}
      etaOk=((globs->wquer*t*normdx)<eta2);
      i++;
      if(i>=10)
	{
	  cerr << "Too many corrections.\n";
	  throw 1;
	}
    }
  return(t);
}
