#include<iostream>
#include<fstream>
#include<math.h>
#include<stdio.h>

#include "def.h"
#include "model.h"
#include "nr.h"

using namespace std;

//#define PRINTINITVALUES

//Begin definition of module prototypes

void intODE(GlobExp *ex,Glob *globs,int doDerivs,int doPlot,long expNr);

double computeRight(Glob *globs,GlobExp *ex);

void solvLin(Glob *globs,GlobExp *ex,int computeCovar);

double dampIt(Glob *globs,GlobExp *ex,double **P0, double ***S0, 
	      double **dP, double ***dS,double *uS);

//End definition of module prototypes




//      Main routine
//      ************

void fitIt(GlobExp ex[],Glob *globs)
{
  long nExp,i,j,k,l,n;
  long tryOnceMore,ok;
  double objF, sum, scale, lambda;
  long maxPoints=0,maxVar=0;
  double minImprovement=globs->minimp;

  double uS=1.;

  //determine maximum number of shooting intervals and variables across experiments
  for (i=1;i<=globs->nrExp;++i) 
    {
      if (ex[i].nPoints>maxPoints) 
	maxPoints=ex[i].nPoints; 
      if (ex[i].nvar>maxVar)       
	maxVar=ex[i].nvar;
    }
  
  double **dPSave=dmatrix(1,globs->nrExp,1,globs->npar);      
  double ***dSSave=d3tensor(1,globs->nrExp,1,maxPoints,1,maxVar);
  double **P0=dmatrix(1,globs->nrExp,1,globs->npar);
  double ***S0=d3tensor(1,globs->nrExp,1,maxPoints,1,maxVar);
  double **dP=dmatrix(1,globs->nrExp,1,globs->npar);      
  double ***dS=d3tensor(1,globs->nrExp,1,maxPoints,1,maxVar);
  
  *dbg << "Minimum norm for Gauss-Newton update: minImprovement = " << minImprovement << "\n\n";
  
  // GGN iteration
  // do the actual fitting
  
  globs->nIter=1; 

  tryOnceMore=TRUE;
  while (globs->nIter < globs->maxit && tryOnceMore) 
    {
      globs->gnuindex=1;
      if(globs->silent!=TRUE) //print current iteration parameters
	{
	  cout << "Integration # " << globs->nIter << '\n';
	  for(j=1;j<=globs->nrExp;j++)
	    {
	      cout << "Parameter value(s) of experiment " << j << ":\n";
	      for (i=1; i<=globs->npar; ++i) 
		{
		  cout << ex[j].par[i] << " ";
		}
	      cout << endl;
	      cout << "Initial value(s) of experiment " << j << ":\n";
	      for (i=1; i<=ex[j].nvar; ++i) 
		{
		  cout << ex[j].yTry[1][i] << " ";
		}
	      cout << endl;
	    }
	  cout << "\n";
	  cout.flush();
	}
      *dbg << "Integration # " << globs->nIter << '\n';
      *dbg << "Parameter value(s):";
      for (i=1; i<=globs->npar; ++i) 
	{
	  if(globs->doP[i]!='L')
	    *dbg << '\t' << ex[1].par[i];
	}
      *dbg << "\n\n";
      
      // all experiments
      sum=0;
      for (nExp=1;nExp<=globs->nrExp;++nExp) 
	{ //for each expt, compute square error in objF, then sum across experiments
#ifdef PRINTINITVALUES
	  *dbg << "Experiment #" << nExp << "\n";
	  *dbg << "Initial values at mesh points:\n";
	  for (i=1; i<ex[nExp].nPoints; ++i) 
	    {
	      *dbg << ex[nExp].mesh[i] << ':';
	      for (j=1; j<=ex[nExp].nvar; ++j) {
		*dbg << '\t' << ex[nExp].yTry[i][j];
	      }
	      *dbg << '\n';
	    }
	  dbg->flush();
#endif
	  //integrate ODE and plot intermediate results (doPlot = TRUE, doDerivs = TRUE)
	  intODE(&ex[nExp],globs,TRUE,TRUE,nExp);
	  
	  // set up internal data (residues, etc.)
	  objF=computeRight(globs,&ex[nExp]);
	  sum += objF;
	  
	  *dbg << "Experiment #" << nExp << ": objective f ="
	       << objF << "\n";
	  *dbg << "Discrepancies: ";
	  for (i=1; i<ex[nExp].nPoints; ++i) 
	    {
	      for (j=1; j<=ex[nExp].nvar; ++j)
		*dbg << ' ' << ex[nExp].h[i][j];
	      *dbg << ",";
	    }
	  *dbg << "\n\n";
	  dbg->flush();
	} //loop over nExp

      if(globs->silent!=TRUE)
	cout << "Objective function f = " << sum << "\n";
      *dbg << "Global objective function f = " << sum << "\n";
      
      for (nExp=1;nExp<=globs->nrExp;++nExp) //save initial guesses at shooting points
	{
	  for (i=1; i<ex[nExp].nPoints; ++i)
	    {
	      for (j=1; j<=ex[nExp].nvar; ++j) 
		{
		  ex[nExp].yTrySave[i][j] = ex[nExp].yTry[i][j];
		}
	    }
	}
      // solve linearized minimisation problem
      solvLin(globs,ex,FALSE);

      //save old values and new updates
      for (nExp=1;nExp<=globs->nrExp;nExp++) 
	{
	  for (i=1; i<ex[nExp].nPoints; i++)
	    {
	      for (j=1; j<=ex[nExp].nvar; j++)
		{
		  S0[nExp][i][j] = ex[nExp].yTry[i][j];
		  dSSave[nExp][i][j]=ex[nExp].dS[i][j];
		  dS[nExp][i][j]=ex[nExp].dS[i][j];
		}
	    }
	  for (i=1; i<=globs->npar; i++)
	    {
	      P0[nExp][i] = ex[nExp].par[i];
	      dPSave[nExp][i]=ex[nExp].dP[i];
	      dP[nExp][i]=ex[nExp].dP[i];
	    }
	}
      
      uS=1.;
      //damping
      if(!globs->nodamp)
	lambda=dampIt(globs,ex,P0,S0,dP,dS,&uS);
      else
	lambda=1.;

      globs->Lambda=lambda;

      globs->chisq=sum;
      *dbg << "Chisq = " << globs->chisq << endl;
      
      // next guess
      *dbg << "using lambda=" << lambda << "\n";
      
      sum=0;
      for (nExp=1;nExp<=globs->nrExp;++nExp) 
	{
	  for (i=1; i<ex[nExp].nPoints; ++i)
	    {
	      for (j=1; j<=ex[nExp].nvar; ++j) 
		{
		  sum += (lambda*uS*dSSave[nExp][i][j])*(lambda*uS*dSSave[nExp][i][j]);
		}
	    }
	  for (i=1; i<=globs->npar; ++i) 
	    {
	      if(globs->doP[i]=='L')
		sum  += (lambda*dPSave[nExp][i])*(lambda*dPSave[nExp][i]);
	    }
	}
      for (i=1; i<=globs->npar; ++i) 
	{
	  if(globs->doP[i]!='L')
	    sum  += (lambda*dPSave[1][i])*(lambda*dPSave[1][i]);
	}
      
      tryOnceMore = (sqrt(sum)/lambda > minImprovement);
      
      if(uS==0)
	tryOnceMore=TRUE;

      if (tryOnceMore) 
	{
	  // set up new initial values for next iteration
	  for (nExp=1;nExp<=globs->nrExp;++nExp) 
	    {
	      for (i=1; i<ex[nExp].nPoints; ++i) //loop only over the starting points of the MS intervals
		{                                //(not the last point as in setInitialValues  
		  for (j=1; j<=ex[nExp].nvar; ++j) 
		    {
		      
		      ex[nExp].yTry[i][j] = S0[nExp][i][j]+lambda*uS*dSSave[nExp][i][j];
		      //122118 commented out setting parameters to zero
		      // positiveness of initial cond --> test
		      //if(ex[nExp].yTry[i][j]<0.0)
		      //ex[nExp].yTry[i][j]=0;
		    }
		}
	      for (i=1; i<=globs->npar; ++i) 
		{
		  ex[nExp].par[i]=P0[nExp][i]+lambda*dPSave[nExp][i];
		}
	    }
	} 

      // else:
      // keep old values; parCovar is covariance matrix at solution point
      // (update is too small to be worth iterating once more)
      
	++globs->nIter;
	if(globs->silent!=TRUE)
	  {
	    cout << "Norm of update vector: " << sqrt(sum) <<"\n\n";
	    cout.flush();
	  }
	*dbg << "Norm of update vector: " << sqrt(sum) <<"\n\n";
	
	dbg->flush();
	if(globs->wait)
	  {
	    cout << "<Press any key to continue>\n";
	    getchar();
	  }
    }
  
  if (globs->nIter>=globs->maxit) 
    {
      if(globs->silent!=TRUE)
	{
	  cout << "fitit: no convergence after " << globs->nIter << " iterations\n";
	  cout << "PRELIMINARY RESULTS! \n";
	}
      globs->fitConverged=FALSE;
      
    }
  else
    globs->fitConverged=TRUE;

  //to obtain the covariance matrix
  *dbg << "last call of solvLin to obtain the covariance matrix\n";
  if(globs->reg=TRUE)
    {
      globs->minimiser=1;
      //globs->reg=FALSE;
      solvLin(globs,ex,TRUE);
    }
  else
    {
      globs->minimiser=1;
      //globs->reg=FALSE;
      solvLin(globs,ex,TRUE);
    }

  if(!globs->wait && !globs->nowait)
    {
      cout << "<Press any key to continue>\n";
      getchar();
    }
  
  
  if(globs->noGnu==FALSE)
    {
      for(j=1;j<=globs->ngnu;j++)
	pclose(globs->gnuFp[j]);
    }

  free_dmatrix(P0,1,globs->nrExp,1,globs->npar);
  free_d3tensor(S0,1,globs->nrExp,1,maxPoints,1,maxVar);
  free_dmatrix(dPSave,1,globs->nrExp,1,globs->npar);
  free_d3tensor(dSSave,1,globs->nrExp,1,maxPoints,1,maxVar);       
  free_dmatrix(dP,1,globs->nrExp,1,globs->npar);
  free_d3tensor(dS,1,globs->nrExp,1,maxPoints,1,maxVar);
}
