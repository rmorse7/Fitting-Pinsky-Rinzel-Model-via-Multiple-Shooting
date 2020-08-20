#include<iostream>
#include<math.h>
#include<stdlib.h>

#include "../def.h"
#include "../model.h"
#include "../nr.h"

using namespace std;

void outFit(GlobExp ex[],Glob *globs)
{
  long i,j;
  long nExp;
  long nvar,npar,nglob=0;
  long ind=1;

  npar=globs->npar;     
  nvar=ex[1].nvar;

  double **errorS=dmatrix(1,globs->nrExp,1,nvar);
  double **errorP=dmatrix(1,globs->nrExp,1,npar);
  

  if(globs->strategy!=2)
    {
      //prepare output
      for(nExp=1;nExp<=globs->nrExp;++nExp) 
	{
	  for(i=1;i<=nvar;i++)
	    errorS[nExp][i]=-1;
	  for(i=1;i<=npar;i++)
	    errorP[nExp][i]=-1;
	  
	  for(i=1;i<=nvar;i++)
	    {
	      if(globs->y0fix[i]!=FALSE)
		{
		  errorS[nExp][i]=sqrt(globs->covar[ind][ind]);
		  ind++;
		}
	    }
	  for(i=1;i<=npar;i++)
	    {
	      if(globs->doP[i]=='L')
		{
		  errorP[nExp][i]=sqrt(globs->covar[ind][ind]);
		  ind++;
		  nglob++;
		}
	    }
	}
      
      for(i=1;i<=npar;i++)
	{
	  if(globs->doP[i]==TRUE)
	    {
	      for(nExp=1;nExp<=globs->nrExp;++nExp) 
		{
		  errorP[nExp][i]=sqrt(globs->covar[ind][ind]);
		}
	      ind++;
	    }
	} 
    }  //end if(globs->strategy..)

  /* print results */
  if(globs->silent!=TRUE)
    {
      cout << "Number of iterations: " << globs->nIter << ", chi^2 = " << globs->chisq << "\n\n";
      cout << "\nBest fit parameters +/- standard errors:\n";
      cout << "----------------------------------------\n\n"; 
      
      ofstream bestfitparamsout;
      bestfitparamsout.open("bestfitparams.dat");

      cout << "Global Parameters:\n";
      for(j=1;j<=npar;j++)
	{ 
	  if(globs->doP[j]!='L')
	    {
	      if(globs->strategy!=2)
		{
		  for(nExp=1;nExp<=globs->nrExp;++nExp) 
		    ex[nExp].errP[j]=errorP[1][j];
		  cout << ParameterNames[j-1];
		  if(errorP[1][j]!=-1)
		    cout << " = " << ex[1].par[j] <<  " +/- " << errorP[1][j] << endl;
		  else
		    cout << " (fixed)\n";
		}
	      else
		cout << ParameterNames[j-1] << " = " << ex[1].par[j] << endl;
		
		
	    }
	}
      cout << endl;

      bestfitparamsout.flush();
      bestfitparamsout.close();

      for(i=1;i<=globs->nrExp;++i) 
	{
	  cout << "Experiment " << i << ":";
	  if(nglob!=0)
	    cout << "\n\nLocal Parameters:\n";
	  else
	    cout << endl;
	  
	  for(j=1;j<=npar;j++)
	    {
	      if(globs->doP[j]=='L')
		{
		  if(globs->strategy!=2)
		    {
		      ex[i].errP[j]=errorP[i][j];
		      cout << ParameterNames[j-1] << " = " << ex[i].par[j];
		      if(errorP[i][j]!=-1)
			cout <<  " +/- " << errorP[i][j] << endl;
		      else
			cout << " (fixed)\n";
		    }
		  else
		    cout << ParameterNames[j-1] << " = " << ex[i].par[j] << endl;
		}
	    }
	  cout << "Initial Values" << endl;
	  for(j=1;j<=nvar;j++)
	    {
	      if(globs->strategy!=2)
		{
		  ex[i].errY0[j]=errorS[i][j];
		  cout << VariableNames[j-1] << " = " << ex[i].yTry[1][j];
		  if(errorS[i][j]!=-1)
		    cout <<  " +/- " << errorS[i][j] << endl;
		  else
		    cout << " (fixed)\n";
		}
	      else
		cout << VariableNames[j-1] << " = " << ex[i].yTry[1][j] << endl;
	    }
	  cout << endl;
	}
    } //end of --> if(silent!=TRUE)

  free_dmatrix(errorS,1,globs->nrExp,1,nvar);
  free_dmatrix(errorP,1,globs->nrExp,1,npar);
}
