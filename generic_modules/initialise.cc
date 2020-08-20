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

void setInitialValues(GlobExp *ex, int nP, double *parameters, double *yValues);

//End definition of module prototypes

//Tolerance of unsatisfied equality constraints
#define EPSTOL 1e-8

//Check if initial values are compatible with constrains
int initValuesOK(Glob *globs,GlobExp *ex)
{
  int i,good=TRUE;
  
  if (ex->me > 0) // me is the number of equality constraints 
    {
      //call to function R2, in model.cc. This automatically generated function computes ex->r2[1..me].
      //The equality constraints are satisfied if each component of ex->r2 is equal to zero
      R2(globs,ex, FALSE);
      for (i=1; i<=ex->me; ++i)
	if (good) 
	  good = fabs(ex->r2[i]) < EPSTOL; //check equality constraints within tolerance, fabs: absolute value, float
    }
  
  if (ex->mg>0) // mg is the number of inequality constraints
    {
      //call to function R3, in model.cc. This automatically generated function computes ex->r3[1..mg].
      //The inequality constraints are statisfied if each component of ex->r3 is non-negative.
      R3(globs,ex, FALSE); 
      for (i=1; i<=ex->mg; ++i) 
	good = good && (ex->r3[i]>=0);
    }
  
  return(good);
}

//Gnuplotting
void setupGnuFp(Glob *globs, GlobExp ex[])
{
  long gnu_n=20,ncomp,k;
  int fullx=700; //max. geometry
  int fully=500;
  int nwi,dx=fullx/2,dy,xpos,ypos;
  char outstr[100];


  ncomp=0;
  for(k=1;k<=globs->nrExp;k++)
    {
      ncomp+=ex[k].nobs;
    }

  // determine nr of windows in y-direction, max of 20 windows
  gnu_n = LMIN (gnu_n, ncomp); //minimum of 2 long integers, see numerical recipes header, nr.h. 
  globs->ngnu=gnu_n;

  // allocate streams
  globs->gnuFp = new FILE *[gnu_n+1];
  
  if((gnu_n % 2)==1)
    nwi=(gnu_n+1)/2;
  else
    nwi=gnu_n/2;
  dy=fully/nwi;
 
  for(k=1;k<=gnu_n;k++) // sets up a 2 columns by gnu_n/2 plotting windows on the screen
    {
      if(k==1)
	{
	  xpos=0;
	  ypos=0;
	}
      else if((k % 2)==0)
	{
	  xpos=dx+10;
	}
      else
	{
	  xpos=0;
	  ypos+=dy+25;
	}
      sprintf(outstr,"gnuplot -noraise -geometry %ix%i+%i+%i -title \"%i\"\0",dx,dy,xpos,ypos,k);
      globs->gnuFp[k]=popen(outstr,"w");  //pipe for kth plot
    }
} // setupGnuFp


//Initialise all vectors, matrices, tensors, etc. 
//used during the iterations
void initialise(GlobExp ex[],Glob *globs,int simit)
{
  long i,j,k,nExp;
  long nrExp=globs->nrExp;
  
  ifstream inSplines;

  if(simit==FALSE) //simit is set to false in the call from diffit
    {
      for (nExp=1;nExp<=nrExp;++nExp) 
	{ 
	  
	  //non-local parameters are equal for each experiment
	  for(i=1;i<=globs->npar;i++)
	    {
	      if(globs->doP[i]==TRUE)
		//note: ex[1].par[i] is defined in parse.cc, where it is set to the default values selected
		//in model.cc by DefParameters
		ex[nExp].par[i]=ex[1].par[i];
	    }
	  
	  // set ex->me and ex->mg from the function setNrConstraints given in model.cc
	  setNrConstraints(&ex[nExp],globs->npar,ex[nExp].par);
	  
	  long nPoints=ex[nExp].nPoints; //number of multiple shooting intervals
	  long nP=globs->npar;
	  long nvar=ex[nExp].nvar, nobs=ex[nExp].nobs;
	  long nMeas=ex[nExp].nMeasure; //number of time points
	  long  me=ex[nExp].me, mg=ex[nExp].mg;
	  
#ifdef PRINTINITVALUES
	  *dbg << "initialize: exp. #" << nExp << ", nPoints=" << nPoints
	       << ", nvar=" << nvar << ", nobs=" << nobs << "\nnMeas="
	       << nMeas << ", nVal=" << nVal << ", me=" << me
	       << ", mg=" << mg << ", nP=" << nP << "\n";
	  dbg->flush();
#endif
          //std. error -> output
          ex[nExp].errP=dvector(1,nP);
          ex[nExp].errY0=dvector(1,nvar);	  
	  // initial guesses at mesh points (initial value for each shooting interval)
	  ex[nExp].yTry=dmatrix(1,nPoints,1,nvar); 
	  ex[nExp].yTrySave=dmatrix(1,nPoints,1,nvar); 
	  // computed values at mesh points   
	  ex[nExp].yComp=dmatrix(1,nPoints,1,nvar);
	  // computed values at measuring points
	  ex[nExp].yPred=dmatrix(1,nMeas,1,nobs);
	  
	  ex[nExp].h=dmatrix(1,nPoints,1,nvar);      // discrepancies at multiple shooting interval boundaries
	  ex[nExp].residues=dmatrix(1,nMeas,1,nobs); // residues
	  
	  // derivatives
	  // d(yComp[i][j]) / d(yTry[i-1][k])
	  ex[nExp].dyds=d3tensor(1,nPoints,1,nvar,1,nvar);
	  // d(yComp[i][j]) / d(p[k])
	  ex[nExp].dydp=d3tensor(1,nPoints,1,nvar,1,nP);
	  // d(yPred[i][j]) / d(yTry[**][k])
	  ex[nExp].dmds=d3tensor(1,nMeas,1,nobs,1,nvar);
	  // d(yPred[i][j]) / d(p[k])
	  ex[nExp].dmdp=d3tensor(1,nMeas,1,nobs,1,nP);
	  if (me>0) 
	    {
	      ex[nExp].r2=dvector(1,me); // equality constraints
	      // d(yR2[i]) / d(yTry[j][k])
	      ex[nExp].dR2ds=d3tensor(1,me,1,nPoints,1,nvar);
	      for(i=1;i<=me;i++)
		for(j=1;j<=nPoints;j++)
		  for(k=1;k<=nvar;k++)
		    ex[nExp].dR2ds[i][j][k]=0.;
	      
	      // d(yR2[i]) / d(p[j])
	      ex[nExp].dR2dp=dmatrix(1,me,1,nP);
	      for(i=1;i<=me;i++)
		for(j=1;j<=nP;j++)
		  ex[nExp].dR2dp[i][j]=0.;
	    } 
	  else 
	    {
	      ex[nExp].r2=NULL; ex[nExp].dR2ds=NULL; ex[nExp].dR2dp=NULL;
	    }
	  if (mg>0) 
	    {
	      ex[nExp].r3=dvector(1,mg); // inequality  constraints
	      // d(yR3[i]) / d(yTry[j][k])
	      ex[nExp].dR3ds=d3tensor(1,mg,1,nPoints,1,nvar);
	      for(i=1;i<=mg;i++)
		for(j=1;j<=nPoints;j++)
		  for(k=1;k<=nvar;k++)
		    ex[nExp].dR3ds[i][j][k]=0.;
	      // d(yR3[i]) / d(p[j])
	      ex[nExp].dR3dp=dmatrix(1,mg,1,nP);
	      for(i=1;i<=mg;i++)
		for(j=1;j<=nP;j++)
		  ex[nExp].dR3dp[i][j]=0.;
	    } 
	  else 
	    {
	      ex[nExp].r3=NULL; ex[nExp].dR3ds=NULL; ex[nExp].dR3dp=NULL;
	    }
	  
	  // compute yTry from ex (all ex read from file still present)
	  setInitialValues(&ex[nExp],globs->npar,ex[nExp].par,ex[nExp].y0);
	  
	  // check if initial values are compatible with constraints
	  // CURRENTLY DISABLED (1. Dez. 2004)
// 	  if (!initValuesOK(globs,&ex[nExp])) 
// 	    {
// 	      cerr << "Experiment #" << nExp 
// 		   << ": initial values are not compatible with constraints\n";
// 	      cerr << "R2:";
// 	      for (i=1;i<=ex[nExp].me;++i) 
// 		cerr << "\t" << ex[nExp].r2[i];
// 	      cerr << "\nR3:";
// 	      for (i=1;i<=ex[nExp].mg;++i) 
// 		cerr << "\t" << ex[nExp].r3[i];
// 	      cerr << "\n";
// 	      exit(1);
// 	    }
	  
#ifdef PRINTEX
	  *dbg << "Ex used for fitting: \n\n";
	  for (i=1; i<=ex[nExp].nMeasure; ++i) 
	    {
	      *dbg << ex[nExp].xMeasure[i];
	      for (j=1; j<=ex[nExp].nobs; ++j)
		*dbg << "\t" << ex[nExp].yMeasure[i][j] << "\t" << ex[nExp].sigma[i][j];
	      *dbg << '\n';
	    }
	  *dbg << '\n';
	  dbg->flush();
#endif
	  //Allocate memory for condensation
	  //least squares
	  ex[nExp].ua=dvector(1,nMeas*nobs);
	  ex[nExp].Ea=dmatrix(1,nMeas*nobs,1,nvar);
	  ex[nExp].Pa=dmatrix(1,nMeas*nobs,1,nP);
	  //equality constraints
	  ex[nExp].ue=dvector(1,me);
	  ex[nExp].Ee=dmatrix(1,me,1,nvar);
	  ex[nExp].Pe=dmatrix(1,me,1,nP);
	  //inequality constraints
	  ex[nExp].ug=dvector(1,mg);
	  ex[nExp].Eg=dmatrix(1,mg,1,nvar);
	  ex[nExp].Pg=dmatrix(1,mg,1,nP);
	  
	  //update steps
	  ex[nExp].dS=dmatrix(1,nPoints,1,nvar);
	  ex[nExp].dP=dvector(1,nP);
	}// matches: for(nExp=1;nExp<=nrExp;++nExp) loop over experiments

      //covar in solvLin allocated
      globs->covar=NULL;
    } // matches with if(simit==FALSE)...

  //Initialise GNUplot
  if(globs->noGnu==FALSE)
    setupGnuFp(globs,ex);
  
  globs->cond=0; //condition number of the linearized system
  globs->Lambda=1; //damping parameter
  
  char line[1000];
  //Initialize Splines
  if(globs->initSpline==TRUE)
    {
      for (nExp=1;nExp<=nrExp;++nExp) 
	{
	  //see def.h for significance 
	  ex[nExp].splineNodes=(double **) malloc((size_t)(NSPLINES+1)*sizeof(double*));
	  ex[nExp].splineY=(double **) malloc((size_t)(NSPLINES+1)*sizeof(double*));
	  ex[nExp].splineGam=(double **) malloc((size_t)(NSPLINES+1)*sizeof(double*));
	  ex[nExp].nNodes=(long unsigned*) malloc((size_t)(NSPLINES+1)*sizeof(long unsigned*));
	  
	  for(i=1;i<=NSPLINES;i++)
	    {
	      j=0;
	      while(ex[nExp].splineFile[i][j]!='\0')
		{
		  line[j]=ex[nExp].splineFile[i][j];
		  j++;
		}
	      line[j]='\0';
	      inSplines.open(line);
	      if(!inSplines.is_open())
		{
		  cerr << "Cannot open " << ex[nExp].splineFile[i] << ".\n";
		  exit(1);
		}
	      //read in number of nodes
	      inSplines >>  ex[nExp].nNodes[i];
	      ex[nExp].splineNodes[i]=(double *) malloc((size_t)(ex[nExp].nNodes[i]+1)*sizeof(double));
	      ex[nExp].splineY[i]=(double *) malloc((size_t)(ex[nExp].nNodes[i]+1)*sizeof(double));
	      ex[nExp].splineGam[i]=(double *) malloc((size_t)(ex[nExp].nNodes[i]+1)*sizeof(double));
	      for(j=1;j<= ex[nExp].nNodes[i];j++)
		{
		  inSplines >>  ex[nExp].splineNodes[i][j];
		  inSplines >>  ex[nExp].splineY[i][j];
		  inSplines >>  ex[nExp].splineGam[i][j];	      
		}
	      inSplines.close();
	    } //loop over splines
	} //loop over experiments
    }//if initSpline == true
}
