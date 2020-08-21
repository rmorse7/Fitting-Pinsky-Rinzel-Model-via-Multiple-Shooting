#include<iostream>
#include<fstream>
#include<stdlib.h>

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
//Module modules/freeMem.cc
void freeMem(GlobExp *ex,Glob *globs,int simit);


//End definition of module prototypes


// DEBUG stream
ofstream *dbg;

main(int argc, char *argv[])
{
  long k,i,j;
  char outstr[100];
  Glob globs;
  GlobExp *ex;
  double *t,time,**y,*state;
  long n;

  ex=parseopts(argc,argv,&globs,outstr);

  //open DEBUG stream
  dbg=new ofstream("simit.dbg");
  if (!dbg) 
    {
      cerr << "Error opening DEBUG file.\n";
      exit(1);
    }
  dbg->precision(4);     // set output format for floating point data
  dbg->unsetf(ios::floatfield);

  //print debugging information
  *dbg << DefModelDescription << "\n";
  *dbg << NPARAMS << " parameter(s):";
  for (k=1; k<=NPARAMS; ++k) 
    {
      *dbg << " " << ex[1].par[k];
    }
  *dbg << "\n";
  dbg->flush();
 
   
  //initial values
  *dbg << "\n\n";
  if (ex[1].y0) 
    {
      *dbg << NEQNS << " starting value(s):";
      for (k=1; k<=NEQNS; ++k) 
	*dbg << ' ' << ex[1].y0[k];
      *dbg << "\n";
    } 
  else 
    {
      *dbg << "No starting values specified\n";
      ex[1].y0=dvector(1,ex[1].nvar);
      for (k=1; k<=NEQNS; ++k) 
	ex[1].y0[k]=DefYValues[k-1];
    }
  *dbg << "\n";
  dbg->flush();

  //setting up time points
  if(globs.dt <= 0.)
    {
      cerr << "illegal integration step: " << globs.dt <<  " \n";
      exit(1);
    }
  *dbg << "Integration step: " << globs.dt << endl;

  time=ex[1].fitstart;
  n=0;
  
  while(time < ex[1].fitend-globs.eps)
    {
      time+=globs.dt;
      n++;
    }
  n++;
  *dbg << "Number of time points : " << n << endl;
  
  t=dvector(1,n);
  t[1]=ex[1].fitstart;
  t[n]=ex[1].fitend;

  for(k=2;k<=n-1;k++)
    t[k]=t[k-1]+globs.dt;
#ifdef PLOTTIMEPOINTS
  *dbg << "Time points : ";
  for(k=1;k<=n;k++)
    {
      *dbg << t[k];
      if(k!=n)
	*dbg << ";";
    }
  *dbg << endl;
#endif
  dbg->flush();

  y=dmatrix(1,n,1,ex[1].nvar);
  state=dvector(1,ex[1].nvar);
  for(k=1;k<=ex[1].nvar;k++)
    state[k]=ex[1].y0[k];

  initialise(ex,&globs,TRUE);

  try
    {
      tabulateValues(&globs,&ex[1],t[1],t[n],t,n,state,y,NULL,NULL,NULL,NULL);
    }
  catch (int i)
    {
      exit(1);
    }

  outSimit(ex,&globs,t,n,y);
  cout << "creating " << ex[1].fileName << endl;

  freeMem(ex,&globs,TRUE);
}
