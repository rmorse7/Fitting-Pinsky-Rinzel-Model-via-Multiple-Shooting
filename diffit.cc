/*! \file diffit.cc \brief Root program for \b diffit (contains function main)

   Diffit is has a modular structure. This file can be considered as
   root module for the entire C stand-alone version. A matlab interface
   is also available under: diffit_mex.cc. The structure of the 
   matlab interface is similar than this file, apart from parsing 
   the data and the program parameters. Note that all global 
   variables are capsuled into two C structures named Glob and GlobExp.
   The first structure Glob stores all experiment unspecific parameters
   such as the used integrator, or the integration accuracy, etc. Experiment
   specific program variables are held in GlobExp, these are, e.g., 
   the number of multiple shooting intervals. Modules called by 
   function \b main is in order of occurrence:

   -# parseopts(): Parsing all program parameters either from the 
                   command line or a file.
   -# readData(): Reading the data.
   -# setMesh(): This module sets up the multiple shooting intervals
                 for each observation and experiment.
   -# initialise(): Allocates memory and initializes vectors and matrices.
   -# simInit(): If some state variables are not directly observed
                 (non trivial observation function), this function 
		 integrates a trajectory using the initial guess to
		 set up the initial guess for the state variables at 
		 multiple shooting intervals. This is activated by the
		 command line argument \b "-siminit". If not used, the 
		 state variables are initialized with 1. Since the 
		 integrated trajectory for the initial guess is continuous,
		 rather than discontinuous the \b "-pert <value>" perturbs
		 this guess using mean zero Gaussian noise. The standard
		 deviation of this noise can be chosen by the value 
		 of the -pert argument.
   -# globOpt(): Experimental !! Global optimizer instead of multiple shooting.
                 \b Currently \b not \b working \b properly.
   -# fitIt(): \b General \b fitting \b function, all the numerics is done in
                here. 
   -# outFit(): This routine writes the output such as parameter estimates
                and its standard error.
   -# freeMem(): Deallocation of vectors and matrices.
*/


#include<iostream>
#include<fstream>
#include<stdlib.h>

#include "def.h"
#include "model.h"
#include "nr.h"

using namespace std;

//Begin definition of module prototypes
//Module modules/parse.cc
GlobExp *parseopts(int argc, char *argv[],Glob *globs,char *outstr);
//Module modules/readData.cc
void readData(GlobExp *ex,Glob *globs,long expNr);
//Module modules/setMesh.cc
void setMesh(GlobExp *ex,Glob *globs,long expNr);
//Module modules/outFit.cc
void outFit(GlobExp ex[],Glob *globs);
//Module modules/freeMem.cc
void freeMem(GlobExp *ex,Glob *globs,int simit);
//Module modules/initialise.cc
void initialise(GlobExp ex[],Glob *globs,int simit);
//Module modules/simInit.cc
void simInit(GlobExp *ex,Glob *globs);

//End definition of module prototypes

//The numerics subprogram
void fitIt(GlobExp ex[],Glob *globs);
//Global optimiser stuff
void globOpt(GlobExp ex[],Glob *globs);

// DEBUG stream
//! debug stream (C++)
ofstream *dbg;

/*! function main  */
main(int argc, char *argv[])
{
  long k,i,j;
  char outstr[100];
  Glob globs; //Glob defined def.h
  GlobExp *ex; //GlobExp defined in def.h
  
  //  parseopt
  ex=parseopts(argc,argv,&globs,outstr);

  //open DEBUG stream
  dbg=new ofstream("diffit.dbg");
  if (!dbg) 
    {
      cerr << "Error opening DEBUG file.\n";
      exit(1);
    }
  dbg->precision(4);     // set output format for floating point data
  dbg->unsetf(ios::floatfield);

  //print debugging information
  *dbg << DefModelDescription << "\n";
  if (globs.nrExp==1) 
    *dbg << "1 experiment\n";
  else
    *dbg << globs.nrExp << " experiments\n";
   *dbg << "\n";

   //writes the list of parameters being optimized 
   for(i=1;i<=globs.nrExp;i++)
     {
       *dbg << "Experiment: " << i << "\n";
       *dbg << NPARAMS << " parameter(s):";
       for (k=1; k<=NPARAMS; ++k) 
	 *dbg << " " << ex[i].par[k];
       *dbg << "\n";
     }
   
  dbg->flush();

  //READ Data for each expt and setup multiple shooting intervals
  for(i=1;i<=globs.nrExp;i++)
    {
      // readData()
      readData(ex,&globs,i);
      //set mesh
      setMesh(ex,&globs,i);

#ifdef PRINTDATA
      *dbg << "\nData: \n";
      for (j=1; j<=ex[i].nMeasure; j++) {
	*dbg << ex[i].xMeasure[j];
	for (k=1; k<=ex[i].nobs; k++)
	  *dbg << "\t" << ex[i].yMeasure[j][k] << "\t" << ex[i].sigma[j][k];
	*dbg << "\n";
      }
      *dbg << "\n";
#endif
      *dbg << "Mesh:";
      for (k=1;k<=ex[i].nPoints;++k) 
	*dbg << " " << ex[i].mesh[k];
      *dbg << "\n";
    } //end readData and setMesh

  //initial values; Note: only those specified at the command line
  for(i=1;i<=globs.nrExp;i++)
    {
      *dbg << "\nExperiment " << i << ":\n";
      if (ex[i].y0) 
	{
	  *dbg << NEQNS << " starting value(s):";
	  for (k=1; k<=NEQNS; ++k) 
	    *dbg << ' ' << ex[i].y0[k];
	  *dbg << "\n";
	} 
      else 
	{
	  *dbg << "No starting values specified\n";
	}
      *dbg << "\n";
      dbg->flush();
    }
  
   initialise(ex,&globs,FALSE);
   
   //simulate initial state
   if(globs.simInit==TRUE)
     {
       for (i=1;i<=globs.nrExp;++i) 
	 simInit(&ex[i],&globs);
     }
   
   // starting of the numerics
   //*************************
      
   if(globs.strategy==2)
     globOpt(ex,&globs);
   else
     {
       try 
	 {	      
	    fitIt(ex,&globs);
	 }
       catch(int i)  //Exception handling
	 {
	   if(globs.noGnu==FALSE)
	     {
	       for(j=1;j<=globs.ngnu;j++)
		 pclose(globs.gnuFp[j]);
	     }       
	   exit(1);
	 }
     }

   
   //*************************

  //Output after convergence
  outFit(ex,&globs);
  if(!globs.noGnu)
    system("rm -f gnuout.dat");

  freeMem(ex,&globs,FALSE);
}