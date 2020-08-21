#include <iostream>
#include <stdlib.h>
#include <fstream>
#include <math.h>
#include<string.h>
#include <getopt.h>
#include "nr.h"

using namespace std;

#define TRUE 1
#define FALSE 0

char usage[]="\n\t    --- createSpline : Creating spline data for diffit ---\n \
 createSpline, v 1.0 25/Jul./2005  Univ. of Freiburg, Author: Martin Peifer\n\n \
      SYNOPSIS \n \t createSpline  [Options]  <File>\n\n \
      DESCRIPTION\n \
      \t -h -? -help \t help\n \
      \t -alpha <val> \t smoothing parameter\n \
      \t -nogcv \t switches cross validation off\n \
      \t -max  <val> \t upper bound for linesearch\n \
      \t -o     <file>\t output file\n"; 

typedef struct 
{
  double *t;
  double *y;
  double *sig;
  double maxLambda;
  double alpha;
  long nKnots;
  int gcv;
  char *outFile;
  char *inFile;
  int outset;
} glob;


void parseTheOpts(int argc, char *argv[],glob *globs)
{
  int longindex,opt;
  long i,j,columns=0,colvor=0;
  double dum;
  char cdum;
  
  ifstream in,inp;
  
  //option definition
  static struct option longopts[]={
    {"help", 0, 0,  'h'},
    {"help", 0, 0,  '?'},
    {"alpha"  , 1, 0,   1 },
    {"nogcv"  , 0, 0,   2 },
    {"max"   , 1, 0,   3 },
    {"o"      , 1, 0,   4 },
    {0, 0, 0, 0}
  };

  //initialise struct glob
  globs->alpha=1;
  globs->gcv=TRUE;
  globs->maxLambda=1e10;
  globs->nKnots=0;
  globs->outset=FALSE;
  
  while ((opt=getopt_long_only(argc,argv,"h?",longopts,&longindex)) != -1)
    {
      switch(opt)
	{
	case 'h':
	  cout << usage << endl;
	  exit(1);
	  break;
	case '?':
	  cerr << usage << endl;
	  exit(1);	
	  break;
	}
    }
  //getFileName
  if((argc-optind)!=1)
    {
      cerr << "No/Too many datafile(s) specified.  \n\n";
      exit(1);
    }
  else
    {
      globs->inFile = new char[strlen(argv[optind])+1];
      strcpy(globs->inFile,argv[optind]);
    }

#ifdef DEBUG
  cout << "Debug ->  Input filename : " << globs->inFile << endl;
#endif 
  optind=0;

  //determine data structure and read data

  in.open(globs->inFile);
  if(!in.is_open())
    {
      cerr << "Input file does not exist.\n";
      exit(1);
    }

  i=0; 
  while(!in.eof())
    {
      cdum=in.peek();
      while((cdum == ' ' || cdum == '\t' || cdum == '\n'))
	{
	  if(cdum=='\n')
	    {
	      if(i==0)
		colvor=columns;
	      i=1;
	      if((columns != 2 && columns != 3)|| colvor != columns)
		{
		  cerr << "Insufficient data format in line " <<  globs->nKnots+1 << ".\n";
		  exit(1);
		}
	      globs->nKnots++;
	      colvor=columns;
	      columns=0;
	    }
	  cdum=in.get();
	  cdum=in.peek();
	}
      columns++;
      in >> dum;    
    }
  columns=colvor;
#ifdef DEBUG
  cout << "Debug -> Columns & Lines : "<< columns << "  " << globs->nKnots << endl;
#endif
  // allocate memory
  in.close();

  globs->t=dvector(1,globs->nKnots);
  globs->y=dvector(1,globs->nKnots);
  globs->sig=dvector(1,globs->nKnots);

  //finally, read the data

  inp.open(globs->inFile);
  for(i=1;i<=globs->nKnots;i++)
    {
      if(columns==2)
	{
	  inp >> globs->t[i] >> globs->y[i];
	  globs->sig[i]=1.;
	}
      else 
	inp>> globs->t[i] >> globs->y[i] >> globs->sig[i];
#ifdef DEBUG
      cout << globs->t[i] << "\t" <<  globs->y[i] << "\t" <<  globs->sig[i] << endl;
#endif
    }
  inp.close();

  //Parse the remaining parameters
  
  while ((opt=getopt_long_only(argc,argv,"h?",longopts,&longindex)) != -1)
    {
      switch(opt)
	{
	case 'h':
	  cout << usage << endl;
	  exit(1);
	  break;
	case '?':
	  cerr << usage << endl;
	  exit(1);	
	  break;
	case 1:
	  globs->alpha=fabs(atof(optarg));
	  break;
	case 2:
	  globs->gcv=FALSE;
	  break;
	case 3:
	  globs->maxLambda=fabs(atof(optarg));
	  break;
	case 4:
	  globs->outFile=new char[strlen(optarg)+1];
	  strcpy(globs->outFile,optarg);
	  globs->outset=TRUE;
	  break;
	default:
	  cerr << endl;
	  cerr << "Parsing parameters produced errors.\n\n";
	  exit(1);
	}
    }

  if(globs->outset==FALSE)
    globs->outFile="spline.dat";
}



main(int argc,char *argv[])
{

  glob globs;
  long i,n;
  double *g,*gam,alpha;

  ofstream out;
  
  parseTheOpts(argc,argv,&globs);
  
  n=globs.nKnots;
  g=dvector(1,n);
  gam=dvector(1,n);
  out.open(globs.outFile);
  alpha=globs.alpha;
  
  if(alpha > globs.maxLambda)
    alpha=globs.maxLambda;

  if(globs.gcv==TRUE)
    splines_gcv(globs.t,globs.y,globs.sig,n,&alpha,g,gam,globs.maxLambda);
  else
    splines(globs.t,globs.y,globs.sig,n,alpha,g,gam);

  //writing output

  out << n <<endl;
  out << globs.t[1] << "\t" << g[1] << "\t" << "0" << endl;
  for(i=2;i<=n-1;i++)
    out << globs.t[i] << "\t" << g[i] << "\t" << gam[i-1] << endl;
  out << globs.t[n] << "\t" << g[n] << "\t" << "0" << endl;

  out.close();
  //free memory
  delete globs.inFile;
  if(globs.outset==TRUE)
    delete globs.outFile;
  free_dvector(globs.y,1,n);
  free_dvector(globs.t,1,n);
  free_dvector(globs.sig,1,n);
  free_dvector(g,1,n);
  free_dvector(gam,1,n);
}
