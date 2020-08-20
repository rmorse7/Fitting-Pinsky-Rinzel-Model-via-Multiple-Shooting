#include<string.h>
#include<stdlib.h>
#include<stddef.h>
#include<iostream>
#include<fstream>
#include<getopt.h>
#include<math.h>

#include "../nr.h"
#include "../def.h"
#include "../model.h"

#define MAXLENGTH 1000

using namespace std;

char usage[]="\n\t --- simit : Simulating differential equations ---\n \
      simit, v 2.3 12/Apr./2006  Univ. of Freiburg, Author: Martin Peifer\n \n \
      SYNOPSIS \n \t simit [Options]  <File>\n\n \
      DESCRIPTION\n \
      \t -h -? -help \t help\n \
      \t -f <File> \t parameter file for simulation \n \
      \t -x0 <val> \t integration starts at time <val>\n \
      \t -x1 <val> \t integration stops at time <val>\n \
      \t -dt <val> \t integration step\n \
      \t -p <list> \t parameters, e.g. -p 0.5,2.3\n \
      \t -y0 <list> \t initial values, e.g. -y0 1,2\n \
      \t -nobs <num> \t number of observed quantities\n \
      \t -eps <val> \t integration accuracy\n \
      \t -sig <val> \t noise level\n \
      \t -int <num> \t integrator: 1-> Runge-Kutta (default)\n \
      \t            \t             2-> CVODES \n \
      \t            \t             3-> ODESSA \n \
      \t -maxstp \t maximal integration steps\n \
      \t -spline <list>  specifies spline data for non-autonomous ode's \n";

long get_list(long n,double *arg,char *string,char *desc)
{
  long k=1,l=0,ind=0;
  char sep=',';
  char *dum=new char[strlen(string)+1];
  while(string[ind]!='\0')
    {
      if(k>n)
	{
	  cerr << "Too many arguments in" << desc << "list." << endl;
       	  exit(1);
	}
      dum[l]=string[ind];
      
      if(dum[l]==sep)
	{
	  dum[l]='\0';
	  arg[k]=atof(dum);
	  l=0;
	  k++;
	  ind++;
	}
      else
	{
	  ind++;
	  l++;
	}
    }
  dum[l]='\0';
  arg[k]=atof(dum);

  delete dum;
  return(k);
} 


GlobExp *parseopts(int argc, char *argv[],Glob *globs,char *outstr)
{
  int longindex,opt;
  int parlistSpecified=FALSE;
  long k,l,nExp,sp;
  char name[MAXLENGTH],line[MAXLENGTH];

  //option definition
  static struct option longopts[]={
    {"help", 0, 0,  'h'},
    {"help", 0, 0,  '?'},
    {"x0"  , 1, 0,   1 },
    {"x1"  , 1, 0,   2 },
    {"dt"  , 1, 0,   3 },
    {"f"    ,1, 0,   4 },
    {"p"    ,1, 0,   5 },
    {"y0"   ,1, 0,   6 },
    {"nobs" ,1, 0,   7 },
    {"eps"  ,1, 0,   8 },
    {"sig"  ,1, 0,   9 },
    {"maxit",1, 0,  12 },
    {"int"  ,1, 0,  15 },
    {"maxstp",1,0,  17 },
    {"spline",1,0,  21 },
    {0, 0, 0, 0}
  };
  //initialise some global parameters in globs
  globs->noGnu=FALSE;
  globs->eps=1e-6;
  globs->npar=NPARAMS;
  globs->noMeasurements=FALSE;
  globs->doP=ivector(1,globs->npar);
  for(k=1;k<=globs->npar;k++)
    globs->doP[k]=TRUE;
  globs->maxit=1000;
  globs->gnuFp=NULL;
  globs->wait=FALSE;
  globs->usesig=FALSE;
  globs->integrator=1;
  globs->stiff=TRUE;
  globs->maxstp=5000;
  globs->minimp=0.05;
  globs->nowait=FALSE;
  globs->elastic=1.;
  globs->reg=FALSE;
  globs->epsilon=1e-10;
  globs->lambda=1e6;
  globs->dt=0.1;
  globs->sig=0.;
  globs->nrExp=1;
  globs->initSpline=TRUE;

  while ((opt=getopt_long_only(argc,argv,"h?",longopts,&longindex)) != -1)
    {
      switch(opt)
	{
	case 4:
	  parlistSpecified=TRUE;
	  strncpy(name,optarg,MAXLENGTH);
	  break;
	}
    }
  optind=0;
  
  GlobExp *ex;
  if(parlistSpecified==FALSE)
    {
      globs->nrExp=1;
      ex=new GlobExp[globs->nrExp+1];  
      //initialise some global parameters
      ex[1].fitstart=0;
      ex[1].fitend=1;
      ex[1].nobs=NOBS;
      ex[1].nvar=NEQNS;
      ex[1].y0=NULL;
      ex[1].par=dvector(1,NPARAMS);
      for(k=1;k<=globs->npar;k++)
	ex[1].par[k]=DefParameters[k-1];
    }
  else //if parlist specified
    {
      ifstream in;     
      ofstream out;
      char str[MAXLENGTH];

      in.open(name);

      if(!in)
	{
	  cerr << "parameter list " << name << "  not found\n";
	  exit(1);
	}
      
      //preprocessing of data list
      globs->nrExp=1;
      ex=new GlobExp[globs->nrExp+1];   
      //initialise some global parameters
      ex[1].fitstart=0;
      ex[1].fitend=1;
      ex[1].nobs=NOBS;
      ex[1].nvar=NEQNS;
      ex[1].par=dvector(1,NPARAMS);
      for(k=1;k<=globs->npar;k++)
	ex[1].par[k]=DefParameters[k-1];

      sprintf(str,"%s.%d",name,globs->nrExp);
      out.open(str);
      while(!in.eof())
	{       
	  in.getline(str,MAXLENGTH,'\n');
	  if(str[0]!='#')
	    {
	      out << str << endl;
	    }
	}
      out.flush();
      out.close();
      in.close();
      
      long _argc;
      char *_argv[MAXLENGTH];

      for(k=0;k<=MAXLENGTH;k++)
	_argv[k]=(char *)malloc((size_t) (500*sizeof(char)));

      
      ifstream inp;
      
      _argc=0;
      sprintf(str,"%s.1",name);
      inp.open(str);
      
      
      while(!inp.eof())
	{
	  inp >> _argv[_argc+1];
	  _argc++;
	}
      inp.close();
      
  
      //PARSING LIST
      optind=0;
      while((opt=getopt_long_only(_argc,_argv,"h?",longopts,&longindex)) != -1)
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
	      ex[1].fitstart=atof(optarg);
	      break;
	    case 2:
	      ex[1].fitend=atof(optarg);
	      break;
	    case 3:
	      globs->dt=fabs(atof(optarg));
	      break;
	    case 4:
	      break;
	    case 5:
	      get_list(globs->npar,ex[1].par,optarg," parameter ");
	      break;
	    case 6:
	      ex[1].y0=dvector(1,ex[1].nvar);
	      for(k=1;k<=ex[1].nvar;k++)
		ex[1].y0[k]=0.;
	      k=get_list(ex[1].nvar,ex[1].y0,optarg," initial values ");
	      if(k!=ex[1].nvar)
		{
		  cerr << ex[1].nvar << " initial values required.\n";
		  exit(1);
		}
	      break;
	    case 7:
	      ex[1].nobs=atol(optarg);
	      break;
	    case 8:
	      globs->eps=fabs(atof(optarg));
	      break;
	    case 9:
	      globs->sig=fabs(atof(optarg));
	      break;
	    case 11:
	      if(strlen(optarg)!=globs->npar)
		{
		  cerr << strlen(optarg) << " parameter(s) specified instead of " <<  globs->npar << ".\n";
		  exit(1);
		}
	      for(k=0;k < globs->npar;k++)
		{
		  if(optarg[k]=='0')
		    globs->doP[k+1]=FALSE;
		  else
		    globs->doP[k+1]=TRUE;
		}
	      break;
	    case 12:
	      globs->maxit=abs(atol(optarg));
	      break;
	    case 15:
	      globs->integrator=abs(atoi(optarg));
	      break;
	    case 17:
	      globs->maxstp=abs(atoi(optarg));
	      break;
	    case 21:
	      k=0;
	      l=1;
	      ex[1].splinesDefined=TRUE;
	      while(optarg[k]!='\0')
		{
		  if(optarg[k]==',')
		    l++;
		  k++;
		}
	      if(l!=NSPLINES)
		{
		  cerr << "Incompatible number of splines.\n";
		  exit(1);
		}
	      ex[1].splineFile=new string[l+1];
	      l=1;
	      sp=0;
	      k=0;
	      while(optarg[k]!='\0')
		{
		  if(optarg[k]==',')
		    {
		      line[sp]='\0';
		      ex[1].splineFile[l]=(string)line;
		      sp=0;
		      l++;
		    }
		  else
		    {
		      line[sp]=optarg[k];
		      sp++;
		    }
		  k++;
		}
	      line[sp]='\0';
	      ex[1].splineFile[l]=(string)line;
	      break;
	    default:
	      cerr << endl;
	      cerr << "Parsing command line options produced errors \n\n";
	      exit(1);
	    }
	}
      nExp=1;
      //Parsing file name
      if((_argc-optind)!=1)
	{
	  cerr << "No/Too many datafile(s) specified for experiment " << nExp << " .  \n\n";
	  exit(1);
	}
      else
	{
	  ex[nExp].fileName=new char[strlen(_argv[optind])+1];
	  strcpy(ex[nExp].fileName,_argv[optind]);
	}
      //checking if x0 < x1
      if(ex[nExp].fitstart >= ex[nExp].fitend)
	{
	  cerr << "x0 must be smaller than x1.\n";
	  exit(1);
	} 
      
      sprintf(str,"rm -f %s.*",name);
      system(str);
    }

  optind=0;
  
  while((opt=getopt_long_only(argc,argv,"h?",longopts,&longindex)) != -1)
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
	  ex[1].fitstart=atof(optarg);
	  break;
	case 2:
	  ex[1].fitend=atof(optarg);
	  break;
	case 3:
	  globs->dt=fabs(atof(optarg));
	  break;
	case 4:
	  break;
	case 5:
	  get_list(globs->npar,ex[1].par,optarg," parameter ");
	  break;
	case 6:
	  ex[1].y0=dvector(1,ex[1].nvar);
	  for(k=1;k<=ex[1].nvar;k++)
	    ex[1].y0[k]=0.;
	  k=get_list(ex[1].nvar,ex[1].y0,optarg," initial values ");
	  if(k!=ex[1].nvar)
	    {
	      cerr << ex[1].nvar << " initial values required.\n";
	      exit(1);
	    }
	  break;
	case 7:
	  ex[1].nobs=atol(optarg);
	  break;
	case 8:
	  globs->eps=fabs(atof(optarg));
	  break;
	case 9:
	  globs->sig=fabs(atof(optarg));
	  break;
	case 11:
	  if(strlen(optarg)!=globs->npar)
	    {
	      cerr << strlen(optarg) << " parameter(s) specified instead of " <<  globs->npar << ".\n";
	      exit(1);
	    }
	  for(k=0;k < globs->npar;k++)
	    {
	      if(optarg[k]=='0')
		globs->doP[k+1]=FALSE;
	      else
		globs->doP[k+1]=TRUE;
	    }
	  break;
	case 12:
	  globs->maxit=abs(atol(optarg));
	  break;
	case 15:
	  globs->integrator=abs(atol(optarg));
	  break;
	case 17:
	  globs->maxstp=abs(atoi(optarg));
	  break;
	case 21:
	  k=0;
	  l=1;
	  ex[1].splinesDefined=TRUE;
	  while(optarg[k]!='\0')
	    {
	      if(optarg[k]==',')
		l++;
	      k++;
	    }
	  if(l!=NSPLINES)
	    {
	      cerr << "Incompatible number of splines.\n";
	      exit(1);
	    }
	  ex[1].splineFile=new string[l+1];
	  l=1;
	  sp=0;
	  k=0;
	  while(optarg[k]!='\0')
	    {
	      if(optarg[k]==',')
		{
		  line[sp]='\0';
		  ex[1].splineFile[l]=(string)line;
		  sp=0;
		  l++;
		}
	      else
		{
		  line[sp]=optarg[k];
		  sp++;
		}
	      k++;
	    }
	  line[sp]='\0';
	  ex[1].splineFile[l]=(string)line;
	  break;
	default:
	  cerr << endl;
	  cerr << "Parsing command line options produced errors \n\n";
	  exit(1);
	}
    }
  //checking if x0 < x1
  if(ex[1].fitstart >= ex[1].fitend)
    {
      cerr << "x0 must be smaller than x1.\n";
      exit(1);
    }
  
  //Parsing file name for a single experiment
  if((argc-optind)!=1 && parlistSpecified==FALSE)
    {
      cerr << "No/Too many datafile(s) specified.  \n\n";
      exit(1);
    }
  else if(parlistSpecified==FALSE)
    {
      ex[1].fileName=new char[strlen(argv[optind])+1];
      strcpy(ex[1].fileName,argv[optind]);
    }

    
  if(NSPLINES!=0)
    {
       if(ex[1].splinesDefined==FALSE)
	    {
	      cerr << "Please specify spline data.\n";
	      exit(1);
	    }
    }
  

  return(ex);
}

