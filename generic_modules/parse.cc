/*!\file parse.cc \brief \b Module for parsing command line arguments

  This module parses the command line arguments. It is based on 
  the standard C routine getopt. The structure \b option further
  below defines the command line arguments and if they require an
  argument or not. A documentation of this structure and further details
  of the getopt routine comes with documentation of the GNU compiler
  (e.g. use the linux/unix command:man 3 getopt).
*/

#include<string.h>
#include<stdlib.h>
#include<stddef.h>
#include<iostream>
#include<fstream>
#include<getopt.h>
#include<string>
#include<math.h>

#include "../nr.h"
#include "../def.h"
#include "../model.h"

//! maximal length of input line
#define MAXLENGTH 1000

using namespace std;

/*! \brief description of command line arguments 
    
    can be displayed by using command line argument -h 
*/
char usage[]="\n\t --- diffit : Fitting parameters in differential equations ---\n \
      diffit, v 2.3 12/Apr./2006  Univ. of Freiburg, \n \
              Authors: Martin Peifer, Kilian Bartholome\n\n				\
      SYNOPSIS \n \t diffit [Options]  <File>\n\n \
      DESCRIPTION\n \
      \t -h -? -help \t help\n \
      \t -x0 <val> \t fit starts at time <val>\n \
      \t -x1 <val> \t fit stops at time <val>\n\
      \t -nms <num> \t number of multiple shooting intervals, default: 1 \n\
      \t -tms <list> \t time boundaries of multiple shooting intervals, including the first and last, e.g. -tms 0.0,5.0,40.0,70.0,80.0\n\
      \t\t              Must come after -nms. If the tms option is absent, the boundaries will be evenly spaced as before.\n\
      \t -f <File> \t option file/multi-experiment definition file\n\
      \t\t Global options first; then @ on one line followed by specific experiments\n\
      \t\t Experiments are separated by @ on one line; comments start with #\n\
      \t\t Experiment specific parameters are : x0,x1,nms,y0,nobs\n\
      \t\t                                      p,spline\n\
      \t -p <list> \t initial parameters, e.g. -p 0.5,2.3\n \
      \t -y0 <list> \t initial values, e.g. -y0 1,2\n\
      \t -nobs <num> \t number of observed quantities\n\
      \t -eps <val> \t integration accuracy\n\
      \t -nognu \t turns gnuplot animation off\n\
       \t -savegnupng <num> \t send gnuplot output in png format to an existing png subdirectory\n\
      \t                   \t             1-> Save each iteration (iter1.png, iter2.png, ...) \n \
      \t                   \t             2-> Save only the last one at the maximum iteration (lastiter.png) \n \
      \t                   \t             3-> Save each iteration to the same file (lastiter.png) \n \
      \t -nomeasure \t using evenly spaced multiple shooting intervals\n\
      \t -doP <list> \t specifies parameters to be fitted, e.g. -doP 110L\n\
      \t -maxit <num> \t maximal number of iterations\n\
      \t -wait \t\t wait after each iteration\n\
      \t -usesig \t use weights from input file\n\
      \t -int <num> \t integrator: 1-> Runge-Kutta (default)\n \
      \t            \t             2-> CVODES \n \
      \t            \t             3-> ODESSA \n \
      \t -maxstp \t maximal integration step\n\
      \t -minimp \t minimal improvement for iteration truncation \n\
      \t -nowait \t do not wait until keypressed after last iteration\n\
      \t -elastic <val>  a factor (0< <val> <= 1) weakens cont. constraints, default enforced constraints: 1.0 \n\
      \t -reg \t\t regularisation for ill posed problems, default: no reg. \n\
      \t -epsilon <val>  singular value threshold, default: 1e-10\n\
      \t -lambda <val>   regularisation parameter, default: 1e6\n\
      \t -spline <list>  specifies spline data for non-autonomous ODEs \n\
      \t -siminit \t simulate initial state\n\
      \t -pert <val>\t perturbate initial try - only in combination with siminit\n\
      \t -y0fix <list>   fixing initial state, default: 11..1 (lenth: # of ODEs), not fixed\n\
      \t -nodamp \t switches damping off, default: on\n\
      \t -strat <var>    optimisation strategy:\n\
      \t                 1: MS 2: SRES (global) \n\
      \t -opt <num> \t selects the minimiser:\n\
      \t                 1: LSEI 2: NAG (default)\n\
      \t -Lconst <list>  set local Parameter Constraints, \n\
      \t                 e.g. -Lconst 2,4,3\n";

// separates variable 'string' in double arguments stored in variable 'arg'
// each floating number arg is separated by a comma ','
// n is the maximal number of arguments expected
// k counts the number of arguments already identified
// ind is index in original string array containing the 
// l is index in buffer variable dum storing the next argument string before conversion to double
// text used in case error messages related to be parsed.
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

/*! \fn GlobExp *parseopts(int argc, char *argv[],Glob *globs,char *outstr)
Main module routine doing all the parsing. 
\param argc The number of command line arguments as passed by \b main.
\param *argv[] String of the arguments.
\param *globs The structure Glob.
\param *outstr Deprecated, not used.
*/
GlobExp *parseopts(int argc, char *argv[],Glob *globs,char *outstr)
{
  int longindex,opt;
  int parlistSpecified=FALSE;
  long k,l,nExp,sp;
  char name[MAXLENGTH],line[MAXLENGTH]; //name of parameter file

  // options definition
  // {name, has_arg, flag, value}
  //has_arg: 0 no arg, 1 required arg, 2 optional arg
  //flag = 0: getopt_long_only returns value as its output integer 
  //value: value returned by getopt_long_only
  static struct option longopts[]={
    {"help", 0, 0,  'h'},
    {"help", 0, 0,  '?'},
    {"x0"  , 1, 0,   1 },
    {"x1"  , 1, 0,   2 },
    {"nms"  ,1, 0,   3 },
    {"f"    ,1, 0,   4 },
    {"p"    ,1, 0,   5 },
    {"y0"   ,1, 0,   6 },
    {"nobs" ,1, 0,   7 },
    {"eps"  ,1, 0,   8 },
    {"nognu",0, 0,   9 },
    {"nomeasure",0,0,10},
    {"doP"  ,1, 0,   11},
    {"maxit",1, 0,   12},
    {"wait" ,0, 0,   13},
    {"usesig",0,0,   14},
    {"int"   ,1,0,   15},
    {"maxstp",1,0,   17},
    {"minimp",1,0,   18},
    {"nowait",0,0,   19},
    {"elastic" ,1,0, 21},
    {"reg"     ,0,0, 22},
    {"epsilon" ,1,0, 23},
    {"lambda"  ,1,0, 24},
    {"spline"  ,1,0, 25},
    {"siminit" ,0,0, 26},
    {"pert"    ,1,0, 27},
    {"y0fix"   ,1,0, 28},
    {"nodamp"  ,0,0, 29},
    {"strat"   ,1,0, 30},
    {"opt"     ,1,0, 31},
    {"Lconst"  ,1,0, 32},
    {"tms"  ,1, 0,   33 },
    {"savegnupng",1, 0,	34},
    {0, 0, 0, 0}
  };
  //initialise some global parameters in globs
  globs->noGnu=FALSE; // by default use gnuplot
  globs->eps=1e-6;  // default integration accuracy
  globs->npar=NPARAMS; //defined in the model.cc file
  globs->noMeasurements=FALSE; //

  //if the parameters will be fitted (by default, yes)
  globs->doP=ivector(1,globs->npar); //allocate an int vector with subscript range from 1 to npar
  for(k=1;k<=globs->npar;k++)
    globs->doP[k]=TRUE;

  //if the initial value of the variables are fixed or not.
  //allocate int vector with subscript from 1 to the number of ODEs in system. NEQNs defined in model.cc
  globs->y0fix=ivector(1,NEQNS); 
  for(k=1;k<=NEQNS;k++)
    globs->y0fix[k]=TRUE; //default value 1, not fixed
  
  globs->maxit=1000; //maximum # of iterations
  globs->gnuFp=NULL;
  globs->wait=FALSE;
  globs->savegnupng=0; 
  globs->usesig=FALSE;
  globs->integrator=1; // Runge-Kutta by default
  globs->stiff=TRUE;
  globs->maxstp=5000; //maximum number of integration steps
  globs->minimp=1e-4; //minimum improvement per iteration
  globs->nowait=FALSE;
  globs->elastic=1.;
  globs->reg=FALSE;
  globs->epsilon=1e-10;
  globs->lambda=1e6;
  globs->simInit=FALSE;
  globs->pert=0.; //perturbate initial state (after simulating the ODE for initialization)
  globs->nodamp=FALSE;
  globs->initSpline=TRUE; //Note: this results in spline initialization in initialise.cc
  globs->wquer=-1;
  globs->silent=FALSE;
  globs->strategy=1;  //multiple shooting
  globs->minimiser=2; //NAG
  globs->faktorL=dvector(1,NPARAMS); //allocate double vector with subscript range 1 to NPARAMS
  globs->faktorLexist=FALSE; // no local parameter constraints

  //initially no parameters are considered local
  for(k=1;k<=NPARAMS;k++)
  {
  	globs->faktorL[k]=-1;
  }

  //long options, start by '-' or '--'
  //"h?" are legitimate single option characters
  // opt: number associated with the options, see value in longopts above
  while ((opt=getopt_long_only(argc,argv,"h?",longopts,&longindex)) != -1)
    {
      //cycle through all the options to see if the arguments are in a file
      switch(opt)
	{
	case 4:  //we have file with the parameters; save its name
	  parlistSpecified=TRUE;
	  // optarg: pointer to text following the option, see getopt documentation in man 3
	  strncpy(name,optarg,MAXLENGTH); //save name of parameter file
	  break;
	}
    }
  //optind, index of the next element to be processed in argv
  //reset option index; next calls to getopt will start from beginning of argv
  optind=0;

  //initializes parameters for specific experiments in ex
  GlobExp *ex;
  if(parlistSpecified==FALSE)
    {
      //no file specifying the parameters, so set the number of experiments to 1
      globs->nrExp=1;
      ex=new GlobExp[globs->nrExp+1];   
      //initialise some global parameters
      ex[1].fitstart=-1e99; //time of fit start and end
      ex[1].fitend=1e99;
      ex[1].nobs=NOBS; //NOBS, NEQNS defined in model.cc
      ex[1].nvar=NEQNS;
      ex[1].nPoints=2; //number of multiple shooting interval boundary points
      ex[1].splinesDefined=FALSE; //no splines by default
      ex[1].y0=NULL; //initial values for the ODE variables
      ex[1].tms=NULL; //time boundaries of multiple shooting intervals
      // initialisation of the parameters
      ex[1].par=dvector(1,NPARAMS);
      for(k=1;k<=globs->npar;k++)
	ex[1].par[k]=DefParameters[k-1]; //DefParameters defined in model.cc
      
      if(NSPLINES!=0)
	{
	  //if there are splines, initialize a long containing the number of nodes for each spline
	  ex[1].nNodes=lvector(1,NSPLINES);
	}
    } // no parameter list file
  else //if parlist specified; open parameter file and parse it
    {
      ifstream in;      ofstream out;
      char str[MAXLENGTH]; //buffer for lines of parameter file, out file name

      in.open(name); //open parameter file

      if(!in)
	{
	  cerr << "parameter list " << name << "  not found\n";
	  exit(1);
	}
      
      //preprocessing of data list
      //essentially read the file and write each expt parameters in a separate file
      //with the same name and a trailing ".#", where # is the experiment index.
      //The first file, with index 0 contains the global parameters.
      globs->nrExp=0;
      sprintf(str,"%s.%d",name,globs->nrExp);
      out.open(str);
      while(!in.eof())
	{
	  in.getline(str,MAXLENGTH,'\n');
	  if(str[0]=='@')  //separator character for experiments; close file and open next one
	    {
	      out.flush();
	      out.close();
	      globs->nrExp++;
	      sprintf(str,"%s.%d",name,globs->nrExp);
	      out.open(str);
	    }
	  else if(str[0]!='#') //not a comment, just copy to parameter file
	    {
	      out << str << endl;
	    }
	}
      out.flush();
      out.close();
      in.close();

      if(globs->nrExp==0)
	{
	  cerr << "insufficient number of experiments, @ missing\n";
	  exit(1);
	}

      //read in the global parameters file as dummy argc, argv list and then parse it
      long _argc;
      char *_argv[MAXLENGTH];

      for(k=0;k<=MAXLENGTH;k++)
	_argv[k]=(char *)malloc((size_t) (150*sizeof(char))); //each argv is at most 150 char long

      {
	ifstream inp;
	
	_argc=0;
	sprintf(str,"%s.0",name); //open file
	inp.open(str);
	
	while(!inp.eof())
	  {
	    inp >> _argv[_argc+1]; //store line in next _argv element
	    _argc++;
	  }
	inp.close();
      }
   
      //PARSING LIST of global parameters
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
	    case 8: //integration accuracy
	      globs->eps=fabs(atof(optarg));
	      break;
	    case 9: // no gnuplotting
	      globs->noGnu=TRUE;
	      break;
	    case 10: //use evenly spaced in time multiple shooting intervals
	      globs->noMeasurements=TRUE;
	      break;
	    case 11: //For each parameter a 0 or a 1 depending on whether to fit or not a given parameter.
	             //Can also be an L for a locak optimalization parameter
	             //optarg defined by getopt, contains the argument string for "doP" option
	      if(strlen(optarg)!=globs->npar)
		{
		  cerr << strlen(optarg) << " parameter(s) specified instead of " <<  globs->npar << ".\n";
		  exit(1);
		}
	      for(k=0;k < globs->npar;k++)
		{
		  if(optarg[k]=='0')
		    globs->doP[k+1]=FALSE;  //no fit
		  else if (optarg[k]=='L')
		    globs->doP[k+1]='L';    //treated as local
		  else
		    globs->doP[k+1]=TRUE;   //do fit
		}
	      break;
	    case 12:  //maximal number of iterations, stored in a long
	      globs->maxit=abs(atol(optarg)); //convert ascii to long and take absolute value
	      break;
	    case 13:
	      globs->wait=TRUE;
	      break;
	    case 14:
	      globs->usesig=TRUE;
	      break;
	    case 15:
	      globs->integrator=abs(atoi(optarg));
	      break;
	    case 17: //max integration step
	      globs->maxstp=abs(atoi(optarg));
	      break;
	    case 18: //minimum improvement
	      globs->minimp=fabs(atof(optarg));
	      break;
	    case 19:
	      globs->nowait=TRUE;
	      break;
	    case 21: //weakening of continuity constraints
	      globs->elastic=fabs(atof(optarg));
	      if(globs->elastic > 1. || globs->elastic == 0.)
		{
		  cerr << "Insufficient range of -elast <val> \n";
		  exit(1);
		}
	      break;
	    case 22: //regularization
	      globs->reg=TRUE;
	      break;
	    case 23: //integration accuracy
	      globs->epsilon=fabs(atof(optarg));
	      break;
	    case 24: //regularization parameter
	      globs->lambda=fabs(atof(optarg));
	      break;
	    case 26: //do initial simulation
	      globs->simInit=TRUE;
	      break;
	    case 27: //random perturbation parameter for initial simulation values
	      globs->pert=fabs(atof(optarg));
	      break;
	    case 28: //list of fixed initial states
	      if(strlen(optarg)!=NEQNS)  //has to be equal to the number of ODEs
		{
		  cerr << strlen(optarg) << " variable(s) specified instead of " <<  NEQNS << ".\n";
		  exit(1);
		}
	      for(k=0;k < NEQNS;k++) //optarg indexed from 0 to NEQNS-1
		{
		  if(optarg[k]=='0') //initial values indexed from 1 to NEQNS
		    globs->y0fix[k+1]=FALSE; //this initial value is fixed
		  else
		    globs->y0fix[k+1]=TRUE; //this initial value is not fixed (default)
		}
	      break;
	    case 29:
	      globs->nodamp=TRUE;
	      break; 
	    case 30:
	      globs->strategy=abs(atoi(optarg));
	      break;
	    case 31:
	      globs->minimiser=abs(atoi(optarg));
	      break;
	    case 32: //store in faktorL the values for the local parameter constraints
	      get_list(globs->npar,globs->faktorL,optarg," local parameter constraints ");
	      globs->faktorLexist=TRUE;
	      break;
	    case 34:
	      globs->savegnupng=abs(atoi(optarg));
	      break;     
	    default:
	      cerr << endl;
	      cerr << "Parsing parameter list produced errors.\n\n";
	      exit(1);
	    }
	} //end parsing first parameter file

      optind=0;
      ex=new GlobExp[globs->nrExp+1]; //array of GlobExp; size is 1 more than nrExp since we index from 1  
  
      //PARSING experiment specific list
      for(nExp=1;nExp<=globs->nrExp;nExp++)
	{
	  //initialise some global parameters
	  ex[nExp].fitstart=-1e99;
	  ex[nExp].fitend=1e99;
	  ex[nExp].nobs=NOBS;
	  ex[nExp].nvar=NEQNS;
	  ex[nExp].nPoints=2; //number of shooting interval boundary points
	  ex[nExp].splinesDefined=FALSE;
	  ex[nExp].y0=NULL;
	  ex[nExp].tms=NULL;
	  // initialisation of the parameters
	  ex[nExp].par=dvector(1,NPARAMS);
	  for(k=1;k<=globs->npar;k++)
	    ex[nExp].par[k]=DefParameters[k-1];
	  
	  if(NSPLINES!=0)
	    {
	      ex[nExp].nNodes=lvector(1,NSPLINES);
	    }

	  _argc=0;

	  //read parameter file for each specific experiment
	  ifstream inp;
	  sprintf(str,"%s.%d",name,nExp);
	  inp.open(str);

	  while(!inp.eof())
	    {
	      inp >> _argv[_argc+1];
	      _argc++;
	    }

	  inp.close();
	  
	  //parse options from parameter file
	  while((opt=getopt_long_only(_argc,_argv,"h?",longopts,&longindex)) != -1)
	     {
	       switch(opt)
		 {
		 case 1:
		   ex[nExp].fitstart=atof(optarg);
		   break;
		 case 2:
		   ex[nExp].fitend=atof(optarg);
		   break;
		 case 3: //number of multiple shooting intervals boundary points
		   ex[nExp].nPoints=abs(atol(optarg)); //user specified number of multiple shooting intervals 
		   ex[nExp].nPoints+=1; //add one extra boundary point
		   break;
		 case 5:  //get list of parameters fit status (0 you fall out, 1 fit, L local)
		   get_list(globs->npar,ex[nExp].par,optarg," parameter ");
		   break;
		 case 6:
		   ex[nExp].y0=dvector(1,ex[nExp].nvar);
		   for(k=1;k<=ex[nExp].nvar;k++)
		     ex[nExp].y0[k]=0.;
		   k=get_list(ex[nExp].nvar,ex[nExp].y0,optarg," initial values ");
		   if(k!=ex[nExp].nvar)
		     {
		       cerr << ex[nExp].nvar << " initial values required.\n";
		       exit(1);
		     }
		   break;
		 case 33:
		   ex[nExp].tms=dvector(1,ex[nExp].nPoints);
		   for(k=1;k<=ex[nExp].nPoints;k++)
		     ex[nExp].tms[k]=0.;
		   k=get_list(ex[nExp].nPoints,ex[nExp].tms,optarg," multiple shooting time boundaries ");
		   if(k!=ex[nExp].nPoints)
		     {
		       cerr << ex[nExp].nPoints << " multiple shooting time boundaries required.\n";
		       exit(1);
		     }
		   break;
		 case 7: //number of observed quantities
		   ex[nExp].nobs=atol(optarg);
		   break;
		 case 25: //look at optarg and determine if there are enough spline file names defined
		   k=0;
		   l=1;
		   ex[nExp].splinesDefined=TRUE;
		   while(optarg[k]!='\0')
		     {
		       if(optarg[k]==',')
			 l++;
		       k++;
		     }
		   if(l!=NSPLINES) //note: NSPLINES is defined in model.cc after parsing of the .mdf file
		     //looking from splines defined as spline1, spline2, etc...
		     {
		       cerr << "Incompatible number of splines.\n";
		       exit(1);
		     }
		   //saves the name of the spline files
		   ex[nExp].splineFile=new string[l+1];
		   l=1;
		   sp=0;
		   k=0;
		   while(optarg[k]!='\0')
		     {
		       if(optarg[k]==',') //finished reading one spline file name; save it
			 {
			   line[sp]='\0';
			   ex[nExp].splineFile[l]=(string)line;
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
		   ex[nExp].splineFile[l]=(string)line;
		   break;
		 default:
		   cerr << endl;
		   cerr << "Parsing parameter list produced errors.\n\n";
		   exit(1);
		 } //end switch(opt)
	     } //end parsing (getopt_long_only)

	   //Parsing file name of data to be fit
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
	   optind=0;
	} //end PARSING experiment specific list 
      
     
      //for(k=0;k<=MAXLENGTH;k++)
      //free(_argv[k]);

      sprintf(str,"rm -f %s.*",name);
      system(str);
    } //end else parlist specified

  //reset index of next element to scan to zero
  optind=0;

  //finally, parse any additional parameters from the command line
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
	case 3: //number of shooting intervals boundary points
	  ex[1].nPoints=abs(atol(optarg)); //user specified number of shooting intervals
	  ex[1].nPoints+=1; //add one extra boundary point
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
	case 33:
          //nExp=1 here 
	  ex[1].tms=dvector(1,ex[1].nPoints);
	  for(k=1;k<=ex[1].nPoints;k++)
	    ex[1].tms[k]=0.;
	  k=get_list(ex[1].nPoints,ex[1].tms,optarg," multiple shooting time boundaries ");
	  if(k!=ex[1].nPoints)
	    {
	      cerr << ex[1].nPoints << " multiple shooting time boundaries required.\n";
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
	  globs->noGnu=TRUE;
	  break;
	case 10:
	  globs->noMeasurements=TRUE;
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
	      else if (optarg[k]=='L')
		globs->doP[k+1]='L'; 
	      else
		globs->doP[k+1]=TRUE;
	    }
	  break;
	case 12:
	  globs->maxit=abs(atol(optarg));
	  break;
	case 13:
	  globs->wait=TRUE;
	  break;
	case 14:
	  globs->usesig=TRUE;
	  break;
	case 15:
	  globs->integrator=abs(atoi(optarg));
	  break;
	case 17:
	  globs->maxstp=abs(atoi(optarg));
	  break;
	case 18:
	  globs->minimp=fabs(atof(optarg));
	  break;
	case 19:
	  globs->nowait=TRUE;
	  break;
	case 21:
	  globs->elastic=fabs(atof(optarg));
	  if(globs->elastic > 1. || globs->elastic == 0.)
	    {
	      cerr << "Insufficient range of -elast <val> \n";
	      exit(1);
	    }
	  break;
	case 22:
	  globs->reg=TRUE;
	  break;
	case 23:
	  globs->epsilon=fabs(atof(optarg));
	  break;
	case 24:
	  globs->lambda=fabs(atof(optarg));
	  break;
	case 25:
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
	case 26:
	  globs->simInit=TRUE;
	  break;
	case 27:
	  globs->pert=fabs(atof(optarg));
	  break;
	case 28:
	  if(strlen(optarg)!=NEQNS)
	    {
	      cerr << strlen(optarg) << " variable(s) specified instead of " <<  NEQNS << ".\n";
	      exit(1);
	    }
	  for(k=0;k < NEQNS;k++)
	    {
	      if(optarg[k]=='0')
		globs->y0fix[k+1]=FALSE;
	      else
		globs->y0fix[k+1]=TRUE;
	    }
	  break;
	case 29:
	  globs->nodamp=TRUE;
	  break;
	case 30:
	  globs->strategy=abs(atoi(optarg));
	  break;
	case 31:
	  globs->minimiser=abs(atoi(optarg));
	  break;
	case 32:
	  get_list(globs->npar,globs->faktorL,optarg," local paramter constraints ");
	  globs->faktorLexist=TRUE;
	  break;
	case 34:
	  globs->savegnupng=abs(atoi(optarg));
	  break;
	default:
	  cerr << endl;
	  cerr << "Parsing command line options produced errors \n\n";
	  exit(1);
	} //switch opt
    } //all command line parameters have been processed
  
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

  //check if spline data has been defined for each experiment
  if(NSPLINES!=0)
    {
      for(k=1;k<=globs->nrExp;k++)
	{
	  if(ex[k].splinesDefined==FALSE)
	    {
	      cerr << "Please specify spline data in experiment " << k << ".\n";
	      exit(1);
	    }
	}
    }
  
  return(ex);
}

