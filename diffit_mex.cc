#include "mex.h"
#include<iostream>
#include<fstream>
#include<stdlib.h>
#include<string.h>
#include<math.h>


#include "def.h"
#include "model.h"
#include "nr.h"

using namespace std;


//Begin definition of module prototypes

void setMesh(GlobExp *ex,Glob *globs,long expNr);
void outFit(GlobExp ex[],Glob *globs);
void freeMem(GlobExp *ex,Glob *globs,int simit);
void initialise(GlobExp ex[],Glob *globs,int simit);
void simInit(GlobExp *ex,Glob *globs);


//End definition of module prototypes

//The numerics subprogram
void fitIt(GlobExp ex[],Glob *globs);

// DEBUG stream
ofstream *dbg;

//input
#define IN prhs[0]
#define INEX prhs[1]
//output
#define OUT plhs[0]
#define OUTEX plhs[1]


void mexFunction( int nlhs,
                  mxArray *plhs[],
                  int nrhs,
                  const mxArray *prhs[])
{
  long k,i,j,nExp;
  long index,nglob=0;
  unsigned long *ndat;
  Glob globs;
  GlobExp *ex;
  mxArray *tmp,*sig;
  double *tmp_,*sig_,*spline_;
  char tmpstr[1000];
  mxArray *splinedat;
  int except=0;
  int init=0;

  //for output
  const char *fieldNamesOUT[]={"wquer","converged","chisq","except","Lambda","cond"};
  const char *fieldNamesOUTEX[]={"p","y0","mesh","errP","errY0"};

  //  int dims[2];
  mwSize dims[2];
  
  //GET DATA FOR GLOBS

  globs.npar=NPARAMS;
  globs.wait=FALSE;
  globs.nowait=TRUE;
  globs.elastic=1.;
  globs.initSpline=FALSE;
  globs.maxstp=5000;
  globs.gnuFp=NULL;
  globs.initSpline=TRUE; 
  globs.wquer=100.;
  globs.silent=FALSE;
  globs.minimiser=2;

  //eps
  tmp=mxGetField(IN,0,"eps");
  tmp_=mxGetPr(tmp);
  globs.eps=fabs(tmp_[0]);
  //eps
  tmp=mxGetField(IN,0,"nognu");
  tmp_=mxGetPr(tmp);
  globs.noGnu=(int)tmp_[0];  
  //nomeasure
  tmp=mxGetField(IN,0,"nomeasure");
  tmp_=mxGetPr(tmp);
  globs.noMeasurements=(int)tmp_[0];   
  //maxit
  tmp=mxGetField(IN,0,"maxit");
  tmp_=mxGetPr(tmp);
  globs.maxit=abs((int)tmp_[0]);  
  //odessa
  tmp=mxGetField(IN,0,"int");
  tmp_=mxGetPr(tmp);
  globs.integrator=abs((int)tmp_[0]);  
  globs.stiff=TRUE;
  //minimp
  tmp=mxGetField(IN,0,"minimp");
  tmp_=mxGetPr(tmp);
  globs.minimp=fabs(tmp_[0]);
  //reg
  tmp=mxGetField(IN,0,"reg");
  tmp_=mxGetPr(tmp);
  globs.reg=(int)tmp_[0];
  //epsilon
  tmp=mxGetField(IN,0,"epsilon");
  tmp_=mxGetPr(tmp);
  globs.epsilon=fabs(tmp_[0]); 
  //lambda
  tmp=mxGetField(IN,0,"lambda");
  tmp_=mxGetPr(tmp);
  globs.lambda=fabs(tmp_[0]); 
  //siminit
  tmp=mxGetField(IN,0,"siminit");
  tmp_=mxGetPr(tmp);
  globs.simInit=(int)tmp_[0];
  //pert
  tmp=mxGetField(IN,0,"pert");
  tmp_=mxGetPr(tmp);
  globs.pert=fabs(tmp_[0]);
  //nodamp
  tmp=mxGetField(IN,0,"nodamp");
  tmp_=mxGetPr(tmp);
  globs.nodamp=(int)tmp_[0];
  //nrExp
  tmp=mxGetField(IN,0,"nexp");
  tmp_=mxGetPr(tmp);
  globs.nrExp=(long)tmp_[0];
  //doP
  globs.doP=ivector(1,globs.npar);
  for(k=1;k<=globs.npar;k++)
    globs.doP[k]=TRUE;
  tmp=mxGetField(IN,0,"doP");
  if(mxGetString(tmp,tmpstr,1000))
    {
      cerr << "Parsing doP failed.\n";
      return;
    }
  if(tmpstr[0]!='N')
    {
      if(strlen(tmpstr)!=globs.npar)
	{
	  cerr << "doP must be of length " << globs.npar << ".\n";
	  return;
	}
      else
	{
	  for(k=0;k < globs.npar;k++)
	    {
	      if(tmpstr[k]=='0')
		globs.doP[k+1]=FALSE;
	      else if (tmpstr[k]=='L')
		globs.doP[k+1]='L'; 
	      else
		globs.doP[k+1]=TRUE;
	    }
	}
    }
  //y0fix
  globs.y0fix=ivector(1,NEQNS);
  for(k=1;k<=NEQNS;k++)
    globs.y0fix[k]=TRUE;
  tmp=mxGetField(IN,0,"y0fix");
  if(mxGetString(tmp,tmpstr,1000))
    {
      cerr << "Parsing y0fix failed.\n";
      return;
    }
  if(tmpstr[0]!='N')
    {
      if(strlen(tmpstr)!=NEQNS)
	{
	  cerr << "y0fix must be of length " << NEQNS << ".\n";
	  return;
	}
      else
	{
	  for(k=0;k < NEQNS ;k++)
	    {
	      if(tmpstr[k]=='0')
		globs.y0fix[k+1]=FALSE;
	      else
		globs.y0fix[k+1]=TRUE;
	    }
	}
    }
  //wquer
  tmp=mxGetField(IN,0,"wquer");
  tmp_=mxGetPr(tmp);
  globs.wquer=tmp_[0];
  //init
  tmp=mxGetField(IN,0,"init");
  tmp_=mxGetPr(tmp);
  init=(int)tmp_[0];
   //silent
  tmp=mxGetField(IN,0,"silent");
  tmp_=mxGetPr(tmp);
  globs.silent=(int)tmp_[0]; 
   //opt
  tmp=mxGetField(IN,0,"opt");
  tmp_=mxGetPr(tmp);
  globs.minimiser=(int)tmp_[0]; 
  //Lconst
  tmp=mxGetField(IN,0,"Lconst");      
  tmp_=mxGetPr(tmp);
  if(mxIsEmpty(tmp))
  {
     globs.faktorLexist=FALSE;
  }
  else
  {
     globs.faktorL=dvector(1,NPARAMS);
     for(k=1;k<=NPARAMS;k++)
  	globs.faktorL[k]=-1;
     for(k=1;k<=mxGetNumberOfElements(tmp) && k<=NPARAMS;k++)
     {
        globs.faktorL[k]=tmp_[k-1];
     }
     globs.faktorLexist=TRUE;
  }     



  // GET DATA FOR EX
  ndat=lvector(1,globs.nrExp);
  ex=new GlobExp[globs.nrExp+1];

  for(nExp=1;nExp<=globs.nrExp;nExp++)
    {
      //general stuff
      ex[nExp].nvar=NEQNS;
      ex[nExp].splinesDefined=FALSE;
      //t0
      tmp=mxGetField(INEX,nExp-1,"t0");
      tmp_=mxGetPr(tmp);
      ex[nExp].fitstart=tmp_[0];
    
      //t1
      tmp=mxGetField(INEX,nExp-1,"t1");
      tmp_=mxGetPr(tmp);
      ex[nExp].fitend=tmp_[0];
      if(ex[nExp].fitstart >= ex[nExp].fitend)
	{
	  cerr << "Fitstart >= fitend in exp. " << nExp << ".\n";
	  return;
	}
      //nobs
      tmp=mxGetField(INEX,nExp-1,"nobs");      
      tmp_=mxGetPr(tmp);
      ex[nExp].nobs=(long)tmp_[0];
      if(ex[nExp].nobs > NOBS)
	{
	  cerr << "More observables than required in exp. " << nExp << ".\n";
	  return;
	}
      //nms
      tmp=mxGetField(INEX,nExp-1,"nms");      
      tmp_=mxGetPr(tmp);
      ex[nExp].nPoints=(long)tmp_[0]+1;
      //p
      ex[nExp].par=dvector(1,NPARAMS);
      tmp=mxGetField(INEX,nExp-1,"p");  
      tmp_=mxGetPr(tmp);
      if(mxIsEmpty(tmp))
	{
	  for(k=1;k<=globs.npar;k++)
	    ex[nExp].par[k]=DefParameters[k-1];
	}
      else
	{
	  if(mxGetNumberOfElements(tmp)!=globs.npar)
	    {
	      cerr << globs.npar << " parameter(s) required in exp. " << nExp <<".\n";
	      return;
	    }
	  else
	    {
	      for(k=1;k<=globs.npar;k++)
		ex[nExp].par[k]= tmp_[k-1];
	    }
	}
      //y0
      tmp=mxGetField(INEX,nExp-1,"y0");      
      tmp_=mxGetPr(tmp);
      if(mxIsEmpty(tmp))
	{
	  ex[nExp].y0=NULL;
	}
      else
	{
	  if(mxGetNumberOfElements(tmp)!=ex[nExp].nvar)
	    {
	      cerr << ex[nExp].nvar << " initial value(s) required in exp. " << nExp <<".\n";
	      return;
	    }
	  else
	    {
	      ex[nExp].y0=dvector(1,ex[nExp].nvar);
	      for(k=1;k<=ex[nExp].nvar;k++)
		ex[nExp].y0[k]= tmp_[k-1];
	    }
	}     
      //n
      tmp=mxGetField(INEX,nExp-1,"n");      
      tmp_=mxGetPr(tmp);
      ndat[nExp]=(long)tmp_[0];
      
    } //end for(nExp=1;nExp<=globs.nrExp;nExp++)

  //open DEBUG stream
  dbg=new ofstream("diffit.dbg");
  if (!dbg) 
    {
      cerr << "Error opening DEBUG file.\n";
      return;
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
  for(i=1;i<=globs.nrExp;i++)
    {
      *dbg << "Experiment: " << i << "\n";
      *dbg << NPARAMS << " parameter(s):";
      for (k=1; k<=NPARAMS; ++k) 
	*dbg << " " << ex[i].par[k];
      *dbg << "\n";
    }
   
  dbg->flush();
  
 
  for(i=1;i<=globs.nrExp;i++)
    {
      //read data and check for for consistency

      //data 
      tmp=mxGetField(INEX,i-1,"data");      
      tmp_=mxGetPr(tmp);
      sig=mxGetField(INEX,i-1,"sig");      
      sig_=mxGetPr(sig);
      //deternine nMeasure
      index=1;
      for(j=1;j<=ndat[i];j++)
	{
	  if(tmp_[j-1]>= ex[i].fitstart && tmp_[j-1]<= ex[i].fitend)
	    index++;
	}
      index--; 
      if(index==0)
	{
	  cerr << "Nothing to fit in exp " << i << ".\n";
	  return;
	}
      ex[i].nMeasure=index;
      //memory allocation
      ex[i].xMeasure=dvector(1,ex[i].nMeasure);
      ex[i].yMeasure=dmatrix(1,ex[i].nMeasure,1,ex[i].nobs);
      ex[i].sigma=dmatrix(1,ex[i].nMeasure,1,ex[i].nobs);
      //copy data
      index=1;
      for(j=1;j<=ndat[i];j++)
	{
	  if(tmp_[j-1]>= ex[i].fitstart && tmp_[j-1]<= ex[i].fitend)
	    {
	      ex[i].xMeasure[index]=tmp_[j-1];
	      for(k=1;k<=ex[i].nobs;k++)
		{
		  ex[i].yMeasure[index][k]=tmp_[ndat[i]*k+j-1];
		  ex[i].sigma[index][k]=sig_[ndat[i]*(k-1)+j-1];
		}
	      index++;
	    }
	}
      //check if time is in ascending order
      for(j=2;j<=ex[i].nMeasure;j++)
	{
	  if(ex[i].xMeasure[j-1] >= ex[i].xMeasure[j])
	    {
	      cerr << "Data not in temporal ascending order in exp " << i << ",\n";
	      return;
	    }
	}
      //setting some variables
      ex[i].firstMeasure=1;
      ex[i].lastMeasure=ex[i].nMeasure;
      
      //set mesh
      tmp=mxGetField(INEX,i-1,"mesh");      
      tmp_=mxGetPr(tmp);
      if(mxIsEmpty(tmp))
	setMesh(ex,&globs,i);
      else
	{
	  ex[i].nPoints=mxGetM(tmp)+2;
	  if (ex[i].nMeasure < ex[i].nPoints) 
	    {
	      cerr << "Too few measurements to construct mesh\n";
	      return;
	    }

	  if(mxGetN(tmp) <= ex[i].nvar+1)
	    {
	      ex[i].mesh=dvector(1,ex[i].nPoints);
	      ex[i].mesh[1]= ex[i].fitstart;
	      ex[i].mesh[ex[i].nPoints]= ex[i].fitend;
	      //reading mesh data
	      for(j=1;j<=ex[i].nPoints-2;j++)
		{
		  ex[i].mesh[j+1]=tmp_[j-1];
		}
	      //check if mesh is in ascending order
	      for(j=2;j<=ex[i].nPoints;j++)
		{
		  if(ex[i].mesh[j-1] >= ex[i].mesh[j])
		    {
		      cerr << "Mesh not in ascending order in exp " << i << ",\n";
		      return;
		    }
		}
	    }
	  else
	    {
	      cerr << "Incorrect format of mesh in exp " << i << ",\n"; 
	      return;
	    }
	}
      
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
    }
  
  //initial values
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
    } // end for(i=1;i<=globs.nrExp;i++)
  // starting of the numerics

  
   // spline section and initialse
  
  if(NSPLINES==0)
    {
      globs.initSpline=TRUE;
      initialise(ex,&globs,FALSE);
    }
  else
    {
      for(nExp=1;nExp<=globs.nrExp;nExp++)
	{
	  globs.initSpline=FALSE;
	  initialise(ex,&globs,FALSE);  
	  ex[nExp].splineNodes=(double **) malloc((size_t)(NSPLINES+1)*sizeof(double*));
	  ex[nExp].splineY=(double **) malloc((size_t)(NSPLINES+1)*sizeof(double*));
	  ex[nExp].splineGam=(double **) malloc((size_t)(NSPLINES+1)*sizeof(double*));
	  ex[nExp].nNodes=(long unsigned*) malloc((size_t)(NSPLINES+1)*sizeof(long unsigned*));
	  tmp=mxGetField(INEX,nExp-1,"spline"); 
	  for(i=1;i<=NSPLINES;i++)
	    {
	      if(mxGetNumberOfElements(tmp) < NSPLINES)
		{
		  cerr << "Please specify " << NSPLINES << " spline(s).\n";
		  return;
		}
	      
	      splinedat=mxGetCell(tmp,i-1);
	      spline_=mxGetPr(splinedat);
	      ex[nExp].nNodes[i]=mxGetM(splinedat);
	      if(mxGetN(splinedat)!=3)
		{
		  cerr << "Invalid data for spline " << i << ".\n";
		  return;
		}
	      ex[nExp].splineNodes[i]=(double *) malloc((size_t)(ex[nExp].nNodes[i]+1)*sizeof(double));
	      ex[nExp].splineY[i]=(double *) malloc((size_t)(ex[nExp].nNodes[i]+1)*sizeof(double));
	      ex[nExp].splineGam[i]=(double *) malloc((size_t)(ex[nExp].nNodes[i]+1)*sizeof(double));	      
	      
	      //copy spline data
	      for(j=1;j<=ex[nExp].nNodes[i];j++)
		{
		  ex[nExp].splineNodes[i][j]=spline_[j-1];
		  ex[nExp].splineY[i][j]=spline_[ex[nExp].nNodes[i]+j-1];
		  ex[nExp].splineGam[i][j]=spline_[2*ex[nExp].nNodes[i]+j-1];	    
		}
	    }
	}
    }
 
  //simulate initial state
  try
    {
      if(globs.simInit==TRUE)
	{
	  for (i=1;i<=globs.nrExp;++i) 
	    simInit(&ex[i],&globs);
	}
    }
  catch (int i)
    {
      cerr << "Cannot integrate trajectory, no siminit!\n";
    }

  //setting initial values from mesh data
  for(nExp=1;nExp<=globs.nrExp;nExp++)
    {
      tmp=mxGetField(INEX,nExp-1,"mesh");      
      tmp_=mxGetPr(tmp);
      if(!mxIsEmpty(tmp))
	{
	  index=mxGetN(tmp);
	  for(i=1;i<=ex[nExp].nPoints-2;i++)
	    {
	      for(j=1;j<=index-1;j++)
		{
		  ex[nExp].yTry[i+1][j]=tmp_[j*(ex[nExp].nPoints-2)+i-1];
		}
	    }
	}
    }
  
  //*************************
   try 
      {

	if(init!=1)
	  fitIt(ex,&globs);
	else
	  globs.fitConverged=0;
      }
   catch(int i)  //Exception handling
     {
       dims[0]=1;
       dims[1]=1;
       OUT = mxCreateStructArray(2,dims,sizeof(fieldNamesOUT)/sizeof(*fieldNamesOUT),fieldNamesOUT);
         //error
       tmp = mxCreateDoubleMatrix(1,1,mxREAL);
       tmp_= mxGetPr(tmp);
       tmp_[0]= 1;
       mxSetField(OUT,0,"except",tmp);    

       dims[0]=1;
       dims[1]=globs.nrExp;
       OUTEX = mxCreateStructArray(2,dims,sizeof(fieldNamesOUTEX)/sizeof(*fieldNamesOUTEX),fieldNamesOUTEX);
       if(globs.noGnu==FALSE)
	 {
	   for(j=1;j<=globs.ngnu;j++)
	     pclose(globs.gnuFp[j]);
	 }
       
       freeMem(ex,&globs,FALSE);
       delete ex;
       return;
     }
  
  //*************************
      
  //Output after convergence
  if(init!=1)
      {    
	outFit(ex,&globs);
	if(!globs.noGnu)
	  system("rm gnuout.dat");
      }
  
  //Mex Output
  
  dims[0]=1;
  dims[1]=1;
  
  
  OUT = mxCreateStructArray(2,dims,sizeof(fieldNamesOUT)/sizeof(*fieldNamesOUT),fieldNamesOUT);
  
  //wquer
  tmp = mxCreateDoubleMatrix(1,1,mxREAL);
  tmp_= mxGetPr(tmp);
  tmp_[0]= globs.wquer;
  mxSetField(OUT,0,"wquer",tmp); 
  //converged
  tmp = mxCreateDoubleMatrix(1,1,mxREAL);
  tmp_= mxGetPr(tmp);
  tmp_[0]= globs.fitConverged;
  mxSetField(OUT,0,"converged",tmp);   
  //chisq
  tmp = mxCreateDoubleMatrix(1,1,mxREAL);
  tmp_= mxGetPr(tmp);
  tmp_[0]= globs.chisq;
  mxSetField(OUT,0,"chisq",tmp); 
  //error
  tmp = mxCreateDoubleMatrix(1,1,mxREAL);
  tmp_= mxGetPr(tmp);
  tmp_[0]= 0;
  mxSetField(OUT,0,"except",tmp); 
  //Lambda
  tmp = mxCreateDoubleMatrix(1,1,mxREAL);
  tmp_= mxGetPr(tmp);
  tmp_[0]= globs.Lambda;
  mxSetField(OUT,0,"Lambda",tmp);  
  //cond
  tmp = mxCreateDoubleMatrix(1,1,mxREAL);
  tmp_= mxGetPr(tmp);
  tmp_[0]= globs.cond;
  mxSetField(OUT,0,"cond",tmp);  

  //Experiment Specific output -> OUTEX

  dims[0]=1;
  dims[1]=globs.nrExp;

  OUTEX = mxCreateStructArray(2,dims,sizeof(fieldNamesOUTEX)/sizeof(*fieldNamesOUTEX),fieldNamesOUTEX);


  for(nExp=1;nExp<=globs.nrExp;nExp++)
    {
      //p
      tmp = mxCreateDoubleMatrix(1,globs.npar,mxREAL);
      tmp_= mxGetPr(tmp);
      for(i=1;i<=globs.npar;i++)
	tmp_[i-1]= ex[nExp].par[i];
      mxSetField(OUTEX,nExp-1,"p",tmp); 
      //errP
      tmp = mxCreateDoubleMatrix(1,globs.npar,mxREAL);
      tmp_= mxGetPr(tmp);
      for(i=1;i<=globs.npar;i++)
	tmp_[i-1]= ex[nExp].errP[i];
      mxSetField(OUTEX,nExp-1,"errP",tmp); 
      //y0
      tmp = mxCreateDoubleMatrix(1,ex[nExp].nvar,mxREAL);
      tmp_= mxGetPr(tmp);
      for(i=1;i<=ex[nExp].nvar;i++)
	tmp_[i-1]= ex[nExp].yTry[1][i];
      mxSetField(OUTEX,nExp-1,"y0",tmp); 
      //errY0
      tmp = mxCreateDoubleMatrix(1,ex[nExp].nvar,mxREAL);
      tmp_= mxGetPr(tmp);
      for(i=1;i<=ex[nExp].nvar;i++)
	tmp_[i-1]= ex[nExp].errY0[i];
      mxSetField(OUTEX,nExp-1,"errY0",tmp); 
      //mesh
      tmp = mxCreateDoubleMatrix(ex[nExp].nPoints-2,ex[nExp].nvar+1,mxREAL);
      tmp_= mxGetPr(tmp);
      for(i=1;i<=ex[nExp].nPoints-2;i++)
	{
	  for(j=1;j<=ex[nExp].nvar+1;j++)
	    {
	      if(j==1)
		tmp_[i-1]=ex[nExp].mesh[i+1];
	      else
		tmp_[(j-1)*(ex[nExp].nPoints-2)+i-1]=ex[nExp].yTry[i+1][j-1];
	    }
	}
      mxSetField(OUTEX,nExp-1,"mesh",tmp); 
    }

  freeMem(ex,&globs,FALSE);
  free_lvector(ndat,1,globs.nrExp);
  delete ex;

}
