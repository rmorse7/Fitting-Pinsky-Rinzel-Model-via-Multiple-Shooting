#include<iostream>
#include<fstream>
#include<math.h>
#include<stdlib.h>
#include<stdint.h>

#include "../def.h"
#include "../model.h"
#include "../nr.h"

#include<nag.h>
#include<nage04.h>

using namespace std;

//#define PRINTDIMENSION
#define PRINTBIGSYS
//#define PRINTSTEP

//Sub-module minimiser

/* Black box routine for linear least squares problem   */
/* with linear equality and inequality constraints      */
/* ACM TOMS 8,3 (1982), 323-333                         */
/* edited for IBM RS/6000 compatibility by TM, 96/02/22 */

// Solves Ex = f, Ax ~= b (least square), Gx >= h, where
// E is me by n, A is ma by n, G is mg by n.
//FG082918 extern "C" int lsei_(double *w, long &mdw, long &me, long &ma, long &mg,
//extern "C" int lsei(double *w, long &mdw, long &me, long &ma, long &mg, 
//		long &n, double *prgopt, double *x, double *rnorme,
//		double *rnorml,
//		long *mode, double *ws, long *ip);
extern "C" void lsei_wrapper(double *w, int32_t &mdw, int32_t &me, int32_t &ma, int32_t &mg, int32_t &n,
			     double *prgopt, double *x, double *rnorme, double *rnorml,
			     int32_t &mode, double *ws, int32_t &wsdim, int32_t *ip, int32_t &ips);

// Solves 0.5*||b - Ax||^2 subject to the constraint l<= {x,Cx} <=u
// Parameters:
// M number of rows of H; N number of variables; NCLIN number of linear constraints
// LDC, LDA first dimension of array C and A
// linear constraint matrix
// BL (BU) vectors of lower (upper) bounds for the variables x and for Cx
// CVEC vector of linear terms for objective function (unused in our case)
// ISTATE
// X initial estimate of the solution
// ITER total number of iterations performed
// CLAMBDA values of the Lagrange multipliers on exit
// IWORK integer array of dim LIWORK
// WORK float array of dim LWORK
// IFAIL on exit = 0 if no errors
// Note: the parameter convention here is valid for Mark 25 
extern "C" int e04ncf_(long &M,long &N,long &NCLIN,long &LDC,long &LDA,
		      double *C,double *BL,double *BU,double *CVEC,
		      long *ISTATE,long *KX,double *X,double *A,double *B,
		      long &ITER,double &OBJ,double *CLAMBDA,long *IWORK,long &LIWORK,
		      double *WORK,long &LWORK,long &IFAIL);

extern "C" int setopte04ncf_();

//***********************************************************************

void condense(Glob *globs,GlobExp *ex)
  //! Condenses each experiment 
{
  long nExp,msP,i,j,k,l;
  long nMeas,nvar,nP;
  long me,mg,nobs,nms;
  
  
  //weakens continuty constraints
  double elastic=globs->elastic;
  double **Ea,**Ee,**Eg; //save the previous E-matrices

  nP=globs->npar;
  for(nExp=1;nExp<=globs->nrExp;++nExp) 
    { 
      nMeas=ex[nExp].nMeasure; //# of measurement time points
      nvar=ex[nExp].nvar; //# of variables
      me=ex[nExp].me; //# of equality constraints
      mg=ex[nExp].mg; //# of inequality constraints
      nobs=ex[nExp].nobs; //# of observations
      nms=ex[nExp].nPoints; //# of shooting intervals
      nms--;
      Ea=dmatrix(1,nMeas*nobs,1,nvar);
      Ee=dmatrix(1,me,1,nvar);
      Eg=dmatrix(1,mg,1,nvar);

      //initialisation of backward recursion
      //------------------------------------
      //initialize least squares:
      for(i=1;i<=nMeas;i++)
	{
	  for(j=1;j<=nobs;j++)
	    {
	      //ua 
	      ex[nExp].ua[(j-1)*nMeas+i]=ex[nExp].residues[i][j]/ex[nExp].sigma[i][j];
	      //Ea
	      for(k=1;k<=nvar;k++)
		{
		  //find out which datapoints are falling in the last MS-interval
		  //the other derivatives are zero
		  if(ex[nExp].xMeasure[i] >= ex[nExp].mesh[nms] && 
		     ex[nExp].xMeasure[i] <= ex[nExp].mesh[nms+1])
		    {
		      ex[nExp].Ea[(j-1)*nMeas+i][k]=ex[nExp].dmds[i][j][k]/ex[nExp].sigma[i][j];
		    }
		  else
		    {
		      ex[nExp].Ea[(j-1)*nMeas+i][k]=0.0;
		    }
		  Ea[(j-1)*nMeas+i][k]=ex[nExp].Ea[(j-1)*nMeas+i][k];
		}
	      //Pa
	      for(k=1;k<=nP;k++)
		{
		  ex[nExp].Pa[(j-1)*nMeas+i][k]=ex[nExp].dmdp[i][j][k]/ex[nExp].sigma[i][j];
		}
	    }
	}// end of least squares initialization
      
      //equality constraints
      for(i=1;i<=me;i++)
	{
	  //ue
	  ex[nExp].ue[i]=ex[nExp].r2[i];
	  //Ee
	  for(j=1;j<=nvar;j++)
	    {
	      ex[nExp].Ee[i][j]=ex[nExp].dR2ds[i][nms][j];
	      Ee[i][j]=ex[nExp].Ee[i][j];
	    }
	  //Pe
	  for(j=1;j<=nP;j++)
	    ex[nExp].Pe[i][j]=ex[nExp].dR2dp[i][j];
	}//end of equality constraints

      //inequality constraints
      for(i=1;i<=mg;i++)
	{
	  //ug
	  ex[nExp].ug[i]=ex[nExp].r3[i];
	  //Ee
	  for(j=1;j<=nvar;j++)
	    {
	      ex[nExp].Eg[i][j]=ex[nExp].dR3ds[i][nms][j];
	      Eg[i][j]=ex[nExp].Eg[i][j];
	    }
	  //Pe
	  for(j=1;j<=nP;j++)
	    ex[nExp].Pg[i][j]=ex[nExp].dR3dp[i][j];
	}//end of inequality constraints
      //initialisation of backward recursion complete
    
      
      //backward recursion 
      //-----------------
      
      //loop over all multiple shooting intervals
      for(msP=nms;msP>=2;msP--)
	{
	  //least squares
	  
	  for(i=1;i<=nMeas;i++)
	    {
	      for(j=1;j<=nobs;j++)
		{
		  //ua 
		  for(k=1;k<=nvar;k++)
		    ex[nExp].ua[(j-1)*nMeas+i]+=ex[nExp].Ea[(j-1)*nMeas+i][k]*ex[nExp].h[msP][k]*elastic;
		  //Ea
		  for(k=1;k<=nvar;k++)
		    {
		      //find out which datapoints are falling in the MS-interval
		      //the other derivatives are zero
		      if(ex[nExp].xMeasure[i] >= ex[nExp].mesh[msP-1] && 
			 ex[nExp].xMeasure[i] <= ex[nExp].mesh[msP])
			{
			  ex[nExp].Ea[(j-1)*nMeas+i][k]=ex[nExp].dmds[i][j][k]/ex[nExp].sigma[i][j];
			}
		      else
			{
			  ex[nExp].Ea[(j-1)*nMeas+i][k]=0.0;
			}
		      for(l=1;l<=nvar;l++)
			ex[nExp].Ea[(j-1)*nMeas+i][k]+=Ea[(j-1)*nMeas+i][l]*ex[nExp].dyds[msP][l][k];
		    }
		  //Pa
		  for(k=1;k<=nP;k++)
		    {
		      for(l=1;l<=nvar;l++)
			{
			  ex[nExp].Pa[(j-1)*nMeas+i][k]+=Ea[(j-1)*nMeas+i][l]*ex[nExp].dydp[msP][l][k];
			}
		    }
		}
	    } //end of least squares

	  //equality constraints
	  for(i=1;i<=me;i++)
	    {
	      //ue
	      for(k=1;k<=nvar;k++)
		ex[nExp].ue[i]+=Ee[i][k]*ex[nExp].h[msP][k]*elastic;
	      //Ee
	      for(j=1;j<=nvar;j++)
		{
		  ex[nExp].Ee[i][j]=ex[nExp].dR2ds[i][msP-1][j];
		    for(l=1;l<=nvar;l++)
		      ex[nExp].Ee[i][j]+=Ee[i][l]*ex[nExp].dyds[msP][l][j];
		}
	      //Pe
	      for(j=1;j<=nP;j++)
		{
		  for(l=1;l<=nvar;l++)
		    ex[nExp].Pe[i][j]+=Ee[i][l]*ex[nExp].dydp[msP][l][j];
		}
	    }//end of equality constraints
	  
	  //inequality constraints
	  for(i=1;i<=mg;i++)
	    {
	      //ue
	      for(k=1;k<=nvar;k++)
		ex[nExp].ug[i]+=Eg[i][k]*ex[nExp].h[msP][k]*elastic;
	      //Ee
	      for(j=1;j<=nvar;j++)
		{
		  ex[nExp].Eg[i][j]=ex[nExp].dR3ds[i][msP-1][j];
		    for(l=1;l<=nvar;l++)
		      ex[nExp].Eg[i][j]+=Eg[i][l]*ex[nExp].dyds[msP][l][j];
		}
	      //Pe
	      for(j=1;j<=nP;j++)
		{
		  for(l=1;l<=nvar;l++)
		    ex[nExp].Pg[i][j]+=Eg[i][l]*ex[nExp].dydp[msP][l][j];
		}
	    }//end of inequality constraints

	  //updating new Ea,Ee,Ea 
	  for(i=1;i<=nvar;i++)
	    {
	      for(j=1;j<=nMeas*nobs;j++)
		Ea[j][i]=ex[nExp].Ea[j][i];
	      for(j=1;j<=me;j++)
		Ee[j][i]=ex[nExp].Ee[j][i];
	      for(j=1;j<=mg;j++)
		Eg[j][i]=ex[nExp].Eg[j][i];	      
	    }

	}//end loop over multiple shooting intervals

      //clearing memory
      free_dmatrix(Ea,1,nMeas*nobs,1,nvar);
      free_dmatrix(Ee,1,me,1,nvar);
      free_dmatrix(Eg,1,mg,1,nvar);
    }//end loop over all experiments
}
 

void decondense(Glob *globs,GlobExp *ex,double *dX)
  //! Decondenses each experiment 
{
  long nExp,msP,i,j,k,l;
  long nvar,npar;
  long nms,ind=0;

  //weakens continuty constraints
  double elastic=globs->elastic;

  npar=globs->npar;     
  nvar=ex[1].nvar;
  
  //initialise
  for(nExp=1;nExp<=globs->nrExp;++nExp) 
    {
      for(i=1;i<=ex[nExp].nPoints;i++)
	for(j=1;j<=nvar;j++)
	  ex[nExp].dS[i][j]=0.;
      for(i=1;i<=npar;i++)
	ex[nExp].dP[i]=0.;

      for(i=1;i<=nvar;i++)
	{
	  if(globs->y0fix[i]!=FALSE)
	    {
	      ex[nExp].dS[1][i]=dX[ind];
	      ind++;
	    }
	}
      for(i=1;i<=npar;i++)
	{
	  if(globs->doP[i]=='L')
	    {
	      ex[nExp].dP[i]=dX[ind];
	      ind++;
	    }
	}
    }

  for(i=1;i<=npar;i++)
    {
      if(globs->doP[i]==TRUE)
	{
	  for(nExp=1;nExp<=globs->nrExp;++nExp) 
	    {
	      ex[nExp].dP[i]=dX[ind];
	    }
	  ind++;
	}
    } 
  
  //end initialise
  
  //forward iteration
  for(nExp=1;nExp<=globs->nrExp;++nExp) 
    { 
      nms=ex[nExp].nPoints;
      nms--;
      
      for(msP=2;msP<=nms;msP++)
	{
	  for(i=1;i<=nvar;i++)
	    {
	      ex[nExp].dS[msP][i]=elastic*ex[nExp].h[msP][i];
	      for(j=1;j<=nvar;j++)
		ex[nExp].dS[msP][i]+=ex[nExp].dyds[msP][i][j]*ex[nExp].dS[msP-1][j];
	      for(j=1;j<=npar;j++)
		ex[nExp].dS[msP][i]+=ex[nExp].dydp[msP][i][j]*ex[nExp].dP[j];
	    }
	}
      
    }//end loop over experiments
}
 
void solvLin(Glob *globs,GlobExp *ex,int computeCovar)
{
  
  long i,j,k,l,xind=0,yind=0,nL;
  long nvarFit=0;  // # der zu fittenden Anfangswerte
  long nLpara=0;   // # lokale Parameter
  long nparL=0;    // nLpara*nrExp
  long npar=globs->npar;
  long nvar=ex[1].nvar; //number of ODEs
  long nrExp=globs->nrExp;
  long nparG=0;
  long fitDim=0;
  long aDim=0,eDim=0,gDim=0;

  //big system (condensed)
  double **Ma,*Ra;
  double **Me,*Re;
  double **Mg,*Rg;

  //solution
  double *dX;
  //condense experiments
  condense(globs,ex);

  //determine dimension of big system
  for(i=1;i<=nvar;i++)
    {
      if(globs->y0fix[i]!=FALSE)
	nvarFit++;
    }
  nvarFit*=nrExp;
  for(i=1;i<=npar;i++)
    {
      if(globs->doP[i]=='L')
	nparL++; //local parameters
      else if(globs->doP[i]!=FALSE)
	nparG++; //global parameters
    }
  nLpara=nparL;
  nparL*=nrExp;
  fitDim=nvarFit+nparL+nparG;
  for(i=1;i<=nrExp;i++)
    {
      aDim+=ex[i].nMeasure*ex[i].nobs;
      eDim+=ex[i].me;
      gDim+=ex[i].mg;
    }
  if(fitDim==0)
    {
      cerr << "Nothing to fit.\n";
      throw 1;
    }
#ifdef PRINTDIMENSION
  //print dimension of bigsys;
  //*dbg << "Dimension of Big-System\n";
  //*dbg << "-----------------------\n";
  //*dbg << "Least-squares:\n";
  //*dbg << "#Rows : " << aDim  << "  #Columns : " << fitDim << endl;
  //*dbg << "Equality constraints:\n";
  //*dbg << "#Rows : " << eDim  << "  #Columns : " << fitDim << endl;  
  //*dbg << "Inequality constraints:\n";
  //*dbg << "#Rows : " << gDim  << "  #Columns : " << fitDim;  
  //*dbg << "\n\n";
  //dbgflush();
#endif


  long aDimold=aDim;
  long eDimold=eDim;
  long fitDimold=fitDim;

  // Beginn Local_Parameter_Constraints
  if(globs->faktorLexist)
  {
  	// Erweiterung der Dimensionen um nLpara*nExp
  	aDim=aDim+nparL;
  	eDim=eDim+nparL;
  	fitDim=fitDim+nparL;
  }
  // Ende Local_Parameter_Constraints
  
  //allocate memory
  Ma=dmatrix(0,aDim,1,fitDim);
  Me=dmatrix(0,eDim,1,fitDim);
  Mg=dmatrix(0,gDim,1,fitDim);
  Ra=dvector(0,aDim);
  Re=dvector(0,eDim);
  Rg=dvector(0,gDim);

  //initialise 
  for(i=1;i<=aDim;i++)
    {
      Ra[i]=0.;
      for(j=1;j<=fitDim;j++)
	Ma[i][j]=0.;
    }
   for(i=1;i<=eDim;i++)
    {
      Re[i]=0.;
      for(j=1;j<=fitDim;j++)
	Me[i][j]=0.;
    } 
   for(i=1;i<=gDim;i++)
     {
      Rg[i]=0.;
      for(j=1;j<=fitDim;j++)
	Mg[i][j]=0.;
    } 

  //least-squares
  for(i=1;i<=nrExp;i++)
    {
      for(k=1;k<=ex[i].nMeasure*ex[i].nobs;k++)    
	Ra[yind+k]=-ex[i].ua[k];
      //initial values
      for(j=1;j<=nvar;j++)
	{
	  if(globs->y0fix[j]!=FALSE)
	    {
	      xind++;
	      for(k=1;k<=ex[i].nMeasure*ex[i].nobs;k++)    
		Ma[yind+k][xind]=ex[i].Ea[k][j];
	    }
	  
	}
      //parameters
      nL=0;
      for(j=1;j<=npar;j++)
	{
	  if(globs->doP[j]=='L')
	    {
	      xind++;
	      for(k=1;k<=ex[i].nMeasure*ex[i].nobs;k++)
		Ma[yind+k][xind]=ex[i].Pa[k][j];
	    }
	  else if(globs->doP[j]!=FALSE)
	    {
	      nL++;
	      for(k=1;k<=ex[i].nMeasure*ex[i].nobs;k++)
		Ma[yind+k][nvarFit+nparL+nL]=ex[i].Pa[k][j];
	    }
	}
      yind=yind+(ex[i].nMeasure*ex[i].nobs);
    }//end loop over experiments
  //-------------> end least squares
  
  //equality constraints
  yind=0;
  xind=0;
  for(i=1;i<=nrExp;i++)
    {
      for(k=1;k<=ex[i].me;k++)    
	Re[yind+k]=-ex[i].ue[k];
      //initial values
      for(j=1;j<=nvar;j++)
	{
	  if(globs->y0fix[j]!=FALSE)
	    {
	      xind++;
	      for(k=1;k<=ex[i].me;k++)    
		Me[yind+k][xind]=ex[i].Ee[k][xind];
	    }
	}
      //parameters
      nL=0;
      for(j=1;j<=npar;j++)
	{
	  if(globs->doP[j]=='L')
	    {
	      xind++;
	      for(k=1;k<=ex[i].me;k++)
		Me[yind+k][xind]=ex[i].Pe[k][j];
	    }
	  else if(globs->doP[j]!=FALSE)
	    {
	      nL++;
	      for(k=1;k<=ex[i].me;k++)
		Me[yind+k][nvarFit+nparL+nL]=ex[i].Pe[k][j];
	    }
	}
      yind+=ex[i].me;
    }//end loop over experiments
  //-------------> end equality  constraints

  //inequality constraints
  yind=0;
  xind=0;
  for(i=1;i<=nrExp;i++)
    {
      for(k=1;k<=ex[i].mg;k++)    
	Rg[yind+k]=-ex[i].ug[k];
      //initial values
      for(j=1;j<=nvar;j++)
	{
	  if(globs->y0fix[j]!=FALSE)
	    {
	      xind++;
	      for(k=1;k<=ex[i].mg;k++)  
		Mg[yind+k][xind]=ex[i].Eg[k][xind];
	    }
	}
      //parameters
      nL=0;
      for(j=1;j<=npar;j++)
	{
	  if(globs->doP[j]=='L')
	    {
	      xind++;
	      for(k=1;k<=ex[i].mg;k++)
		Mg[yind+k][xind]=ex[i].Pg[k][j];
	    }
	  else if(globs->doP[j]!=FALSE)
	    {
	      nL++;
	      for(k=1;k<=ex[i].mg;k++)
		Mg[yind+k][nvarFit+nparL+nL]=ex[i].Pg[k][j];
	    }
	}
      yind+=ex[i].mg;
    }//end loop over experiments
  //-------------> end equality  constraints
  
  //Big System Ready to Solve
  //-------------------------
  //
  //  Ma*dX  ~  Ra
  //  Me*dX  =  Re
  //  Mg*dX >=  Rg
  //
  //-------------------------
  
  //****************************************************************************
  // Begin Local_Parameter_Constraints
  // Erweiterung der Systems auf Local_Parameter_Constraints:
  //
  // pLmean und pLsig bestimmen
  if(globs->faktorLexist)
  {
  	double *pLmean, *pLsig;
  	int index;
  	int index_local_parameter;
	pLmean=dvector(1,nLpara);
  	pLsig=dvector(1,nLpara);
  
//   	for(i=1;i<=nLpara;i++)
//   	{
//   	   cerr<<"i="<< i<<" faktorL[i]="<<globs->faktorL[i]<<endl;
//   	   if(globs->faktorL[i]<0)
//   	   {
//   	      //i;
//   	      break;
//   	   }
//   	}
//   	cerr<<"i="<< i<<" faktorL[i]="<<globs->faktorL[i]<<endl;
//   	if(i<=nLpara)
//   	{
//   	   cerr<<"number of local parameter constraints too small"<<endl;
//   	   throw(1);
//   	}
  	for(i=1;i<=nLpara;i++)
  	{
  		pLmean[i]=0;
  		index=0;
  		index_local_parameter=0;
  		while(index_local_parameter!=i)
  		{
  			index++;
//   			cerr<<"doP[index]="<<globs->doP[index]<<" index="<<index<<" i="<<i<<endl;
	  		if(globs->doP[index]=='L')
  			{
  				index_local_parameter++;
//   				cerr<<"index="<<index<<" i="<<i<<endl;
  			}
  			
  		}  		
  		for(j=1;j<=nrExp;j++)
  		{
	  		pLmean[i]=pLmean[i]+ex[j].par[index]; 
// 	  		cerr<<"ex[j].par[index]="<<ex[j].par[index]<<endl;
	  	}
	  	pLmean[i]=pLmean[i]/double(nrExp);
	  //	cerr<<"pLmean["<<i<<"]="<<pLmean[i]<<" globs->faktorL["<<i<<"]"<<globs->faktorL[i]<<endl;
  	}
  	// Berechnung pLsig
  	for(i=1;i<=nLpara;i++)
  	{
	  	pLsig[i]=(globs->faktorL[i]*pLmean[i]-pLmean[i])/(1+globs->faktorL[i]);
	  	if(pLsig[i]<0)
	  	{
	  	cerr<<"ERROR - Negative varianz bei den Lokalen Parameter Constraints..."<<endl;
	  	}
  	}  	
  	// Ra,Re
  	k=1;
  	for(i=1;i<=nrExp;i++)
  	{
	  	for(j=1;j<=nLpara;j++)
	  	{	
  			index=0;
  			index_local_parameter=0;
  			while(index_local_parameter!=j)
  			{
	  			index++;
  				if(globs->doP[index]=='L')
  				{
	  				index_local_parameter++;
  				}
  			
  			}  
  			Ra[aDimold+k]=(ex[i].par[index]-pLmean[j])/pLsig[j];
  			Re[eDimold+k]=0;
  			double pLsig_real=0;
  			for(int z=1;z<=nrExp;z++)
  			{
  				pLsig_real=pLsig_real+(ex[z].par[index]-pLmean[j])*(ex[z].par[index]-pLmean[j]);
  			}
  			pLsig_real=sqrt(pLsig_real/nrExp);
  			////*dbg<<"ex["<<i<<"].par["<<index<<"]="<<ex[i].par[index]<<" pLmean["<<j<<"]="<<pLmean[j]<<" pLsig["<<j<<"]="<<pLsig[j]<<" pLsig_real="<<pLsig_real<<endl;
  			// Ma
  			for(int l=1;l<=nrExp;l++)
  			{
	  			Ma[aDimold+(i-1)*nLpara+j][fitDimold+(l-1)*nLpara+j]=1/(pLsig[j]*nrExp);
  				if(i == l)
  				{
	  				Ma[aDimold+(i-1)*nLpara+j][fitDimold+(l-1)*nLpara+j]=Ma[aDimold+(i-1)*nLpara+j][fitDimold+(l-1)*nLpara+j]-1/pLsig[j];
  				}
  			}
  			k++;
  		}
  	}
  	// Me
  	for(i=eDimold+1;i<=eDim;i++)
  	{
	  	for(j=nvarFit+1;j<=fitDim;j++)   
	  	{
  			if(j==fitDimold+i)
  			{
	  			Me[i][j]=Me[i][j]+1;
  			}
  			if(j==nvarFit+i)
  			{
	  			Me[i][j]=Me[i][j]-1;
  			}
  		}
  	}
  }
  // Ende Local_Parameter_Constraints
  //****************************************************************************
  
  //regularisation; Note: default is no regularisation
  if(globs->reg==TRUE)
    {
      double *diag,**V,**MMa;
      double cond,maxx;
      double reg=0.;

      V=dmatrix(1,fitDim,1,fitDim);
      MMa=dmatrix(1,aDim,1,fitDim);
      diag=dvector(1,fitDim);

      for(i=1;i<=aDim;i++)
	{
	  for(j=1;j<=fitDim;j++)
	    {
	      MMa[i][j]=Ma[i][j];
	    }     
	}
      
      //SVD
      svdcmp(MMa,aDim,fitDim,diag,V);
      
      //determine condition number
      //*dbg << "\ncondition number : ";
      maxx=diag[1];
      cond=diag[1];
      for(i=1;i<=fitDim;i++)
	{
	  if(diag[i]>maxx)
	    maxx=diag[i];
	  if(diag[i]<cond)
	    cond=diag[i];
	}

      //Is the matrix approximately singular?
      for(i=1;i<=fitDim;i++)
	{
	  if(fabs(diag[i]/maxx)<=globs->epsilon)
	    {
	      diag[i]=0.;
	      reg=globs->lambda*maxx;
	    }
	}
      
      cond*=1./maxx;
      //*dbg << cond << endl;
      globs->cond=cond;
      
      //writing regularised system back
      for(i=1;i<=aDim;i++)
	{
	  for(j=1;j<=fitDim;j++)
	    {
	      Ma[i][j]=0.;
	      for(k=1;k<=fitDim;k++)
		{
		  Ma[i][j]+=MMa[i][k]*diag[k]*V[j][k];
		}
	      if(i <=fitDim)
		{
		  if(diag[i]==0.)
		    Ma[i][j]+=reg*V[j][i];
		}
	    }
	}
     /* if(reg!=0.)
	{
	  //*dbg << "regularisation active, linear dependency of parameters:\n";
	  for(i=1;i<=fitDim;i++)
	      {
	        if(diag[i]==0.)
		  //*dbg << "line " << i << ": ";
		for(j=1;j<=fitDim;j++)
		  {
		    if(diag[i]==0.)
		      //*dbg << V[j][i] << " ";
		  }
		if(diag[i]==0.)
		  //*dbg << endl;
	      }
	}*/

      free_dmatrix(MMa,1,aDim,1,fitDim);
      free_dmatrix(V,1,fitDim,1,fitDim);
      free_dvector(diag,1,fitDim);
    }//end if(globs->reg==TRUE)

#ifdef PRINTBIGSYS
  //*dbg << "\nLeast-squares :\n";
   for(i=1;i<=aDim;i++)
    {
      for(j=1;j<=fitDim;j++)
	{
	  //*dbg << Ma[i][j] << " ";
	}
      //*dbg << "| " << Ra[i] << endl;
    }  
   //*dbg << "\nEquality constraints :\n";
   for(i=1;i<=eDim;i++)
    {
      for(j=1;j<=fitDim;j++)
	{
	  //*dbg << Me[i][j] << " ";
	}
      //*dbg << "| " << Re[i] << endl;
    }
   //*dbg << "\nInequality constraints :\n";
   for(i=1;i<=gDim;i++)
    {
      for(j=1;j<=fitDim;j++)
	{
	  //*dbg << Mg[i][j] << " ";
	}
      //*dbg << "| " << Rg[i] << endl;
    }
   //dbgflush();
#endif
  //allocate solution vector
  dX=(double*)malloc((fitDim+10)*sizeof(double));
  //allocate covariance matrix
  globs->fitdim=fitDim;
  if(globs->covar==NULL)
    globs->covar=dmatrix(1,fitDim,1,fitDim);

  if(globs->minimiser==1)
    {
      //*dbg << "right before LSEI \n";
      //dbgflush();
      // ###############
      // #     LSEI    #
      // ###############
      
      //allocation/initialisation for LSEI
      double *w;
      double *ws;
      double *prgopt;
      double rnorme = 0,rnorml = 0;
      int32_t mdw= (int32_t) (aDim+eDim+gDim); //long mdw=aDim+eDim+gDim;
      int32_t modeLSEI = 0; //long modeLSEI = 0;
      int32_t ipstor = gDim + 2* ((int32_t) fitDim) + 2; //long ipstor=mdw+npar+2*fitDim+2;
      int32_t wsstor=2*(eDim+ (int32_t) fitDim)+(aDim+ (int32_t) gDim)+(gDim + (int32_t) 2)*((int32_t) fitDim + 7);
      //long wsstor=3*(mdw+npar)+(mdw+npar+2)*(fitDim+8)+fitDim+npar;
      int32_t *ip; //long *ip;

      //recast explicitly as 32 bit integers for passing to lsei_wrapper
      int32_t aDim32 = (int32_t) aDim;
      int32_t eDim32 = (int32_t) eDim;
      int32_t gDim32 = (int32_t) gDim;
      int32_t fitDim32 = (int32_t) fitDim;
      
      long wdim = (fitDim+1)*mdw;
      int32_t wsdim = wsstor+1; //long wsdim = wsstor+1;
      w=(double*)malloc(((fitDim+1)*mdw)*sizeof(double));
      ws=(double*)malloc((wsstor+1)*sizeof(double));
      //ip=(long*)malloc((ipstor+10)*sizeof(long));
      ip=(int32_t*)malloc((ipstor+10)*sizeof(int32_t));
      prgopt=(double*)malloc(30*sizeof(double));

      //*dbg << "initialization \n";
      //dbgflush();


      //initialization helps when passing to lsei_wrapper
      for (i=0;i<30;i++)
	prgopt[i]=0.0;

      for (i=0;i<wdim;i++)
	w[i] = 0.0;

      for (i=0;i<wsdim;i++)
	ws[i] = 0.0;
      
      xind=0;
      yind=0;
      
      for(i=1;i<=fitDim;i++)
	{
	  for(j=1;j<=eDim;j++)
	    {
	      w[xind*mdw+yind]=Me[j][i];
	      yind++;
	    }
	  for(j=1;j<=aDim;j++)
	    {
	      w[xind*mdw+yind]=Ma[j][i];
	      yind++;
	    }
	  for(j=1;j<=gDim;j++)
	    {
	      w[xind*mdw+yind]=Mg[j][i];
	      yind++;
	    } 
	  yind=0;
	  xind++;
	}
 
      for(j=1;j<=eDim;j++)
	{
	  w[(fitDim)*mdw+yind]=Re[j];
	  yind++;
	}
      for(j=1;j<=aDim;j++)
	{
	  w[(fitDim)*mdw+yind]=Ra[j];
	  yind++;
	}
      for(j=1;j<=gDim;j++)
	{
	  w[(fitDim)*mdw+yind]=Rg[j];
	  yind++;
	}
      
      ip[0]=wsstor;
      ip[1]=ipstor;
      ip[2] = 0;
      
      if(computeCovar==TRUE)
	{
	  prgopt[0]=4.0; //link to next option index 3 since we start at zero
	  prgopt[1]=1.0; //covariance matrix key; compute covariance matrix
	  prgopt[2]=1.0; //option data, must be 1 when covariance matrix is requested
	  prgopt[3]=7.0; //link to next option, index 6
	  prgopt[4]=10.0; 
	  prgopt[5]=1.0;

	  //prgopt[3] = 1.0; //FG tentative, new end of list
	}
      else
	{
	  prgopt[0]=4.0; 
	  prgopt[1]=0.0; 
	  prgopt[2]=1.0;
	  prgopt[3]=7.0;
	  prgopt[4]=10.0;
	  prgopt[5]=0.0;

	  //prgopt[3] = 1.0; //FG tentative, new end of list
	}
      prgopt[6]=1.0; //end of list

      *dbg << "at LSEI wrapper \n";
      *dbg << sizeof(*w)/sizeof(w[0]);
      *dbg << " :w \n";
      *dbg << sizeof(*ws)/sizeof(ws[0]);
      *dbg << " :ws \n";
      *dbg << sizeof(*ip)/sizeof(ip[0]);
      *dbg << " :ip \n";
      dbg->flush();



      //NOte use converted data to 32 bit integers in the call to the wrapper
      
      //call LSEI
      //FG082918 lsei_(w,mdw,eDim,aDim,gDim,fitDim,prgopt,dX,&rnorme,&rnorml,&mode,ws,ip);
      //lsei(w,mdw,eDim,aDim,gDim,fitDim,prgopt,dX,&rnorme,&rnorml,&mode,ws,ip);
      //lsei_wrapper(w,mdw,eDim32,aDim32,gDim32,fitDim32,prgopt,dX,&rnorme,&rnorml,modeLSEI,ws,wsdim,ip,ipstor);

      *dbg << "after LSEI wrapper \n";
      
      if(modeLSEI>=4)
	{
	  cerr << "Minimisation (LSEI) failed.\n";
	  throw 1;
	}
      
      if(computeCovar==TRUE)
	{
	  for (i=1; i<=fitDim; ++i)
	    {
	      for (j=1;j<=fitDim;++j)
		globs->covar[i][j] = 0.1; //w[(j-1)*mdw+i-1];
	    } 
	}
      
      free(w);
      free(ws);
      free(ip);
      free(prgopt);
      
      // #######  END OF LSEI ######
    }
  else if(globs->minimiser==2)
    {
      //to solve ----->
      //-------------------------
      //
      //  Ma*dX  ~  Ra
      //  Me*dX  =  Re
      //  Mg*dX >=  Rg
      //
      //-------------------------

      // ###############
      // #     NAG     #
      // ###############
        long m = aDim, n = fitDim, nclin = eDim + gDim, nbnd = n + nclin, tda = n, tdh = n;
	long *kx = (long*)malloc(n*sizeof(long));
  	double *a = (double*)malloc(nclin*tda*sizeof(double));
	double *bl = (double*)malloc((nclin+n)*sizeof(double));
	double *bu = (double*)malloc((nclin+n)*sizeof(double));
       	double *cvec = (double*)malloc(n*sizeof(double));
        double *b = (double*)malloc(m*sizeof(double));
        double *h = (double*)malloc(m*tdh*sizeof(double));
	double *x = (double*)malloc(n*sizeof(double));
        double objf;
	
  	Nag_Comm comm;
  	NagError fail;
	INIT_FAIL(fail);
	Nag_E04_Opt options; e04xxc(&options);
	options.list = Nag_FALSE;
      	options.print_level = Nag_NoPrint;
	
   
      //right-hand side
     	for(i=0;i<m;++i)
	  {
		b[i]=Ra[i+1];
	  }

       	for(i=0;i<m;++i)
	  {
	  for(j=0;j<n;++j)
	    {
		h[i*tdh + j]=Ma[i+1][j+1];
	    }
	  }

       	for(i=0;i<nclin;++i)
      	  {
	    for(j=0;j<n;++j)
	      {
	        if(i<eDim)
		  a[i*tda + j]=Me[i+1][j+1];
	        else
		  a[i*tda + j]=Mg[i-eDim+1][j+1];
	      }
	  }
      //bounds
      	for(i=0;i<n;++i)
	  {
	    bl[i]=-1e25;
	    bu[i]=1e25;
	    x[i]=0.;  //initialisation of X
	  }

      	for(i=n;i<n+eDim;++i)
	  {
	    bl[i]=Re[i-n+1]; //equality constraints
	    bu[i]=Re[i-n+1];
	  }
      	for(i=n+eDim;i<n+eDim+gDim;++i)
	  {
	    bl[i]=Rg[i-n-eDim+1]; //inequality constraints
	    bu[i]=1e25;
	  }
	
      e04ncc(m,n,nclin,a,tda,bl,bu,cvec,b,h,tdh,kx,x,&objf,&options,&comm,&fail);

      if (fail.code != NE_NOERROR) {
      	cout << "Error from nag_opt_lin_lsq (e04ncc).\n%s\n" << fail.message;
	}

      for(i=0;i<n;++i)
      {
	dX[i]=x[i];
      }
      
      //free memory
      free(a);
      free(bl);
      free(bu);
      free(cvec);
      free(kx);
      free(x);
      free(h);
      free(b);
  
    } // #######  END OF NAG ######

  //decondense experiments
  decondense(globs,ex,dX);

#ifdef PRINTSTEP
  for(i=1;i<=nrExp;i++)
    {
      //*dbg << "\nExperiment #" << i << ":\n\n";
      for(k=1;k<=ex[i].nPoints-1;k++)
	{
	  //*dbg << "dS(msP=" << k << ") = ";
	  for(j=1;j<=nvar;j++)
	    //*dbg << ex[i].dS[k][j] << " ";
	  //*dbg << endl;
	}
      //*dbg << "dP = ";
      for(j=1;j<=npar;j++)
	//*dbg << ex[i].dP[j] << " ";
    }
  //*dbg << endl;
  //dbgflush();
#endif
  

  //free memory
  free_dmatrix(Ma,0,aDim,1,fitDim);
  free_dmatrix(Me,0,eDim,1,fitDim);
  free_dmatrix(Mg,0,gDim,1,fitDim);
  free_dvector(Ra,0,aDim);
  free_dvector(Re,0,eDim);
  free_dvector(Rg,0,gDim);

  dX[0];
  free(dX);
}

/******************************************************************************
 * error routines called from Fortran modules constr.f and odessa.f
 */

// for constr.f
extern "C" void
myerr_ (int *messg, int &nmessg)
{
  char *s = (char *) messg;
  cerr << endl;
  for (int i = 0; i < nmessg; i++)
    cerr << s[i];
  cerr << endl;
}

/* for odessa.f
 * nerr ignored
 * iert==2 checked, otherwise ignored
 * e.g.: CALL XERR ('AT T=R1, MXSTEP=I1 STEPS WERE TAKEN... ',
     1   201, 2, 1, MXSTEP, 0, 1, TN, ZERO)
 */
extern "C" void
xerr_ (char *msg, int &nerr, int &iert, int &ni, int &i1,
       int &i2, int &nr, double &r1, double &r2)
{
  cerr << msg << endl;
  if (ni > 0)
    cerr << i1 << " ";
  if (ni > 1)
    cerr << i2 << " ";
  if (nr > 0)
    cerr << r1 << " ";
  if (nr > 1)
    cerr << r2 << " ";
  if (ni > 0 || nr > 0)
    cerr << endl;
  if (iert >= 2)
    cerr << endl;
}
