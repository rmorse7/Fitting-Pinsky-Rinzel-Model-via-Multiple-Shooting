#include <math.h>
#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <string.h>

#include "nr.h"
#include "def.h"

using namespace std;

extern const char *DefModelDescription = 
"Implements the Morris Lecar model used in the NIPS18 submission.";

extern const int NEQNS   = 2;
extern const int NOBS    = 1;
extern const int NPARAMS = 4;
extern const int NSPLINES  = 0;

//Spline data
double **splineNodes;
double **splineY;
double **splineGam;
unsigned long *nNodes;


const int	iu	=1;
const int	iw	=2;
#define dudt	f[1]
#define dwdt	f[2]
#define u	y[1]
#define w	y[2]
#define gCa	p[1]
#define gK	p[2]
#define gL	p[3]
#define twMax	p[4]


extern const double DefParameters[] = {4.0,8.0,2.0,15.0};
extern const double DefYValues[]    = {-30.0,0.0};
extern const string VariableNames[] = {"u","w"};
extern const string ParameterNames[] = {"gCa","gK","gL","twMax"};

void ode(double t, double *y, double *f, double *p)
{
  dudt=-gCa*(0.5*(1+tanh((u-(-1.2))/(18))))*(u-(120))-gK*w*(u-(-84))-gL*(u-(-60))+42;
  dwdt=((0.5*(1+tanh((u-(12))/(17.4))))-w)/(twMax/cosh((u-(12))/(2*(17.4))));
} // ode

void unobsinit(double *y, double *p)
{
  w = 0.5*(1+tanh((u-(12))/(17.4)));
}

/* set ex->r3[1..data->mg]; (ineq. constr., components must be >= 0) */
void R3(Glob *globs, GlobExp *ex, int computeDerivs)
{
  long i,j;
  // parameter constraints
  // parameter: gCa
  ex->r3[1]=ex->par[1]-(3.6);
  if(computeDerivs) {
    for (i=1;i<=globs->npar;i++) {
      if (i == 1)
        ex->dR3dp[1][i]=1.0;
     }
  }
  ex->r3[2]=-ex->par[1]+(4.4);
  if(computeDerivs) {
    for (i=1;i<=globs->npar;i++) {
      if (i == 1)
        ex->dR3dp[2][i]=-1.0;
    }
  }
  // parameter: gK
  ex->r3[3]=ex->par[2]-(1.0);
  if(computeDerivs) {
    for (i=1;i<=globs->npar;i++) {
      if (i == 2)
        ex->dR3dp[3][i]=1.0;
     }
  }
  ex->r3[4]=-ex->par[2]+(25.0);
  if(computeDerivs) {
    for (i=1;i<=globs->npar;i++) {
      if (i == 2)
        ex->dR3dp[4][i]=-1.0;
    }
  }
  // parameter: gL
  ex->r3[5]=ex->par[3]-(0.1);
  if(computeDerivs) {
    for (i=1;i<=globs->npar;i++) {
      if (i == 3)
        ex->dR3dp[5][i]=1.0;
     }
  }
  ex->r3[6]=-ex->par[3]+(10.0);
  if(computeDerivs) {
    for (i=1;i<=globs->npar;i++) {
      if (i == 3)
        ex->dR3dp[6][i]=-1.0;
    }
  }
  // parameter: twMax
  ex->r3[7]=ex->par[4]-(5.0);
  if(computeDerivs) {
    for (i=1;i<=globs->npar;i++) {
      if (i == 4)
        ex->dR3dp[7][i]=1.0;
     }
  }
  ex->r3[8]=-ex->par[4]+(25.0);
  if(computeDerivs) {
    for (i=1;i<=globs->npar;i++) {
      if (i == 4)
        ex->dR3dp[8][i]=-1.0;
    }
  }
  // variable constraints
  // variable: u
  ex->r3[9]=ex->yTry[1][1]-(-85.0);
  if(computeDerivs) {
    for (i=1;i<=ex->nvar;i++) {
      if (i == 1)
        ex->dR3ds[9][1][i]=1.0;
     }
  }
  ex->r3[10]=-ex->yTry[1][1]+(125.0);
  if(computeDerivs) {
    for (i=1;i<=ex->nvar;i++) {
      if (i == 1)
        ex->dR3ds[10][1][i]=-1.0;
    }
  }
  // variable: w
  ex->r3[11]=ex->yTry[1][2]-(0.0);
  if(computeDerivs) {
    for (i=1;i<=ex->nvar;i++) {
      if (i == 2)
        ex->dR3ds[11][1][i]=1.0;
     }
  }
  ex->r3[12]=-ex->yTry[1][2]+(1.0);
  if(computeDerivs) {
    for (i=1;i<=ex->nvar;i++) {
      if (i == 2)
        ex->dR3ds[12][1][i]=-1.0;
    }
  }
}
void setNrConstraints(GlobExp *ex, int nP, double *parameters)
{
  ex->me=0; ex->mg=12;
}
int observation (Glob *globs,GlobExp *ex,double t, double *y, double *gy,double *p, double **dgdy,double **dgdp)
{
  int generic=FALSE; // user defined observation function
  long k;
  gy[1]=u;
  
  if(dgdy) {
    dgdy[1][1]=1.0;
    dgdy[1][2]=0.0;
  }

  if(dgdp) {
     dgdp[1][1]=0.0;
     dgdp[1][2]=0.0;
     dgdp[1][3]=0.0;
     dgdp[1][4]=0.0;
  }

  return(generic);
}
void initInt(Glob *globs,GlobExp *ex)
{
 splineNodes=ex->splineNodes;
 splineY=ex->splineY;
 splineGam=ex->splineGam;
 nNodes=ex->nNodes;
}

/* set ex->r2[1..data->me]; (equal. constr., components must be == 0) */
void R2(Glob *globs,GlobExp *ex,int computeDerivs)
{

}

void jacobi(double t, double *y, double **jac, double *p) {
  jac[iu][iu]=-gK*w-( 5.0000000e-01*tanh( u/18.0+6.6666667e-02)+5.0000000e-01)*gCa-gL-( u-120.0)*( -2.7777778e-02*pow(tanh( u/18.0+6.6666667e-02),2.0)+2.7777778e-02)*gCa;
  jac[iw][iu]=-gK*( u+84.0);
  jac[iu][iw]= 2.8735632e-02*1.0/twMax*( 5.0000000e-01*tanh( 5.7471264e-02*u-6.8965517e-01)-w+5.0000000e-01)*sinh( 2.8735632e-02*u-3.4482759e-01)+cosh( 2.8735632e-02*u-3.4482759e-01)/twMax*( -2.8735632e-02*pow(tanh( 5.7471264e-02*u-6.8965517e-01),2.0)+2.8735632e-02);
  jac[iw][iw]=-cosh( 2.8735632e-02*u-3.4482759e-01)/twMax;

} // automatically generated elements

void inhomo1(double t, double *y, double *inh1, double *p, int &jp)
{
if(jp==1){
   inh1[iu]=-( u-120.0)*( 5.0000000e-01*tanh( u/18.0+6.6666667e-02)+5.0000000e-01);
   inh1[iw]=0.0;
}
if(jp==2){
   inh1[iu]=-w*( u+84.0);
   inh1[iw]=0.0;
}
if(jp==3){
   inh1[iu]=-u-60.0;
   inh1[iw]=0.0;
}
if(jp==4){
   inh1[iu]=0.0;
   inh1[iw]=-cosh( 2.8735632e-02*u-3.4482759e-01)/(twMax*twMax)*( 5.0000000e-01*tanh( 5.7471264e-02*u-6.8965517e-01)-w+5.0000000e-01);
}
} // automatically generated elements
void inhomo(double t, double *y, double **inh, double *p)
{
  int jp=1;

  for (int jp=1; jp<=NPARAMS; jp++)
      inhomo1(t,y,inh[jp],p,jp);
}

	
