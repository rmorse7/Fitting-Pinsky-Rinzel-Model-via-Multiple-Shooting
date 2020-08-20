#include <math.h>
#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <string.h>

#include "nr.h"
#include "def.h"

using namespace std;

extern const char *DefModelDescription = 
"dlfb";

extern const int NEQNS   = 4;
extern const int NOBS    = 2;
extern const int NPARAMS = 1;
extern const int NSPLINES  = 0;

//Spline data
double **splineNodes;
double **splineY;
double **splineGam;
unsigned long *nNodes;


const int	ip00	=1;
const int	ip10	=2;
const int	ip01	=3;
const int	ip11	=4;
#define dp00dt	f[1]
#define dp10dt	f[2]
#define dp01dt	f[3]
#define dp11dt	f[4]
#define p00	y[1]
#define p10	y[2]
#define p01	y[3]
#define p11	y[4]
#define gam	p[1]


extern const double DefParameters[] = {40};
extern const double DefYValues[]    = {1,0,0,0};
extern const string VariableNames[] = {"p00","p10","p01","p11"};
extern const string ParameterNames[] = {"gam"};

void ode(double t, double *y, double *f, double *p)
{
  dp00dt=-(1+gam)*p00+p10-p00+p01;
  dp10dt=(1+gam)*p00-p10-p10+p11;
  dp01dt=-p01+p11+p00-p01;
  dp11dt=p01-p11+p10-p11;
} // ode

/* set ex->r3[1..data->mg]; (ineq. constr., components must be >= 0) */
void R3(Glob *globs, GlobExp *ex, int computeDerivs)
{
  long i,j;
  // parameter constraints
  // parameter: gam
  // variable constraints
  // variable: p00
  // variable: p10
  // variable: p01
  // variable: p11
}
void setNrConstraints(GlobExp *ex, int nP, double *parameters)
{
  ex->me=0; ex->mg=0;
}
int observation (Glob *globs,GlobExp *ex,double t, double *y, double *gy,double *p, double **dgdy,double **dgdp)
{
  int generic=FALSE; // user defined observation function
  long k;
  gy[1]=p00+p01;
  gy[2]=p10+p11;
  
  if(dgdy) {
    dgdy[1][1]=1.0;
    dgdy[1][2]=0.0;
    dgdy[1][3]=1.0;
    dgdy[1][4]=0.0;
    dgdy[2][1]=0.0;
    dgdy[2][2]=1.0;
    dgdy[2][3]=0.0;
    dgdy[2][4]=1.0;
  }

  if(dgdp) {
     dgdp[1][1]=0.0;
     dgdp[2][1]=0.0;
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
  jac[ip00][ip00]=-gam-2.0;
  jac[ip10][ip00]=1.0;
  jac[ip01][ip00]=1.0;
  jac[ip11][ip00]=0.0;
  jac[ip00][ip10]=gam+1.0;
  jac[ip10][ip10]=-2.0;
  jac[ip01][ip10]=0.0;
  jac[ip11][ip10]=1.0;
  jac[ip00][ip01]=1.0;
  jac[ip10][ip01]=0.0;
  jac[ip01][ip01]=-2.0;
  jac[ip11][ip01]=1.0;
  jac[ip00][ip11]=0.0;
  jac[ip10][ip11]=1.0;
  jac[ip01][ip11]=1.0;
  jac[ip11][ip11]=-2.0;

} // automatically generated elements

void inhomo1(double t, double *y, double *inh1, double *p, int &jp)
{
if(jp==1){
   inh1[ip00]=-1.0*p00;
   inh1[ip10]=p00;
   inh1[ip01]=0.0;
   inh1[ip11]=0.0;
}
} // automatically generated elements
void inhomo(double t, double *y, double **inh, double *p)
{
  int jp=1;

  for (int jp=1; jp<=NPARAMS; jp++)
      inhomo1(t,y,inh[jp],p,jp);
}

