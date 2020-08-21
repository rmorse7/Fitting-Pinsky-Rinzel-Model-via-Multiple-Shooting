#include <math.h>
#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <string.h>

#include "nr.h"
#include "def.h"

using namespace std;

extern const char *DefModelDescription = 
"Einfaches modell";

extern const int NEQNS   = 1;
extern const int NOBS    = 1;
extern const int NPARAMS = 1;
extern const int NSPLINES  = 0;

//Spline data
double **splineNodes;
double **splineY;
double **splineGam;
unsigned long *nNodes;


const int	iX	=1;
#define dXdt	f[1]
#define X	y[1]
#define phi	p[1]


extern const double DefParameters[] = {1};
extern const double DefYValues[]    = {1};
extern const string VariableNames[] = {"X"};
extern const string ParameterNames[] = {"phi"};

void ode(double t, double *y, double *f, double *p)
{
  dXdt=-phi*X;
} // ode

/* set ex->r3[1..data->mg]; (ineq. constr., components must be >= 0) */
void R3(Glob *globs, GlobExp *ex, int computeDerivs)
{
  long i,j;
  // parameter constraints
  // parameter: phi
  ex->r3[1]=ex->par[1]-(0);
  if(computeDerivs) {
    for (i=1;i<=globs->npar;i++) {
      if (i == 1)
        ex->dR3dp[1][i]=1.0;
     }
  }
  ex->r3[2]=-ex->par[1]+(10);
  if(computeDerivs) {
    for (i=1;i<=globs->npar;i++) {
      if (i == 1)
        ex->dR3dp[2][i]=-1.0;
    }
  }
  // variable constraints
  // variable: X
  ex->r3[3]=ex->yTry[1][1]-(0);
  if(computeDerivs) {
    for (i=1;i<=ex->nvar;i++) {
      if (i == 1)
        ex->dR3ds[3][1][i]=1.0;
     }
  }
  ex->r3[4]=-ex->yTry[1][1]+(10);
  if(computeDerivs) {
    for (i=1;i<=ex->nvar;i++) {
      if (i == 1)
        ex->dR3ds[4][1][i]=-1.0;
    }
  }
}
void setNrConstraints(GlobExp *ex, int nP, double *parameters)
{
  ex->me=0; ex->mg=4;
}
int observation (Glob *globs,GlobExp *ex,double t, double *y, double *gy,double *p, double **dgdy,double **dgdp)
{
  int generic=TRUE; // generic observation function
  long k;

  for(k=1;k<=ex->nobs;k++)
    gy[k]=y[k];

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
  jac[iX][iX]=-phi;

} // automatically generated elements

void inhomo1(double t, double *y, double *inh1, double *p, int &jp)
{
if(jp==1){
   inh1[iX]=-X;
}
} // automatically generated elements
void inhomo(double t, double *y, double **inh, double *p)
{
  int jp=1;

  for (int jp=1; jp<=NPARAMS; jp++)
      inhomo1(t,y,inh[jp],p,jp);
}

