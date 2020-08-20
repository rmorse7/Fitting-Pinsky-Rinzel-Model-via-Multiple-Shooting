#include <math.h>
#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <string.h>

#include "nr.h"
#include "def.h"

using namespace std;

extern const char *DefModelDescription = 
"Lotka-Volterra-System";

extern const int NEQNS   = 2;
extern const int NOBS    = 2;
extern const int NPARAMS = 4;
extern const int NSPLINES  = 0;

//Spline data
double **splineNodes;
double **splineY;
double **splineGam;
unsigned long *nNodes;


const int	iy1	=1;
const int	iy2	=2;
#define dy1dt	f[1]
#define dy2dt	f[2]
#define y1	y[1]
#define y2	y[2]
#define p1	p[1]
#define p2	p[2]
#define p3	p[3]
#define p4	p[4]


extern const double DefParameters[] = {1.0,1.0,1.0,1.0};
extern const double DefYValues[]    = {0.4,1.0};
extern const string VariableNames[] = {"y1","y2"};
extern const string ParameterNames[] = {"p1","p2","p3","p4"};

void ode(double t, double *y, double *f, double *p)
{
  dy1dt=p1*y1-p2*y1*y2;
  dy2dt=-p3*y2+p4*y1*y2;
} // ode

/* set ex->r3[1..data->mg]; (ineq. constr., components must be >= 0) */
void R3(Glob *globs, GlobExp *ex, int computeDerivs)
{
  long i,j;
  // parameter constraints
  // parameter: p1
  ex->r3[1]=ex->par[1]-(0.0);
  if(computeDerivs) {
    for (i=1;i<=globs->npar;i++) {
      if (i == 1)
        ex->dR3dp[1][i]=1.0;
     }
  }
  ex->r3[2]=-ex->par[1]+(10.0);
  if(computeDerivs) {
    for (i=1;i<=globs->npar;i++) {
      if (i == 1)
        ex->dR3dp[2][i]=-1.0;
    }
  }
  // parameter: p2
  ex->r3[3]=ex->par[2]-(0.0);
  if(computeDerivs) {
    for (i=1;i<=globs->npar;i++) {
      if (i == 2)
        ex->dR3dp[3][i]=1.0;
     }
  }
  ex->r3[4]=-ex->par[2]+(10.0);
  if(computeDerivs) {
    for (i=1;i<=globs->npar;i++) {
      if (i == 2)
        ex->dR3dp[4][i]=-1.0;
    }
  }
  // parameter: p3
  ex->r3[5]=ex->par[3]-(0.0);
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
  // parameter: p4
  ex->r3[7]=ex->par[4]-(0.0);
  if(computeDerivs) {
    for (i=1;i<=globs->npar;i++) {
      if (i == 4)
        ex->dR3dp[7][i]=1.0;
     }
  }
  ex->r3[8]=-ex->par[4]+(10.0);
  if(computeDerivs) {
    for (i=1;i<=globs->npar;i++) {
      if (i == 4)
        ex->dR3dp[8][i]=-1.0;
    }
  }
  // variable constraints
  // variable: y1
  // variable: y2
}
void setNrConstraints(GlobExp *ex, int nP, double *parameters)
{
  ex->me=0; ex->mg=8;
}
int observation (Glob *globs,GlobExp *ex,double t, double *y, double *gy,double *p, double **dgdy,double **dgdp)
{
  int generic=TRUE; // generic obsevation function
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
  jac[iy1][iy1]=p1-y2*p2;
  jac[iy2][iy1]=-1.0*y1*p2;
  jac[iy1][iy2]=y2*p4;
  jac[iy2][iy2]=-p3+p4*y1;

} // automatically generated elements

void inhomo1(double t, double *y, double *inh1, double *p, int &jp)
{
if(jp==1){
   inh1[iy1]=y1;
   inh1[iy2]=0.0;
}
if(jp==2){
   inh1[iy1]=-1.0*y2*y1;
   inh1[iy2]=0.0;
}
if(jp==3){
   inh1[iy1]=0.0;
   inh1[iy2]=-1.0*y2;
}
if(jp==4){
   inh1[iy1]=0.0;
   inh1[iy2]=y2*y1;
}
} // automatically generated elements
void inhomo(double t, double *y, double **inh, double *p)
{
  int jp=1;

  for (int jp=1; jp<=NPARAMS; jp++)
      inhomo1(t,y,inh[jp],p,jp);
}

