#include <math.h>
#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <string.h>

#include "nr.h"
#include "def.h"

using namespace std;

extern const char *DefModelDescription = 
"Ein Modell der Repeats;A->C    T->G    p1;A->T    T->A    p2;C->G    G->C    p3;C->A    G->T    p4;A->G    T->C    p5;G->A    C->T    p6";

extern const int NEQNS   = 4;
extern const int NOBS    = 4;
extern const int NPARAMS = 6;
extern const int NSPLINES  = 0;

//Spline data
double **splineNodes;
double **splineY;
double **splineGam;
unsigned long *nNodes;


const int	iy1	=1;
const int	iy2	=2;
const int	iy3	=3;
const int	iy4	=4;
#define dy1dt	f[1]
#define dy2dt	f[2]
#define dy3dt	f[3]
#define dy4dt	f[4]
#define y1	y[1]
#define y2	y[2]
#define y3	y[3]
#define y4	y[4]
#define p01	p[1]
#define p02	p[2]
#define p03	p[3]
#define p04	p[4]
#define p05	p[5]
#define p06	p[6]


extern const double DefParameters[] = {1.0,1.0,1.0,1.0,1.0,1.0};
extern const double DefYValues[]    = {1.0,0.0,0.0,0.0};
extern const string VariableNames[] = {"y1","y2","y3","y4"};
extern const string ParameterNames[] = {"p01","p02","p03","p04","p05","p06"};

void ode(double t, double *y, double *f, double *p)
{
  dy1dt=-(p01+p05+p02)*y1+(p04)*y2+(p06)*y3+(p02)*y4;
  dy2dt=(p01)*y1-(p04+p03+p06)*y2+(p03)*y3+(p05)*y4;
  dy3dt=(p05)*y1+(p03)*y2-(p06+p03+p04)*y3+(p01)*y4;
  dy4dt=(p02)*y1+(p06)*y2+(p04)*y3-(p02+p05+p01)*y4;
} // ode

/* set ex->r3[1..data->mg]; (ineq. constr., components must be >= 0) */
void R3(Glob *globs, GlobExp *ex, int computeDerivs)
{
  long i,j;
  // parameter constraints
  // parameter: p01
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
  // parameter: p02
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
  // parameter: p03
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
  // parameter: p04
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
  // parameter: p05
  ex->r3[9]=ex->par[5]-(0.0);
  if(computeDerivs) {
    for (i=1;i<=globs->npar;i++) {
      if (i == 5)
        ex->dR3dp[9][i]=1.0;
     }
  }
  ex->r3[10]=-ex->par[5]+(10.0);
  if(computeDerivs) {
    for (i=1;i<=globs->npar;i++) {
      if (i == 5)
        ex->dR3dp[10][i]=-1.0;
    }
  }
  // parameter: p06
  ex->r3[11]=ex->par[6]-(0.0);
  if(computeDerivs) {
    for (i=1;i<=globs->npar;i++) {
      if (i == 6)
        ex->dR3dp[11][i]=1.0;
     }
  }
  ex->r3[12]=-ex->par[6]+(10.0);
  if(computeDerivs) {
    for (i=1;i<=globs->npar;i++) {
      if (i == 6)
        ex->dR3dp[12][i]=-1.0;
    }
  }
  // variable constraints
  // variable: y1
  // variable: y2
  // variable: y3
  // variable: y4
}
void setNrConstraints(GlobExp *ex, int nP, double *parameters)
{
  ex->me=0; ex->mg=12;
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
  jac[iy1][iy1]=-p02-p01-p05;
  jac[iy2][iy1]=p04;
  jac[iy3][iy1]=p06;
  jac[iy4][iy1]=p02;
  jac[iy1][iy2]=p01;
  jac[iy2][iy2]=-p03-p06-p04;
  jac[iy3][iy2]=p03;
  jac[iy4][iy2]=p05;
  jac[iy1][iy3]=p05;
  jac[iy2][iy3]=p03;
  jac[iy3][iy3]=-p03-p06-p04;
  jac[iy4][iy3]=p01;
  jac[iy1][iy4]=p02;
  jac[iy2][iy4]=p06;
  jac[iy3][iy4]=p04;
  jac[iy4][iy4]=-p02-p01-p05;

} // automatically generated elements

void inhomo1(double t, double *y, double *inh1, double *p, int &jp)
{
if(jp==1){
   inh1[iy1]=-1.0*y1;
   inh1[iy2]=y1;
   inh1[iy3]=y4;
   inh1[iy4]=-1.0*y4;
}
if(jp==2){
   inh1[iy1]=-y1+y4;
   inh1[iy2]=0.0;
   inh1[iy3]=0.0;
   inh1[iy4]=y1-y4;
}
if(jp==3){
   inh1[iy1]=0.0;
   inh1[iy2]=y3-y2;
   inh1[iy3]=-y3+y2;
   inh1[iy4]=0.0;
}
if(jp==4){
   inh1[iy1]=y2;
   inh1[iy2]=-1.0*y2;
   inh1[iy3]=-1.0*y3;
   inh1[iy4]=y3;
}
if(jp==5){
   inh1[iy1]=-1.0*y1;
   inh1[iy2]=y4;
   inh1[iy3]=y1;
   inh1[iy4]=-1.0*y4;
}
if(jp==6){
   inh1[iy1]=y3;
   inh1[iy2]=-1.0*y2;
   inh1[iy3]=-1.0*y3;
   inh1[iy4]=y2;
}
} // automatically generated elements
void inhomo(double t, double *y, double **inh, double *p)
{
  int jp=1;

  for (int jp=1; jp<=NPARAMS; jp++)
      inhomo1(t,y,inh[jp],p,jp);
}

