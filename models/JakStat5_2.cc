#include <math.h>
#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <string.h>

#include "nr.h"
#include "def.h"

using namespace std;

extern const char *DefModelDescription = 
"Thorstens JakStat model";

extern const int NEQNS   = 12;
extern const int NOBS    = 2;
extern const int NPARAMS = 5;
extern const int NSPLINES  = 1;

//Spline data
double **splineNodes;
double **splineY;
double **splineGam;
unsigned long *nNodes;


const int	ix1	=1;
const int	ix2	=2;
const int	ix3	=3;
const int	ix4	=4;
const int	iq1	=5;
const int	iq2	=6;
const int	iq3	=7;
const int	iq4	=8;
const int	iq5	=9;
const int	iq6	=10;
const int	iq7	=11;
const int	iq8	=12;
#define dx1dt	f[1]
#define dx2dt	f[2]
#define dx3dt	f[3]
#define dx4dt	f[4]
#define dq1dt	f[5]
#define dq2dt	f[6]
#define dq3dt	f[7]
#define dq4dt	f[8]
#define dq5dt	f[9]
#define dq6dt	f[10]
#define dq7dt	f[11]
#define dq8dt	f[12]
#define x1	y[1]
#define x2	y[2]
#define x3	y[3]
#define x4	y[4]
#define q1	y[5]
#define q2	y[6]
#define q3	y[7]
#define q4	y[8]
#define q5	y[9]
#define q6	y[10]
#define q7	y[11]
#define q8	y[12]
#define k1	p[1]
#define k3	p[2]
#define k4	p[3]
#define k5	p[4]
#define k6	p[5]
#define spline1	spline(splineNodes[1],splineY[1],splineGam[1],nNodes[1],t)


extern const double DefParameters[] = {1.4634,0.24,0.3,1,1};
extern const double DefYValues[]    = {1,0,0,0,0,0,0,0,0,0,0,0};
extern const string VariableNames[] = {"x1","x2","x3","x4","q1","q2","q3","q4","q5","q6","q7","q8"};
extern const string ParameterNames[] = {"k1","k3","k4","k5","k6"};

void ode(double t, double *y, double *f, double *p)
{
  dx1dt=-k1*x1*spline1+k3*q8;
  dx2dt=-x2*x2+k1*x1*spline1;
  dx3dt=-k3*x3+x2*x2;
  dx4dt=k3*x3-k3*q8;
  dq1dt=k4*x3-k4*q1;
  dq2dt=k4*q1-k4*q2;
  dq3dt=k4*q2-k4*q3;
  dq4dt=k4*q3-k4*q4;
  dq5dt=k4*q4-k4*q5;
  dq6dt=k4*q5-k4*q6;
  dq7dt=k4*q6-k4*q7;
  dq8dt=k4*q7-k4*q8;
} // ode

/* set ex->r3[1..data->mg]; (ineq. constr., components must be >= 0) */
void R3(Glob *globs, GlobExp *ex, int computeDerivs)
{
  long i,j;
  // parameter constraints
  // parameter: k1
  ex->r3[1]=ex->par[1]-(0);
  if(computeDerivs) {
    for (i=1;i<=globs->npar;i++) {
      if (i == 1)
        ex->dR3dp[1][i]=1.0;
     }
  }
  // parameter: k3
  ex->r3[2]=ex->par[2]-(0);
  if(computeDerivs) {
    for (i=1;i<=globs->npar;i++) {
      if (i == 2)
        ex->dR3dp[2][i]=1.0;
     }
  }
  // parameter: k4
  ex->r3[3]=ex->par[3]-(0.1);
  if(computeDerivs) {
    for (i=1;i<=globs->npar;i++) {
      if (i == 3)
        ex->dR3dp[3][i]=1.0;
     }
  }
  ex->r3[4]=-ex->par[3]+(10);
  if(computeDerivs) {
    for (i=1;i<=globs->npar;i++) {
      if (i == 3)
        ex->dR3dp[4][i]=-1.0;
    }
  }
  // parameter: k5
  ex->r3[5]=ex->par[4]-(0);
  if(computeDerivs) {
    for (i=1;i<=globs->npar;i++) {
      if (i == 4)
        ex->dR3dp[5][i]=1.0;
     }
  }
  // parameter: k6
  ex->r3[6]=ex->par[5]-(0);
  if(computeDerivs) {
    for (i=1;i<=globs->npar;i++) {
      if (i == 5)
        ex->dR3dp[6][i]=1.0;
     }
  }
  // variable constraints
  // variable: x1
  // variable: x2
  // variable: x3
  // variable: x4
  // variable: q1
  // variable: q2
  // variable: q3
  // variable: q4
  // variable: q5
  // variable: q6
  // variable: q7
  // variable: q8
}
void setNrConstraints(GlobExp *ex, int nP, double *parameters)
{
  ex->me=0; ex->mg=6;
}
int observation (Glob *globs,GlobExp *ex,double t, double *y, double *gy,double *p, double **dgdy,double **dgdp)
{
  int generic=FALSE; // user defined observation function
  long k;
  gy[1]=k5*(x2+x3);
  gy[2]=k6*(x1+x2+x3);
  
  if(dgdy) {
    dgdy[1][1]=0.0;
    dgdy[1][2]=k5;
    dgdy[1][3]=k5;
    dgdy[1][4]=0.0;
    dgdy[1][5]=0.0;
    dgdy[1][6]=0.0;
    dgdy[1][7]=0.0;
    dgdy[1][8]=0.0;
    dgdy[1][9]=0.0;
    dgdy[1][10]=0.0;
    dgdy[1][11]=0.0;
    dgdy[1][12]=0.0;
    dgdy[2][1]=k6;
    dgdy[2][2]=k6;
    dgdy[2][3]=k6;
    dgdy[2][4]=0.0;
    dgdy[2][5]=0.0;
    dgdy[2][6]=0.0;
    dgdy[2][7]=0.0;
    dgdy[2][8]=0.0;
    dgdy[2][9]=0.0;
    dgdy[2][10]=0.0;
    dgdy[2][11]=0.0;
    dgdy[2][12]=0.0;
  }

  if(dgdp) {
     dgdp[1][1]=0.0;
     dgdp[1][2]=0.0;
     dgdp[1][3]=0.0;
     dgdp[1][4]=x3+x2;
     dgdp[1][5]=0.0;
     dgdp[2][1]=0.0;
     dgdp[2][2]=0.0;
     dgdp[2][3]=0.0;
     dgdp[2][4]=0.0;
     dgdp[2][5]=x3+x2+x1;
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
  jac[ix1][ix1]=-1.0*spline1*k1;
  jac[ix2][ix1]=0.0;
  jac[ix3][ix1]=0.0;
  jac[ix4][ix1]=0.0;
  jac[iq1][ix1]=0.0;
  jac[iq2][ix1]=0.0;
  jac[iq3][ix1]=0.0;
  jac[iq4][ix1]=0.0;
  jac[iq5][ix1]=0.0;
  jac[iq6][ix1]=0.0;
  jac[iq7][ix1]=0.0;
  jac[iq8][ix1]=k3;
  jac[ix1][ix2]=spline1*k1;
  jac[ix2][ix2]=-2.0*x2;
  jac[ix3][ix2]=0.0;
  jac[ix4][ix2]=0.0;
  jac[iq1][ix2]=0.0;
  jac[iq2][ix2]=0.0;
  jac[iq3][ix2]=0.0;
  jac[iq4][ix2]=0.0;
  jac[iq5][ix2]=0.0;
  jac[iq6][ix2]=0.0;
  jac[iq7][ix2]=0.0;
  jac[iq8][ix2]=0.0;
  jac[ix1][ix3]=0.0;
  jac[ix2][ix3]=2.0*x2;
  jac[ix3][ix3]=-1.0*k3;
  jac[ix4][ix3]=0.0;
  jac[iq1][ix3]=0.0;
  jac[iq2][ix3]=0.0;
  jac[iq3][ix3]=0.0;
  jac[iq4][ix3]=0.0;
  jac[iq5][ix3]=0.0;
  jac[iq6][ix3]=0.0;
  jac[iq7][ix3]=0.0;
  jac[iq8][ix3]=0.0;
  jac[ix1][ix4]=0.0;
  jac[ix2][ix4]=0.0;
  jac[ix3][ix4]=k3;
  jac[ix4][ix4]=0.0;
  jac[iq1][ix4]=0.0;
  jac[iq2][ix4]=0.0;
  jac[iq3][ix4]=0.0;
  jac[iq4][ix4]=0.0;
  jac[iq5][ix4]=0.0;
  jac[iq6][ix4]=0.0;
  jac[iq7][ix4]=0.0;
  jac[iq8][ix4]=-1.0*k3;
  jac[ix1][iq1]=0.0;
  jac[ix2][iq1]=0.0;
  jac[ix3][iq1]=k4;
  jac[ix4][iq1]=0.0;
  jac[iq1][iq1]=-1.0*k4;
  jac[iq2][iq1]=0.0;
  jac[iq3][iq1]=0.0;
  jac[iq4][iq1]=0.0;
  jac[iq5][iq1]=0.0;
  jac[iq6][iq1]=0.0;
  jac[iq7][iq1]=0.0;
  jac[iq8][iq1]=0.0;
  jac[ix1][iq2]=0.0;
  jac[ix2][iq2]=0.0;
  jac[ix3][iq2]=0.0;
  jac[ix4][iq2]=0.0;
  jac[iq1][iq2]=k4;
  jac[iq2][iq2]=-1.0*k4;
  jac[iq3][iq2]=0.0;
  jac[iq4][iq2]=0.0;
  jac[iq5][iq2]=0.0;
  jac[iq6][iq2]=0.0;
  jac[iq7][iq2]=0.0;
  jac[iq8][iq2]=0.0;
  jac[ix1][iq3]=0.0;
  jac[ix2][iq3]=0.0;
  jac[ix3][iq3]=0.0;
  jac[ix4][iq3]=0.0;
  jac[iq1][iq3]=0.0;
  jac[iq2][iq3]=k4;
  jac[iq3][iq3]=-1.0*k4;
  jac[iq4][iq3]=0.0;
  jac[iq5][iq3]=0.0;
  jac[iq6][iq3]=0.0;
  jac[iq7][iq3]=0.0;
  jac[iq8][iq3]=0.0;
  jac[ix1][iq4]=0.0;
  jac[ix2][iq4]=0.0;
  jac[ix3][iq4]=0.0;
  jac[ix4][iq4]=0.0;
  jac[iq1][iq4]=0.0;
  jac[iq2][iq4]=0.0;
  jac[iq3][iq4]=k4;
  jac[iq4][iq4]=-1.0*k4;
  jac[iq5][iq4]=0.0;
  jac[iq6][iq4]=0.0;
  jac[iq7][iq4]=0.0;
  jac[iq8][iq4]=0.0;
  jac[ix1][iq5]=0.0;
  jac[ix2][iq5]=0.0;
  jac[ix3][iq5]=0.0;
  jac[ix4][iq5]=0.0;
  jac[iq1][iq5]=0.0;
  jac[iq2][iq5]=0.0;
  jac[iq3][iq5]=0.0;
  jac[iq4][iq5]=k4;
  jac[iq5][iq5]=-1.0*k4;
  jac[iq6][iq5]=0.0;
  jac[iq7][iq5]=0.0;
  jac[iq8][iq5]=0.0;
  jac[ix1][iq6]=0.0;
  jac[ix2][iq6]=0.0;
  jac[ix3][iq6]=0.0;
  jac[ix4][iq6]=0.0;
  jac[iq1][iq6]=0.0;
  jac[iq2][iq6]=0.0;
  jac[iq3][iq6]=0.0;
  jac[iq4][iq6]=0.0;
  jac[iq5][iq6]=k4;
  jac[iq6][iq6]=-1.0*k4;
  jac[iq7][iq6]=0.0;
  jac[iq8][iq6]=0.0;
  jac[ix1][iq7]=0.0;
  jac[ix2][iq7]=0.0;
  jac[ix3][iq7]=0.0;
  jac[ix4][iq7]=0.0;
  jac[iq1][iq7]=0.0;
  jac[iq2][iq7]=0.0;
  jac[iq3][iq7]=0.0;
  jac[iq4][iq7]=0.0;
  jac[iq5][iq7]=0.0;
  jac[iq6][iq7]=k4;
  jac[iq7][iq7]=-1.0*k4;
  jac[iq8][iq7]=0.0;
  jac[ix1][iq8]=0.0;
  jac[ix2][iq8]=0.0;
  jac[ix3][iq8]=0.0;
  jac[ix4][iq8]=0.0;
  jac[iq1][iq8]=0.0;
  jac[iq2][iq8]=0.0;
  jac[iq3][iq8]=0.0;
  jac[iq4][iq8]=0.0;
  jac[iq5][iq8]=0.0;
  jac[iq6][iq8]=0.0;
  jac[iq7][iq8]=k4;
  jac[iq8][iq8]=-1.0*k4;

} // automatically generated elements

void inhomo1(double t, double *y, double *inh1, double *p, int &jp)
{
if(jp==1){
   inh1[ix1]=-1.0*spline1*x1;
   inh1[ix2]=spline1*x1;
   inh1[ix3]=0.0;
   inh1[ix4]=0.0;
   inh1[iq1]=0.0;
   inh1[iq2]=0.0;
   inh1[iq3]=0.0;
   inh1[iq4]=0.0;
   inh1[iq5]=0.0;
   inh1[iq6]=0.0;
   inh1[iq7]=0.0;
   inh1[iq8]=0.0;
}
if(jp==2){
   inh1[ix1]=q8;
   inh1[ix2]=0.0;
   inh1[ix3]=-1.0*x3;
   inh1[ix4]=x3-q8;
   inh1[iq1]=0.0;
   inh1[iq2]=0.0;
   inh1[iq3]=0.0;
   inh1[iq4]=0.0;
   inh1[iq5]=0.0;
   inh1[iq6]=0.0;
   inh1[iq7]=0.0;
   inh1[iq8]=0.0;
}
if(jp==3){
   inh1[ix1]=0.0;
   inh1[ix2]=0.0;
   inh1[ix3]=0.0;
   inh1[ix4]=0.0;
   inh1[iq1]=x3-q1;
   inh1[iq2]=-q2+q1;
   inh1[iq3]=-q3+q2;
   inh1[iq4]=q3-q4;
   inh1[iq5]=-q5+q4;
   inh1[iq6]=-q6+q5;
   inh1[iq7]=-q7+q6;
   inh1[iq8]=-q8+q7;
}
if(jp==4){
   inh1[ix1]=0.0;
   inh1[ix2]=0.0;
   inh1[ix3]=0.0;
   inh1[ix4]=0.0;
   inh1[iq1]=0.0;
   inh1[iq2]=0.0;
   inh1[iq3]=0.0;
   inh1[iq4]=0.0;
   inh1[iq5]=0.0;
   inh1[iq6]=0.0;
   inh1[iq7]=0.0;
   inh1[iq8]=0.0;
}
if(jp==5){
   inh1[ix1]=0.0;
   inh1[ix2]=0.0;
   inh1[ix3]=0.0;
   inh1[ix4]=0.0;
   inh1[iq1]=0.0;
   inh1[iq2]=0.0;
   inh1[iq3]=0.0;
   inh1[iq4]=0.0;
   inh1[iq5]=0.0;
   inh1[iq6]=0.0;
   inh1[iq7]=0.0;
   inh1[iq8]=0.0;
}
} // automatically generated elements
void inhomo(double t, double *y, double **inh, double *p)
{
  int jp=1;

  for (int jp=1; jp<=NPARAMS; jp++)
      inhomo1(t,y,inh[jp],p,jp);
}

