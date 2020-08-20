#include <math.h>
#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <string.h>

#include "nr.h"
#include "def.h"

using namespace std;

extern const char *DefModelDescription = 
"Ein Modell der Repeats mit CpG;A->C    T->G    p1;A->T    T->A    p2;C->G    G->C    p3;C->A    G->T    p4;A->G    T->C    p5;G->A    C->T    p6";

extern const int NEQNS   = 14;
extern const int NOBS    = 4;
extern const int NPARAMS = 8;
extern const int NSPLINES  = 0;

//Spline data
double **splineNodes;
double **splineY;
double **splineGam;
unsigned long *nNodes;


const int	ipM1	=1;
const int	ipM2	=2;
const int	ipM3	=3;
const int	ipM4	=4;
const int	ipL1	=5;
const int	ipL2	=6;
const int	ipL3	=7;
const int	ipL4	=8;
const int	ipR1	=9;
const int	ipR2	=10;
const int	ipR3	=11;
const int	ipR4	=12;
const int	ipCGL	=13;
const int	ipCGR	=14;
#define dpM1dt	f[1]
#define dpM2dt	f[2]
#define dpM3dt	f[3]
#define dpM4dt	f[4]
#define dpL1dt	f[5]
#define dpL2dt	f[6]
#define dpL3dt	f[7]
#define dpL4dt	f[8]
#define dpR1dt	f[9]
#define dpR2dt	f[10]
#define dpR3dt	f[11]
#define dpR4dt	f[12]
#define dpCGLdt	f[13]
#define dpCGRdt	f[14]
#define pM1	y[1]
#define pM2	y[2]
#define pM3	y[3]
#define pM4	y[4]
#define pL1	y[5]
#define pL2	y[6]
#define pL3	y[7]
#define pL4	y[8]
#define pR1	y[9]
#define pR2	y[10]
#define pR3	y[11]
#define pR4	y[12]
#define pCGL	y[13]
#define pCGR	y[14]
#define q01	p[1]
#define q02	p[2]
#define q03	p[3]
#define q04	p[4]
#define q05	p[5]
#define q06	p[6]
#define qQuer	p[7]
#define qCpG	p[8]


extern const double DefParameters[] = {1.0,1.0,1.0,1.0,1.0,1.0,1.0,40};
extern const double DefYValues[]    = {1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0};
extern const string VariableNames[] = {"pM1","pM2","pM3","pM4","pL1","pL2","pL3","pL4","pR1","pR2","pR3","pR4","pCGL","pCGR"};
extern const string ParameterNames[] = {"q01","q02","q03","q04","q05","q06","qQuer","qCpG"};

void ode(double t, double *y, double *f, double *p)
{
  dpM1dt=-(q01+q05+q02)*pM1+(q04)*pM2+(q06)*pM3+(q02)*pM4+qCpG*pCGL;
  dpM2dt=(q01)*pM1-(q04+q03+q06)*pM2+(q03)*pM3+(q05)*pM4-qCpG*pCGR;
  dpM3dt=(q05)*pM1+(q03)*pM2-(q06+q03+q04)*pM3+(q01)*pM4-qCpG*pCGL;
  dpM4dt=(q02)*pM1+(q06)*pM2+(q04)*pM3-(q02+q05+q01)*pM4+qCpG*pCGR;
  dpL1dt=-(q01+q05+q02)*pL1+(q04)*pL2+(q06)*pL3+(q02)*pL4;
  dpL2dt=(q01)*pL1-(q04+q03+q06)*pL2+(q03)*pL3+(q05)*pL4-qCpG*pCGL;
  dpL3dt=(q05)*pL1+(q03)*pL2-(q06+q03+q04)*pL3+(q01)*pL4;
  dpL4dt=(q02)*pL1+(q06)*pL2+(q04)*pL3-(q02+q05+q01)*pL4+qCpG*pCGL;
  dpR1dt=-(q01+q05+q02)*pR1+(q04)*pR2+(q06)*pR3+(q02)*pR4+qCpG*pCGR;
  dpR2dt=(q01)*pR1-(q04+q03+q06)*pR2+(q03)*pR3+(q05)*pR4;
  dpR3dt=(q05)*pR1+(q03)*pR2-(q06+q03+q04)*pR3+(q01)*pR4-qCpG*pCGR;
  dpR4dt=(q02)*pR1+(q06)*pR2+(q04)*pR3-(q02+q05+q01)*pR4;
  dpCGLdt=-2*(qCpG+4*qQuer+q05+q06-q01-q05)*pCGL+qQuer*(pM3+pL2);
  dpCGRdt=-2*(qCpG+4*qQuer+q05+q06-q01-q05)*pCGR+qQuer*(pM2+pL3);
} // ode

/* set ex->r3[1..data->mg]; (ineq. constr., components must be >= 0) */
void R3(Glob *globs, GlobExp *ex, int computeDerivs)
{
  long i,j;
  // parameter constraints
  // parameter: q01
  ex->r3[1]=ex->par[1]-(0.0);
  if(computeDerivs) {
    for (i=1;i<=globs->npar;i++) {
      if (i == 1)
        ex->dR3dp[1][i]=1.0;
     }
  }
  // parameter: q02
  ex->r3[2]=ex->par[2]-(0.0);
  if(computeDerivs) {
    for (i=1;i<=globs->npar;i++) {
      if (i == 2)
        ex->dR3dp[2][i]=1.0;
     }
  }
  // parameter: q03
  ex->r3[3]=ex->par[3]-(0.0);
  if(computeDerivs) {
    for (i=1;i<=globs->npar;i++) {
      if (i == 3)
        ex->dR3dp[3][i]=1.0;
     }
  }
  // parameter: q04
  ex->r3[4]=ex->par[4]-(0.0);
  if(computeDerivs) {
    for (i=1;i<=globs->npar;i++) {
      if (i == 4)
        ex->dR3dp[4][i]=1.0;
     }
  }
  // parameter: q05
  ex->r3[5]=ex->par[5]-(0.0);
  if(computeDerivs) {
    for (i=1;i<=globs->npar;i++) {
      if (i == 5)
        ex->dR3dp[5][i]=1.0;
     }
  }
  // parameter: q06
  ex->r3[6]=ex->par[6]-(0.0);
  if(computeDerivs) {
    for (i=1;i<=globs->npar;i++) {
      if (i == 6)
        ex->dR3dp[6][i]=1.0;
     }
  }
  // parameter: qQuer
  ex->r3[7]=ex->par[7]-(0.0);
  if(computeDerivs) {
    for (i=1;i<=globs->npar;i++) {
      if (i == 7)
        ex->dR3dp[7][i]=1.0;
     }
  }
  // parameter: qCpG
  ex->r3[8]=ex->par[8]-(0.0);
  if(computeDerivs) {
    for (i=1;i<=globs->npar;i++) {
      if (i == 8)
        ex->dR3dp[8][i]=1.0;
     }
  }
  // variable constraints
  // variable: pM1
  // variable: pM2
  // variable: pM3
  // variable: pM4
  // variable: pL1
  // variable: pL2
  // variable: pL3
  // variable: pL4
  // variable: pR1
  // variable: pR2
  // variable: pR3
  // variable: pR4
  // variable: pCGL
  // variable: pCGR
}
void setNrConstraints(GlobExp *ex, int nP, double *parameters)
{
  ex->me=0; ex->mg=8;
}
int observation (Glob *globs,GlobExp *ex,double t, double *y, double *gy,double *p, double **dgdy,double **dgdp)
{
  int generic=FALSE; // user defined observation function
  long k;
  gy[1]=pM1;
  gy[2]=pM2;
  gy[3]=pM3;
  gy[4]=pM4;
  
  if(dgdy) {
    dgdy[1][1]=1.0;
    dgdy[1][2]=0.0;
    dgdy[1][3]=0.0;
    dgdy[1][4]=0.0;
    dgdy[1][5]=0.0;
    dgdy[1][6]=0.0;
    dgdy[1][7]=0.0;
    dgdy[1][8]=0.0;
    dgdy[1][9]=0.0;
    dgdy[1][10]=0.0;
    dgdy[1][11]=0.0;
    dgdy[1][12]=0.0;
    dgdy[1][13]=0.0;
    dgdy[1][14]=0.0;
    dgdy[2][1]=0.0;
    dgdy[2][2]=1.0;
    dgdy[2][3]=0.0;
    dgdy[2][4]=0.0;
    dgdy[2][5]=0.0;
    dgdy[2][6]=0.0;
    dgdy[2][7]=0.0;
    dgdy[2][8]=0.0;
    dgdy[2][9]=0.0;
    dgdy[2][10]=0.0;
    dgdy[2][11]=0.0;
    dgdy[2][12]=0.0;
    dgdy[2][13]=0.0;
    dgdy[2][14]=0.0;
    dgdy[3][1]=0.0;
    dgdy[3][2]=0.0;
    dgdy[3][3]=1.0;
    dgdy[3][4]=0.0;
    dgdy[3][5]=0.0;
    dgdy[3][6]=0.0;
    dgdy[3][7]=0.0;
    dgdy[3][8]=0.0;
    dgdy[3][9]=0.0;
    dgdy[3][10]=0.0;
    dgdy[3][11]=0.0;
    dgdy[3][12]=0.0;
    dgdy[3][13]=0.0;
    dgdy[3][14]=0.0;
    dgdy[4][1]=0.0;
    dgdy[4][2]=0.0;
    dgdy[4][3]=0.0;
    dgdy[4][4]=1.0;
    dgdy[4][5]=0.0;
    dgdy[4][6]=0.0;
    dgdy[4][7]=0.0;
    dgdy[4][8]=0.0;
    dgdy[4][9]=0.0;
    dgdy[4][10]=0.0;
    dgdy[4][11]=0.0;
    dgdy[4][12]=0.0;
    dgdy[4][13]=0.0;
    dgdy[4][14]=0.0;
  }

  if(dgdp) {
     dgdp[1][1]=0.0;
     dgdp[1][2]=0.0;
     dgdp[1][3]=0.0;
     dgdp[1][4]=0.0;
     dgdp[1][5]=0.0;
     dgdp[1][6]=0.0;
     dgdp[1][7]=0.0;
     dgdp[1][8]=0.0;
     dgdp[2][1]=0.0;
     dgdp[2][2]=0.0;
     dgdp[2][3]=0.0;
     dgdp[2][4]=0.0;
     dgdp[2][5]=0.0;
     dgdp[2][6]=0.0;
     dgdp[2][7]=0.0;
     dgdp[2][8]=0.0;
     dgdp[3][1]=0.0;
     dgdp[3][2]=0.0;
     dgdp[3][3]=0.0;
     dgdp[3][4]=0.0;
     dgdp[3][5]=0.0;
     dgdp[3][6]=0.0;
     dgdp[3][7]=0.0;
     dgdp[3][8]=0.0;
     dgdp[4][1]=0.0;
     dgdp[4][2]=0.0;
     dgdp[4][3]=0.0;
     dgdp[4][4]=0.0;
     dgdp[4][5]=0.0;
     dgdp[4][6]=0.0;
     dgdp[4][7]=0.0;
     dgdp[4][8]=0.0;
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
  jac[ipM1][ipM1]=-q02-q01-q05;
  jac[ipM2][ipM1]=q04;
  jac[ipM3][ipM1]=q06;
  jac[ipM4][ipM1]=q02;
  jac[ipL1][ipM1]=0.0;
  jac[ipL2][ipM1]=0.0;
  jac[ipL3][ipM1]=0.0;
  jac[ipL4][ipM1]=0.0;
  jac[ipR1][ipM1]=0.0;
  jac[ipR2][ipM1]=0.0;
  jac[ipR3][ipM1]=0.0;
  jac[ipR4][ipM1]=0.0;
  jac[ipCGL][ipM1]=qCpG;
  jac[ipCGR][ipM1]=0.0;
  jac[ipM1][ipM2]=q01;
  jac[ipM2][ipM2]=-q04-q03-q06;
  jac[ipM3][ipM2]=q03;
  jac[ipM4][ipM2]=q05;
  jac[ipL1][ipM2]=0.0;
  jac[ipL2][ipM2]=0.0;
  jac[ipL3][ipM2]=0.0;
  jac[ipL4][ipM2]=0.0;
  jac[ipR1][ipM2]=0.0;
  jac[ipR2][ipM2]=0.0;
  jac[ipR3][ipM2]=0.0;
  jac[ipR4][ipM2]=0.0;
  jac[ipCGL][ipM2]=0.0;
  jac[ipCGR][ipM2]=-1.0*qCpG;
  jac[ipM1][ipM3]=q05;
  jac[ipM2][ipM3]=q03;
  jac[ipM3][ipM3]=-q04-q03-q06;
  jac[ipM4][ipM3]=q01;
  jac[ipL1][ipM3]=0.0;
  jac[ipL2][ipM3]=0.0;
  jac[ipL3][ipM3]=0.0;
  jac[ipL4][ipM3]=0.0;
  jac[ipR1][ipM3]=0.0;
  jac[ipR2][ipM3]=0.0;
  jac[ipR3][ipM3]=0.0;
  jac[ipR4][ipM3]=0.0;
  jac[ipCGL][ipM3]=-1.0*qCpG;
  jac[ipCGR][ipM3]=0.0;
  jac[ipM1][ipM4]=q02;
  jac[ipM2][ipM4]=q06;
  jac[ipM3][ipM4]=q04;
  jac[ipM4][ipM4]=-q02-q01-q05;
  jac[ipL1][ipM4]=0.0;
  jac[ipL2][ipM4]=0.0;
  jac[ipL3][ipM4]=0.0;
  jac[ipL4][ipM4]=0.0;
  jac[ipR1][ipM4]=0.0;
  jac[ipR2][ipM4]=0.0;
  jac[ipR3][ipM4]=0.0;
  jac[ipR4][ipM4]=0.0;
  jac[ipCGL][ipM4]=0.0;
  jac[ipCGR][ipM4]=qCpG;
  jac[ipM1][ipL1]=0.0;
  jac[ipM2][ipL1]=0.0;
  jac[ipM3][ipL1]=0.0;
  jac[ipM4][ipL1]=0.0;
  jac[ipL1][ipL1]=-q02-q01-q05;
  jac[ipL2][ipL1]=q04;
  jac[ipL3][ipL1]=q06;
  jac[ipL4][ipL1]=q02;
  jac[ipR1][ipL1]=0.0;
  jac[ipR2][ipL1]=0.0;
  jac[ipR3][ipL1]=0.0;
  jac[ipR4][ipL1]=0.0;
  jac[ipCGL][ipL1]=0.0;
  jac[ipCGR][ipL1]=0.0;
  jac[ipM1][ipL2]=0.0;
  jac[ipM2][ipL2]=0.0;
  jac[ipM3][ipL2]=0.0;
  jac[ipM4][ipL2]=0.0;
  jac[ipL1][ipL2]=q01;
  jac[ipL2][ipL2]=-q04-q03-q06;
  jac[ipL3][ipL2]=q03;
  jac[ipL4][ipL2]=q05;
  jac[ipR1][ipL2]=0.0;
  jac[ipR2][ipL2]=0.0;
  jac[ipR3][ipL2]=0.0;
  jac[ipR4][ipL2]=0.0;
  jac[ipCGL][ipL2]=-1.0*qCpG;
  jac[ipCGR][ipL2]=0.0;
  jac[ipM1][ipL3]=0.0;
  jac[ipM2][ipL3]=0.0;
  jac[ipM3][ipL3]=0.0;
  jac[ipM4][ipL3]=0.0;
  jac[ipL1][ipL3]=q05;
  jac[ipL2][ipL3]=q03;
  jac[ipL3][ipL3]=-q04-q03-q06;
  jac[ipL4][ipL3]=q01;
  jac[ipR1][ipL3]=0.0;
  jac[ipR2][ipL3]=0.0;
  jac[ipR3][ipL3]=0.0;
  jac[ipR4][ipL3]=0.0;
  jac[ipCGL][ipL3]=0.0;
  jac[ipCGR][ipL3]=0.0;
  jac[ipM1][ipL4]=0.0;
  jac[ipM2][ipL4]=0.0;
  jac[ipM3][ipL4]=0.0;
  jac[ipM4][ipL4]=0.0;
  jac[ipL1][ipL4]=q02;
  jac[ipL2][ipL4]=q06;
  jac[ipL3][ipL4]=q04;
  jac[ipL4][ipL4]=-q02-q01-q05;
  jac[ipR1][ipL4]=0.0;
  jac[ipR2][ipL4]=0.0;
  jac[ipR3][ipL4]=0.0;
  jac[ipR4][ipL4]=0.0;
  jac[ipCGL][ipL4]=qCpG;
  jac[ipCGR][ipL4]=0.0;
  jac[ipM1][ipR1]=0.0;
  jac[ipM2][ipR1]=0.0;
  jac[ipM3][ipR1]=0.0;
  jac[ipM4][ipR1]=0.0;
  jac[ipL1][ipR1]=0.0;
  jac[ipL2][ipR1]=0.0;
  jac[ipL3][ipR1]=0.0;
  jac[ipL4][ipR1]=0.0;
  jac[ipR1][ipR1]=-q02-q01-q05;
  jac[ipR2][ipR1]=q04;
  jac[ipR3][ipR1]=q06;
  jac[ipR4][ipR1]=q02;
  jac[ipCGL][ipR1]=0.0;
  jac[ipCGR][ipR1]=qCpG;
  jac[ipM1][ipR2]=0.0;
  jac[ipM2][ipR2]=0.0;
  jac[ipM3][ipR2]=0.0;
  jac[ipM4][ipR2]=0.0;
  jac[ipL1][ipR2]=0.0;
  jac[ipL2][ipR2]=0.0;
  jac[ipL3][ipR2]=0.0;
  jac[ipL4][ipR2]=0.0;
  jac[ipR1][ipR2]=q01;
  jac[ipR2][ipR2]=-q04-q03-q06;
  jac[ipR3][ipR2]=q03;
  jac[ipR4][ipR2]=q05;
  jac[ipCGL][ipR2]=0.0;
  jac[ipCGR][ipR2]=0.0;
  jac[ipM1][ipR3]=0.0;
  jac[ipM2][ipR3]=0.0;
  jac[ipM3][ipR3]=0.0;
  jac[ipM4][ipR3]=0.0;
  jac[ipL1][ipR3]=0.0;
  jac[ipL2][ipR3]=0.0;
  jac[ipL3][ipR3]=0.0;
  jac[ipL4][ipR3]=0.0;
  jac[ipR1][ipR3]=q05;
  jac[ipR2][ipR3]=q03;
  jac[ipR3][ipR3]=-q04-q03-q06;
  jac[ipR4][ipR3]=q01;
  jac[ipCGL][ipR3]=0.0;
  jac[ipCGR][ipR3]=-1.0*qCpG;
  jac[ipM1][ipR4]=0.0;
  jac[ipM2][ipR4]=0.0;
  jac[ipM3][ipR4]=0.0;
  jac[ipM4][ipR4]=0.0;
  jac[ipL1][ipR4]=0.0;
  jac[ipL2][ipR4]=0.0;
  jac[ipL3][ipR4]=0.0;
  jac[ipL4][ipR4]=0.0;
  jac[ipR1][ipR4]=q02;
  jac[ipR2][ipR4]=q06;
  jac[ipR3][ipR4]=q04;
  jac[ipR4][ipR4]=-q02-q01-q05;
  jac[ipCGL][ipR4]=0.0;
  jac[ipCGR][ipR4]=0.0;
  jac[ipM1][ipCGL]=0.0;
  jac[ipM2][ipCGL]=0.0;
  jac[ipM3][ipCGL]=qQuer;
  jac[ipM4][ipCGL]=0.0;
  jac[ipL1][ipCGL]=0.0;
  jac[ipL2][ipCGL]=qQuer;
  jac[ipL3][ipCGL]=0.0;
  jac[ipL4][ipCGL]=0.0;
  jac[ipR1][ipCGL]=0.0;
  jac[ipR2][ipCGL]=0.0;
  jac[ipR3][ipCGL]=0.0;
  jac[ipR4][ipCGL]=0.0;
  jac[ipCGL][ipCGL]=-2.0*qCpG-8.0*qQuer-2.0*q06+2.0*q01;
  jac[ipCGR][ipCGL]=0.0;
  jac[ipM1][ipCGR]=0.0;
  jac[ipM2][ipCGR]=qQuer;
  jac[ipM3][ipCGR]=0.0;
  jac[ipM4][ipCGR]=0.0;
  jac[ipL1][ipCGR]=0.0;
  jac[ipL2][ipCGR]=0.0;
  jac[ipL3][ipCGR]=qQuer;
  jac[ipL4][ipCGR]=0.0;
  jac[ipR1][ipCGR]=0.0;
  jac[ipR2][ipCGR]=0.0;
  jac[ipR3][ipCGR]=0.0;
  jac[ipR4][ipCGR]=0.0;
  jac[ipCGL][ipCGR]=0.0;
  jac[ipCGR][ipCGR]=-2.0*qCpG-8.0*qQuer-2.0*q06+2.0*q01;

} // automatically generated elements

void inhomo1(double t, double *y, double *inh1, double *p, int &jp)
{
if(jp==1){
   inh1[ipM1]=-1.0*pM1;
   inh1[ipM2]=pM1;
   inh1[ipM3]=pM4;
   inh1[ipM4]=-1.0*pM4;
   inh1[ipL1]=-1.0*pL1;
   inh1[ipL2]=pL1;
   inh1[ipL3]=pL4;
   inh1[ipL4]=-1.0*pL4;
   inh1[ipR1]=-1.0*pR1;
   inh1[ipR2]=pR1;
   inh1[ipR3]=pR4;
   inh1[ipR4]=-1.0*pR4;
   inh1[ipCGL]=2.0*pCGL;
   inh1[ipCGR]=2.0*pCGR;
}
if(jp==2){
   inh1[ipM1]=-pM1+pM4;
   inh1[ipM2]=0.0;
   inh1[ipM3]=0.0;
   inh1[ipM4]=pM1-pM4;
   inh1[ipL1]=-pL1+pL4;
   inh1[ipL2]=0.0;
   inh1[ipL3]=0.0;
   inh1[ipL4]=pL1-pL4;
   inh1[ipR1]=pR4-pR1;
   inh1[ipR2]=0.0;
   inh1[ipR3]=0.0;
   inh1[ipR4]=-pR4+pR1;
   inh1[ipCGL]=0.0;
   inh1[ipCGR]=0.0;
}
if(jp==3){
   inh1[ipM1]=0.0;
   inh1[ipM2]=pM3-pM2;
   inh1[ipM3]=-pM3+pM2;
   inh1[ipM4]=0.0;
   inh1[ipL1]=0.0;
   inh1[ipL2]=pL3-pL2;
   inh1[ipL3]=-pL3+pL2;
   inh1[ipL4]=0.0;
   inh1[ipR1]=0.0;
   inh1[ipR2]=pR3-pR2;
   inh1[ipR3]=-pR3+pR2;
   inh1[ipR4]=0.0;
   inh1[ipCGL]=0.0;
   inh1[ipCGR]=0.0;
}
if(jp==4){
   inh1[ipM1]=pM2;
   inh1[ipM2]=-1.0*pM2;
   inh1[ipM3]=-1.0*pM3;
   inh1[ipM4]=pM3;
   inh1[ipL1]=pL2;
   inh1[ipL2]=-1.0*pL2;
   inh1[ipL3]=-1.0*pL3;
   inh1[ipL4]=pL3;
   inh1[ipR1]=pR2;
   inh1[ipR2]=-1.0*pR2;
   inh1[ipR3]=-1.0*pR3;
   inh1[ipR4]=pR3;
   inh1[ipCGL]=0.0;
   inh1[ipCGR]=0.0;
}
if(jp==5){
   inh1[ipM1]=-1.0*pM1;
   inh1[ipM2]=pM4;
   inh1[ipM3]=pM1;
   inh1[ipM4]=-1.0*pM4;
   inh1[ipL1]=-1.0*pL1;
   inh1[ipL2]=pL4;
   inh1[ipL3]=pL1;
   inh1[ipL4]=-1.0*pL4;
   inh1[ipR1]=-1.0*pR1;
   inh1[ipR2]=pR4;
   inh1[ipR3]=pR1;
   inh1[ipR4]=-1.0*pR4;
   inh1[ipCGL]=0.0;
   inh1[ipCGR]=0.0;
}
if(jp==6){
   inh1[ipM1]=pM3;
   inh1[ipM2]=-1.0*pM2;
   inh1[ipM3]=-1.0*pM3;
   inh1[ipM4]=pM2;
   inh1[ipL1]=pL3;
   inh1[ipL2]=-1.0*pL2;
   inh1[ipL3]=-1.0*pL3;
   inh1[ipL4]=pL2;
   inh1[ipR1]=pR3;
   inh1[ipR2]=-1.0*pR2;
   inh1[ipR3]=-1.0*pR3;
   inh1[ipR4]=pR2;
   inh1[ipCGL]=-2.0*pCGL;
   inh1[ipCGR]=-2.0*pCGR;
}
if(jp==7){
   inh1[ipM1]=0.0;
   inh1[ipM2]=0.0;
   inh1[ipM3]=0.0;
   inh1[ipM4]=0.0;
   inh1[ipL1]=0.0;
   inh1[ipL2]=0.0;
   inh1[ipL3]=0.0;
   inh1[ipL4]=0.0;
   inh1[ipR1]=0.0;
   inh1[ipR2]=0.0;
   inh1[ipR3]=0.0;
   inh1[ipR4]=0.0;
   inh1[ipCGL]=pM3-8.0*pCGL+pL2;
   inh1[ipCGR]=pM2+pL3-8.0*pCGR;
}
if(jp==8){
   inh1[ipM1]=pCGL;
   inh1[ipM2]=-1.0*pCGR;
   inh1[ipM3]=-1.0*pCGL;
   inh1[ipM4]=pCGR;
   inh1[ipL1]=0.0;
   inh1[ipL2]=-1.0*pCGL;
   inh1[ipL3]=0.0;
   inh1[ipL4]=pCGL;
   inh1[ipR1]=pCGR;
   inh1[ipR2]=0.0;
   inh1[ipR3]=-1.0*pCGR;
   inh1[ipR4]=0.0;
   inh1[ipCGL]=-2.0*pCGL;
   inh1[ipCGR]=-2.0*pCGR;
}
} // automatically generated elements
void inhomo(double t, double *y, double **inh, double *p)
{
  int jp=1;

  for (int jp=1; jp<=NPARAMS; jp++)
      inhomo1(t,y,inh[jp],p,jp);
}

