// file modified my TM Nov 20, 1995
// modified by w.h. 5/2000: added rkqs_ign: 
// number of state vector components that are not relevant in error control
// if >0, only the leading n-rkqs_ign components are compared
// rkqs_ign=0 is the std behaviour

extern int rkqs_ign;

#include <math.h>
#define NRANSI
#include "../nr.h"
#define SAFETY 0.9
#define PGROW -0.2
#define PSHRNK -0.25
#define ERRCON 1.89e-4

typedef void (*derivType)(double, double [], double [], double[]);

void rkqs(double y[], double dydx[], double p[], long n, double *x,
	  double htry, double eps, double yscal[], double *hdid, double *hnext,long rkqs_ign,
	  derivType derivs) 
{
  void rkck(double y[], double dydx[], double p[], long n, double x,
	    double h, double yout[], double yerr[], derivType derivs);
  int i;
  double errmax,h,htemp,xnew,*yerr,*ytemp;
  
  yerr=dvector(1,n);
  ytemp=dvector(1,n);
  h=htry;
  for (;;) {
    rkck(y,dydx,p,n,*x,h,ytemp,yerr,derivs);
    errmax=0.0;
    for (i=1;i<=n-rkqs_ign;i++) errmax=FMAX(errmax,fabs(yerr[i]/yscal[i]));
    errmax /= eps;
    if (errmax > 1.0) {
      htemp=SAFETY*h*pow(errmax,PSHRNK);
      h=(h >= 0.0 ? FMAX(htemp,0.1*h) : FMIN(htemp,0.1*h));
      xnew=(*x)+h;
      if (xnew == *x) nrerror("stepsize underflow in rkqs");
      continue;
    } else {
      if (errmax > ERRCON) *hnext=SAFETY*h*pow(errmax,PGROW);
      else *hnext=5.0*h;
      *x += (*hdid=h);
      for (i=1;i<=n;i++) y[i]=ytemp[i];
      break;
    }
  }
  free_dvector(ytemp,1,n);
  free_dvector(yerr,1,n);
}
#undef SAFETY
#undef PGROW
#undef PSHRNK
#undef ERRCON
#undef NRANSI
/* (C) Copr. 1986-92 Numerical Recipes Software !"v1`2+!-0. */

