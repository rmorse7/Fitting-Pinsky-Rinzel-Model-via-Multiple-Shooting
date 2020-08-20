#include <math.h>
#include <iostream>
#include "../nr.h"
#include "../def.h"

using namespace std;

#define TINY 1.0e-30

/* extension by w.h. 1/98: f means full
if ode_fxp!=NULL, ALL intermediate x-values will be stored there.
Then we can calculate derivatives w.resp. to parameters or initial values
by integrating the curve with the same steps and a light routine like rk4.
We need no intermediate y-values for this purpose => We don't want to set
ode_dxsav=0 and use ode_xp. 
If desired, allocate ode_fxp[1..ode_fkmax] in the calling program. */


void integrateRK(Glob *globs,double ystart[], double p[], long nvar, double x1, double x2,double eps,
		 double h1, double hmin, long *nok, long *nbad,
		 void (*derivs) (double, double[], double[], double[]))
{
  long ode_fkmax = 0, ode_fkount = 0, ode_kmax = 0, ode_kount = 0;
  double ode_dxsav = 0, *ode_fxp = 0, *ode_xp = 0, **ode_yp = 0;
  long rkqs_ign = 0;
  long nstp,i;
  double xsav,x,hnext,hdid,h;
  double *yscal,*y,*dydx;
  long odeint_MAXSTP=globs->maxstp;

  yscal=dvector(1,nvar);
  y=dvector(1,nvar);
  dydx=dvector(1,nvar);
  x=x1;
  h=SIGN(h1,x2-x1);
  *nok = (*nbad) = ode_kount = 0;
  for (i=1;i<=nvar;i++) y[i]=ystart[i];
  if (ode_dxsav<0) nrerror("odeint:dxsav<0");
  if (ode_kmax>0 && ode_kount<ode_kmax-1) {
     ode_xp[++ode_kount]=x;
     for (i=1;i<=nvar;i++) ode_yp[i][ode_kount]=y[i];
     xsav=x;
  } /* save first point, w.h. 19.03.99 */
  for (nstp=1;nstp<=odeint_MAXSTP;nstp++) {
    (*derivs)(x,y,dydx,p);
    for (i=1;i<=nvar;i++)
      yscal[i]=fabs(y[i])+fabs(dydx[i]*h)+TINY;
    if (ode_fkmax) {
      if (ode_fkount>=ode_fkmax) nrerror("ode_fkount exceeds ode_fkmax");
      ode_fxp[++ode_fkount]=x;
    }
    if (ode_kmax>0 && ode_kount<ode_kmax-1 && fabs(x-xsav)>ode_dxsav) {
      ode_xp[++ode_kount]=x;
      for (i=1;i<=nvar;i++) ode_yp[i][ode_kount]=y[i];
      xsav=x;
    }
    if ((x+h-x2)*(x+h-x1) > 0.0) h=x2-x;
    rkqs(y,dydx,p,nvar,&x,h,eps,yscal,&hdid,&hnext,globs->rkqs_ign,derivs);
    if (hdid == h) ++(*nok); 
    else ++(*nbad);
    if ((x-x2)*(x2-x1) >= 0.0) {
      for (i=1;i<=nvar;i++) ystart[i]=y[i];
      if (ode_kmax>0 && ode_kount<ode_kmax-1) {
	ode_xp[++ode_kount]=x;
	for (i=1;i<=nvar;i++) ode_yp[i][ode_kount]=y[i];
      }
      free_dvector(dydx,1,nvar);
      free_dvector(y,1,nvar);
      free_dvector(yscal,1,nvar);
      return;
    }
    if (fabs(hnext) <= hmin)
      {
	cerr << "Step size too small in odeint\n";
	throw 1;
      }
    h=hnext;
  }
  cerr << "Too many steps in routine odeint; try -maxstp <n> greater n\n";
  throw 1;
}
#undef TINY
