#ifndef __MODELPUB_H__
#define __MODELPUB_H__


/* global variables describing the model */

extern int NEQNS;       // number of ODEs
extern int NOBS;        // number of observed variables
extern int NPARAMS;     // number of parameters
extern int NCONDS;      // number of experimental conditions
extern int NSPLINES;    // number of user defined splines

/* default values */

extern char *DefModelDescription;
extern double DefParameters[];
extern double DefConditions[];
extern double DefYValues[];
extern string ParameterNames[];
extern string VariableNames[];


extern long ode_fkmax,ode_fkount,ode_kmax,ode_kount,rkqs_ign;
extern double ode_dxsav,*ode_fxp,*ode_xp,**ode_yp;

/* model-specific functions */

void ode (double t, double *y, double *f, double *p);
void jacobi (double t, double *y, double **jac, double *p);
void inhomo (double t, double *y, double **inh, double *p);
void inhomo1 (double t, double *y, double *inh1, double *p, int &jp);
int observation (Glob *globs,GlobExp *ex,double t, double *y,double *gy, double *p, double **dgdy,double **dgdp);
void unobsinit (double *y, double *p);

/* constraints */

/* set data->me and data->mg according to number of constraints */
extern void setNrConstraints(GlobExp *data, int nP, double *parameters);

/* set data->r2[1..data->me]; (equal. constr., components must be == 0) */
extern void R2(Glob *globs,GlobExp *ex,int computeDerivs);

/* set data->r3[1..data->mg]; (ineq. constr., components must be >= 0) */
extern void R3(Glob *globs,GlobExp *ex,int computeDerivs);

void initInt(Glob *globs,GlobExp *ex);
#endif

