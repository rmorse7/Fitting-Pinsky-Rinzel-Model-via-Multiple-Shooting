#include <math.h>
#include <iostream>
#include "../nr.h"
#include "../def.h"
#include "../model.h"
#include <string.h>
#include "../libCVODES/sundialstypes.h"   
#include "../libCVODES/cvodes.h"          
#include "../libCVODES/cvdense.h"         
#include "../libCVODES/nvector_serial.h"  
#include "../libCVODES/dense.h"           
#include "../libCVODES/sundialsmath.h"

using namespace std;

/* Accessor macros */

#define Ith(v,i)    NV_Ith_S(v,i-1)       // i-th vector component
#define IJth(A,i,j) DENSE_ELEM(A,i-1,j-1) // (i,j)-th matrix component

#define RTOL  RCONST(1e-8) 
#define ATOL RCONST(1e-8) 


//user data
typedef struct {
  double *p;
  long nvar;
  int sensi;
} *UserData;


//prototypes 
void f(realtype t, N_Vector y, N_Vector ydot, void *f_data);

static int ewt(N_Vector y, N_Vector w, void *e_data);

void sensDerivs (double t, double *Yt, double *Ft, double *p);

void integrateCVODES(Glob *globs,double ystart[], double p[], long nvar, double x1, 
		     double x2,double eps,int sensi)
{
  
  void *cvode_mem;
  UserData data;
  realtype t,t0,t1;
  N_Vector y;
  long i;

  cvode_mem = NULL;
  data      = NULL;
  y         = NULL;

 
  data = (UserData) malloc(sizeof *data);
  
  //user data
  data->p=p;
  data->nvar=nvar;
  data->sensi=sensi;
  
  //  Initial conditions 
  y = N_VNew_Serial(nvar);
  for(i=1;i<=nvar;i++)
    Ith(y,i)=ystart[i];

  t0= RCONST(x1);
  t1= RCONST(x2);
  t = RCONST(t1);

  /* Create CVODES object */
  cvode_mem = CVodeCreate(CV_BDF, CV_NEWTON);
  /* Allocate space for CVODES */
  realtype atol=ATOL;
  CVodeMalloc(cvode_mem, f, t0, y, CV_SS, RTOL, &atol);
  /* Attach user data */
  CVodeSetFdata(cvode_mem, data);
  /* Attach linear solver */
  CVDense(cvode_mem, nvar);
  
  //integrate
  CVode(cvode_mem, t1, y, &t, CV_NORMAL);

  //copy output back
  for(i=1;i<=nvar;i++)
    ystart[i]=Ith(y,i);

  //free memory
  N_VDestroy_Serial(y);
  free(data);
  CVodeFree(cvode_mem);
}

/*
 * f routine. Compute f(t,y). 
 */

void f(realtype t, N_Vector y, N_Vector ydot, void *f_data)
{
  UserData data;
  data = (UserData) f_data;
  long i,nvar=data->nvar;
  double *ydot_nr=(double*)malloc((nvar+1)*sizeof(double));
  double *y_nr=(double*)malloc((nvar+1)*sizeof(double));

  for(i=1;i<=nvar;i++)
    y_nr[i]=Ith(y,i);
  
  if(data->sensi==TRUE)
    sensDerivs((double)t,y_nr,ydot_nr,data->p);
  else
    ode((double)t,y_nr,ydot_nr,data->p);
  
  for(i=1;i<=nvar;i++)
      Ith(ydot,i)=ydot_nr[i];

  free(y_nr);
  free(ydot_nr);

}

 
