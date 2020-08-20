#include<iostream>
#include<math.h>
#include<stdlib.h>
#include<time.h>

#include "../def.h"
#include "../model.h"
#include "../nr.h"

using namespace std;

void outSimit(GlobExp ex[],Glob *globs,double *t,long n,double **y)
{
  long k,l;
  long idum=-time(NULL);
  long nobs=ex[1].nobs;
  ofstream out;
  
  out.open(ex[1].fileName);

  if(globs->sig==0)
    {
      for(k=1;k<=n;k++)
	{
	  out << t[k];
	  for(l=1;l<=nobs;l++)
	    out << " " << y[k][l] << " 1";
	  out << endl;
	}

    }
  else
    {
      double *max=dvector(1,nobs),*min=dvector(1,nobs);
      for(l=1;l<=nobs;l++)
	{
	  max[l]=y[1][l];
	  min[l]=y[1][l];
	}
      for(k=2;k<=n;k++)
	{
	  for(l=1;l<=nobs;l++)
	    {
	      if(max[l]<y[k][l])
		max[l]=y[k][l];
	      if(min[l]>y[k][l])
		min[l]=y[k][l];
	    }
	}
      for(k=1;k<=n;k++)
	{
	  out << t[k] ;
	  for(l=1;l<=nobs;l++)
	    out << " " << y[k][l]+globs->sig*fabs(max[l]-min[l])*gasdev(&idum) << " " << globs->sig*fabs(max[l]-min[l]);
	  out << endl;
	}
    }

}
