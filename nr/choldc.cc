#include <math.h>
#include <iostream>
#include "../nr.h"
#include <stdlib.h>

using namespace std;

long geo(long in)
//geo: greater/equal one
{
  if(in>1)
    return(in);
  else
    return(1);
}

long len(long in,long n)
//gen: less/equal n
{
  if(in<n)
    return(in);
  else
    return(n);
}

void Choldc(Bandmatrix *a, int n, double p[])
{
  void nrerror(char error_text[]);
  long k,l,m;
  double sum;

  for(k=1;k<=n;k++)
    {
      for(l=k;l<=len(k+2,n);l++)
	{
	  for(sum=a->r(k,l),m=k-1;m>=1;m--)
	    {
	      if(labs(k-m)<=2 && labs(l-m)<=2)
		sum-=(a->r(k,m)*a->r(l,m));
	    }
	  if(k==l)
	    {
	      if(sum<=0.0)
		{
		  cerr << "Sum : " << sum << endl;
		  nrerror("choldc failed");
		}
	      p[k]=sqrt(sum);
	    }
	  else
	    a->w(sum/p[k],l,k);
	}
    }
}

