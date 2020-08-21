#include <math.h>
#include <stdlib.h>
#include "../nr.h"

using namespace std;

void Cholsl(Bandmatrix *a, int n, double p[], double b[], double x[])
{
	int i,k;
	double sum;

	for (i=1;i<=n;i++) {
		for (sum=b[i],k=i-1;k>=i-2;k--) 
		  {
		      if(k>=1 && abs(i-k)<=2)
		      sum -= a->r(i,k)*x[k];
		  }
		x[i]=sum/p[i];
	}
	for (i=n;i>=1;i--) {
		for (sum=x[i],k=i+1;k<=i+2;k++) 
		  {
		    if(k<=n && abs(i-k)<=2)
		      sum -= a->r(k,i)*x[k];
		  }
		x[i]=sum/p[i];
	}
}
