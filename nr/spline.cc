
double spline(double *x,double *g,double *gam,long n,double t)
// Auswerten der Splines
// double *x : n-Vektor Design-Points
// double *g : n-Vektor Stuetzstellen
// double *gam : n-Vektor zweite Ableitung Vorsicht: Andere Notation gam[1]=0 und gam[n]=0
// long n : #Punkte
// double t : Auswertstelle
{
  double out,dum;
  long k;

  if(t<=x[1])
    {
      dum=(g[2]-g[1])/(x[2]-x[1])-(1./6.)*gam[2]*(x[2]-x[1]);
      out=g[1]-(x[1]-t)*dum;
    }
  else if(t>=x[n])
    {
      dum=(g[n]-g[n-1])/(x[n]-x[n-1])+(1./6.)*gam[n-1]*(x[n]-x[n-1]);
      out=g[n]+(t-x[n])*dum;
    }
  else
    {
      k=1;
      while(t>x[k])
	{
	  k++;
	}
      k--;
      dum=x[k+1]-x[k];
      out=((t-x[k])*g[k+1]+(x[k+1]-t)*g[k])/dum-(1./6.)*(t-x[k])*(x[k+1]-t)*((1+(t-x[k])/dum)*gam[k+1]+(1+(x[k+1]-t)/dum)*gam[k]);
    }
  return(out);
}
