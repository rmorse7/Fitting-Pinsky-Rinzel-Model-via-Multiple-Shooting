#include<fstream>
#include<iostream>
#include<stdlib.h>
#include<math.h>
#include "../nr.h"

using namespace std;

//----------------------------------------------------------------------
//                  SPEZIELLE KLASSE FUER BANDMATRIZEN
//----------------------------------------------------------------------


//Konstruktor
Bandmatrix::Bandmatrix(long ndata,long bandbreite) : n(ndata), bw(bandbreite)
{
  long k,l;
  if(bw<0)
    {
      cerr << "Bandwidth insufficient !!\n";
      throw 1;
    }
  band=dmatrix(1,n,0,2*bw+1);
  //die Bandmatrix mit Null initialisieren
  for(k=1;k<=n;k++)
    {
      for(l=0;l<=2*bw+1;l++)
	{
	  band[k][l]=0.;
	}
    }
}

//Lesen der Eintraege
double Bandmatrix::r(long i,long j)
{
  if(abs(i-j)>bw)
    {
//      cerr << "*";
      return(0.);
    }
  else
    {
      return(band[i][bw-(i-j)]);
    }
}

//Schreiben der Eintraege
void Bandmatrix::w(double value,long i,long j)
{
  if(abs(i-j)>bw);
  else
    {
      band[i][bw-(i-j)]=value;
    } 
}

//Destruktor
Bandmatrix::~Bandmatrix()
{
  free_dmatrix(band,1,n,0,bw);
}


//------------------------------------------------------------

void splines(double *x,double *y,double *w,long n,double alpha,double *g,double *gam)
//Berechnet smoothing splines (natural cubic splines) nach dem
//Reinsch-Algorithmus
// double *x : n-Vektor, design points
// double *y : n-Vektor, Daten
// double *w : n-Vektor, Gewichte
// long n : #Datenpunkte
// double alpha: smoothing-Parameter
// double *g : n-Vektor, Ausgabewerte
// double *gam: n-2 Vector, Ausgabe zweite Ableitung
{
  double *h,*b,*p,dummy;
  long k,l,m;

  if(n<=2)
    {
      cerr << "Error : n < 3 !!! \n";
      throw 1;
    }
  if(alpha<0)
    alpha=alpha*(-1);

  h=dvector(1,n-1);
  Bandmatrix Q(n,2),A(n-2,2);
   
  b=dvector(1,n-2);
  p=dvector(1,n-2);

  //Initialisierung von h,b
  //h
  for(k=1;k<=n-1;k++)
    h[k]=x[k+1]-x[k];
  //b
  for(k=1;k<=n-2;k++)
    b[k]=(y[k+2]-y[k+1])/h[k+1]-(y[k+1]-y[k])/h[k];
  // Q & A

  for(k=1;k<=n-2;k++)
    {
      Q.w(1./h[k],k,k);
      Q.w(-1./h[k]-1./h[k+1],k+1,k);
      Q.w(1./h[k+1],k+2,k);
    }

  // Ausnutzen der Bandstruktur fuer schnelle Berechnung von Q^t*Q

  for(l=-2;l<=2;l++)
    {
      for(m=-2;m<=2;m++)
	{
	  for(k=1;k<=n;k++)
	    {
	      if(k+l >=1 && k+m >=1 && k+l <=n-2 && k+l <=n-2)
		A.w(A.r(k+l,k+m)+Q.r(k,k+l)*(1./w[k])*Q.r(k,k+m),k+l,k+m);
	    }
	}
    }

  for(k=1;k<=n-2;k++)
    {
      for(l=-2;l<=2;l++)
	{
	  if(k+l >=1 && k+l <=n-2)
	  A.w(alpha*A.r(k,k+l),k,k+l);
	}
    }

  for(k=1;k<=n-2;k++)
    {
      A.w(A.r(k,k)+(1./3.)*(h[k]+h[k+1]),k,k);
      if(k<=n-3)
	{
	  A.w(A.r(k,k+1)+h[k+1]/6.,k,k+1);
	  A.w(A.r(k+1,k)+h[k+1]/6.,k+1,k);
	}
    }
  //Loesen der Gleichungen mittels Cholesky-Zerlegung
  Choldc(&A,n-2,p);
  Cholsl(&A,n-2,p,b,gam);
   
   //Berechen der Ausgabewerte g

  for(k=1;k<=n;k++)
    {
      dummy=0.;
      for(l=1;l<=n-2;l++)
	{
	  dummy+=Q.r(k,l)*gam[l];
	}
      g[k]=y[k]-alpha*(1./w[k])*dummy;
    }
  
   
   
   free_dvector(h,1,n-1);
   free_dvector(b,1,n-2);
   free_dvector(p,1,n-2);
}


double Score(double alpha,double *x,double *y,double *w,double *b,double *g,double *gam,Bandmatrix *Q,Bandmatrix *R,long n)
  //Berechnet wird der GCV-Score
{
  double score=0.,dummy;
  double spur=0.; // tr(S)
  long k,l,m;
  double *d;

  Bandmatrix A(n-2,2),B(n-2,2);
  d=dvector(1,n-2);

  //Spline berechnen (Reinsch Algorithmus)

  //Uebertragen der Elemente von R auf A
  for(k=1;k<=n-2;k++)
    {
      for(l=-2;l<=2;l++)
	{
	  if(k+l>=1 && k+l <=n-2)
	    A.w(R->r(k,k+l),k,k+l);
	}
    }
  
  for(l=-2;l<=2;l++)
    {
      for(m=-2;m<=2;m++)
	{
	  for(k=1;k<=n;k++)
	    {
	      if(k+l >=1 && k+m >=1 && k+l <=n-2 && k+l <=n-2)
		A.w(A.r(k+l,k+m)+alpha*Q->r(k,k+l)*(1./w[k])*Q->r(k,k+m),k+l,k+m);
	    }
	}
    }

  Choldc(&A,n-2,d);
  Cholsl(&A,n-2,d,b,gam);  

  
  for(k=1;k<=n;k++)
    {
      g[k]=y[k];
      for(l=k-2;l<=k+2;l++)
	{
	  if(l>=1 && l<=n-2)
	    g[k]-=alpha*(1./w[k])*Q->r(k,l)*gam[l];
	}
    }
  //Ende Spline berechnen
  cerr << "*";
  //Spur berechnen nach nach Hutchison & de Hoog
  B.w(1./(d[n-2]*d[n-2]),n-2,n-2);
  B.w((-A.r(n-2,n-3)*B.r(n-2,n-2))/d[n-3],n-3,n-2);
  B.w(B.r(n-3,n-2),n-2,n-3);
  B.w(((1./d[n-3])-A.r(n-2,n-3)*B.r(n-3,n-2))/d[n-3],n-3,n-3);

  for(k=n-4;k>=1;k--)
    {
      B.w((-A.r(k+1,k)*B.r(k+1,k+2)-A.r(k+2,k)*B.r(k+2,k+2))/d[k],k,k+2);
      B.w(B.r(k,k+2),k+2,k);
      B.w((-A.r(k+1,k)*B.r(k+1,k+1)-A.r(k+2,k)*B.r(k+1,k+2))/d[k],k,k+1);
      B.w(B.r(k,k+1),k+1,k);
      B.w(((1./d[k])-A.r(k+1,k)*B.r(k,k+1)-A.r(k+2,k)*B.r(k,k+2))/d[k],k,k);
    }


  
  for(k=1;k<=n;k++)
    {
      dummy=0.;
      for(l=-2;l<=2;l++)
	{
	  for(m=-2;m<=2;m++)
	    {
	      if(k+l >=1 && k+m >=1 && k+l <=n-2 && k+l <=n-2)
		dummy+=Q->r(k,k+l)*B.r(k+l,k+m)*Q->r(k,k+m);
	    }
	}
      spur+=1.+alpha*(1./w[k])*dummy;
    }

  free_dvector(d,1,n-2);
 
  //Ende Spur berechnen

  //Score berechnen
  for(k=1;k<=n;k++)
    {
      score+=w[k]*pow(y[k]-g[k],2);
    }
  score*=1./((double)n*pow(1.-spur/(double)n,2));
  return(score);
  
}

//Golden Section Minimierungsroutiene -> Recipes

#define RR 0.61803399
#define C (1.0-RR)
#define SHFT2(a,b,c) (a)=(b);(b)=(c);
#define SHFT3(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);

double golden(double ax, double bx, double cx, double tol,double *x,double *y,double *w,double *b,double *g,double *gam,Bandmatrix *Q,Bandmatrix *R,long n)
{
  double f1,f2,x0,x1,x2,x3;
  
  x0=ax;
  x3=cx;
  if (fabs(cx-bx) > fabs(bx-ax)) {
    x1=bx;
    x2=bx+C*(cx-bx);
  } else {
    x2=bx;
    x1=bx-C*(bx-ax);
  }
  f1=log(Score(x1,x,y,w,b,g,gam,Q,R,n));
  f2=log(Score(x2,x,y,w,b,g,gam,Q,R,n));
  while (fabs(x3-x0) > tol*(fabs(x1)+fabs(x2))) {
    if (f2 < f1) {
      SHFT3(x0,x1,x2,RR*x1+C*x3)
	SHFT2(f1,f2,log(Score(x2,x,y,w,b,g,gam,Q,R,n)))
	} else {
	  SHFT3(x3,x2,x1,RR*x2+C*x0)
	    SHFT2(f2,f1,log(Score(x1,x,y,w,b,g,gam,Q,R,n)))
	    }
  }
  if (f1 < f2) {
    
    return x1;
  } else {
    
    return x2;
  }
}

#undef C
#undef RR
#undef SHFT2
#undef SHFT3


void splines_gcv(double *x,double *y,double *w,long n,double *alpha,double *g,double *gam,double maxa)
//Berechnet smoothing splines (natural cubic splines) nach dem
//Reinsch-Algorithmus der Glaettungsparameter wird durch Generalized Crossvalidation
//geschaetzt.
// double *x : n-Vektor, design points
// double *y : n-Vektor, Daten
// double *w : n-Vektor, Gewichte
// long n : #Datenpunkte
// double alpha: smoothing-Parameter
// double *g : n-Vektor, Ausgabewerte
// double *gam: n-2 Vector, Ausgabe zweite Ableitung
{
  double *h,*b,dummy;


  long k,l,m;

  if(n<=2)
    {
      cerr << "Error : n < 3 !!! \n";
      throw 1;
    }
  h=dvector(1,n-1);
  b=dvector(1,n-2);
  Bandmatrix Q(n,2),R(n,2);
  //Initialisierung von h,b
  //h
  for(k=1;k<=n-1;k++)
    h[k]=x[k+1]-x[k];
  //b
  for(k=1;k<=n-2;k++)
    b[k]=(y[k+2]-y[k+1])/h[k+1]-(y[k+1]-y[k])/h[k];

  //Q
  for(k=1;k<=n-2;k++)
    {
      Q.w(1./h[k],k,k);
      Q.w(-1./h[k]-1./h[k+1],k+1,k);
      Q.w(1./h[k+1],k+2,k);
    }
  //R
  for(k=1;k<=n-2;k++)
    {
      R.w((1./3.)*(h[k]+h[k+1]),k,k);
      if(k<=n-3)
	{
	  R.w((1./6.)*h[k+1],k,k+1);
	  R.w((1./6.)*h[k+1],k+1,k);
	}
    }
  
  
  //Minimum vom Score berechnen
  double ax=1e-10;
  double cx=maxa;

  cerr << "GCV : ";
  
  (*alpha)=golden(ax,(*alpha),cx,1e-4,x,y,w,b,g,gam,&Q,&R,n);
  
  //Golden Section Minimierungsroutine

  cerr << "\n Alpha = " << *alpha << endl;
  free_dvector(h,1,n-1);
  free_dvector(b,1,n-2);
}

