/****************************************************************/
/*                                                              */
/* Programm zum Aufruf von ODESSA                               */
/* Aufruf-Syntax ist diesselbe wie bei CALL_DOP853              */
/* Ziel: Sensitivitätsanalyse der Differentialgleichung ODE     */
/* bereitgestellte Routinen:                                    */
/*                          ODE                                 */
/*                          AD_ODE, OD_ODE                      */
/*                          JAC                                 */
/* Thorsten Mueller, Apr, 21, 1999                              */
/* und September, 18, 2001 (Krakau, Polen)                      */
/*                                                              */
/****************************************************************/
/* modified for inclusion in odin by w.h. 11/01
 */

#include <math.h>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../def.h"
#include "../model.h"
#include "../nr.h"

using namespace std;


inline void _dcopy (double * dest, double * src, long n)
{
  memcpy (dest + 1, src + 1, n * sizeof (double));
}

inline void _dfill0 (double * x, long n)
{
  memset (x + 1, 0, n * sizeof (double));
}

void _initYt (double *Yt, double *state)
{
  // initialise Y
  if (state)
    _dcopy (Yt, state, NEQNS);    // original state vector
  _dfill0 (Yt + NEQNS, (NPARAMS + NEQNS) * NEQNS);    // 0 in middle and rear part
  for (long i = 1; i <= NEQNS; i++)       // identity matrix in rear part dy/dy0
    Yt[(NPARAMS + i) * NEQNS + i] = 1;
}                       


extern "C" int odessa_ (int NEQ[],
                        double Y[],
                        double PAR[],
                        double *T,
                        double *TOUT,
                        int *ITOL,
                        double *RTOL,
                        double *ATOL,
                        int *ITASK,
                        int *ISTATE,
                        int IOPT[],
                        double RWORK[],
                        int *LRW, int IWORK[], int *LIW, int *MF);

extern "C" void dop_komponente_ (int NEQ[], int &J_OPT, int &J_GES);

void call_odessa(Glob *globs,GlobExp *ex,int N, int M_PAR_GES, char *doP,
   int D_FLAG, int MAX_NR_STEP, int STIFF, int inhomogen, int nmesh,
   double *tmesh, double t0, double t1, double *zustand, double **y,
   double ***dmds, double ***dmdp, double **dyds, double **dydp,
   double *parameter, double rtol, double atol,
   double min_stepsize, double max_stepsize, double initial_stepsize)
{

/* from DOKUMENTATION/call_odessa.doc
obsolet: N	nVar=anzahl der gleichungen
obsolet: M_PAR_GES, nPExp=gesamt_anzahl der parameter
char *doP	[1..nPExp] parameter fixen ("0") oder optimieren (sonst)
hieraus berechnet:
    M_PAR_SENS	anzahl zu variierender parameter
    M_GES = 	nVar+M_PAR_SENS; anzahl fit-variablen
int D_FLAG	verwende OD_ODE (1) oder AD_ODE (0)
int MAX_NR_STEP	maximale anzahl der integrationsschritte
int STIFF	steifes / nicht-steifes problem (1 / 0)
int inhomogen	liegt Inhomogenitätsmatrix vor (1)? bestimme die ableitungen
		ansonsten über finite differenzen (0)
int nmesh	anzahl der punkte, auf denen der zustand berechnet werden soll
double *tmesh		zeitpunkte, an denen der zustand berechnet werden soll
obsolet: sens_mesh; 	Sensitivitäten werden entweder gar nicht berechnet
			(dmds==NULL) oder am Endpunkt UND an den Datenpunkten
			(dmds,dmdp,dyds,dydp!=NULL)
double t0		Startpunkt, kann < tmesh[1] sein
double t1		Endpunkt, kann > tmesh[nmesh] sein
double **y		[k=1..nmesh][1..nVar] output: y(tmesh[k])
double *zustand		[1..nVar] input: y(t0), output: y(t1)
wie in tabVal:
double ***dmds		[k=1..nmesh][i=1..nVar][j=1..nVar] dy[i]/dy0[j](jacobi)
double ***dmdp		[k=1..nmesh][i=1..nVar][j=1..nPExp] (hiess sensitiv)
			output: Sensitivitäten dy[i]/dp[j] bei t=tmesh[k]
double **dyds		[i=1..nVar][j=1..nVar]         dy[i]/dy0[j] bei t=t1
double **dydp		[i=1..nVar][j=1..nPExp] dy[i]/dp[j]  bei t=t1
double *parameter	[1..nPExp] parameter-vektor
double rtol		error test passes if yerr<rtol*|y|+atol
double atol
double min_stepsize	minimale schrittweite
double max_stepsize	maximale schrittweite
double initial_stepsize
 */

  /*==================*/
  /* lokale variablen */
  /*==================*/

  int i = 0, j = 0, k = 0, M_GES, M_PAR_SENS = 0;
  long nPExp=globs->npar;
  long nVar=ex->nvar;
  double * gy=dvector(1,NOBS);
  int generic=TRUE;
  
  /*===============================*/
  /* berechne M_GES und M_PAR_SENS */
  /*===============================*/

  for (i = 1; i <= nPExp; i++)
    if (doP[i] != '0')
      M_PAR_SENS++;
  M_GES = nVar + M_PAR_SENS;

  /*=======================*/
  /* stelle NEQ() zusammen */
  /*=======================*/

  int IPAR_offset = 10;
  // wenn diese größe verändert wird, muß berechne_doP_komponente.f
  // ebenfalls angepaßt werden !!

  int *NEQ;
  NEQ = ivector (1, IPAR_offset + nPExp);

  // stecke den gesamten doP-vektor (in int umgewandelt) in NEQ,
  // sodass dies global zugänglich ist!
  for (i = 1; i <= nPExp; i++)
    NEQ[IPAR_offset + i] = (doP[i] != '0' ? 1 : 0);

  // alle NEQ-einträge die startwerte betreffend werden
  // in berechne_komponente auf null gesetzt, sodass keine
  // inhomogenitäten für diese "Parameter" berechnet werden

  NEQ[1] = nVar;
  NEQ[2] = M_GES;               // = nVar + M_PAR_SENS
  NEQ[3] = M_PAR_SENS;          // Anzahl Fitparameter, <=nPExp
  NEQ[4] = D_FLAG;
  NEQ[5] = nPExp;

  /*===================*/
  /* stelle Y zusammen */
  /*===================*/

  double **Yt = dmatrix (0, M_GES, 1, nVar);
  //Y = dvector(1,nVar + nVar*M_GES);
  for (i = 1; i <= nVar; i++)
    Yt[0][i] = zustand[i];

  for (j = 1; j <= M_PAR_SENS + nVar; j++)
    {

      // initialisiere: Sij = 0 für dy/dp

      if (j <= M_PAR_SENS)
        for (i = 1; i <= nVar; i++)
          Yt[j][i] = 0.0;

      // initialisiere: Sii = 1 für dyi/dyi°

      else
        for (i = 1; i <= nVar; i++)
          if ((j - M_PAR_SENS) == i)
            Yt[j][i] = 1.0;
          else
            Yt[j][i] = 0.0;
    }                           //nächster zustand i, der abgeleitet werden soll
  //Yt[2][2]=0;//!!

  /*=====================*/
  /* stelle PAR zusammen */
  /*=====================*/

  double *PARAMETER;
  PARAMETER = dvector (1, nPExp + nVar * nVar);
  //(1,nPExp) wird von ODESSA erwartet, der rest ist zur vereinfachung

  //übergebe an PARAMETER alle Parameter!
  int zaehler = 1;
  for (i = 1; i <= nPExp; i++)
    PARAMETER[i] = parameter[i];

  /*=================================*/
  /* stelle fehlerparameter zusammen */
  /*=================================*/

  int SqualCheck=FALSE;
  int ITOL = SqualCheck ? 1 : 4;        // erlaubt: 1 or 4
  int nTOL = (ITOL == 4) ? nVar * (M_GES + 1) : 1;

  double *RTOL = dvector (1, nTOL);
  double *ATOL = dvector (1, nTOL);
  if (ITOL == 4)
    {
      for (i = 1; i <= nVar; i++)
        {
          RTOL[i] = rtol;       // 1..nVar: rtol
          ATOL[i] = atol;
        }
      for (; i <= nTOL; i++)
        {
          RTOL[i] = rtol * 1e8; // nVar+1..: fast ignoriert
          ATOL[i] = atol * 1e8;
        }
    }
  else
    {
      RTOL[1] = rtol;           // alle: rtol
      ATOL[1] = atol;
    }

  /*===================================*/
  /* stelle kontrollparameter zusammen */
  /*===================================*/

  int ITASK = 1;                //siehe ODESSA-beschreibung -> normale integration von T bis TOUT
  // variable der status-abfrage
  int ISTATE = 1;               //-> erster aufruf von ODESSA!

  int *IOPT;
  IOPT = ivector (1, 3);
  IOPT[1] = 1;                  //zusätzliche angaben über RWORK-feld
  IOPT[2] = 1;                  //sensitivitäten werden berechnet !
  IOPT[3] = inhomogen;          //verwende inhomogenitäten, sonst finite differenzen

  int MITER = 1;                // -> volle Jacobi-Matrix vorhanden
  int METH;                     // normaler / steifer Integrator
  int MAXORD;                   // maximale ordnung der verwendeten approximationen

  if (STIFF == 0)
    {
      METH = 1;
      MAXORD = 12;
    }
  else
    {
      METH = 2;
      MAXORD = 5;
    }
  int MF = 10 * METH + MITER;
  //cout<<"mf = "<<MF<<LF;

  /*=================================*/
  /* stelle arbeitsspeicher zusammen */
  /*=================================*/

  int LRW;
  int LIW;

  if (IOPT[1] == 0)
    {                           //ohne sensitivitäten
      LRW = 20 + nVar * (MAXORD + 1) + 3 * nVar + nVar * nVar + 2 + 1;
      LIW = 20 + nVar;
    }
  if (IOPT[1] == 1)
    {                           //mit sensitivitäten
      LRW = 20 + nVar * (M_GES + 1) * (MAXORD + 1) + 2 * nVar * (M_GES + 1) +
        nVar * nVar + 2 + nVar + 1;
      LIW = 21 + nVar + M_GES;
    }

  double *RWORK;
  RWORK = dvector (1, LRW);

  int *IWORK;
  IWORK = ivector (1, LIW);

  /*==================================*/
  /* stelle optionalen input zusammen */
  /*==================================*/

  // odeint detects stepsize underflow independently of hmin
  // odessa relies on min_stepsize being large enough to avoid prevent it
  if (t0 + min_stepsize == t0)
    {
      cerr<< t0 << "\t" << min_stepsize << endl;
      throw 1;
    }
  if (t1 + min_stepsize == t0)
    {
      cerr<< t0 << "\t" << min_stepsize << endl;
      throw 1;
    }

  RWORK[5] = initial_stepsize;
  RWORK[6] = max_stepsize;
  RWORK[7] = min_stepsize;

  IWORK[5] = MAXORD;
  IWORK[6] = MAX_NR_STEP;
  IWORK[7] = 3;                 // not more than 3 warnings of type T=T+H

  /*==================================*/
  /* bestimme T und TOUT mittels tmesh */
  /*==================================*/


  double T = t0;
  double TOUT;                  //wird mittels gitter vorgegeben


  if (dmds)
    {                           // Sensitivitäten speichern
      double **dgdy = NULL, **dgdp = NULL;
      dgdy = dmatrix (1, NOBS, 1, nVar);
      dgdp = dmatrix (1, NOBS, 1, nPExp);
      _dfill0 (dgdy[1], NOBS * nVar);
      _dfill0 (dgdp[1], NOBS * nPExp);


      for (k = 1; k <= nmesh + 1; k++)
        {
          if (k > nmesh)
            TOUT = t1;          // Endpunkt
          else
            TOUT = tmesh[k];    // Datenpunkt

          /*=================*/
          /* rufe ODESSA auf */
          /*=================*/

          odessa_ (&NEQ[1],
                   &Yt[0][1],
                   &PARAMETER[1],
                   &T, &TOUT, &ITOL, &RTOL[1], &ATOL[1],
                   &ITASK, &ISTATE, &IOPT[1], &RWORK[1], &LRW, &IWORK[1],
                   &LIW, &MF);

          if (ISTATE < 0)
            break;              // skip remaining k

          // Zustand speichern
          int jlong;            // 1..nPExp
          if (k > nmesh)
            {                   // Endpunkt
              _dcopy (zustand, Yt[0], nVar);

              for (j = 1; j <= M_PAR_SENS; j++)
                {
                  dop_komponente_ (&NEQ[1], j, jlong);

                  //übernehme die ableitung nach den parametern
                  for (i = 1; i <= nVar; i++)
                    dydp[i][jlong] = Yt[j][i];
                  // jlong and j were exchanged until 3/02
                }

              for (j = 1; j <= nVar; j++)
                //übernehme die ableitung nach den startwerten
                for (i = 1; i <= nVar; i++)
		  {
		    dyds[i][j] = Yt[M_PAR_SENS + j][i];
		  }
            }
          else
            {                   // k<=nmesh: Datenpunkt
              _dcopy (y[k], Yt[0], nVar);

              // write sensitivities to dm/dp and dm/ds
              // support for extra observables as in odeExpment::tabVal
              //if (nVarObs > nVar)
              //  observation (globs,ex,TOUT, y[k], globs->par, dgdy, dgdp);
	      generic=observation(globs,ex,TOUT, y[k],gy, ex->par, dgdy, dgdp);
	      for(i=1;i<=NOBS;i++)
		y[k][i]=gy[i];
	      
              for (i = 1; i <= NOBS; i++)
                {
                  int ivar = i;
                  for (j = 1; j <= M_PAR_SENS; j++)
                    {
		      dop_komponente_ (&NEQ[1], j, jlong);
                      double &dest = dmdp[k][i][jlong];
		      if (!generic)
                        {
                          dest = dgdp[i][j];
                          for (long l = 1; l <= nVar; l++)
                            dest += dgdy[i][l] * Yt[j][l];
                        }
                      else
                        dest = Yt[j][i];
                    }           // for j
                  for (j = 1; j <= nVar; j++)
                    {
                      double &dest = dmds[k][i][j];
                      if (!generic)
                        {
                          dest = 0;
                          for (long l = 1; l <= nVar; l++)
                            dest +=  dgdy[ivar][l] * Yt[M_PAR_SENS + j][l];
                        }
                      else
                        dest = Yt[M_PAR_SENS + j][ivar];
                    }           // for j
                }               // for i
            }                   // if k>nmesh else

          T = TOUT;             // für nächstes k
        }                       //nächster zeitpunkt k
      //BEGIN eingefuegt von Felix
      if (dgdy)
        free_dmatrix (dgdy, 1, NOBS, 1, nVar);
      if (dgdp)
        free_dmatrix (dgdp, 1, NOBS, 1, nPExp);
      //END
    }
  else
    {                           // !dmds

      for (k = 1; k <= nmesh + 1; k++)
        {
          if (k > nmesh)
            TOUT = t1;          // Endpunkt
          else
            TOUT = tmesh[k];    // Datenpunkt

          odessa_ (&NEQ[1],
                   &Yt[0][1],
                   &PARAMETER[1],
                   &T, &TOUT, &ITOL, &RTOL[1], &ATOL[1],
                   &ITASK, &ISTATE, &IOPT[1], &RWORK[1], &LRW, &IWORK[1],
                   &LIW, &MF);

          if (ISTATE < 0)
            break;              // skip remaining k

          // Zustand speichern
          if (k > nmesh)        // Endpunkt
            _dcopy (zustand, Yt[0], nVar);
          else
            {                   // k<=nmesh: Datenpunkt
              _dcopy (y[k], Yt[0], nVar);
	      
              // compute extra observables
              generic=observation (globs,ex,TOUT, y[k],gy,ex->par, NULL, NULL);
	      for(i=1;i<=NOBS;i++)
		y[k][i]=gy[i];
	    }
	  T = TOUT;             // für nächstes k
        }                       //nächster zeitpunkt k

    }                           // if dmds else


  //////////// speicherfreigabe ////////////////

  // dvector

  free_dmatrix (Yt, 0, M_GES, 1, nVar);
  free_dvector (PARAMETER, 1, nPExp + nVar * nVar);

  free_dvector (RTOL, 1, nTOL);
  free_dvector (ATOL, 1, nTOL);
  free_dvector (RWORK, 1, LRW);
  free_dvector (gy,1,NOBS);

  // ivector

  free_ivector (NEQ, 1, IPAR_offset + nPExp);
  free_ivector (IWORK, 1, LIW);
  free_ivector (IOPT, 1, 3);

  /*===================================================*/
  /* abfang eventueller fehler / auswertung von ISTATE */
  /*===================================================*/

  char *odessaErrStr[6] = {
    "excess work done on this call (perhaps wrong MF)",
    "excess accuracy requested (tolerances too small)",
    "illegal input detected (see printed message)",
    "repeated error test failures (check all inputs)",
    "repeated convergence failures\n"
      "(perhaps bad jacobian supplied or wrong choice of MF or tolerances)",
    "error weight became zero during problem\n"
      "(solution component I,J vanished, and ATOL or ATOL(I,J) = 0.0)"
  };

  if (ISTATE < 0)
    {
      cerr << odessaErrStr[-ISTATE - 1] <<endl;
      throw 1;
    }
}                              // odeExpment::call_odessa


/* wrapper routines between odessa and the dynamical equations generated
 * from the .mdf
 */



extern "C" void ode_ (int &n, int &m, double &x, double *y, double *dydx, double *p)
{
  ode (x, y - 1, dydx - 1, p - 1);
}                               // ode_

extern "C" void jac_ (int *NEQ, double &t, double *y, int &n, double *p, int &M, double *pd)
{
  double **dfdy = new double *[NEQNS] - 1;
  for (int i = 1; i <= NEQNS; i++)
    dfdy[i] = pd - 1 + (i - 1) * NEQNS;
  jacobi (t, y - 1, dfdy, p - 1);
  delete[]++ dfdy;
}

extern "C" void od_ode_ (int *NEQ, double &t, double *y, int &N, double *p,
         int &M, double *dfdp, int &j)
{
  if (!j)
    return;
  inhomo1 (t, y - 1, dfdp - 1, p - 1, j);
}                               // od_ode_

extern "C" void ad_ode_ (int &, int &, int &, double &, double *, double *, double *,
         int &, double *, double *, int &)
{
  cerr << "unexpected call to AD_ODE" <<endl;
}
