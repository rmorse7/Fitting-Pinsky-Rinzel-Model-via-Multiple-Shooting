#include<fstream>
#include<stdio.h>
#include<string>

using namespace std;


#ifndef TRUE
#define TRUE 1
#define FALSE 0
#endif

typedef struct 
{
  // experimental data
  //! File name of experiment
  char *fileName;
  //! Number of observed quantities
  long nobs;
  /*! \brief Number of variables/equations 

    \f[ {\bf nobs} \le {\bf nvar}\f]
  */
  long nvar;    

  // the actual data
  //! Number of measured time points
  long nMeasure;
  /*! \brief Temporal data points

      Vector of dimension: \f[ {\bf nMeasure} \f]
  */
  double *xMeasure;
  /*! \brief Data of the observations

      Matrix of dimension: \f[ {\bf nMeasure}\times{\bf nobs}\f]
  */
  
  double **yMeasure;
  /*! \brief Standard deviation of data

      Matrix of dimension:  \f[ {\bf nMeasure}\times{\bf nobs}\f]
  */
  double **sigma;
  //! Fitting range
  double fitstart, fitend; 
  //! First point in fitting range
  long firstMeasure;
  //! Last point in fitting range
  long lastMeasure;
  // multiple shooting data
  /*! \brief Number of multiple shooting intervals

    Command line argument is \b -nms and not \b -nPoints !!
  */
  long nPoints;
   /*! \brief Multiple shooting number of intervals
     
       Not command line
   */
  double *tms;
   /*! \brief Multiple shooting time boundaries
     
       Vector  of dimension:  \f[ {\bf nPoints} \f]
   */

  double *mesh;   
  /*! \brief Initial guesses at mesh points
     
       Matrix of dimension:  \f[ {\bf nPoints} \times {\bf nvar} \f]
   */
  double **yTry,**yTrySave;
  //! Number of equality constraints
  long me;
  //! Number of inequality constraints
  long mg;
  
  // data obtained by integration
  //! Current value of the objective function
  double objF;
   /*! \brief Predicted trajectory (intergation)
    
      Matrix of dimension:  \f[ {\bf nMeasure}\times{\bf nvar}\f]
  */ 
  double **yComp;
  /*! \brief Predicted observables (intergation)
    
      Matrix of dimension:  \f[ {\bf nMeasure}\times{\bf nobs}\f]
  */
  double **yPred;
  /*! \brief  Discrepancies/deviation from continuous trajectories
    
      \f[{\bf h}_{i,j}={\bf yComp}_{i,j}-{\bf yTry}_{i,j} \f]
      Matrix of dimension:  \f[ {\bf nPoints}\times{\bf nvar}\f]
  */
  double **h;
  /*! \brief Current residues

     \f[ residues_{i,j}=yPred_{i,j}-yMeasure_{i,j}\f]
     Therefore matrix of dimension:  \f[ {\bf nMeasure}\times{\bf nobs}\f]
  */
  double **residues;
  /*! \brief Right-hand side of equality constraints

      Vector of dimension:  \f[ {\bf me}\f]
  */
  double *r2;
  /*! \brief Right-hand side of inequality constraints

      Vector of dimension:  \f[ {\bf mg}\f]
  */  
  double *r3;
  /*! \brief Derivative of the trajectory with respect to 
    intial values inside the multiple shooting interval
    
    3-tensor of dimension:  \f[ {\bf nPoints} \times {\bf nvar} \times  {\bf nvar} \f]
  */
  double ***dyds;
  /*! \brief Derivative of the trajectory with respect to 
    parameters inside the multiple shooting interval
    
    3-tensor of dimension:  \f[ {\bf nPoints} \times {\bf nvar} \times  {\bf npar} \f]
  */
  double ***dydp;
  /*! \brief Derivative of the observed quantities with respect to 
    intial values inside the multiple shooting interval
    
    3-tensor of dimension:  \f[ {\bf nMeasure} \times {\bf nobs} \times  {\bf nvar} \f]
  */
  double ***dmds;
  /*! \brief Derivative of the observed quantities with respect to 
    parameters inside the multiple shooting interval
    
    3-tensor of dimension:  \f[ {\bf nMeasure} \times {\bf nobs} \times  {\bf npar} \f]
  */  
  double ***dmdp;
   /*! \brief Derivative of the equality constraints with respect to 
    intial values inside the multiple shooting interval
    
    3-tensor of dimension:  \f[ {\bf me} \times {\bf nPoints} \times  {\bf nvar} \f]
  */ 
  double ***dR2ds;
   /*! \brief Derivative of the equality constraints with respect to 
    parameters
    
    Matrix of dimension:  \f[ {\bf me} \times {\bf npar} \f]
  */   
  double **dR2dp;
   /*! \brief Derivative of the inequality constraints with respect to 
    intial values inside the multiple shooting interval
    
    3-tensor of dimension:  \f[ {\bf mg} \times {\bf nPoints} \times  {\bf nvar} \f]
  */   
  double ***dR3ds;
   /*! \brief Derivative of the inequality constraints with respect to 
    parameters
    
    Matrix of dimension:  \f[ {\bf mg} \times {\bf npar} \f]
  */ 
  double **dR3dp;
  /*! \brief Condensation of the least squares

      \f[||u_a+E_a \Delta s_0 +P_a \Delta p ||={\bf min} \f]
  */
  double *ua,**Ea,**Pa;
  /*! \brief Condensation of the equality constraints

      \f[u_e+E_e \Delta s_0 +P_e \Delta p=0 \f]
  */  
  double *ue,**Ee,**Pe;
  /*! \brief Condensation of the inequality constraints

      \f[u_g+E_g \Delta s_0 +P_g \Delta p \ge 0 \f]
  */ 
  double *ug,**Eg,**Pg;
  /*! \brief Initial value
    
     Vector of dimension:  \f[ {\bf nvar}\f]
  */
  double *y0;
  /*! \brief Parameters

      Vector of dimension:  \f[ {\bf npar}\;,\f]
      where \b npar is defined in structure \b Glob.
  */
  double *par;
  double **dS;
  double *dP;
  //! Name of spline data-file
  string *splineFile;
  //! Spline nodes
  double **splineNodes;
  //! Value at the spline nodes
  double **splineY;
  //! Second derivative at the spline nodes
  double **splineGam;
  //! Number of spline nodes
  unsigned long *nNodes;
  //! Flag: TRUE if splines are defined
  int splinesDefined;
  double *errP;
  double *errY0;
} GlobExp;

typedef struct
{
  //flags
  //! Flag: TRUE turns GnuPlot animation off
  int noGnu;
  int savegnupng;
  int noMeasurements;
  //! Flag: TRUE for waiting after each iteration 
  int wait;
  //! Flag: TRUE for using variances from datafile
  int usesig;
  //! Integrator selection
  int integrator;
  int stiff;
  //! Flag: TRUE to prevent waiting after last iteration
  int nowait;
  //! Flag: TRUE for turning automatic regularisation on
  int reg;
  //! Flag: TRUE for using simulation based data for initialisation
  int simInit;
  //! Flag: TRUE for suppress damping
  int nodamp;
  //all the rest
  //! Number of parameters
  long npar;
  //! Maximal number of iterations
  long maxit;
  //! Iteration counter
  long nIter;
  //! Number of experiments
  int nrExp;
  /*! \brief Fix and local parameters definition vector
      
      TRUE  : Fit parameter\n
      FALSE : Fix parameter\n
      'L'   : Local Parameter\n

      Vector of dimension:  \f[ {\bf npar}\f]
  */
  int *doP;
  /*! \brief Covariance matrix of parameters
    
      Matrix of dimension: \f[ {\bf npar} \times  {\bf npar}\f]
   */
  double **covar;
  //! Current chi-square
  double chisq;
  //! Integration accuracy
  double eps;
  //! Controls the integration of the variational differential equations
  long rkqs_ign;
  //! Maximal number of integration steps (Runge-Kutta)
  int maxstp;
  //! Minimal improvement of parameters
  double minimp;
  /*! \brief Factor for weakening the continuity constraints
    
      \f[{\bf elastic}\in (0,1] \f]
  */
  double elastic;
  double epsilon;
  double lambda;
  //! Sampling for simit
  double dt;
  //! Relative Gaussian for simulations
  double sig;
  //! Perturbation of the simulated initial guess in \%
  double pert;
  /*! \brief Fix initial values
      
      TRUE  : Fit initial value\n
      FALSE : Fix initial value\n

      Vector of dimension:  \f[ {\bf nvar}\f]
  */
  int *y0fix;
  //damping section
  //! Controls damping
  double wquer;

  //GNUplot
  //! File pointer vector for GnuPlot animation
  FILE **gnuFp;
  //! Number of GnuPlot windows
  long ngnu;
  //! Index of the current GnuPlot window
  long gnuindex;
  /*! Total number of fitted variables,
       needed for the covariance matrix
  */
  long fitdim;
  //! Flag: TRUE if splines should be initialised
  int initSpline;
  //! 1 if fit is converged, 0 else
  int fitConverged;
  /*! Condition number of the linearised
      system to be solved
  */
  double cond;
  /*! Actual damping parameter
      (not to be confused with lambda)
  */
  double Lambda;
  /*!  Relaxation Parameter
  */
  int silent;
  /*!  if silent=1 no shell output is printed
  */
 
  int strategy;

  int minimiser;

  int faktorLexist;
  double *faktorL;
} Glob;

/* DEBUG stream */
extern ofstream *dbg;

