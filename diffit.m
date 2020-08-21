function [out,outex] = diffit(in,inex)
%Matlab interface for diffit, v 1.1 06/Sep./2005
%-----------------------------------------------
% 
%SYNOPSIS;
%  [out outex]= diffit(in,inex)
% 
%DESCRIPTION
%     in: Structure containing experiment unspecific options.
%   inex: Structure array containing experiment specific options.
%    out: Experiment unspecific output structure.
%  outex: Experiment specific output structure array.
% 
%FIELDS OF STRUCTURE in:
%  in.eps        <double>   integration accuracy
%  in.nognu      <0/1>      turns gnuplot animation off 
%  in.nomeasure  <0/1>      use evently spaced m.s. intervals
%  in.doP        <string>   specifies parameters to be fitted,
%                           e.g. 110L
%  in.maxit      <integer>  maximal number of iterations
%  in.int        <1/2/3>    0: RK (default), 1: CVODES, 3: ODESSA
%  in.minimp     <double>   stopping criterion, def.: 1e-4
%  in.reg        <0/1>      switchs regularisation on
%  in.epsilon    <double>   singular val. threshold, def.: 1e-10
%  in.lambda     <double>   reg. parameter, def.: 1e6
%  in.siminit    <0/1>      simulate initial state
%  in.pert       <double>   pertubate simulated init. state
%  in.y0fix      <string>   fixing initial values, e.g. 1100
%  in.wquer      <double>   controls damping,
%                           set negative for initialisation
%  in.nodamp     <0/1>      switches damping off
%  in.init       <0/1>      1: only initialisation, no fitting
%  in.silent     <0/1>      1: silent mode
%  in.opt        <1/2>      optimiser 1: LSEI 2: NAG
%  in.Lconst     <vector>   set local Parameter Constraints
% 
%FIELDS OF STRUCTURE ARRAY inex:
%  inex(exp).data    <matrix>    data to be fitted
%  inex(exp).sig     <matrix>    variance of data
%  inex(exp).t0      <double>    starting time
%  inex(exp).t1      <double>    stopping time
%  inex(exp).nms     <integer>   number of m.s. intervals
%  inex(exp).p       <vector>    parameters
%  inex(exp).y0      <vector>    initial values
%  inex(exp).spline  <cell mat.> spline data
%  inex(exp).mesh    <matrix>    mesh data
%   
%  The contents of in and inex and its default values are copied
%  to the output variables out and outex respectively.
%   
%  ADDITIONAL OUTPUT FIELDS OF STRUCTURE ARRAY out:
%    out.wquer     <double>   actual value of wquer (see above)
%    out.converged <integer>  1 if fit is converged,0 else
%    out.chisq     <double>   actual chi^2 value
%    out.except    <integer>  1 if exception occurred,0 else
%    out.Lambda    <double>   actual damping parameter
%    out.cond      <double>   condition number << 1 for
%                             badly conditioned systems
%                             (only disp. if in.reg=1)
%   
%  ADDITIONAL OUTPUT FIELDS OF STRUCTURE ARRAY outex:
%    outex(exp).p      <vector>    estimated parameters
%    outex(exp).errP   <vector>    error of parameters,
%                                  -1 if parameter is fixed;
%    outex(exp).y0     <vector>    estimated initial values
%    outex(exp).errY0  <vector>    error of initial values,
%                                  -1 if initial value is fixed
%    outex(exp).mesh   <matrix>    the actual mesh data

if nargin < 2
  disp('Matlab interface for diffit, v 1.1 06/Sep./2005')
  disp('-----------------------------------------------')
  disp(' ')
  disp('SYNOPSIS')
  disp('  [out outex]= diffit(in,inex)')
  disp(' ')
  disp('DESCRIPTION')
  disp('     in: Structure containing experiment unspecific options.')
  disp('   inex: Structure array containing experiment specific options.')
  disp('    out: Experiment unspecific output structure.')
  disp('  outex: Experiment specific output structure array.')
  disp(' ')
  disp('FIELDS OF STRUCTURE in:')
  disp('  in.eps        <double>   integration accuracy')
  disp('  in.nognu      <0/1>      turns gnuplot animation off') 
  disp('  in.nomeasure  <0/1>      use evently spaced m.s. intervals')
  disp('  in.doP        <string>   specifies parameters to be fitted,')
  disp('                           e.g. 110L')
  disp('  in.maxit      <integer>  maximal number of iterations')
  disp('  in.int        <1/2/3>    0: RK (default), 1: CVODES, 3: ODESSA')
  disp('  in.minimp     <double>   stopping criterion, def.: 1e-4')
  disp('  in.reg        <0/1>      switchs regularisation on')
  disp('  in.epsilon    <double>   singular val. threshold, def.: 1e-10')
  disp('  in.lambda     <double>   reg. parameter, def.: 1e6')
  disp('  in.siminit    <0/1>      simulate initial state')
  disp('  in.pert       <double>   pertubate simulated init. state')
  disp('  in.y0fix      <string>   fixing initial values, e.g. 1100')
  disp('  in.wquer      <double>   controls damping,')
  disp('                           set negative for initialisation')
  disp('  in.nodamp     <0/1>      switches damping off')
  disp('  in.init       <0/1>      1: only initialisation, no fitting')
  disp('  in.silent     <0/1>      1: silent mode')
  disp('  in.opt        <1/2>      optimiser 1: LSEI 2: NAG')
  disp('  in.Lconst     <vector>   set local Parameter Constraints')
  disp(' ')
  disp('FIELDS OF STRUCTURE ARRAY inex:')
  disp('  inex(exp).data    <matrix>    data to be fitted')
  disp('  inex(exp).sig     <matrix>    variance of data')
  disp('  inex(exp).t0      <double>    starting time')
  disp('  inex(exp).t1      <double>    stopping time')
  disp('  inex(exp).nms     <integer>   number of m.s. intervals')
  disp('  inex(exp).p       <vector>    parameters')
  disp('  inex(exp).y0      <vector>    initial values')
  disp('  inex(exp).spline  <cell mat.> spline data')
  disp('  inex(exp).mesh    <matrix>    mesh data')
  disp(' ')
  disp('The contents of in and inex and its default values are copied')
  disp('to the output variables out and outex respectively.')
  disp(' ')
  disp('ADDITIONAL OUTPUT FIELDS OF STRUCTURE ARRAY out:')
  disp('  out.wquer     <double>   actual value of wquer (see above)')
  disp('  out.converged <integer>  1 if fit is converged,0 else')
  disp('  out.chisq     <double>   actual chi^2 value')
  disp('  out.except    <integer>  1 if exception occurred,0 else')
  disp('  out.Lambda    <double>   actual damping parameter')
  disp('  out.cond      <double>   condition number << 1 for')
  disp('                           badly conditioned systems')
  disp('                           (only disp. if in.reg=1)')
  disp(' ')
  disp('ADDITIONAL OUTPUT FIELDS OF STRUCTURE ARRAY outex:')
  disp('  outex(exp).p      <vector>    estimated parameters')
  disp('  outex(exp).errP   <vector>    error of parameters,')
  disp('                                -1 if parameter is fixed');
  disp('  outex(exp).y0     <vector>    estimated initial values')
  disp('  outex(exp).errY0  <vector>    error of initial values,')
  disp('                                -1 if initial value is fixed')
  disp('  outex(exp).mesh   <matrix>    the actual mesh data')
  return
end

% defaults of structure in

if isfield(in,'eps')==0 | isempty(in.eps)
  in.eps=1e-6;
end
if isfield(in,'nognu')==0 | isempty(in.nognu)
  in.nognu=0;
end
if isfield(in,'nomeasure')==0 | isempty(in.nomeasure)
  in.nomeasure=0;
end
if isfield(in,'doP')==0 | isempty(in.doP)
  in.doP='N'; % label for: not defined
end
if isfield(in,'maxit')==0 | isempty(in.maxit)
  in.maxit=10000;
end
if isfield(in,'int')==0 | isempty(in.int)
  in.int=1;
end
if isfield(in,'minimp')==0 | isempty(in.minimp)
  in.minimp=1e-4;
end
if isfield(in,'reg')==0 | isempty(in.reg)
  in.reg=0;
end
if isfield(in,'epsilon')==0 | isempty(in.epsilon)
  in.epsilon=1e-10;
end
if isfield(in,'lambda')==0 | isempty(in.lambda)
  in.lambda=1e6;
end
if isfield(in,'siminit')==0 | isempty(in.siminit)
  in.siminit=0;
end
if isfield(in,'pert')==0 | isempty(in.pert)
  in.pert=0;
end
if isfield(in,'y0fix')==0 | isempty(in.y0fix)
  in.y0fix='N';
end
if isfield(in,'wquer')==0 | isempty(in.wquer)
  in.wquer=-1;
end
if isfield(in,'nodamp')==0 | isempty(in.nodamp)
  in.nodamp=0;
end
if isfield(in,'init')==0 | isempty(in.init)
  in.init=0;
end
if isfield(in,'silent')==0 | isempty(in.silent)
  in.silent=0;
end
if isfield(in,'opt')==0 | isempty(in.opt)
  in.opt=2;
end
if isfield(in,'Lconst')==0 | isempty(in.Lconst)
  in.Lconst=[];
end

% defaults of structure inex
%data
if isfield(inex,'data')==0
  error('Please specify inex(exp).data!');
end
in.nexp=size(inex,2);
for i=1:in.nexp
  inex(i).nobs=size(inex(i).data,2)-1;
  if inex(i).nobs < 1
    error(strcat('No obserables for experiment=',num2str(i)));
  end
  inex(i).n=size(inex(i).data,1);
end
%sig
if isfield(inex,'sig')==0
  for  i=1:in.nexp
    inex(i).sig=ones(inex(i).n,inex(i).nobs);
  end
else
  for  i=1:in.nexp
    if isempty(inex(i).sig)
      inex(i).sig=ones(inex(i).n,inex(i).nobs);
    elseif size(inex(i).sig,1)~=inex(i).n | size(inex(i).sig,2)~= inex(i).nobs
      error(strcat('Incorrect dimension of inex(',num2str(i),').sig.'));
    end
  end
end
%t0
if isfield(inex,'t0')==0
  for i=1:in.nexp
    inex(i).t0=inex(i).data(1,1);
  end
else
  for i=1:in.nexp
    if isempty(inex(i).t0)
      inex(i).t0=inex(i).data(1,1);
    end
  end
end
%t1
if isfield(inex,'t1')==0
  for i=1:in.nexp
    inex(i).t1=inex(i).data(inex(i).n,1);
  end
else
  for i=1:in.nexp
    if isempty(inex(i).t1)
      inex(i).t1=inex(i).data(inex(i).n,1);
    end
  end
end
%p
if isfield(inex,'p')==0
  for i=1:in.nexp
    inex(i).p=[];
  end
else
  for i=1:in.nexp
    if isempty(inex(i).p)
      inex(i).p=[];
    end
  end
end
%y0
if isfield(inex,'y0')==0
  for i=1:in.nexp
    inex(i).y0=[];
  end
else
  for i=1:in.nexp
    if isempty(inex(i).y0)
      inex(i).y0=[];
    end
  end
end
%spline
if isfield(inex,'spline')==0
  for i=1:in.nexp
    inex(i).spline= {};
  end
else
  for i=1:in.nexp
    if isempty(inex(i).spline)
      inex(i).spline={};
    end
  end
end
%mesh
if isfield(inex,'mesh')==0
  for i=1:in.nexp
    inex(i).mesh=[];
  end
else
  for i=1:in.nexp
    if isempty(inex(i).mesh)
      inex(i).mesh=[];
    end
  end
end
%nms
if isfield(inex,'nms')==0
  for  i=1:in.nexp
    inex(i).nms=1;
  end
else
  for  i=1:in.nexp
    if isempty(inex(i).nms)
      inex(i).nms=1;
    end
  end
end
out=in;
outex=inex;

[o oex]=diffit_mex(in,inex);

if o.except==1
  out.except=1;
  out.converged=0;
else
  %OUTPUT -> out
  out.wquer=o.wquer;
  out.converged=o.converged;
  out.chisq=o.chisq;
  out.except=o.except;
  out.Lambda=o.Lambda;
  if in.reg==1
    out.cond=o.cond;
  end
  %OUTPUT -> outex
  for i=1:in.nexp
    outex(i).p=oex(i).p;
    outex(i).errP=oex(i).errP;
    outex(i).y0=oex(i).y0;
    outex(i).errY0=oex(i).errY0;
    outex(i).mesh=oex(i).mesh;
  end
end

