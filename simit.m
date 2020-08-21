function varargout = simit(in)
%Matlab interface for simit, v 1.0 06/Sep./2005
%-----------------------------------------------
% 
%SYNOPSIS;
%  [yout (dyds dydp)]= simit(in)
% 
%  DESCRIPTION
%       in: Structure containing the input parameters
%   
%  FIELDS OF STRUCTURE in: 
%    in.t       <vector>      time points      
%    in.spline  <cell mat.>   spline data
%    in.t0      <double>      first time point
%    in.y0      <vector>      initial values
%    in.p       <vector>      parameters
%    in.eps     <double>      integration accuracy
%    in.int     <integer>     1: Runge Kutta
%                             2: CVODES (default)
%                             3: ODESSA
%    in.hidden  <integer>     1: calculate hidden states
%   
%  OUTPUT (empty if exception occured)
%    yout  <matrix>   trajectory
%    dyds  <3tensor>  sensitvities resp. initial val.
%    dydp  <3tensor>  sensitvities resp. parameters

if nargin < 1
  help simit
  return
end

% defaults of structure in

if isfield(in,'t')==0 | isempty(in.t)
  error('Please specify in.t!');
end
if isfield(in,'spline')==0 | isempty(in.spline)
  in.spline ={};
end
if isfield(in,'t0')==0 | isempty(in.t0)
  in.t0=0;
end
if isfield(in,'y0')==0 | isempty(in.y0)
  in.y0=[];
end
if isfield(in,'p')==0 | isempty(in.p)
  in.p=[];
end
if isfield(in,'eps')==0 | isempty(in.eps)
  in.eps= 10e-6;
end
if isfield(in,'int')==0 | isempty(in.int)
  in.int=2;
end
if isfield(in,'hidden')==0 | isempty(in.hidden)
  in.hidden=0;
end

% simit_mexaufruf
if nargout > 1
  [N yout dydy0 dydp] = simit_mex(in);
  varargout{1} = yout;
  % reshape
  
  NEQNS=N(1);
  NOBS=N(2);
  NPAR=N(3);

  dydy0 = reshape(dydy0,NOBS,NEQNS,length(in.t));
  dydp = reshape(dydp,NOBS,NPAR,length(in.t));
  varargout{2} = dydy0;
  varargout{3} = dydp;
else
  [N yout] = simit_mex(in);
  varargout{1} = yout;
end
