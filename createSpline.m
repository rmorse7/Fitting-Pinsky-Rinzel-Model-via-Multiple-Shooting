function varargout = createSpline(varargin)
% createSpline wrapper

if nargin < 2
  disp('Usage: splinedat =createSpline(''i'', <input data matrix>)')
  disp('options: ''alpha'',  <double>')
  disp('         ''max '', <double>')
  disp('         ''gcv'', 0/1')
  return
end

% defaults
alpha=1;
maxa=1e10;
gcv=1;

for i=1:2:nargin
  switch varargin{i},
   case 'alpha' 
    alpha = varargin{i+1};
   case 'i' 
    inp = varargin{i+1};
   case 'max' 
    maxa = varargin{i+1};
   case 'gcv' 
    gcv = varargin{i+1};
  end
end

if isempty(inp)
  error('No data!')
end


% createSpline_mex aufruf

[splinedat,alpha_out] = createSpline_mex(inp,alpha,maxa,gcv);
varargout{1} = splinedat;
varargout{2} = alpha_out;