%
%  Gabbiani
%

%
%  mlsolve.m, modified from prsolve.m, which solves the Pinsky-Rinzel model
%
%  solves the Morris-Lecar single compartment model
%
%  the state is y = [u w]
%  with: 
% 
%  u = membrane potential, w = recovery variable.
%
%  usage   ml_solve(T,I)     e.g.,  ml_solve(80,42)
%     or   [t,y] = ml_solve(T,I)
%
%  where   T = duration of simulation (ms)
%          I = current injection (uA/cm^2)
%          
%  Alternatively, ml_solve(T,I,p) or [t,y] = ml_solve(T,I,p), where
%  
%  p.gCa, p.gK, p.gL, p.VCa, p.VK, p.VL, p.Cm, p.v1, p.v2, p.v3, p.v4, p.tau_w_max  
%
%  are the model parameters.

%  figure produces V vs time
%
% References: NIPS manuscript and 
% http://www.math.pitt.edu/~bard/bardware/meth3/node18.html
%

function [t,y] = ml_solve(T,I,varargin)

if ( isempty(varargin) )
    %maximal conductances in mS/cm^2
    p.gCa = 4.0; %calcium current
    p.gK = 8.0; %potassium current
    p.gL = 2.0; %leak
    
    %reversal potentials (in mV)
    p.VCa = 120;
    p.VK = -84;
    p.VL = -60;
    
    %capacitance uF/cm^2
    p.Cm = 1;
    
    %activation function parameters for Ca and K currents
    p.v1 = -1.2; %Ca
    p.v2 = 18.0; %Ca
    p.v3 = 12.0; %K
    p.v4 = 17.4; %K
    p.tau_w_max = 15; %ms 
else
    p = varargin{1,1};
end

%initial value of the state variables
y0 = [-30.0 0.0];

dt = 0.1;
tspan = 0:dt:T;

%  ode23 used here is a Runge-Kutta (2,3) algorithm with
%  3rd order accuracy and an embedded 2nd order method allowing
%  for variable step size. See Num Recipes in C, chapt 16, sec 2.
[t,y] = ode23(@red,tspan,y0,[],I,p);

if ( nargout == 0 )
    figure; 
    subplot(2,1,1);
    plot(t,y(:,1))
    legend('V')
    xlabel('time (ms)','fontsize',12);
    ylabel('Membrane Potential  (mV)','fontsize',12);
    subplot(2,1,2);
    plot(t,y(:,2));
    legend('w');
    xlabel('time (ms)','fontsize',12);
    ylabel('Recovery variable  (dimensionless)','fontsize',12);
end
    
return

%
% the reduced traub model of Pinsky & Rinzel
%

function dy = red(t,y,I,p)

dy = zeros(2,1);

%leak current
IL = p.gL*(y(1)-p.VL);

%steady-state calcium activation (instantaneous)
minf = m_infty(y(1),p.v1,p.v2);

%sodium current (y(3) is h, inactivation of sodium current)
ICa = p.gCa*minf.*(y(1)-p.VCa);

%potassium current (y(2) is w, activation of K)
IK = p.gK*y(2)*(y(1)-p.VK);

%derivative update of membrane potential 
dy(1) = (-IL - ICa - IK + I)/p.Cm;

%steady-state K activation
winf = w_infty(y(1),p.v3,p.v4);

%K time constant
tw = tau_w(y(1),p.v3,p.v4,p.tau_w_max);

%derivative update of w
dy(2) = (winf-y(2))/tw;

return

%
%For following rate constants, see eq. 6 of [PR94] and erratum
%

%steady-state value of gCa
function val = m_infty(v,v1,v2)
val = 0.5*( 1 + tanh( (v-v1)/v2 ) );

%steady-state value of gK
function val = w_infty(v,v3,v4)
val = 0.5*( 1 + tanh( (v-v3)/v4 ) );

%time constant for gK
function val = tau_w(v,v3,v4,tau_w_max)
val = tau_w_max/( cosh((v-v3)/(2*v4) ) );
