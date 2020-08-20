%  Adapted from prsolve, but uses the 4th order classic Runge-Kutta method
%
%  prsolve_rk4.m
%
%  solve the pinsky-rinzel 2 compartment Ca3 model
%
%  the state is y = [Vs Vd h n s c q Ca]
%  with: 
% 
%  Vs = soma potential, Vd = dendrite potential,
%  h = fast sodium inactivation, n = delayed rectifier activation,
%  s = Ca current activation, c = fast Ca and voltage dependent K current
%  activation, q = slow Ca dependent K current activation
%size(
%  usage   prsolve(T,[Is Id gc])     e.g.,  prsolve(50,[1 0 2.1])
%
%  where   T = duration of simulation (ms)
%          Is = somatic current injection (microA/cm^2)
%          Id = dendritic current injection (microA/cm^2)
%          gc = coupling between the two compartments

%  figure produces Vs and Vd vs time
%
% Reference: [PR94] Pinsky PF, Rinzel J (1994) Intrinsic and Network
% Rhythmogenesis in a Reduced Traub Model for CA3 Neurons. J Comp Neurosci
% 1:39-60. Erratum in J Comp Neurosci 2:275, 1995. 
%


function [tspan, y] = prsolve_rk4_sigmoid(T,params,gs)

%initial value of the state variables
y0 = [-4.6 -4.5 0.999 0.001 0.009 0.007 0.01 0.2];
n_y0 = length(y0);

Is = params(1);
Id = params(2);
gc = params(3);

dt = 0.05;
tspan = 0:dt:T;
n_tspan = length(tspan);

y = zeros(n_tspan,n_y0);
y(1,:) = y0;
yc = y0;

%[t,y] = ode23(@red,tspan,y0,[],Is,Id,gc);

dt2 = dt/2;
for i = 2:n_tspan
    dy1 = dt2*red(tspan(i-1),yc,Is,Id,gc,gs);
    dy2 = dt2*red(tspan(i-1)+dt2,yc+dy1,Is,Id,gc,gs);
    dy3 = dt2*red(tspan(i-1)+dt2,yc+dy2,Is,Id,gc,gs);
    dy4 = dt2*red(tspan(i-1)+dt,yc+2*dy3,Is,Id,gc,gs);
    
    yc = yc + (1/3)*(dy1 + 2*dy2 + 2*dy3 + dy4);
    y(i,:) = yc;
end

% plot(tspan,y(:,1),tspan,y(:,2))
% legend('V_s','V_d')
% xlabel('time  (ms)','fontsize',16)
% ylabel('Membrane Potential  (mV)','fontsize',16)

return

%
% the reduced traub model of Pinsky & Rinzel
%

function dy = red(t,y,Is,Id,gc,gs)

%maximal conductances in mS/cm^2
gL = gs(1); %leak
gNa = gs(2); %fast sodium
gKDR = gs(3); %delayed rectifier
gCa = gs(4); %fast calcium current
gKAHP = gs(5); %calcium dependent potassium current (slow)
gKC = gs(6); %voltage and calcium dependent potassium current (fast)

%reversal potentials (in mV)
VNa = 120; VCa = 140; VK = -15; VL = 0;

% Is = -0.5; passed as parameter
Isyn = 0; %used in [PR94] but not here
%Id = 0.5;

%conductance coupling coefficient between soma and dendrite compartment (in
%mS/cm^2)
%gc = 2.1; passed as parameter
%gc = 0;

%fraction of cable length assigned to soma (1-p for dendrite)
p = 0.5;

%capacitance microF/cm^2
Cm = 3;

dy = zeros(1,8);

%somatic leak current
Ils = gL*(y(1)-VL);

%steady-state sodium activation (instantaneous)
minf = am(y(1))/(am(y(1))+bm(y(1)));

%sodium current (y(3) is h, inactivation of sodium current)
INa = gNa*minf.^2.*y(3).*(y(1)-VNa);

%delayed rectifier current (y(4) is n, activation of DR)
IKDR = gKDR*y(4)*(y(1)-VK);

%derivative update of somatic membrane potential, eq. 1 of [PR94] 
dy(1) = (-Ils - INa - IKDR + (gc/p)*(y(2)-y(1)) + Is/p)/Cm;

%dendritic leak current
Ild = gL*(y(2)-VL);

%dendritic calcium current (y(5) is s activation variable)
ICa = gCa*y(5).^2.*(y(2)-VCa);

%voltage and calcium dependent K current (y(6) is c activation variable,
%y(8) is Ca)
%IKC = gKC*y(6)*min(y(8)/250,1)*(y(2)-VK);
IKC = gKC*y(6)*( (1-1/(1+exp(-(y(8)-250))))*y(8)/250.0+(1/(1+exp(-(y(8)-250)))) )*(y(2)-VK);

%calcium dependent K current (y(7) is q activation variable)
IKAHP = gKAHP*y(7).*(y(2)-VK);

%derivative update of dendritic membrane potential eq. 1 of [PR94]
dy(2) = (-Ild-ICa-IKAHP-IKC-Isyn/(1-p)+(gc/(1-p))*(y(1)-y(2))+Id/(1-p))/Cm;

%derivative update of h, sodium inactivation
dy(3) = ah(y(1))*(1-y(3)) - bh(y(1))*y(3);

%derivative update of n, DR activation
dy(4) = an(y(1))*(1-y(4)) - bn(y(1))*y(4);

%derivative update of s, Ca activation
dy(5) = as(y(2))*(1-y(5)) - bs(y(2))*y(5);

%derivative update of c, fast Ca, V dependent K current
dy(6) = ac(y(2))*(1-y(6)) - bc(y(2))*y(6);

%derivative update of q, slow Ca dependent K current
dy(7) = aq(y(8))*(1-y(7)) - bq(y(8))*y(7);

%derivative update of calcium, first term is calcium influx through Ca
%channels, second term is exponential decay eq. 3 of [PR94]
dy(8) = -0.13*ICa - 0.075*y(8);

return

%
%For following rate constants, see eq. 6 of [PR94] and erratum
%

%forward rate constant for fast sodium
function val = am(v)
val = 0.32*(13.1-v)/(exp((13.1-v)/4)-1);

%backward rate constant for fast sodium
function val = bm(v)
val = 0.28*(v-40.1)/(exp((v-40.1)/5)-1);

%forward rate constant for DR activation
function val = an(v)
val = 0.016*(35.1-v)/(exp((35.1-v)/5)-1);

%backward rate constant for DR activation
function val = bn(v)
val = 0.25*exp(0.5-0.025*v);

%forward rate constant for sodium inactivation
function val = ah(v)
val = 0.128*exp((17-v)/18);

%backward rate constant for sodium inactivation
function val = bh(v)
val = 4./(1+exp((40-v)/5));

%forward rate constant for Ca current activation
function val = as(v)
val = 1.6./(1+exp(-0.072*(v-65)));

%backward rate constant for Ca current activation
function val = bs(v)
val = 0.02*(v-51.1)./(exp((v-51.1)/5)-1);

%forward rate constant for fast Ca and V dependent K current
function val = ac(v)
%if v <= 50
%   %see erratum 
%   val = exp((v-10)/11-(v-6.5)/27)/18.975;
%else
%   val = 2*exp((6.5-v)/27);
%end
val=((1-1/(1+exp(-(v-50))))*((exp((v-10)/11-(v-6.5)/27))/18.975)+(1/(1+exp(-(v-50))))*(2*exp((6.5-v)/27)));

%backward rate constant for fast Ca and V dependent K current
function val = bc(v)
%if v <= 50
%   val = 2*exp((6.5-v)/27) - exp((v-10)/11-(v-6.5)/27)/18.975;
%else
%   val = 0;
%end
val=(1-1/(1+exp(-(v-50))))*(2*exp((6.5-v)/27)-(exp((v-10)/11-(v-6.5)/27))/18.975);

%forward rate constant for slow Ca dependent K current
function val = aq(v)
%val = min((0.00002)*v,0.01);
val=((1-1/(1+exp(-(0.00002*v-0.01)/0.00002)))*(0.00002*v)+(1/(1+exp(-(0.00002*v-0.01)/0.00002)))*0.01);

%backward rate constant for slow Ca dependent K current
function val = bq(v)
val = 0.001;
