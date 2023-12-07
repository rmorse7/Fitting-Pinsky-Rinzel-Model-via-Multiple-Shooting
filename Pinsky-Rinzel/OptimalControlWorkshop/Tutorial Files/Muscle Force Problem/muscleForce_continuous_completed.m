function output = muscleForce_continuous_completed(input)
%% Unpack auxdata
auxdata = input.auxdata;
scale = auxdata.scale;
spline = auxdata.spline;

tact = auxdata.tact;
tdeact = auxdata.tdeact;
muscParams = auxdata.muscParams;
FaCoefs = auxdata.FaCoefs;
FvCoefs = auxdata.FvCoefs;
FpCoefs = auxdata.FpCoefs;
kT = auxdata.kT;


%% Unscale states, time, and controls
t = input.phase.time/scale.time;
a = input.phase.state(:,1:5);
lMtilda = input.phase.state(:,6:10);
e = input.phase.control(:,1:5);
vMtilda = input.phase.control(:,6:10);
resAct = input.phase.control(:,11:12)/scale.torque;


%% Sample Splines
R1A = ppval(spline.R1A,t);
R2A = ppval(spline.R2A,t);
R5A = ppval(spline.R5A,t);
R3K = ppval(spline.R3K,t);
R4K = ppval(spline.R4K,t);
R5K = ppval(spline.R5K,t);

TA_ref = ppval(spline.TA_ref,t);
TK_ref = ppval(spline.TK_ref,t);

lMT = ppval(spline.lMT_ref,t)';


%% Activation Dynamics
ad = activationDynamics(e,a,tact,tdeact);


%% Contraction Dynamics
[Ferror,FT] = forceEquilibrium(a,lMtilda,vMtilda,lMT,muscParams,FaCoefs,FvCoefs,FpCoefs,kT);


%% Calculate Torques
TA = resAct(:,1) + FT(:,1).*R1A - FT(:,2).*R2A + FT(:,5).*R5A;
TK = resAct(:,2) + -FT(:,3).*R3K + FT(:,4).*R4K - FT(:,5).*R5K;

torqueError = [TA - TA_ref,TK - TK_ref];


%% Output
npts = length(t);
lMo = repmat(muscParams(2,:),npts,1);
vMmax = repmat(muscParams(5,:),npts,1);
lMtildad = vMtilda.*vMmax./lMo;

output.dynamics = [ad/scale.time,lMtildad/scale.time];
output.integrand = (sum(e.^2,2) + 1e5*sum((resAct*1e-3).^2,2))*scale.cost;
output.path = [torqueError*scale.torque,Ferror*scale.force];

end