function [] = muscleForce_postAnalysis(auxdata,output)
close all
solution = output.result.solution;

scale = auxdata.scale;
spline = auxdata.spline;

muscParams = auxdata.muscParams;
cCoefs = auxdata.FaCoefs;
dCoefs = auxdata.FvCoefs;
eCoefs = auxdata.FpCoefs;

tact = auxdata.tact;
tdeact = auxdata.tdeact;
kT = auxdata.kT;

% Unscale states, time, and controls
t = solution.phase.time/scale.time;
a = solution.phase.state(:,1:5);
lMtilda = solution.phase.state(:,6:10);
e = solution.phase.control(:,1:5);
vMtilda = solution.phase.control(:,6:10);
resAct = solution.phase.control(:,11:12)/scale.angAccel;
npts = length(t);

lMo = repmat(muscParams(2,:),npts,1);
vMmax = repmat(muscParams(5,:),npts,1);
lMtildad = vMtilda.*vMmax./lMo;

% Resample persistent variables if time vector has changed
TA_ref = ppval(spline.TA_ref,t);
TK_ref = ppval(spline.TK_ref,t);

R1A = ppval(spline.R1A,t);
R2A = ppval(spline.R2A,t);
R5A = ppval(spline.R5A,t);
R3K = ppval(spline.R3K,t);
R4K = ppval(spline.R4K,t);
R5K = ppval(spline.R5K,t);

lMT = ppval(spline.lMT_ref,t)';


% Run activation dynamics and calculate muscle forces
ad = activationDynamics(e,a,tact,tdeact);
[Ferror,FT] = forceEquilibrium(a,lMtilda,vMtilda,lMT,muscParams,cCoefs,dCoefs,eCoefs,kT);


TA = FT(:,1).*R1A - FT(:,2).*R2A + FT(:,5).*R5A;
TK = -FT(:,3).*R3K + FT(:,4).*R4K - FT(:,5).*R5K;

torqueError = [TA - TA_ref,TK - TK_ref];

output.dynamics = [ad/scale.time,lMtildad/scale.time];
output.integrand = (sum(e.^2,2) + 1e5*sum((resAct*1e-3).^2,2))*scale.cost;
output.path = [torqueError*scale.torque,Ferror*scale.force];


load('referenceData.mat');
Q_ref = state_ref(:,1:2);
Qd_ref = state_ref(:,3:4);
a_ref = state_ref(:,5:9);


figure(1);
e1 = plot(t,e(:,1),'r*',t_ref,e_ref(:,1),'r-');
hold all;
e2 = plot(t,e(:,2),'g*',t_ref,e_ref(:,2),'g-');
e3 = plot(t,e(:,3),'b*',t_ref,e_ref(:,3),'b-');
e4 = plot(t,e(:,4),'k*',t_ref,e_ref(:,4),'k-');
e5 = plot(t,e(:,5),'m*',t_ref,e_ref(:,5),'m-');
xlabel('Time (s)');
ylabel('Excitations');
legend([e4],'output','reference');

figure(2);
a1 = plot(t,a(:,1),'r*',t_ref,a_ref(:,1),'r-');
hold all;
a2 = plot(t,a(:,2),'g*',t_ref,a_ref(:,2),'g-');
a3 = plot(t,a(:,3),'b*',t_ref,a_ref(:,3),'b-');
a4 = plot(t,a(:,4),'k*',t_ref,a_ref(:,4),'k-');
a5 = plot(t,a(:,5),'m*',t_ref,a_ref(:,5),'m-');
xlabel('Time (s)');
ylabel('Activations');
legend([a4],'output','reference');

figure(3)
plot(t,TA,'r*',t,TA_ref,'r-');
hold all;
plot(t,TK,'b*',t,TK_ref,'b-');
xlabel('Time (s)');
ylabel('Joint Torques');
legend('TA output','TA ref','TK output','TK ref');

figure(4)
F1 = plot(t,FT(:,1),'r*',t_ref,F_ref(:,1),'r-');
hold all;
F2 = plot(t,FT(:,2),'g*',t_ref,F_ref(:,2),'g-');
F3 = plot(t,FT(:,3),'b*',t_ref,F_ref(:,3),'b-');
F4 = plot(t,FT(:,4),'k*',t_ref,F_ref(:,4),'k-');
F5 = plot(t,FT(:,5),'m*',t_ref,F_ref(:,5),'m-');
xlabel('Time (s)');
ylabel('Muscle Forces');
legend([F4],'output','reference');

figure(5)
plot(t,resAct);
xlabel('Time (s)');
ylabel('Residual Torques');

end