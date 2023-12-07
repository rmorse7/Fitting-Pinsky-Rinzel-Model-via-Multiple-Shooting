%-------------------------------------------------------------------------%
% Optimal squatting problem
%-------------------------------------------------------------------------%
clear all; close all; clc

% Scaling
scale.time = 30;
scale.mass = 1;
scale.length = 1;
scale.vel = scale.length/scale.time;
scale.accel = scale.vel/scale.time;
scale.force = scale.mass*scale.accel;
scale.torque = scale.force*scale.length;
scale.angle = 1;
scale.angVel = scale.angle/scale.time;
scale.angAccel = scale.angVel/scale.time;
scale.inertia = scale.torque/scale.angAccel;

scale.cost = 1e1/scale.time;
auxdata.scale = scale;


% Intial and intermediate joint angles
initialQ1 = 0*scale.angle;
initialQ2 = 0*scale.angle;
Q1via = 30*scale.angle; % target angle (deg)
Q2via = 60*scale.angle;

% Constant parameters
auxdata.TA =  0.0;                  % Nm
auxdata.TK =  0.0;                  % Nm
auxdata.g =  9.81;                  % m/s^2
auxdata.r =  0.1;                   % m
auxdata.IC =  0.065;                % Kg m^2
auxdata.ID =  0.126;                % Kg m^2
auxdata.LC =  0.435;                % m
auxdata.LD =  0.4;                  % m
auxdata.mC =  7.5;                  % Kg
auxdata.mD =  15.15;                % Kg
auxdata.mH =  51.22;                % Kg
auxdata.rhoC =  0.2;                % m
auxdata.rhoD =  0.2;                % m
auxdata.IB = 0.008;                 % kg m^2
auxdata.mB = 2.2;                   % kg
auxdata.HB = 0.1;                   % m
auxdata.rhoB = 0.1;                 % m
auxdata.LB = 0.15;                   % m

auxdata.tact = 0.01;
auxdata.tdeact = 0.04;
auxdata.muscParams = [[18000,3000,8700,18900,7700];
    [0.128652740182846,0.140039645883349,0.169934531781695,0.188994102302640,0.353638533815950];
    [0.223606797749979,0.173205080756888,0.140121893133717,0.239387238590531,0.567940357784160]/2;
    pi/180*[28.3,9.6,11.6,29.6,9.9];
    [1.28652740182846,1.40039645883349,1.69934531781695,1.88994102302640,3.53638533815950]];
auxdata.cCoefs = 6.52843814543065;
auxdata.dCoefs = [0.705320356648889;6.50627486128284];
auxdata.eCoefs = [2.33892909857083;1.51871677360826;-2.50603008180404];


% Set minimum and maximum values for states and controls
xmin = [0, 0]*scale.angle;
xmax  = pi/180*[50,80]*scale.angle;
vmin = 1/1*[-10, -10]*scale.angVel;
vmax  = 1/1*[10,10]*scale.angVel;
amin = zeros(1,5);
amax = ones(1,5);

emin = 0*ones(1,5);
emax = 1*ones(1,5);

% Define initial, phase 1 final, and phase 2 final times
t0 = 0*scale.time; 
tf1 = 0.5*scale.time;
tf2 = 1.0*scale.time;

% initial states for simulation
x0 = pi/180*[initialQ1, initialQ2];
v0 = [0,0]*scale.angVel;
a0 = zeros(1,5);


% final states for phase 1
xfmin1 = pi/180*([Q1via, Q2via]-0.1*scale.angle); % convert target angle (deg) to radians.
xfmax1 = pi/180*([Q1via, Q2via]+0.1*scale.angle);
vfmin1 = ([0, 0] - pi/180)*scale.angVel;
vfmax1 = ([0, 0] + pi/180)*scale.angVel ;
afmin1 = zeros(1,5);
afmax1 = ones(1,5);


% final states for phase 2
xfmin2 = pi/180*([initialQ1, initialQ2]); % convert target angle (deg) to radians.
xfmax2 = pi/180*([initialQ1, initialQ2]);
vfmin2 = v0;
vfmax2 = v0;
afmin2 = zeros(1,5);
afmax2 = ones(1,5);



%-------------------------------------------------------------------------%
% Provide Bounds and Guess in Each Phase of Problem
%-------------------------------------------------------------------------%
% Phase 1
iphase = 1;
bounds.phase(iphase).initialtime.lower = t0; 
bounds.phase(iphase).initialtime.upper = t0;
bounds.phase(iphase).finaltime.lower = tf1; 
bounds.phase(iphase).finaltime.upper = tf1;
bounds.phase(iphase).initialstate.lower = [x0,v0,amin];
bounds.phase(iphase).initialstate.upper = [x0,v0,amax];
bounds.phase(iphase).state.lower = [xmin,vmin,amin];
bounds.phase(iphase).state.upper = [xmax,vmax,amax];
bounds.phase(iphase).finalstate.lower = [xfmin1,vfmin1,afmin1];
bounds.phase(iphase).finalstate.upper = [xfmax1,vfmax1,afmax1];
bounds.phase(iphase).control.lower = emin;
bounds.phase(iphase).control.upper = emax;
bounds.phase(iphase).integral.lower = 0; 
bounds.phase(iphase).integral.upper = 100*scale.cost;


guess.phase(iphase).time    = [t0; tf1]; 
guess.phase(iphase).state   = [x0,v0,a0; xfmin1,vfmin1,afmin1];
guess.phase(iphase).control = zeros(2,5);
guess.phase(iphase).integral = 0.1;

% Phase 2
iphase = 2;
bounds.phase(iphase).initialtime.lower = tf1; 
bounds.phase(iphase).initialtime.upper = tf1;
bounds.phase(iphase).finaltime.lower = tf2; 
bounds.phase(iphase).finaltime.upper = tf2;
bounds.phase(iphase).initialstate.lower = [xfmin1,vfmin1,afmin1];
bounds.phase(iphase).initialstate.upper = [xfmax1,vfmax1,afmax1];
bounds.phase(iphase).state.lower = [xmin,vmin,amin];
bounds.phase(iphase).state.upper = [xmax,vmax,amax];
bounds.phase(iphase).finalstate.lower = [xfmin2,vfmin2,afmin2];
bounds.phase(iphase).finalstate.upper = [xfmax2,vfmax2,afmax2];
bounds.phase(iphase).control.lower = emin;
bounds.phase(iphase).control.upper = emax;
bounds.phase(iphase).integral.lower = 0; 
bounds.phase(iphase).integral.upper = 100*scale.cost;


guess.phase(iphase).time    = [tf1; tf2]; 
guess.phase(iphase).state   = [xfmin1,vfmin1,afmin1; xfmin2,vfmin2,afmin2];
guess.phase(iphase).control = zeros(2,5);
guess.phase(iphase).integral = 0.1;

% Event constraints
bounds.eventgroup.lower = zeros(1,10);
bounds.eventgroup.upper = zeros(1,10);

%-------------------------------------------------------------------------%
%-------------- Provide an Initial Mesh in Each Phase --------------------%
%-------------------------------------------------------------------------%
for i=1:2
meshphase(i).colpoints = 4*ones(1,20);
meshphase(i).fraction = 1/20*ones(1,20);
end

%-------------------------------------------------------------------------%
%------------- Assemble Information into Problem Structure ---------------%        
%-------------------------------------------------------------------------%
setup.name = 'Squat';
setup.functions.continuous = @squat_continuous;
setup.functions.endpoint = @squat_endpoint;
setup.auxdata = auxdata;
setup.bounds = bounds;
setup.guess = guess;

setup.nlp.solver = 'ipopt';
setup.nlp.ipoptoptions.tolerance = 1e-7;
setup.derivatives.supplier = 'sparseCD';
setup.derivatives.derivativelevel = 'first';
setup.derivatives.dependencies = 'sparse';

setup.mesh.phase = meshphase;
setup.mesh.method = 'hp-LiuRao';
setup.mesh.tolerance = 1e-5;
setup.mesh.maxiterations = 10;
setup.mesh.colpointsmin = 4;
setup.mesh.colpointsmax = 10;

setup.method = 'RPM-Integration';


%-------------------------------------------------------------------------%
% Solve Problem Using GPOP2
%-------------------------------------------------------------------------%
tic
output = gpops2(setup);
solution = output.result.solution;
toc

save('squat_output','output');


%-------------------------------------------------------------------------%
% Plot Solution
%-------------------------------------------------------------------------%

TA = auxdata.TA;
TK = auxdata.TK;
g = auxdata.g;
r = auxdata.r;
IC = auxdata.IC;
ID = auxdata.ID;
LC = auxdata.LC;
LD = auxdata.LD;
mC = auxdata.mC;
mD = auxdata.mD;
mH = auxdata.mH;
rhoC = auxdata.rhoC;
rhoD = auxdata.rhoD;


tact = auxdata.tact;
tdeact = auxdata.tdeact;
muscParams = auxdata.muscParams;
cCoefs = auxdata.cCoefs;
dCoefs = auxdata.dCoefs;
eCoefs = auxdata.eCoefs;


scale = auxdata.scale;


%-------------------------------------------------------------------------%
% Phase 1
%-------------------------------------------------------------------------%
t1 = solution.phase(1).time;
x1 = solution.phase(1).state(:,1:4);
a1 = solution.phase(1).state(:,5:9);
e1 = solution.phase(1).control(:,1:5);
numPts1 = size(x1,1);

Q1 = x1(:,1)/scale.angle;
Q2 = x1(:,2)/scale.angle;
Q1d = x1(:,3)/scale.angVel;
Q2d = x1(:,4)/scale.angVel;


c1 = cos(Q1);
s1 = sin(Q1);
c2 = cos(Q2);
s2 = sin(Q2);
c1_2 = cos((Q1-Q2));
s1_2 = sin((Q1-Q2));


lMT1(:,1) = (r^2+rhoC^2+2*r*rhoC*s1).^0.5;
lMT1(:,2) = (r^2+rhoC^2-2*r*rhoC*s1).^0.5;
lMT1(:,3) = (r^2+rhoD^2-2*r*rhoD*s2).^0.5;
lMT1(:,4) = (r^2+rhoD^2+2*r*rhoD*s2).^0.5;
lMT1(:,5) = (r^2+LC^2+rhoD^2+2*r*LC*s1+2*LC*rhoD*c2+2*r*rhoD*s1_2).^0.5;
vMT1(:,1) = r*rhoC*c1.*Q1d./(r^2+rhoC^2+2*r*rhoC*s1).^0.5;
vMT1(:,2) = -r*rhoC*c1.*Q1d./(r^2+rhoC^2-2*r*rhoC*s1).^0.5;
vMT1(:,3) = -r*rhoD*c2.*Q2d./(r^2+rhoD^2-2*r*rhoD*s2).^0.5;
vMT1(:,4) = r*rhoD*c2.*Q2d./(r^2+rhoD^2+2*r*rhoD*s2).^0.5;
vMT1(:,5) = -(LC*rhoD*s2.*Q2d-r*LC*c1.*Q1d-r*rhoD*c1_2.*(Q1d-Q2d))./(r^2+LC^2+rhoD^2+2*r*LC*s1+2*LC*rhoD*c2+2*r*rhoD*s1_2).^0.5;


a1d = activationDynamics(e1,a1,tact,tdeact);
[F1,lMtilda1,FMltilda1,vMtilda1,FMvtilda1] = rigidTendonModel_simple(a1(:,1),lMT1(:,1),vMT1(:,1),muscParams(:,1),cCoefs,dCoefs,eCoefs);
[F2,lMtilda2,FMltilda2,vMtilda2,FMvtilda2] = rigidTendonModel_simple(a1(:,2),lMT1(:,2),vMT1(:,2),muscParams(:,2),cCoefs,dCoefs,eCoefs);
[F3,lMtilda3,FMltilda3,vMtilda3,FMvtilda3] = rigidTendonModel_simple(a1(:,3),lMT1(:,3),vMT1(:,3),muscParams(:,3),cCoefs,dCoefs,eCoefs);
[F4,lMtilda4,FMltilda4,vMtilda4,FMvtilda4] = rigidTendonModel_simple(a1(:,4),lMT1(:,4),vMT1(:,4),muscParams(:,4),cCoefs,dCoefs,eCoefs);
[F5,lMtilda5,FMltilda5,vMtilda5,FMvtilda5] = rigidTendonModel_simple(a1(:,5),lMT1(:,5),vMT1(:,5),muscParams(:,5),cCoefs,dCoefs,eCoefs);

    
Q1d2 = ((ID+LD*mH*(LD+LC*c2)+mD*rhoD*(rhoD+LC*c2)).*(TK+g*LD*mH*s1_2 ... 
    +rhoD*(g*mD*s1_2-r*c2.*(F3./(r^2+rhoD^2-2*r*rhoD*s2).^0.5-F4./(r^2+rhoD^2+2*r*rhoD*s2).^0.5) ... 
    -F5.*(LC*s2+r*c1_2)./(r^2+LC^2+rhoD^2+2*r*LC*s1+2*LC*rhoD*c2 ... 
    +2*r*rhoD*s1_2).^0.5)+LC*(LD*mH+mD*rhoD)*s2.*Q1d.^2)+(ID+mD*rhoD^2+mH*LD^2)*(TA+r*F5.*(LC*c1 ... 
    +rhoD*c1_2)./(r^2+LC^2+rhoD^2+2*r*LC*s1+2*LC*rhoD*c2+2*r*rhoD*s1_2).^0.5 ... 
    -g*LC*mD*s1-g*mD*rhoD*s1_2-g*mH*(LC*s1+LD*s1_2)-rhoC*(g*mC*s1 ... 
    -r*c1.*(F1./(r^2+rhoC^2+2*r*rhoC*s1).^0.5-F2./(r^2+rhoC^2-2*r*rhoC*s1).^0.5)) ... 
    -LC*(LD*mH+mD*rhoD)*s2.*(Q1d.^2-(Q1d-Q2d).^2)))./((ID+LD*mH*(LD+LC*c2)+mD*rhoD*(rhoD+LC*c2)).^2 ... 
    -(ID+mD*rhoD^2+mH*LD^2)*(IC+ID+mC*rhoC^2+mD*(LC^2+rhoD^2+2*LC*rhoD*c2)+mH*(LC^2+LD^2+2*LC*LD*c2)));

Q2d2 = ((IC+ID+mC*rhoC^2+mD*(LC^2+rhoD^2+2*LC*rhoD*c2)+mH*(LC^2+LD^2+2*LC*LD*c2)).*(TK+g*LD*mH*s1_2 ... 
    +rhoD*(g*mD*s1_2-r*c2.*(F3./(r^2+rhoD^2-2*r*rhoD*s2).^0.5-F4./(r^2+rhoD^2+2*r*rhoD*s2).^0.5) ... 
    -F5.*(LC*s2+r*c1_2)./(r^2+LC^2+rhoD^2+2*r*LC*s1+2*LC*rhoD*c2+2*r*rhoD*s1_2).^0.5) ... 
    +LC*(LD*mH+mD*rhoD)*s2.*Q1d.^2)+(ID+LD*mH*(LD+LC*c2)+mD*rhoD*(rhoD+LC*c2)).*(TA+r*F5.*(LC*c1 ... 
    +rhoD*c1_2)./(r^2+LC^2+rhoD^2+2*r*LC*s1+2*LC*rhoD*c2+2*r*rhoD*s1_2).^0.5-g*LC*mD*s1 ... 
    -g*mD*rhoD*s1_2-g*mH*(LC*s1+LD*s1_2)-rhoC*(g*mC*s1-r*c1.*(F1./(r^2+rhoC^2+2*r*rhoC*s1).^0.5 ... 
    -F2./(r^2+rhoC^2-2*r*rhoC*s1).^0.5))-LC*(LD*mH+mD*rhoD)*s2.*(Q1d.^2-(Q1d-Q2d).^2)))./((ID+LD*mH*(LD+LC*c2) ... 
    +mD*rhoD*(rhoD+LC*c2)).^2-(ID+mD*rhoD^2+mH*LD^2)*(IC+ID+mC*rhoC^2+mD*(LC^2+rhoD^2+2*LC*rhoD*c2)+mH*(LC^2+LD^2+2*LC*LD*c2)));



phaseout(1).dynamics = [Q1d*scale.angVel,Q2d*scale.angVel,Q1d2*scale.angAccel,Q2d2*scale.angAccel,a1d/scale.time]; %[U1,U2,U1d,U2d,a1d,vMT1];
phaseout(1).integrand = scale.cost*sum(e1.^2,2);


%-------------------------------------------------------------------------%
% Phase 2
%-------------------------------------------------------------------------%
t2 = solution.phase(2).time;
x2 = solution.phase(2).state(:,1:4);
a2 = solution.phase(2).state(:,5:9);
e2 = solution.phase(2).control(:,1:5);
numPts2 = size(x2,1);

Q1 = x2(:,1)/scale.angle;
Q2 = x2(:,2)/scale.angle;
Q1d = x2(:,3)/scale.angVel;
Q2d = x2(:,4)/scale.angVel;


c1 = cos(Q1);
s1 = sin(Q1);
c2 = cos(Q2);
s2 = sin(Q2);
c1_2 = cos((Q1-Q2));
s1_2 = sin((Q1-Q2));


lMT2(:,1) = (r^2+rhoC^2+2*r*rhoC*s1).^0.5;
lMT2(:,2) = (r^2+rhoC^2-2*r*rhoC*s1).^0.5;
lMT2(:,3) = (r^2+rhoD^2-2*r*rhoD*s2).^0.5;
lMT2(:,4) = (r^2+rhoD^2+2*r*rhoD*s2).^0.5;
lMT2(:,5) = (r^2+LC^2+rhoD^2+2*r*LC*s1+2*LC*rhoD*c2+2*r*rhoD*s1_2).^0.5;
vMT2(:,1) = r*rhoC*c1.*Q1d./(r^2+rhoC^2+2*r*rhoC*s1).^0.5;
vMT2(:,2) = -r*rhoC*c1.*Q1d./(r^2+rhoC^2-2*r*rhoC*s1).^0.5;
vMT2(:,3) = -r*rhoD*c2.*Q2d./(r^2+rhoD^2-2*r*rhoD*s2).^0.5;
vMT2(:,4) = r*rhoD*c2.*Q2d./(r^2+rhoD^2+2*r*rhoD*s2).^0.5;
vMT2(:,5) = -(LC*rhoD*s2.*Q2d-r*LC*c1.*Q1d-r*rhoD*c1_2.*(Q1d-Q2d))./(r^2+LC^2+rhoD^2+2*r*LC*s1+2*LC*rhoD*c2+2*r*rhoD*s1_2).^0.5;


a2d = activationDynamics(e2,a2,tact,tdeact);
[F1,lMtilda1,FMltilda1,vMtilda1,FMvtilda1] = rigidTendonModel_simple(a2(:,1),lMT2(:,1),vMT2(:,1),muscParams(:,1),cCoefs,dCoefs,eCoefs);
[F2,lMtilda2,FMltilda2,vMtilda2,FMvtilda2] = rigidTendonModel_simple(a2(:,2),lMT2(:,2),vMT2(:,2),muscParams(:,2),cCoefs,dCoefs,eCoefs);
[F3,lMtilda3,FMltilda3,vMtilda3,FMvtilda3] = rigidTendonModel_simple(a2(:,3),lMT2(:,3),vMT2(:,3),muscParams(:,3),cCoefs,dCoefs,eCoefs);
[F4,lMtilda4,FMltilda4,vMtilda4,FMvtilda4] = rigidTendonModel_simple(a2(:,4),lMT2(:,4),vMT2(:,4),muscParams(:,4),cCoefs,dCoefs,eCoefs);
[F5,lMtilda5,FMltilda5,vMtilda5,FMvtilda5] = rigidTendonModel_simple(a2(:,5),lMT2(:,5),vMT2(:,5),muscParams(:,5),cCoefs,dCoefs,eCoefs);

    
Q1d2 = ((ID+LD*mH*(LD+LC*c2)+mD*rhoD*(rhoD+LC*c2)).*(TK+g*LD*mH*s1_2 ... 
    +rhoD*(g*mD*s1_2-r*c2.*(F3./(r^2+rhoD^2-2*r*rhoD*s2).^0.5-F4./(r^2+rhoD^2+2*r*rhoD*s2).^0.5) ... 
    -F5.*(LC*s2+r*c1_2)./(r^2+LC^2+rhoD^2+2*r*LC*s1+2*LC*rhoD*c2 ... 
    +2*r*rhoD*s1_2).^0.5)+LC*(LD*mH+mD*rhoD)*s2.*Q1d.^2)+(ID+mD*rhoD^2+mH*LD^2)*(TA+r*F5.*(LC*c1 ... 
    +rhoD*c1_2)./(r^2+LC^2+rhoD^2+2*r*LC*s1+2*LC*rhoD*c2+2*r*rhoD*s1_2).^0.5 ... 
    -g*LC*mD*s1-g*mD*rhoD*s1_2-g*mH*(LC*s1+LD*s1_2)-rhoC*(g*mC*s1 ... 
    -r*c1.*(F1./(r^2+rhoC^2+2*r*rhoC*s1).^0.5-F2./(r^2+rhoC^2-2*r*rhoC*s1).^0.5)) ... 
    -LC*(LD*mH+mD*rhoD)*s2.*(Q1d.^2-(Q1d-Q2d).^2)))./((ID+LD*mH*(LD+LC*c2)+mD*rhoD*(rhoD+LC*c2)).^2 ... 
    -(ID+mD*rhoD^2+mH*LD^2)*(IC+ID+mC*rhoC^2+mD*(LC^2+rhoD^2+2*LC*rhoD*c2)+mH*(LC^2+LD^2+2*LC*LD*c2)));

Q2d2 = ((IC+ID+mC*rhoC^2+mD*(LC^2+rhoD^2+2*LC*rhoD*c2)+mH*(LC^2+LD^2+2*LC*LD*c2)).*(TK+g*LD*mH*s1_2 ... 
    +rhoD*(g*mD*s1_2-r*c2.*(F3./(r^2+rhoD^2-2*r*rhoD*s2).^0.5-F4./(r^2+rhoD^2+2*r*rhoD*s2).^0.5) ... 
    -F5.*(LC*s2+r*c1_2)./(r^2+LC^2+rhoD^2+2*r*LC*s1+2*LC*rhoD*c2+2*r*rhoD*s1_2).^0.5) ... 
    +LC*(LD*mH+mD*rhoD)*s2.*Q1d.^2)+(ID+LD*mH*(LD+LC*c2)+mD*rhoD*(rhoD+LC*c2)).*(TA+r*F5.*(LC*c1 ... 
    +rhoD*c1_2)./(r^2+LC^2+rhoD^2+2*r*LC*s1+2*LC*rhoD*c2+2*r*rhoD*s1_2).^0.5-g*LC*mD*s1 ... 
    -g*mD*rhoD*s1_2-g*mH*(LC*s1+LD*s1_2)-rhoC*(g*mC*s1-r*c1.*(F1./(r^2+rhoC^2+2*r*rhoC*s1).^0.5 ... 
    -F2./(r^2+rhoC^2-2*r*rhoC*s1).^0.5))-LC*(LD*mH+mD*rhoD)*s2.*(Q1d.^2-(Q1d-Q2d).^2)))./((ID+LD*mH*(LD+LC*c2) ... 
    +mD*rhoD*(rhoD+LC*c2)).^2-(ID+mD*rhoD^2+mH*LD^2)*(IC+ID+mC*rhoC^2+mD*(LC^2+rhoD^2+2*LC*rhoD*c2)+mH*(LC^2+LD^2+2*LC*LD*c2)));


phaseout(2).dynamics = [Q1d*scale.angVel,Q2d*scale.angVel,Q1d2*scale.angAccel,Q2d2*scale.angAccel,a2d/scale.time]; %[U1,U2,U1d,U2d,a2d,vMT2];
phaseout(2).integrand = scale.cost*sum(e2.^2,2);



figure(1)
subplot(1,2,1);
plot(t1,x1(:,1),'b-o',t2,x2(:,1),'b-*')
xlabel('Time (s)');
ylabel('q1');

subplot(1,2,2);
plot(t1,x1(:,2),'b-o',t2,x2(:,2),'b-*')
xlabel('Time (s)');
ylabel('q2');

figure(2)
plot(t1,a1(:,1),'r-o',t2,a2(:,1),'r-*');
hold on
plot(t1,a1(:,2),'g-o',t2,a2(:,2),'g-*');
plot(t1,a1(:,3),'b-o',t2,a2(:,3),'b-*');
plot(t1,a1(:,4),'k-o',t2,a2(:,4),'k-*');
plot(t1,a1(:,5),'m-o',t2,a2(:,5),'m-*');
xlabel('Time (s)');
ylabel('Activation');

figure(3)
plot(t1,e1(:,1),'r-o',t2,e2(:,1),'r-*');
hold on
plot(t1,e1(:,2),'g-o',t2,e2(:,2),'g-*');
plot(t1,e1(:,3),'b-o',t2,e2(:,3),'b-*');
plot(t1,e1(:,4),'k-o',t2,e2(:,4),'k-*');
plot(t1,e1(:,5),'m-o',t2,e2(:,5),'m-*');
xlabel('Time (s)');
ylabel('Excitation');



figure(4);
q1 = [solution.phase(1).state(:,1);solution.phase(2).state(:,1)];
q2 = [solution.phase(1).state(:,2);solution.phase(2).state(:,2)];
N = length(q1);
duration = tf2;
LC = auxdata.LC;
LD = auxdata.LD;

% make movie of the solution
fps = N/duration;
avi = VideoWriter('squat.avi','Uncompressed AVI');
open(avi);

set(gcf, 'color', 'white');
s = 2*LC;
for i=1:N
	plot([-s s],[0 0],'k','LineWidth',2);
    hold on
    plot([0 LC*sin(q1(i))], [0 LC*cos(q1(i))],'b-o','LineWidth',2);
    plot([LC*sin(q1(i)) LC*sin(q1(i))+LD*sin(q1(i)-q2(i))],[LC*cos(q1(i)), LC*cos(q1(i))+LD*cos(q1(i)-q2(i))],'k-o','LineWidth',4);
    axis('equal');
    axis('square');
    axis([-s s -0.5*s 1.5*s]);
    if (i==1)
        F = getframe(gca);
        frame = [1 1 size(F.cdata,2) size(F.cdata,1)];
    else
        F = getframe(gca,frame);
    end

    writeVideo(avi,F);
    if i ~= N
        cla;
    end
end

close(avi);

