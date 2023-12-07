function phaseout = squat_continuous(input)

% Unpack auxdata
TA = input.auxdata.TA;
TK = input.auxdata.TK;
g = input.auxdata.g;
r = input.auxdata.r;
IC = input.auxdata.IC;
ID = input.auxdata.ID;
LC = input.auxdata.LC;
LD = input.auxdata.LD;
mC = input.auxdata.mC;
mD = input.auxdata.mD;
mH = input.auxdata.mH;
rhoC = input.auxdata.rhoC;
rhoD = input.auxdata.rhoD;


tact = input.auxdata.tact;
tdeact = input.auxdata.tdeact;
muscParams = input.auxdata.muscParams;
cCoefs = input.auxdata.cCoefs;
dCoefs = input.auxdata.dCoefs;
eCoefs = input.auxdata.eCoefs;


scale = input.auxdata.scale;


%-------------------------------------------------------------------------%
%% Phase 1
%-------------------------------------------------------------------------%

% Unpack time, states, and controls
t1 = input.phase(1).time;
x1 = input.phase(1).state(:,1:4);
a1 = input.phase(1).state(:,5:9);
e1 = input.phase(1).control(:,1:5);
numPts1 = size(x1,1);

Q1 = x1(:,1)/scale.angle;
Q2 = x1(:,2)/scale.angle;
Q1d = x1(:,3)/scale.angVel;
Q2d = x1(:,4)/scale.angVel;

% Precalculate sines and cosines for speed
c1 = cos(Q1);
s1 = sin(Q1);
c2 = cos(Q2);
s2 = sin(Q2);
c1_2 = cos((Q1-Q2));
s1_2 = sin((Q1-Q2));

% Calculate muscle-tendon lengths and velocities
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

% Muscle dynamics
a1d = activationDynamics(e1,a1,tact,tdeact);
[F1,lMtilda1,FMltilda1,vMtilda1,FMvtilda1] = rigidTendonModel_simple(a1(:,1),lMT1(:,1),vMT1(:,1),muscParams(:,1),cCoefs,dCoefs,eCoefs);
[F2,lMtilda2,FMltilda2,vMtilda2,FMvtilda2] = rigidTendonModel_simple(a1(:,2),lMT1(:,2),vMT1(:,2),muscParams(:,2),cCoefs,dCoefs,eCoefs);
[F3,lMtilda3,FMltilda3,vMtilda3,FMvtilda3] = rigidTendonModel_simple(a1(:,3),lMT1(:,3),vMT1(:,3),muscParams(:,3),cCoefs,dCoefs,eCoefs);
[F4,lMtilda4,FMltilda4,vMtilda4,FMvtilda4] = rigidTendonModel_simple(a1(:,4),lMT1(:,4),vMT1(:,4),muscParams(:,4),cCoefs,dCoefs,eCoefs);
[F5,lMtilda5,FMltilda5,vMtilda5,FMvtilda5] = rigidTendonModel_simple(a1(:,5),lMT1(:,5),vMT1(:,5),muscParams(:,5),cCoefs,dCoefs,eCoefs);

% Skeletal dynamics
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


% Phase 1 output
phaseout(1).dynamics = [Q1d*scale.angVel,Q2d*scale.angVel,Q1d2*scale.angAccel,Q2d2*scale.angAccel,a1d/scale.time]; %[U1,U2,U1d,U2d,a1d,vMT1];
phaseout(1).integrand = scale.cost*sum(e1.^2,2);


%-------------------------------------------------------------------------%
%% Phase 2
%-------------------------------------------------------------------------%

% Unpack time, states, and contols
t2 = input.phase(2).time;
x2 = input.phase(2).state(:,1:4);
a2 = input.phase(2).state(:,5:9);
e2 = input.phase(2).control(:,1:5);
numPts2 = size(x2,1);

Q1 = x2(:,1)/scale.angle;
Q2 = x2(:,2)/scale.angle;
Q1d = x2(:,3)/scale.angVel;
Q2d = x2(:,4)/scale.angVel;

% Precalculate sines and cosines for speed
c1 = cos(Q1);
s1 = sin(Q1);
c2 = cos(Q2);
s2 = sin(Q2);
c1_2 = cos((Q1-Q2));
s1_2 = sin((Q1-Q2));

% Calculate muscle-tendon lengths and velocities
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

% Muscle dynamics
a2d = activationDynamics(e2,a2,tact,tdeact);
[F1,lMtilda1,FMltilda1,vMtilda1,FMvtilda1] = rigidTendonModel_simple(a2(:,1),lMT2(:,1),vMT2(:,1),muscParams(:,1),cCoefs,dCoefs,eCoefs);
[F2,lMtilda2,FMltilda2,vMtilda2,FMvtilda2] = rigidTendonModel_simple(a2(:,2),lMT2(:,2),vMT2(:,2),muscParams(:,2),cCoefs,dCoefs,eCoefs);
[F3,lMtilda3,FMltilda3,vMtilda3,FMvtilda3] = rigidTendonModel_simple(a2(:,3),lMT2(:,3),vMT2(:,3),muscParams(:,3),cCoefs,dCoefs,eCoefs);
[F4,lMtilda4,FMltilda4,vMtilda4,FMvtilda4] = rigidTendonModel_simple(a2(:,4),lMT2(:,4),vMT2(:,4),muscParams(:,4),cCoefs,dCoefs,eCoefs);
[F5,lMtilda5,FMltilda5,vMtilda5,FMvtilda5] = rigidTendonModel_simple(a2(:,5),lMT2(:,5),vMT2(:,5),muscParams(:,5),cCoefs,dCoefs,eCoefs);

% Skeletal dynamics
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



% Phase 2 output
phaseout(2).dynamics = [Q1d*scale.angVel,Q2d*scale.angVel,Q1d2*scale.angAccel,Q2d2*scale.angAccel,a2d/scale.time]; %[U1,U2,U1d,U2d,a2d,vMT2];
phaseout(2).integrand = scale.cost*sum(e2.^2,2);
