%=========================================================================%
% Squatting problem: Optimize muscle excitations to match joint
% accelerations from previous rigid-tendon optimization
% States: a,lMtilda
% Controls: e,vMtilda,resAct
%=========================================================================%
clear all; close all; clc;

%% Scaling
scale.time = 30;
scale.mass = 1;
scale.length = 1;
scale.vel = scale.length/scale.time;
scale.accel = scale.vel/scale.time;
scale.force = scale.mass*scale.accel;
scale.torque = scale.force*scale.length;
scale.angVel = 1/scale.time;
scale.angAccel = scale.angVel/scale.time;
scale.inertia = scale.torque/scale.angAccel;

scale.cost = 1e1/scale.time;

%% Load Data
load('muscleForce_data.mat');
auxdata.scale = scale;

auxdata.tact = 0.01;
auxdata.tdeact = 0.04;
auxdata.kT = 1e9;


%% Bounds
a_min = zeros(1,5);
a_max = ones(1,5);
lMtilda_min = zeros(1,5);
lMtilda_max = 3*ones(1,5);

e_min = zeros(1,5);
e_max = ones(1,5);
vMtilda_min = -ones(1,5);
vMtilda_max = ones(1,5);
resAct_min = -200*ones(1,2);
resAct_max = 200*ones(1,2);

path_min = zeros(1,7);
path_max = zeros(1,7);

t0 = 0;
tf = 1;
a0_min = auxdata.a0_ref;
a0_max = auxdata.a0_ref;
af_min = auxdata.af_ref;
af_max = auxdata.af_ref;
lMtilda0_min = auxdata.lMtilda0_ref;
lMtilda0_max = auxdata.lMtilda0_ref;
lMtildaf_min = auxdata.lMtildaf_ref;
lMtildaf_max = auxdata.lMtildaf_ref;

bounds.phase.initialtime.lower = t0*scale.time; 
bounds.phase.initialtime.upper = t0*scale.time;
bounds.phase.finaltime.lower = tf*scale.time; 
bounds.phase.finaltime.upper = tf*scale.time;
bounds.phase.initialstate.lower = [a0_min,lMtilda_min];
bounds.phase.initialstate.upper = [a0_max,lMtilda_max];
bounds.phase.state.lower = [a_min,lMtilda_min];
bounds.phase.state.upper = [a_max,lMtilda_max];
bounds.phase.finalstate.lower = [af_min,lMtildaf_min];
bounds.phase.finalstate.upper = [af_max,lMtildaf_max]; 
bounds.phase.control.lower = [e_min,vMtilda_min,resAct_min*scale.torque];
bounds.phase.control.upper = [e_max,vMtilda_max,resAct_max*scale.torque];
bounds.phase.integral.lower = 0; 
bounds.phase.integral.upper = 1e5;
bounds.phase.path.lower = path_min;
bounds.phase.path.upper = path_max;

%-------------------------------------------------------------------------%
%% Guess
%-------------------------------------------------------------------------%

guess.phase.time    = [t0*scale.time; tf*scale.time]; 
guess.phase.state   = [a_min,ones(1,5); a_min,ones(1,5)];
guess.phase.control = zeros(2,12);
guess.phase.integral = 0.1*scale.cost;


%-------------------------------------------------------------------------%
%% Problem Structure
%-------------------------------------------------------------------------%
setup.name = 'MuscleOptimization_rigidTendon';
setup.functions.continuous = @muscleForce_continuous_completed;
setup.functions.endpoint = @muscleForce_endpoint_completed;
setup.auxdata = auxdata;
setup.bounds = bounds;
setup.guess = guess;

setup.nlp.solver = 'ipopt';
setup.nlp.ipoptoptions.tolerance = 1e-7;
setup.derivatives.supplier = 'sparseCD';
setup.derivatives.derivativelevel = 'second';
setup.derivatives.dependencies = 'sparse';

meshphase.colpoints = 4*ones(1,20);
meshphase.fraction = 1/20*ones(1,20);
setup.mesh.phase = meshphase;
setup.mesh.method = 'hp-PattersonRao';
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
toc

save('muscleForce_output','output');



%-------------------------------------------------------------------------%
%% Post Analysis
%-------------------------------------------------------------------------%
muscleForce_postAnalysis(auxdata,output);
