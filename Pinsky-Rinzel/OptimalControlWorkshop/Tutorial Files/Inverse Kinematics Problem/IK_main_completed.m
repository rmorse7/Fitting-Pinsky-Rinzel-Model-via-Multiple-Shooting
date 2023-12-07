%=========================================================================%
% Marker Tracking Problem for Optimal Control Workshop
% Objective: Control joint torques to minimize marker error while
% maintaining dynamic consitency
% States: q1 = ankle angle, q2 = knee angle, q3 = angle of thigh relative
% to world, q4 = hip x position in world, q5 = hip y position in world, q1d
% = d(q1)/dt, q2d = d(q2)/dt, q3d = d(q3)/dt, q4d = d(q4)/dy, q5d =
% d(q5)/dt
% Controls: TA = ankle torque, TK = knee torque
% External Inputs: Ground reaction loads and reference marker trajectories
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

scale.cost = 1e2/scale.time;

%% Load Data
load('IK_data.mat')
auxdata.scale = scale;


%% Define Bounds
t0 = 0;
tf = 1;

q1_min = -pi/3;
q1_max = pi/3;
q2_min = -pi/12;
q2_max = pi/2;
q3_min = -pi/2;
q3_max = pi/2;
q4_min = -0.3;
q4_max = 0.3;
q5_min = 0.75;
q5_max = 1;

u1_min = -3*pi;
u1_max = 3*pi;
u2_min = -3*pi;
u2_max = 3*pi;
u3_min = -3*pi;
u3_max = 3*pi;
u4_min = -0.4;
u4_max = 0.4;
u5_min = -1.5;
u5_max = 1.5;

TA_min = -500;
TA_max = 500;
TK_min = -500;
TK_max = 500;

cost_min = 0;
cost_max = 10;

bounds.phase.initialtime.lower = t0*scale.time; 
bounds.phase.initialtime.upper = t0*scale.time;
bounds.phase.finaltime.lower = tf*scale.time; 
bounds.phase.finaltime.upper = tf*scale.time;
bounds.phase.initialstate.lower = [q1_min,q2_min,q3_min,q4_min*scale.length,q5_min*scale.length, ... 
    u1_min*scale.angVel,u2_min*scale.angVel,u3_min*scale.angVel, ... 
    u4_min*scale.vel,u5_min*scale.vel];
bounds.phase.initialstate.upper = [q1_max,q2_max,q3_max,q4_max*scale.length,q5_max*scale.length, ... 
    u1_max*scale.angVel,u2_max*scale.angVel,u3_max*scale.angVel, ... 
    u4_max*scale.vel,u5_max*scale.vel];
bounds.phase.finalstate.lower = [q1_min,q2_min,q3_min,q4_min*scale.length,q5_min*scale.length, ... 
    u1_min*scale.angVel,u2_min*scale.angVel,u3_min*scale.angVel, ...
    u4_min*scale.vel,u5_min*scale.vel];
bounds.phase.finalstate.upper = [q1_max,q2_max,q3_max,q4_max*scale.length,q5_max*scale.length, ... 
    u1_max*scale.angVel,u2_max*scale.angVel,u3_max*scale.angVel, ... 
    u4_max*scale.vel,u5_max*scale.vel];
bounds.phase.state.lower = [q1_min,q2_min,q3_min,q4_min*scale.length,q5_min*scale.length, ... 
    u1_min*scale.angVel,u2_min*scale.angVel,u3_min*scale.angVel, ... 
    u4_min*scale.vel,u5_min*scale.vel];
bounds.phase.state.upper = [q1_max,q2_max,q3_max,q4_max*scale.length,q5_max*scale.length, ... 
    u1_max*scale.angVel,u2_max*scale.angVel,u3_max*scale.angVel, ... 
    u4_max*scale.vel,u5_max*scale.vel];
bounds.phase.control.lower = [TA_min*scale.torque,TK_min*scale.torque];
bounds.phase.control.upper = [TA_max*scale.torque,TK_max*scale.torque];
bounds.phase.integral.lower = cost_min;
bounds.phase.integral.upper = cost_max*scale.cost;


%% Provide Guess
guess.phase.time = [t0*scale.time; tf*scale.time]; 
guess.phase.state = [zeros(1,3),0.2*scale.length,0.935*scale.length,zeros(1,5); 
    zeros(1,3),0.2*scale.length,0.935*scale.length,zeros(1,5)];
guess.phase.control = zeros(2,2);
guess.phase.integral = 0.001*scale.cost;


%-------------------------------------------------------------------------%
%% Assemble Problem Structure
%-------------------------------------------------------------------------%
setup.name = 'Inverse_Kinematics';
setup.functions.continuous = @IK_continuous_completed;
setup.functions.endpoint = @IK_endpoint_completed;
setup.auxdata = auxdata;
setup.bounds = bounds;
setup.guess = guess;

setup.nlp.solver = 'ipopt';
setup.nlp.ipoptoptions.tolerance = 1e-7;
setup.derivatives.supplier = 'sparseCD';
setup.derivatives.derivativelevel = 'second';
setup.derivatives.dependencies = 'sparse';

setup.mesh.phase.colpoints = 4*ones(1,20);
setup.mesh.phase.fraction = 1/20*ones(1,20);
setup.mesh.method = 'hp-LiuRao';
setup.mesh.tolerance = 1e-4;
setup.mesh.maxiterations = 10;
setup.mesh.colpointsmin = 4;
setup.mesh.colpointsmax = 10;

setup.method = 'RPM-Integration';


%-------------------------------------------------------------------------%
%% Solve Problem Using GPOP2
%-------------------------------------------------------------------------%
tic
output = gpops2(setup);
toc

save('IK_working_output','output');


%-------------------------------------------------------------------------%
%% Post Analysis
%-------------------------------------------------------------------------%
IK_postAnalysis(auxdata,output);
