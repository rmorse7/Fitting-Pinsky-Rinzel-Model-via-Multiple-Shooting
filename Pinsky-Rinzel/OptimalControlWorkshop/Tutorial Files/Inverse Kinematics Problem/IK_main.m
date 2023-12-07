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

%-------------------------------------------------------------------------%
%% Scaling
%-------------------------------------------------------------------------%



%-------------------------------------------------------------------------%
%% Load Data
%-------------------------------------------------------------------------%



%-------------------------------------------------------------------------%
%% Define Bounds
%-------------------------------------------------------------------------%



%-------------------------------------------------------------------------%
%% Provide Guess
%-------------------------------------------------------------------------%



%-------------------------------------------------------------------------%
%% Assemble Problem Structure
%-------------------------------------------------------------------------%



%-------------------------------------------------------------------------%
%% Solve Problem Using GPOP2
%-------------------------------------------------------------------------%



%-------------------------------------------------------------------------%
%% Post Analysis
%-------------------------------------------------------------------------%
