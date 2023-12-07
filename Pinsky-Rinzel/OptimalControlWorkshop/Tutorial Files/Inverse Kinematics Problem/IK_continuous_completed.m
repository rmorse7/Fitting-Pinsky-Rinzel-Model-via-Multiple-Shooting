%=========================================================================%
% Continuous function for the inverse kinematics problem
%=========================================================================%

function output = IK_continuous_completed(input)
%% Unpack Auxdata
scale = input.auxdata.scale;
spline = input.auxdata.spline;


%% Unscale Time, States, and Controls
t = input.phase.time/scale.time;
q1 = input.phase.state(:,1);
q2 = input.phase.state(:,2);
q3 = input.phase.state(:,3);
q4 = input.phase.state(:,4)/scale.length;
q5 = input.phase.state(:,5)/scale.length;
u1 = input.phase.state(:,6)/scale.angVel;
u2 = input.phase.state(:,7)/scale.angVel;
u3 = input.phase.state(:,8)/scale.angVel;
u4 = input.phase.state(:,9)/scale.vel;
u5 = input.phase.state(:,10)/scale.vel;

TA = input.phase.control(:,1)/scale.torque;
TK = input.phase.control(:,2)/scale.torque;


%% Sample Splines
M_ref.xB1 = ppval(spline.markers.xB1,t);
M_ref.yB1 = ppval(spline.markers.yB1,t);
M_ref.xB2 = ppval(spline.markers.xB2,t);
M_ref.yB2 = ppval(spline.markers.yB2,t);
M_ref.xC1 = ppval(spline.markers.xC1,t);
M_ref.yC1 = ppval(spline.markers.yC1,t);
M_ref.xC2 = ppval(spline.markers.xC2,t);
M_ref.yC2 = ppval(spline.markers.yC2,t);
M_ref.xD1 = ppval(spline.markers.xD1,t);
M_ref.yD1 = ppval(spline.markers.yD1,t);
M_ref.xD2 = ppval(spline.markers.xD2,t);
M_ref.yD2 = ppval(spline.markers.yD2,t);


%% Skeletal Dynamics and Marker Positions
[u1d,u2d,u3d,u4d,u5d,markers] = skeletalDynamics(input.auxdata,t,q1,q2,q3,q4,q5,u1,u2,u3,u4,u5,TA,TK);


%% Unpack Marker Positions
xB1 = markers.xB1;
yB1 = markers.yB1;
xB2 = markers.xB2;
yB2 = markers.yB2;
xC1 = markers.xC1;
yC1 = markers.yC1;
xC2 = markers.xC2;
yC2 = markers.yC2;
xD1 = markers.xD1;
yD1 = markers.yD1;
xD2 = markers.xD2;
yD2 = markers.yD2;


%% Output
cost = ((xB1-M_ref.xB1).^2 + (yB1-M_ref.yB1).^2 + (xB2-M_ref.xB2).^2 + ... 
    (yB2-M_ref.yB2).^2 + (xC1-M_ref.xC1).^2 + (yC1-M_ref.yC1).^2 + ... 
    (xC2-M_ref.xC2).^2 + (yC2-M_ref.yC2).^2 + (xD1-M_ref.xD1).^2 + ... 
    (yD1-M_ref.yD1).^2 + (xD2-M_ref.xD2).^2 + (yD2-M_ref.yD2).^2);

output.integrand = cost*scale.length^2*scale.cost;
output.dynamics = [u1*scale.angVel,u2*scale.angVel,u3*scale.angVel,u4*scale.vel, ... 
    u5*scale.vel,u1d*scale.angAccel,u2d*scale.angAccel, ... 
    u3d*scale.angAccel,u4d*scale.accel,u5d*scale.accel];


end