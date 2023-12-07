function [] = IK_postAnalysis(auxdata,output)
close all

solution = output.result.solution;

F1 = 0; F2 = 0; F3 = 0; F4 = 0; F5 = 0;
FR1 = 0; FR2 = 0; TR3 = 0;

% Unpack auxdata
scale = auxdata.scale;

g = auxdata.g;
r = auxdata.r;
IB = auxdata.IB;
IC = auxdata.IC;
ID = auxdata.ID;
LB = auxdata.LB;
LC = auxdata.LC;
LD = auxdata.LD;
mB = auxdata.mB;
mC = auxdata.mC;
mD = auxdata.mD;
mH = auxdata.mH;
rhoB = auxdata.rhoB;
rhoC = auxdata.rhoC;
rhoD = auxdata.rhoD;
HB = auxdata.HB;


markers = auxdata.markers;
px1 = markers.px1; % Relative to A in B coord sys
py1 = markers.py1;
px2 = markers.px2;
py2 = markers.py2;
px3 = markers.px3; % Relative to K in C coord sys
py3 = markers.py3;
px4 = markers.px4;
py4 = markers.py4;
px5 = markers.px5; % Relative to H in D coord sys
py5 = markers.py5;
px6 = markers.px6;
py6 = markers.py6;


% Unpack and unscale time, states, and controls
t = solution.phase.time/scale.time;
npts = length(t);
Q1 = solution.phase.state(:,1);
Q2 = solution.phase.state(:,2);
Q3 = solution.phase.state(:,3);
Q4 = solution.phase.state(:,4)/scale.length;
Q5 = solution.phase.state(:,5)/scale.length;
Q1d = solution.phase.state(:,6)/scale.angVel;
Q2d = solution.phase.state(:,7)/scale.angVel;
Q3d = solution.phase.state(:,8)/scale.angVel;
Q4d = solution.phase.state(:,9)/scale.vel;
Q5d = solution.phase.state(:,10)/scale.vel;


TA = solution.phase.control(:,1)/scale.torque;
TK = solution.phase.control(:,2)/scale.torque;


% Get values at correct time points for ground reaction loads and marker
% position from splines.
FG1 = ppval(auxdata.spline.GRFx,t);
FG2 = ppval(auxdata.spline.GRFy,t);
TG3 = ppval(auxdata.spline.GRTz,t);

% M.rel = auxdata.markers.rel;
M_ref.B1x = ppval(auxdata.spline.markers.xB1,t);
M_ref.B1y = ppval(auxdata.spline.markers.yB1,t);
M_ref.B2x = ppval(auxdata.spline.markers.xB2,t);
M_ref.B2y = ppval(auxdata.spline.markers.yB2,t);
M_ref.C1x = ppval(auxdata.spline.markers.xC1,t);
M_ref.C1y = ppval(auxdata.spline.markers.yC1,t);
M_ref.C2x = ppval(auxdata.spline.markers.xC2,t);
M_ref.C2y = ppval(auxdata.spline.markers.yC2,t);
M_ref.D1x = ppval(auxdata.spline.markers.xD1,t);
M_ref.D1y = ppval(auxdata.spline.markers.yD1,t);
M_ref.D2x = ppval(auxdata.spline.markers.xD2,t);
M_ref.D2y = ppval(auxdata.spline.markers.yD2,t);


% Solve dynamics equations for ud

% Initialize z matrix
z = zeros(npts,240);

% Evaluate constants
z(:,1) = cos(Q1);
z(:,2) = sin(Q1);
z(:,3) = cos(Q2);
z(:,4) = sin(Q2);
z(:,5) = cos(Q3);
z(:,6) = sin(Q3);
z(:,7) = z(:,3).*z(:,5) + z(:,4).*z(:,6);
z(:,8) = z(:,4).*z(:,5) - z(:,3).*z(:,6);
z(:,9) = z(:,3).*z(:,6) - z(:,4).*z(:,5);
z(:,10) = z(:,1).*z(:,7) + z(:,2).*z(:,8);
z(:,11) = z(:,1).*z(:,8) - z(:,2).*z(:,7);
z(:,12) = z(:,1).*z(:,9) + z(:,2).*z(:,7);
z(:,13) = z(:,1).*z(:,7) - z(:,2).*z(:,9);
z(:,14) = rhoD - LD;
z(:,15) = rhoC - LC;
z(:,16) = LD.*z(:,3);
z(:,17) = LD.*z(:,4);
z(:,18) = z(:,16) - z(:,15);
z(:,19) = LC + z(:,16);
z(:,20) = z(:,17) - r;
z(:,21) = r + z(:,17);
z(:,22) = LB - rhoB;
z(:,23) = 0.25.*rhoB - HB;
z(:,24) = z(:,1).*z(:,19) + z(:,2).*z(:,17);
z(:,25) = LC.*z(:,1);
z(:,26) = z(:,1).*z(:,17) - z(:,2).*z(:,19);
z(:,27) = LC.*z(:,2);
z(:,28) = z(:,24) - z(:,23);
z(:,29) = z(:,23) - z(:,25);
z(:,30) = z(:,22) + z(:,26);
z(:,31) = z(:,27) - z(:,22);
z(:,32) = HB + z(:,24);
z(:,33) = -HB - z(:,25);
z(:,34) = LB + z(:,26);
z(:,35) = z(:,27) - LB;
z(:,36) = z(:,26) - r;
z(:,37) = r + z(:,27);
z(:,38) = r + z(:,26);
z(:,39) = z(:,27) - r;
z(:,40) = Q3d - Q2d;
z(:,41) = Q1d + Q3d - Q2d;
z(:,42) = z(:,5).*Q4d + z(:,6).*Q5d - z(:,14).*Q3d;
z(:,43) = z(:,5).*Q5d - z(:,6).*Q4d;
z(:,44) = Q3d.*(z(:,5).*Q5d-z(:,6).*Q4d);
z(:,45) = Q3d.*(z(:,5).*Q4d+z(:,6).*Q5d);
z(:,46) = z(:,44) - Q3d.*z(:,43);
z(:,47) = Q3d.*z(:,42) - z(:,45);
z(:,48) = z(:,15).*Q2d + z(:,7).*Q4d + z(:,9).*Q5d + z(:,18).*Q3d;
z(:,49) = z(:,7).*Q5d + z(:,8).*Q4d + z(:,17).*Q3d;
z(:,50) = z(:,3).*z(:,6).*Q2d + z(:,4).*z(:,5).*Q3d - z(:,3).*z(:,6).*Q3d - z(:,4).*z(:,5).*Q2d;
z(:,51) = z(:,3).*z(:,5).*Q3d + z(:,4).*z(:,6).*Q3d - z(:,3).*z(:,5).*Q2d - z(:,4).*z(:,6).*Q2d;
z(:,52) = LD.*z(:,4).*Q2d;
z(:,53) = Q4d.*z(:,50) + Q5d.*z(:,51) - Q3d.*z(:,52);
z(:,54) = z(:,3).*z(:,5).*Q2d + z(:,4).*z(:,6).*Q2d - z(:,3).*z(:,5).*Q3d - z(:,4).*z(:,6).*Q3d;
z(:,55) = LD.*z(:,3).*Q2d;
z(:,56) = Q3d.*z(:,55) + Q4d.*z(:,54) + Q5d.*z(:,50);
z(:,57) = z(:,53) - z(:,40).*z(:,49);
z(:,58) = z(:,56) + z(:,40).*z(:,48);
z(:,59) = z(:,10).*Q4d + z(:,12).*Q5d + z(:,28).*Q3d + z(:,29).*Q2d - z(:,23).*Q1d;
z(:,60) = z(:,22).*Q1d + z(:,11).*Q4d + z(:,13).*Q5d + z(:,30).*Q3d + z(:,31).*Q2d;
z(:,61) = z(:,1).*z(:,8).*Q1d + z(:,1).*z(:,50) + z(:,2).*z(:,54) - z(:,2).*z(:,7).*Q1d;
z(:,62) = z(:,1).*z(:,7).*Q1d + z(:,1).*z(:,51) + z(:,2).*z(:,50) - z(:,2).*z(:,9).*Q1d;
z(:,63) = z(:,1).*z(:,17).*Q1d + z(:,2).*z(:,55) - z(:,2).*z(:,19).*Q1d - z(:,1).*z(:,52);
z(:,64) = LC.*z(:,2).*Q1d;
z(:,65) = Q2d.*z(:,64) + Q3d.*z(:,63) + Q4d.*z(:,61) + Q5d.*z(:,62);
z(:,66) = z(:,1).*z(:,54) - z(:,1).*z(:,7).*Q1d - z(:,2).*z(:,8).*Q1d - z(:,2).*z(:,50);
z(:,67) = z(:,1).*z(:,50) - z(:,1).*z(:,9).*Q1d - z(:,2).*z(:,7).*Q1d - z(:,2).*z(:,51);
z(:,68) = z(:,1).*z(:,55) + z(:,2).*z(:,52) - z(:,1).*z(:,19).*Q1d - z(:,2).*z(:,17).*Q1d;
z(:,69) = LC.*z(:,1).*Q1d;
z(:,70) = Q2d.*z(:,69) + Q3d.*z(:,68) + Q4d.*z(:,66) + Q5d.*z(:,67);
z(:,71) = z(:,65) - z(:,41).*z(:,60);
z(:,72) = z(:,70) + z(:,41).*z(:,59);
z(:,73) = -LC - z(:,15);
z(:,74) = -LD - z(:,14);
z(:,75) = r.^2 + z(:,73).^2;
z(:,76) = r.*z(:,73);
z(:,77) = z(:,75) - 2.*z(:,76).*z(:,2);
z(:,78) = z(:,77).^0.5;
z(:,79) = r./z(:,78);
z(:,80) = z(:,73)./z(:,78);
z(:,81) = z(:,75) + 2.*z(:,76).*z(:,2);
z(:,82) = z(:,81).^0.5;
z(:,83) = r./z(:,82);
z(:,84) = z(:,73)./z(:,82);
z(:,85) = r.^2 + z(:,74).^2;
z(:,86) = r.*z(:,74);
z(:,87) = z(:,85) + 2.*z(:,86).*z(:,4);
z(:,88) = z(:,87).^0.5;
z(:,89) = r./z(:,88);
z(:,90) = z(:,74)./z(:,88);
z(:,91) = z(:,85) - 2.*z(:,86).*z(:,4);
z(:,92) = z(:,91).^0.5;
z(:,93) = r./z(:,92);
z(:,94) = z(:,74)./z(:,92);
z(:,95) = z(:,1).*z(:,3) + z(:,2).*z(:,4);
z(:,96) = z(:,2).*z(:,3) - z(:,1).*z(:,4);
z(:,98) = LC.^2 + r.^2 + z(:,74).^2;
z(:,99) = LC.*r;
z(:,100) = LC.*z(:,74);
z(:,101) = z(:,98) + 2.*z(:,99).*z(:,2) - 2.*z(:,86).*z(:,96) - 2.*z(:,100).*z(:,3);
z(:,102) = z(:,101).^0.5;
z(:,103) = r./z(:,102);
z(:,104) = LC./z(:,102);
z(:,105) = z(:,74)./z(:,102);
z(:,106) = g.*mB;
z(:,107) = g.*mC;
z(:,108) = g.*mD;
z(:,109) = g.*mH;
z(:,110) = FR2 - z(:,109);
z(:,111) = Q5 + LB.*z(:,12) - HB.*z(:,13) - LC.*z(:,7) - LD.*z(:,5);
z(:,112) = HB.*z(:,11) + LC.*z(:,8) - Q4 - LB.*z(:,10) - LD.*z(:,6);
z(:,113) = z(:,111) + HB.*z(:,10) + LB.*z(:,11);
z(:,114) = z(:,112) + HB.*z(:,12) + LB.*z(:,13);
z(:,115) = r.*(z(:,1).*z(:,104)-z(:,95).*z(:,105));
z(:,116) = r.*z(:,1).*z(:,80);
z(:,117) = r.*z(:,1).*z(:,84);
z(:,118) = z(:,106).*(z(:,22).*z(:,13)-z(:,23).*z(:,12));
z(:,119) = z(:,10).*z(:,33) + z(:,11).*z(:,35) - z(:,111);
z(:,120) = z(:,12).*z(:,33) + z(:,13).*z(:,35) - z(:,112);
z(:,121) = z(:,2).*z(:,25).*z(:,80) - z(:,25).*z(:,79) - z(:,15).*z(:,1).*z(:,79) - z(:,1).*z(:,37).*z(:,80);
z(:,122) = z(:,25).*z(:,83) + z(:,15).*z(:,1).*z(:,83) + z(:,2).*z(:,25).*z(:,84) - z(:,1).*z(:,39).*z(:,84);
z(:,123) = z(:,1).*z(:,37).*z(:,104) + z(:,25).*z(:,96).*z(:,105) - z(:,25).*z(:,103) - z(:,2).*z(:,25).*z(:,104) - z(:,37).*z(:,95).*z(:,105);
z(:,124) = r.*z(:,3).*z(:,90);
z(:,125) = r.*z(:,3).*z(:,94);
z(:,126) = z(:,15).*z(:,107);
z(:,127) = -z(:,126).*z(:,9) - z(:,106).*(z(:,12).*z(:,29)+z(:,13).*z(:,31));
z(:,128) = z(:,111) + z(:,10).*z(:,32) + z(:,11).*z(:,34);
z(:,129) = z(:,112) + z(:,12).*z(:,32) + z(:,13).*z(:,34);
z(:,130) = z(:,17).*z(:,80) + z(:,24).*z(:,79) - z(:,1).*z(:,18).*z(:,79) - z(:,1).*z(:,36).*z(:,80) - z(:,2).*z(:,17).*z(:,79) - z(:,2).*z(:,24).*z(:,80);
z(:,131) = z(:,17).*z(:,84) + z(:,1).*z(:,18).*z(:,83) + z(:,2).*z(:,17).*z(:,83) - z(:,24).*z(:,83) - z(:,1).*z(:,38).*z(:,84) - z(:,2).*z(:,24).*z(:,84);
z(:,132) = z(:,24).*z(:,103) + z(:,1).*z(:,36).*z(:,104) + z(:,2).*z(:,24).*z(:,104) + z(:,14).*(z(:,4).*z(:,104)+z(:,95).*z(:,103)) - z(:,24).*z(:,96).*z(:,105) - z(:,36).*z(:,95).*z(:,105);
z(:,133) = z(:,16).*z(:,89) + z(:,14).*z(:,3).*z(:,89) + z(:,4).*z(:,16).*z(:,90) - z(:,3).*z(:,20).*z(:,90);
z(:,134) = z(:,4).*z(:,16).*z(:,94) - z(:,16).*z(:,93) - z(:,14).*z(:,3).*z(:,93) - z(:,3).*z(:,21).*z(:,94);
z(:,135) = z(:,14).*z(:,108);
z(:,136) = TR3 + z(:,135).*z(:,6) - z(:,107).*z(:,7).*z(:,17) - z(:,107).*z(:,9).*z(:,18) - z(:,106).*(z(:,12).*z(:,28)+z(:,13).*z(:,30));
z(:,137) = z(:,8).*z(:,80) + z(:,10).*z(:,79) - z(:,1).*z(:,7).*z(:,79) - z(:,1).*z(:,11).*z(:,80) - z(:,2).*z(:,8).*z(:,79) - z(:,2).*z(:,10).*z(:,80);
z(:,138) = z(:,8).*z(:,84) + z(:,1).*z(:,7).*z(:,83) + z(:,2).*z(:,8).*z(:,83) - z(:,10).*z(:,83) - z(:,1).*z(:,11).*z(:,84) - z(:,2).*z(:,10).*z(:,84);
z(:,139) = z(:,10).*z(:,103) + z(:,1).*z(:,11).*z(:,104) + z(:,2).*z(:,10).*z(:,104) + z(:,3).*z(:,6).*z(:,104) + z(:,6).*z(:,96).*z(:,103) - z(:,6).*z(:,105) - z(:,4).*z(:,5).*z(:,104) - z(:,5).*z(:,95).*z(:,103) - z(:,10).*z(:,96).*z(:,105) - z(:,11).*z(:,95).*z(:,105);
z(:,140) = z(:,7).*z(:,89) + z(:,4).*z(:,7).*z(:,90) - z(:,6).*z(:,90) - z(:,3).*z(:,5).*z(:,89) - z(:,3).*z(:,8).*z(:,90) - z(:,4).*z(:,6).*z(:,89);
z(:,141) = z(:,3).*z(:,5).*z(:,93) + z(:,4).*z(:,6).*z(:,93) + z(:,4).*z(:,7).*z(:,94) - z(:,6).*z(:,94) - z(:,7).*z(:,93) - z(:,3).*z(:,8).*z(:,94);
z(:,142) = z(:,10).^2 + z(:,11).^2;
z(:,143) = z(:,10).*z(:,12) + z(:,11).*z(:,13);
z(:,144) = FR1 - z(:,107).*z(:,7).*z(:,8) - z(:,107).*z(:,7).*z(:,9) - z(:,106).*(z(:,10).*z(:,12)+z(:,11).*z(:,13));
z(:,145) = z(:,7).*z(:,80) + z(:,12).*z(:,79) - z(:,1).*z(:,9).*z(:,79) - z(:,1).*z(:,13).*z(:,80) - z(:,2).*z(:,7).*z(:,79) - z(:,2).*z(:,12).*z(:,80);
z(:,146) = z(:,7).*z(:,84) + z(:,1).*z(:,9).*z(:,83) + z(:,2).*z(:,7).*z(:,83) - z(:,12).*z(:,83) - z(:,1).*z(:,13).*z(:,84) - z(:,2).*z(:,12).*z(:,84);
z(:,147) = z(:,5).*z(:,90) + z(:,9).*z(:,89) + z(:,4).*z(:,5).*z(:,89) + z(:,4).*z(:,9).*z(:,90) - z(:,3).*z(:,6).*z(:,89) - z(:,3).*z(:,7).*z(:,90);
z(:,148) = z(:,5).*z(:,94) + z(:,3).*z(:,6).*z(:,93) + z(:,4).*z(:,9).*z(:,94) - z(:,9).*z(:,93) - z(:,3).*z(:,7).*z(:,94) - z(:,4).*z(:,5).*z(:,93);
z(:,149) = z(:,5).*z(:,105) + z(:,12).*z(:,103) + z(:,1).*z(:,13).*z(:,104) + z(:,2).*z(:,12).*z(:,104) - z(:,3).*z(:,5).*z(:,104) - z(:,4).*z(:,6).*z(:,104) - z(:,5).*z(:,96).*z(:,103) - z(:,6).*z(:,95).*z(:,103) - z(:,12).*z(:,96).*z(:,105) - z(:,13).*z(:,95).*z(:,105);
z(:,150) = z(:,12).^2 + z(:,13).^2;
z(:,151) = z(:,110) - z(:,107).*z(:,7).^2 - z(:,107).*z(:,9).^2 - z(:,108).*z(:,5).^2 - z(:,108).*z(:,6).^2 - z(:,106).*(z(:,12).^2+z(:,13).^2);
z(:,152) = mB.*(z(:,22).*z(:,31)-z(:,23).*z(:,29)) - IB;
z(:,153) = IB + mB.*(z(:,22).^2+z(:,23).^2);
z(:,154) = IB + mB.*(z(:,22).*z(:,30)-z(:,23).*z(:,28));
z(:,155) = mB.*(z(:,22).*z(:,11)-z(:,23).*z(:,10));
z(:,156) = mB.*(z(:,22).*z(:,13)-z(:,23).*z(:,12));
z(:,157) = mB.*(z(:,22).*z(:,72)-z(:,23).*z(:,71));
z(:,158) = IB + IC + mC.*z(:,15).^2;
z(:,159) = z(:,158) + mB.*(z(:,29).^2+z(:,31).^2);
z(:,160) = mC.*z(:,15);
z(:,161) = z(:,160).*z(:,18) + mB.*(z(:,28).*z(:,29)+z(:,30).*z(:,31)) - IB - IC;
z(:,162) = z(:,160).*z(:,7) + mB.*(z(:,10).*z(:,29)+z(:,11).*z(:,31));
z(:,163) = z(:,160).*z(:,9) + mB.*(z(:,12).*z(:,29)+z(:,13).*z(:,31));
z(:,164) = z(:,160).*z(:,57) + mB.*(z(:,29).*z(:,71)+z(:,31).*z(:,72));
z(:,165) = IB + IC + ID + mD.*z(:,14).^2;
z(:,166) = z(:,165) + mB.*(z(:,28).^2+z(:,30).^2) + mC.*(z(:,17).^2+z(:,18).^2);
z(:,167) = mD.*z(:,14);
z(:,168) = mB.*(z(:,10).*z(:,28)+z(:,11).*z(:,30)) + mC.*(z(:,7).*z(:,18)+z(:,8).*z(:,17)) - z(:,167).*z(:,5);
z(:,169) = mB.*(z(:,12).*z(:,28)+z(:,13).*z(:,30)) + mC.*(z(:,7).*z(:,17)+z(:,9).*z(:,18)) - z(:,167).*z(:,6);
z(:,170) = mB.*(z(:,28).*z(:,71)+z(:,30).*z(:,72)) + mC.*(z(:,17).*z(:,58)+z(:,18).*z(:,57)) - z(:,167).*z(:,46);
z(:,171) = mH + mB.*(z(:,10).^2+z(:,11).^2) + mC.*(z(:,7).^2+z(:,8).^2) + mD.*(z(:,5).^2+z(:,6).^2);
z(:,172) = mC.*z(:,7).*(z(:,8)+z(:,9)) + mB.*(z(:,10).*z(:,12)+z(:,11).*z(:,13));
z(:,173) = mB.*(z(:,10).*z(:,71)+z(:,11).*z(:,72)) + mC.*(z(:,7).*z(:,57)+z(:,8).*z(:,58)) + mD.*(z(:,5).*z(:,46)-z(:,6).*z(:,47));
z(:,174) = mH + mB.*(z(:,12).^2+z(:,13).^2) + mC.*(z(:,7).^2+z(:,9).^2) + mD.*(z(:,5).^2+z(:,6).^2);
z(:,175) = mB.*(z(:,12).*z(:,71)+z(:,13).*z(:,72)) + mC.*(z(:,7).*z(:,58)+z(:,9).*z(:,57)) + mD.*(z(:,5).*z(:,47)+z(:,6).*z(:,46));


% Construct mass matrix and right hand side
Amat11 = z(:,153);
Amat12 = z(:,152);
Amat13 = z(:,154);
Amat14 = z(:,155);
Amat15 = z(:,156);
Amat21 = z(:,152);
Amat22 = z(:,159);
Amat23 = z(:,161);
Amat24 = z(:,162);
Amat25 = z(:,163);
Amat31 = z(:,154);
Amat32 = z(:,161);
Amat33 = z(:,166);
Amat34 = z(:,168);
Amat35 = z(:,169);
Amat41 = z(:,155);
Amat42 = z(:,162);
Amat43 = z(:,168);
Amat44 = z(:,171);
Amat45 = z(:,172);
Amat51 = z(:,156);
Amat52 = z(:,163);
Amat53 = z(:,169);
Amat54 = z(:,172);
Amat55 = z(:,174);
bvec1 = TG3 + F1.*z(:,116) + FG1.*z(:,113) + FG2.*z(:,114) - TA - z(:,118) - F2.*z(:,117) - F5.*z(:,115) - z(:,157);
bvec2 = z(:,127) + F1.*z(:,121) + F2.*z(:,122) + F4.*z(:,125) + F5.*z(:,123) + FG1.*z(:,119) + FG2.*z(:,120) - TG3 - TK - F3.*z(:,124) - z(:,164);
bvec3 = TG3 + z(:,136) + F1.*z(:,130) + F2.*z(:,131) + F3.*z(:,133) + F4.*z(:,134) + F5.*z(:,132) + FG1.*z(:,128) + FG2.*z(:,129) - z(:,170);
bvec4 = z(:,144) + F1.*z(:,137) + F2.*z(:,138) + F3.*z(:,140) + F4.*z(:,141) + F5.*z(:,139) + FG1.*z(:,142) + FG2.*z(:,143) - z(:,173);
bvec5 = z(:,151) + F1.*z(:,145) + F2.*z(:,146) + F3.*z(:,147) + F4.*z(:,148) + F5.*z(:,149) + FG1.*z(:,143) + FG2.*z(:,150) - z(:,175);

Qd2 = zeros(npts,5);
for i = 1:npts
    Amat = [Amat11(i),Amat12(i),Amat13(i),Amat14(i),Amat15(i);
        Amat21(i),Amat22(i),Amat23(i),Amat24(i),Amat25(i);
        Amat31(i),Amat32(i),Amat33(i),Amat34(i),Amat35(i);
        Amat41(i),Amat42(i),Amat43(i),Amat44(i),Amat45(i);
        Amat51(i),Amat52(i),Amat53(i),Amat54(i),Amat55(i)];
    bvec = [bvec1(i);bvec2(i);bvec3(i);bvec4(i);bvec5(i)];
    Qd2(i,:) = (Amat\bvec)';
end
Q1d2 = Qd2(:,1);
Q2d2 = Qd2(:,2);
Q3d2 = Qd2(:,3);
Q4d2 = Qd2(:,4);
Q5d2 = Qd2(:,5);

% Calculate marker positions and accelerations
xB1 = Q4 + LD.*z(:,6) + px1.*z(:,10) + py1.*z(:,11) - LC.*z(:,8);
yB1 = Q5 + px1.*z(:,12) + py1.*z(:,13) - LC.*z(:,7) - LD.*z(:,5);
xB2 = Q4 + LD.*z(:,6) + px2.*z(:,10) + py2.*z(:,11) - LC.*z(:,8);
yB2 = Q5 + px2.*z(:,12) + py2.*z(:,13) - LC.*z(:,7) - LD.*z(:,5);
xC1 = Q4 + LD.*z(:,6) + px3.*z(:,7) + py3.*z(:,8);
yC1 = Q5 + px3.*z(:,9) + py3.*z(:,7) - LD.*z(:,5);
xC2 = Q4 + LD.*z(:,6) + px4.*z(:,7) + py4.*z(:,8);
yC2 = Q5 + px4.*z(:,9) + py4.*z(:,7) - LD.*z(:,5);
xD1 = Q4 + px5.*z(:,5) - py5.*z(:,6);
yD1 = Q5 + px5.*z(:,6) + py5.*z(:,5);
xD2 = Q4 + px6.*z(:,5) - py6.*z(:,6);
yD2 = Q5 + px6.*z(:,6) + py6.*z(:,5);



cost_pos = ((xB1-M_ref.B1x).^2 + (yB1-M_ref.B1y).^2 + (xB2-M_ref.B2x).^2 + (yB2-M_ref.B2y).^2 + ...
    (xC1-M_ref.C1x).^2 + (yC1-M_ref.C1y).^2 + (xC2-M_ref.C2x).^2 + (yC2-M_ref.C2y).^2 + ...
    (xD1-M_ref.D1x).^2 + (yD1-M_ref.D1y).^2 + (xD2-M_ref.D2x).^2 + (yD2-M_ref.D2y).^2)*scale.length^2;
output.integrand = (cost_pos)*scale.cost;
output.dynamics = [Q1d*scale.angVel,Q2d*scale.angVel,Q3d*scale.angVel,Q4d*scale.vel,Q5d*scale.vel,...
    Q1d2*scale.angAccel,Q2d2*scale.angAccel,Q3d2*scale.angAccel,Q4d2*scale.accel,Q5d2*scale.accel];


load('referenceData.mat');
Q_ref = state_ref(:,1:2);
Qd_ref = state_ref(:,3:4);
Qd2_ref = ud_ref;
a_ref = state_ref(:,5:9);

Q1_ref = Q_ref(:,1);
Q2_ref = Q_ref(:,2);
Q1d_ref = Qd_ref(:,1);
Q2d_ref = Qd_ref(:,2);
Q1d2_ref = Qd2_ref(:,1);
Q2d2_ref = Qd2_ref(:,2);
Q3_ref = Q2_ref - Q1_ref;
Q3d_ref = Q2d_ref - Q1d_ref;
Q3d2_ref = Q2d2_ref - Q1d2_ref;
Q4_ref = LC*sin(Q1_ref) + LD*sin(Q1_ref-Q2_ref) - LB;
Q4d_ref = LC*cos(Q1_ref).*Q1d_ref + LD*cos(Q1_ref-Q2_ref).*(Q1d_ref-Q2d_ref);
Q4d2_ref = LC*cos(Q1_ref).*Q1d2_ref + LD*cos(Q1_ref-Q2_ref).*(Q1d2_ref-Q2d2_ref) - LC*sin(Q1_ref).*Q1d_ref.^2 - LD*sin(Q1_ref-Q2_ref).*(Q1d_ref-Q2d_ref).^2;
Q5_ref = HB + LC*cos(Q1_ref) + LD*cos(Q1_ref-Q2_ref);
Q5d_ref = -LC*sin(Q1_ref).*Q1d_ref - LD*sin(Q1_ref-Q2_ref).*(Q1d_ref-Q2d_ref);
Q5d2_ref = -LC*cos(Q1_ref).*Q1d_ref.^2 - LD*cos(Q1_ref-Q2_ref).*(Q1d_ref-Q2d_ref).^2 - LC*sin(Q1_ref).*Q1d2_ref - LD*sin(Q1_ref-Q2_ref).*(Q1d2_ref-Q2d2_ref);

r = auxdata.r;
rhoC = auxdata.rhoC;
rhoD = auxdata.rhoD;
R1A = abs(r*rhoC*cos(Q1_ref)./(r^2+rhoC^2+2*r*rhoC*sin(Q1_ref)).^0.5);
R2A = abs(r*rhoC*cos(Q1_ref)./(r^2+rhoC^2-2*r*rhoC*sin(Q1_ref)).^0.5);
R5A = abs(r*(LC*cos(Q1_ref)+rhoD*cos(Q1_ref-Q2_ref))./(LC^2+r^2+rhoD^2+2*LC*r*sin(Q1_ref)+2*LC*rhoD*cos(Q2_ref)+2*r*rhoD*sin(Q1_ref-Q2_ref)).^0.5);
R3K = abs(r*rhoD*cos(Q2_ref)./(r^2+rhoD^2-2*r*rhoD*sin(Q2_ref)).^0.5);
R4K = abs(r*rhoD*cos(Q2_ref)./(r^2+rhoD^2+2*r*rhoD*sin(Q2_ref)).^0.5);
R5K = abs(rhoD*(LC*sin(Q2_ref)+r*cos(Q1_ref-Q2_ref))./(LC^2+r^2+rhoD^2+2*LC*r*sin(Q1_ref)+2*LC*rhoD*cos(Q2_ref)+2*r*rhoD*sin(Q1_ref-Q2_ref)).^0.5);

TA_ref = F_ref(:,1).*R1A - F_ref(:,2).*R2A + F_ref(:,5).*R5A;
TK_ref = -F_ref(:,3).*R3K + F_ref(:,4).*R4K - F_ref(:,5).*R5K;


figure(1)
subplot(2,3,1);
plot(t,Q1,'b*',t_ref,Q_ref(:,1),'k-');
xlabel('Time (s)');
ylabel('q1 (rad)');
legend('output','reference','location','south');
subplot(2,3,2);
plot(t,Q2,'b*',t_ref,Q_ref(:,2),'k-');
xlabel('Time (s)');
ylabel('q2 (rad)');
subplot(2,3,3);
plot(t,Q3,'b*',t_ref,Q3_ref,'k');
xlabel('Time (s)');
ylabel('q3 (rad)');
subplot(2,3,4);
plot(t,Q4,'b*',t_ref,Q4_ref,'k');
xlabel('Time (s)');
ylabel('q4 (m)');
subplot(2,3,5);
plot(t,Q5,'b*',t_ref,Q5_ref,'k');
xlabel('Time (s)');
ylabel('q5 (m)');

figure(2)
subplot(2,3,1);
plot(t,Q1d,'b*',t_ref,Q1d_ref,'k-');
xlabel('Time (s)');
ylabel('u1 (rad/s)');
legend('output','reference','location','south');
subplot(2,3,2);
plot(t,Q2d,'b*',t_ref,Q2d_ref,'k-');
xlabel('Time (s)');
ylabel('u2 (rad/s)');
subplot(2,3,3);
plot(t,Q3d,'b*',t_ref,Q3d_ref,'k');
xlabel('Time (s)');
ylabel('u3 (rad/s)');
subplot(2,3,4);
plot(t,Q4d,'b*',t_ref,Q4d_ref,'k');
xlabel('Time (s)');
ylabel('u4 (m/s)');
subplot(2,3,5);
plot(t,Q5d,'b*',t_ref,Q5d_ref,'k');
xlabel('Time (s)');
ylabel('u5 (m/s)');

figure(3)
subplot(2,3,1);
plot(t,Q1d2,'b*',t_ref,Q1d2_ref,'k-');
xlabel('Time (s)');
ylabel('u1d (m/s)');
legend('output','reference','location','south');
subplot(2,3,2);
plot(t,Q2d2,'b*',t_ref,Q2d2_ref,'k-');
xlabel('Time (s)');
ylabel('u2d (m/s)');
subplot(2,3,3);
plot(t,Q3d2,'b*',t_ref,Q3d2_ref,'k');
xlabel('Time (s)');
ylabel('u3d (m/s)');
subplot(2,3,4);
plot(t,Q4d2,'b*',t_ref,Q4d2_ref,'k');
xlabel('Time (s)');
ylabel('u4d (m/s)');
subplot(2,3,5);
plot(t,Q5d2,'b*',t_ref,Q5d2_ref,'k');
xlabel('Time (s)');
ylabel('u5d (m/s)');

figure(4)
plot(t,TA,'r*',t_ref,TA_ref,'r',t,TK,'b*',t_ref,TK_ref,'b')
title('Joint Torques');
xlabel('Time (s)');
ylabel('Torque (Nm)');
legend('TA output','TA reference','TK output','TK refence');

figure(5)
subplot(1,2,1); hold off;
out = plot(t,xB1,'b*');
hold all;
plot(t,xB2,'b*',t,xC1,'b*',t,xC2,'b*',t,xD1,'b*',t,xD2,'b*');
ref = plot(t,M_ref.B1x,'k');
plot(t,M_ref.B2x,'k',t,M_ref.C1x,'k',t,M_ref.C2x,'k',t,M_ref.D1x,'k',t,M_ref.D2x,'k');
xlabel('Time (s)');
ylabel('Marker x positions (m)');
legend([out,ref],'output','reference');

subplot(1,2,2); hold off;
plot(t,yB1,'b*',t,yB2,'b*',t,yC1,'b*',t,yC2,'b*',t,yD1,'b*',t,yD2,'b*');
hold on;
plot(t,M_ref.B1y,'k',t,M_ref.B2y,'k',t,M_ref.C1y,'k',t,M_ref.C2y,'k',t,M_ref.D1y,'k',t,M_ref.D2y,'k');
xlabel('Time (s)');
ylabel('Marker y Positions (m)');
end