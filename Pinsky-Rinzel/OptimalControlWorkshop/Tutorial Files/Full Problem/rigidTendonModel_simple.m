function [FT,lMtilda,FMltilda,vMtilda,FMvtilda] = rigidTendonModel_simple(a,lMT,vMT,params,cCoefs,dCoefs,eCoefs)

% Extract muscle model parameter values
FMo = params(1,:);
lMo = params(2,:);
lTs = params(3,:);
alphao = params(4,:);
vMmax = params(5,:);

nMusc = size(FMo,2);

% Extract normalized curve coefficient values
% cCoefs = [c1];
c1 = cCoefs;

% dCoefs = [d1; d2];
d1 = dCoefs(1,1);
d2 = dCoefs(2,1);

% eCoefs = [e1; e2; e3];
e1 = eCoefs(1,1);
e2 = eCoefs(2,1);
e3 = eCoefs(3,1);

% Calculate kinematic and geometric inputs
for m = 1:nMusc
    lM = sqrt((lMo*sin(alphao)).^2+(lMT-lTs).^2);
%     lM(:,m) = lMT(:,m)-lTs(m);
%     alpha(:,m) = asin(lMo(m)*sin(alphao(m))./lM(:,m));
    alpha(:,m) = acos((lMT(:,m)-lTs(m))./lM(:,m));
    vM(:,m) = vMT(:,m).*cos(alpha(:,m));

    % Calculate normalized kinematic inputs
    lMtilda(:,m) = lM(:,m)/lMo(m);
    vMtilda(:,m) = vM(:,m)/vMmax(m);

    % Calculate locations on normalized curves
    FMltilda(:,m) = exp(-c1*(lMtilda(:,m)-1).^2); % Active force-length
    FMvtilda(:,m)  = d1*atan(d2*vMtilda(:,m))+1; % Force-velocity
    FMpetilda(:,m) = e1*(exp(e2*(lMtilda(:,m)+e3).^3)); % Passive force-length

    % Calculate normalized forces
    FMcetilda(:,m) = a(:,m).*FMltilda(:,m).*FMvtilda(:,m);
    FMtilda(:,m) = FMcetilda(:,m)+FMpetilda(:,m);
    FTtilda(:,m) = FMtilda(:,m).*cos(alpha(:,m));

    % Calculate tendon force output
    FT(:,m) = FMo(m)*FTtilda(:,m);
end

return