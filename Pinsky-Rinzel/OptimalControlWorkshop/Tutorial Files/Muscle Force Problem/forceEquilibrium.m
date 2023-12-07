function [error, FT] = forceEquilibrium(a,lMtilda,vMtilda,lMT,muscParams,FaParam,FvParam,FpParam,kT)

% Unpack and vectorize muscle parameters
npts = size(a,1);
FMo = repmat(muscParams(1,:),npts,1);
lMo = repmat(muscParams(2,:),npts,1);
lTs = repmat(muscParams(3,:),npts,1);
alphao = repmat(muscParams(4,:),npts,1);
% vMmax = repmat(muscParams(5,:),npts,1);

% Calculate muscle and tendon lengths
lM = lMtilda.*lMo;
w = lMo.*sin(alphao);

lT = lMT - sqrt(abs(lM.^2 - w.^2));

% Calculate tendon force
FT = kT*(lT-lTs);

% Extract normalized curve coefficient values
c1 = FaParam;

d1 = FvParam(1,1);
d2 = FvParam(2,1);

e1 = FpParam(1,1);
e2 = FpParam(2,1);
e3 = FpParam(3,1);

% Calculate muscle force
FMltilda = exp(-c1*(lMtilda-1).^2); % Active force-length
FMvtilda  = d1*atan(d2*vMtilda)+1; % Force-velocity

FMcetilda = a.*FMltilda.*FMvtilda; % Active normalized muscle force
FMpetilda = e1*(exp(e2*(lMtilda+e3).^3)); % Passive normalized muscle force

FM = FMo.*(FMcetilda+FMpetilda); % Muscle force

% Calculate error between muscle and tendon force
error =  FM.*(lMT-lT)./lM-FT;

end