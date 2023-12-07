function dadt = activationDynamics(e,a,tact,tdeact)

% Winters_continuous
d1 = 1./(tact*(0.5+1.5*a));
d2 = (0.5+1.5*a)/tdeact;
b = 0.1;
f = 0.5*tanh(b*(e-a));
dadt = (d1.*(f+0.5) + d2.*(-f+0.5)).*(e-a);
end