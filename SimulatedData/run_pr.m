T = 2000; Is = 0.1; Id = 0; gc = 1.85;
params = [Is,Id,gc];
gs = [0.1, 30, 15, 10, 0.8, 15];
gtest = [0.0999335 29.9721 15.0022 9.9992 0.800086 14.9954];

params4 = [0.4, Id, gc];
gs4 = [0.110372 26.6248 20 13.3 0.735293 14.7863]; %500 iters --> 117k error
params5 = [0.4, Id, gc];
gs5 = [0.0979703 20 20 13.3 0.750437 14.6061]; %500 iters --> 197k error
%Is=0.6 blew up so not even worth considering. Above ones did not fit
%burst.

[tspan,y] = prsolve_rk4_sigmoid(T,params,gs);
%[~,ytest] = prsolve_rk4_sigmoid(T,params,gtest);
[~,y4] = prsolve_rk4_sigmoid(T,params4,gs4);
[~,y5] = prsolve_rk4_sigmoid(T,params4,gs5);
%y = readmatrix('prmod_Is0.1.txt');
%Vs = y(:,2);
%tspan = y(:,1);
Vs = y(:,1);
tspan = tspan';
intervals = get_intervals(tspan,Vs);
%intervals = [0, intervals, T];
%writematrix(intervals(:)','intervalsmore.txt','Delimiter',',')
%disp(intervals(2,:)-intervals(1,:))

frange = 1:10;
freqs = zeros(10,3);
for i = frange
    params(1) = i/10;
    [~,y] = prsolve_rk4_sigmoid(T,params,gs);
    intervals = get_intervals(tspan,y(:,1));
    freq = (1000/(intervals(end)-intervals(end-1)));
    freqs(i,1) = freq;
    
    [~,y] = prsolve_rk4_sigmoid(T,params,gs4);
    intervals = get_intervals(tspan,y(:,1));
    freq = (1000/(intervals(end)-intervals(end-1)));
    freqs(i,2) = freq;
    
    [~,y] = prsolve_rk4_sigmoid(T,params,gs5);
    intervals = get_intervals(tspan,y(:,1));
    freq = (1000/(intervals(end)-intervals(end-1)));
    freqs(i,3) = freq;
end
%out = [tspan, Vs, ones(length(tspan),1)]; %padding data
%writematrix(out,'prmod_Is0.1.txt','Delimiter','space');

%testing sigmoid vs piecewise interpolations for differences:
%differences seem to be very small if existent at all
%[~,ynotsig] = prsolve_rk4(T,params);
%Vtest = ynotsig(:,1);

set(0, 'defaultaxesfontsize',16,'defaultaxeslinewidth',2,...
       'defaultlinelinewidth',2,'defaultpatchlinewidth',2,...
       'defaulttextfontsize',16);

%plot(tspan,y(:,1),'b',tspan,y4(:,1),'r',tspan,y5(:,1),'g')
plot(freqs(:,1),'rx')
hold on
plot(freqs(:,2),'bx')
plot(freqs(:,3),'gx')
xlabel('Is (\umA/cm^2)')
ylabel('Frequency (Hz)')
title('Spike/Burst Frequency of PR Models at Various Is')
legend('Is=0.1','Is=0.4','Is=0.5')
