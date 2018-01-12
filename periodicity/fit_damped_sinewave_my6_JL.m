function [xr resr]=fit_damped_sinewave_my6_JL(sig)
% v4. just pure fit.
% v6. 2.23.2016. use the JL's GoF concept.

%%%resr = [Ar alphar wr phr]
%resr = [Ar alphar wr phr]
%%%% MINLI_XU ADDED at 6/5/2015
%%%% http://en.wikipedia.org/wiki/Damped_sine_wave
%[A alpha w ph]
% A :  initial amplitude
% alpha : (negative)  decay constant. (also -alpha = lambda)
% w (omega): Angular frequency.
% ph : phase angle
% FREQUENCY  : f = w/(2*pi);
% PERIOD = 1/f
% LAMBDA = -alpha
% HALF-LIFE = log(2)/ LAMBDA



%clear all;close all;
%sig = [0.6801,0.3466, -0.0038, -0.1561, -0.2340, -0.2815, -0.3185, -0.2832, -0.1552,0.0556,0.2278,0.3044,0.2898,0.2227,0.1249,0.0036, -0.1202, -0.2148, -0.2766, -0.2982, -0.2455, -0.1630, -0.0452,0.0187,0.0566,0.0541,0.0483,0.0566,0.0651,0.0324, -0.0463, -0.1710, -0.2263, -0.2204, -0.1532, -0.0940, -0.0540, -0.0158,0.0670,0.1491,0.2227,0.2205,0.1865,0.1432,0.1040,0.0467, -0.0335, -0.1086, -0.1146, -0.0572,0.0259,0.1072,0.1508,0.1862,0.2075,0.2101,0.1753,0.1173,0.0338, -0.0438, -0.1401, -0.2013, -0.2222, -0.1763, -0.0928, -0.0124,0.0658,0.1381,0.1898,0.2024,0.1308,0.0242, -0.0796, -0.1423, -0.1623, -0.1950, -0.2048, -0.1664, -0.0574,0.1020,0.2111,0.2122,0.1251,0.0119, -0.0390, -0.0524, -0.0774, -0.1300, -0.1831, -0.1805, -0.1107, -0.0252,0.0529,0.0899,0.0967,0.0988,0.0690,0.0466, -0.0234];
x = sig;
t = 1:length(x);

% Define curve model functions
%expsin = @(a, omega, phi, lambda, t)a * sin(omega * t + phi) .* exp(-lambda * t);
%expsin2 = @(a, omega, phi, lambda)a * sin(omega * t + phi) .* exp(-lambda * t);
%lsqexpsin = @(p, t)expsin(p(1), p(2), p(3), p(4), t);
%GoFfunc = @(a, omega, phi, lambda)sum((x-expsin(a, omega, phi, lambda,t)).^2./var(expsin(a, omega, phi, lambda,t)));
%GoFfunc2 = @(a, omega, phi, lambda)sum((x-expsin(a, omega, phi, lambda,t)).^2);


expsin3 = @(a, prd, phi, hl)a * sin((2*pi/prd) * t + phi) .* exp(-log(2)/hl * t);
%GoFfunc3 = @(a, prd, phi, hl)sum((x-expsin3(a, prd, phi, hl)).^2);
%GoFfunc3 = @(a, prd, phi, hl)  abs(  sum(   ((x-expsin3(a, prd, phi, hl)).^2)   ./ abs(x)  )  ); % 2.26.2016 GoF change JL
GoFfunc3 = @(a, prd, phi, hl)  abs(  sum(   ((x-expsin3(a, prd, phi, hl)).^2)   ./ abs(expsin3(a, prd, phi, hl))  )  ); % 2.26.2016 GoF change JL


%x


% Setup data params
A = max(sig);          % gain . [0.8-1.2] * height of first peak
P =  10;		 % period [ 7 -15]
f = 1/P;         % frequency
phi = pi;     % phase angle l [0, 2*pi]
hl = 0;			% half life  [0, 4]
lambda = 0.01;% time constant
omega = 2 * pi * f; % angular freq


%n = 200000;
%n1 = 20;
%A = 1;
%lb = [0.8*A, 7,  0,   0.01];  % A, P, phi, hl
%ub = [1.2*A, 15, 2*pi,4   ];  % A, P, phi, hl

rf1 = 1; %resolution factor
wf1 = 0; %search width factor
As   = linspace(0.1*A*(1-wf1),	1.5*A*(1+wf1), 14*rf1);
Ps   = linspace(0,				25*(1+wf1),    70*rf1);
Phis = linspace(0,				2*pi*(1+wf1),  12*rf1);
hls  = linspace(5*(1-wf1),		40*(1+wf1),     8*rf1);

omegas  = 2*pi ./ Ps;
lambdas = log(2) ./ hls;

%[x1,x2,x3,x4] = ndgrid(As, omegas,Phis,lambdas);
[x1,x2,x3,x4] = ndgrid(As, Ps,Phis,hls);
%tic
DDD = arrayfun(GoFfunc3,x1,x2,x3,x4);
%toc
[min_v,min_i] = min(DDD(:));
%min_v
[mi1,mi2,mi3,mi4] = ind2sub(size(DDD),min_i);
min_sse = DDD(mi1,mi2,mi3,mi4);
min_para = [As(mi1), Ps(mi2),Phis(mi3),hls(mi4)];

As2   = linspace(0.9*As(mi1), 1.1*As(mi1), 11);
Ps2   = linspace(Ps(mi2)-0.5, Ps(mi2)+0.5, 11);
Phis2 = linspace(Phis(mi3)-pi/10, Phis(mi3)+pi/10, 11);
hls2  = linspace(0.8*hls(mi4), 1.2*hls(mi4), 9);

[xx1,xx2,xx3,xx4] = ndgrid(As2, Ps2,Phis2,hls2);
%tic
DDD2 = arrayfun(GoFfunc3,xx1,xx2,xx3,xx4);
%toc
[min_vv,min_ii] = min(DDD2(:));
[mii1,mii2,mii3,mii4] = ind2sub(size(DDD2),min_ii);
min_sse2 = DDD2(mii1,mii2,mii3,mii4);
min_para2 = [As2(mii1), Ps2(mii2),Phis2(mii3),hls2(mii4)];

x = sig;
t = 1:length(x);
expsin55 = @(a, prd, phi, hl)a * sin((2*pi/prd) * t + phi) .* exp(-log(2)/hl * t);
xr = expsin55(min_para2(1),min_para2(2),min_para2(3),min_para2(4));
%xr = expsin55(min_para(1),min_para(2),min_para(3),min_para(4));%%experiment 2.23.2016
SSD = sum((x-xr).^2);
ChiSq = sum((x-xr).^2./xr);
GoF = sum((x-xr).^2./var(xr));
%GoF2 = sum(   ((x-xr).^2 ./ abs(xr))    ) / var(x); % 2.26.2016

resr = [Ps2(mii2) Phis2(mii3) As2(mii1) hls2(mii4) SSD ChiSq GoF];

return

% ffgggg
% 
% %tic
% %AAA = arrayfun(GoFfunc,As,omegas,Phis,lambdas);
% %toc
% 
% tic
% BBB = arrayfun(GoFfunc2,As,omegas,Phis,lambdas);
% toc
% 
% tic
% CCC = arrayfun(GoFfunc2,x1,x2,x3,x4);
% toc
% 
% [min_v,min_i] = min(CCC(:));
% [mi1,mi2,mi3,mi4] = ind2sub(size(CCC),min_i);
% min_sse = CCC(mi1,mi2,mi3,mi4)
% 
% min_para = [As(mi1), omegas(mi2),Phis(mi3),lambdas(mi4)]
% 
% 
% %  1.0263    0.4756    1.6535    0.1733 =>  1.6976
% %  0.3523    0.4546    1.7080    0.0125 =>  0.6626
% 
%  arrayfun(GoFfunc2,1.0263,0.4756,1.6535,0.1733);
% 
% arrayfun(GoFfunc2,0.3523,0.4546,1.7080,0.0125);

% for i1 = 1:n
% 	i1
% 	for i2 = 1:n
% 		i2
% 		for i3 = 1:n
% 			for i4 = 1:n
% 				xr = expsin(As(i1),omegas(i2),Phis(i3),lambdas(i4),t);
% 				%A
% 				GoF = sum((x-xr).^2./var(xr));
% 				if (GoF < val_min)
% 					val_min = GoF;
% 					ind_min = [i1,i2,i3,i4];
% 				end
% 				%size(A)
% 				%sum(A)
% 				%GoFfunc(As(i1),omegas(i2),Phis(i3),lambdas(i4))
% 			end
% 		end
% 	end
% end

%close all; figure; hold on;
%plot(t, x, 'k-', 'LineWidth', 2);

% spd=fft(x);
% [mv,mind]=max(abs(spd(1:floor(end/2))));
% fEstimate = mind/ length(t);
% omegaEstimate = 2 * pi * fEstimate;
%plot(t./length(x),abs(spd));

fff
% Fit model to data
init = [A, omega, phi, lambda];
lb = [  0*max(sig), 2*pi*(1/17), 0   , 0 ];
ub = [ 1.5*max(sig), 2*pi*(1/5) , 4*pi, log(2)/4];
opts = optimset('Display','off','TolFun',1e-5);
[newparams, err] = lsqcurvefit(lsqexpsin, init, t, x, lb, ub, opts)

%plot(t, lsqexpsin(newparams, t));



xr = lsqexpsin(newparams, t);
SSD = sum((x-xr).^2);
ChiSq = sum((x-xr).^2./xr);
%GoF = sum((x-xr).^2./var(xr))
GoF = sum(   ((x-xr).^2 ./ x)    ) / var(x); % 2.26.2016


A_r = newparams(1);
w_r = newparams(2);
phi_r = newparams(3);
lambda_r = newparams(4);
f_r = w_r / (2*pi);
period_r = 1 / f_r;
hl_r = log(2) / (lambda_r);

resr = [period_r phi_r A_r hl_r SSD ChiSq GoF];


return

