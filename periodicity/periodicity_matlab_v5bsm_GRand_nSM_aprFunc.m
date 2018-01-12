function SameDMoutput  = periodicity_matlab_v5bsm_GRand_nSM_aprFunc (c_GRand_allDist,GenomeLen,nMotif,binCount,binSize)
%tic %4_8_2016
%global c_GRand_allDist;

startlag = 2;

allDist = c_GRand_allDist;

allDistInRange = allDist(allDist<=(binCount+0)*binSize);

 % toc %4_8_2016

[histHeight,histCenter] = hist(allDistInRange,binCount);

binCount = binCount - 1; %%%% WHAT IS THIS !!!!????? %%%%%%%%%%%%%
histHeight = histHeight(2:end); %%%% WHAT IS THIS !!!!????? %%%%%%%%%%%%%
%histHeight2 = histHeight2(2:end); %%%% WHAT IS THIS !!!!????? %%%%%%%%%%%%%


 [ACF,lags,bounds] = autocorr(histHeight,floor(binCount/2));

  % toc %4_8_2016
 % plot(lags(2:end),ACF(2:end));
 x = lags(startlag:end);
 y = ACF(startlag:end);
 p = polyfit(x,y,2);
 y1 = polyval(p,x);

y2 = y - y1;
y2 = y2 - mean(y2);
%[sigr resr]=fit_damped_sinewave_my5b_JL_apr(y2);



%% copy from script function 'fit_damped_sinewave_my5b_JL_apr.m'
sig = y2;
x = sig;
t = 1:length(x);
expsin3 = @(a, prd, phi, hl)a * sin((2*pi/prd) * t + phi) .* exp(-log(2)/hl * t);
GoFfunc3 = @(a, prd, phi, hl)sum((x-expsin3(a, prd, phi, hl)).^2);

A = max(sig);          % gain . [0.8-1.2] * height of first peak
P =  10;		 % period [ 7 -15]
f = 1/P;         % frequency
phi = pi;     % phase angle l [0, 2*pi]
hl = 0;			% half life  [0, 4]
lambda = 0.01;% time constant
omega = 2 * pi * f; % angular freq

As   = linspace(0.1*A, 1.5*A, 14);
Ps   = linspace(5,     17,    17);
Phis = linspace(0,     2*pi,  12);
hls  = linspace(5,  40,     8);

omegas  = 2*pi ./ Ps;
lambdas = log(2) ./ hls;

%[x1,x2,x3,x4] = ndgrid(As, omegas,Phis,lambdas);
[x1,x2,x3,x4] = ndgrid(As, Ps,Phis,hls);
%tic
DDD = arrayfun(GoFfunc3,x1,x2,x3,x4);
%toc
[min_v,min_i] = min(DDD(:));
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
SSD = sum((x-xr).^2);
ChiSq = sum((x-xr).^2./xr);
GoF = sum((x-xr).^2./var(xr));

resr = [Ps2(mii2) Phis2(mii3) As2(mii1) hls2(mii4) SSD ChiSq GoF];
%%% end of 'copy from script function'

SameDMoutput = [resr nMotif binSize*resr(1,1)];

return


 


