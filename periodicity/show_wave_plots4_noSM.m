function SameDMoutput  = show_wave_plots4 (file_cord,binCount,binSize,motifSize,flagC2,flagCR)
%function [P,phase,amp,hl,ssd,chi2,gof,count,RealP]  = periodicity_matlab_v4_func (file_cord,binCount,binSize,motifSize,flagC2,flagCR)
% 8.13.2015 build from 3. xxjr. allow rand (green)

% Readme:
% file_cord: The file name for the motif coordinates (First line is the genome size)
% binCount	(default 200)
% binSize	(default 500)
% motifSize	(default 8)
% flagC2	(default 0. c2 control motif sets, coordinates in the input file)
% flagCR	(default 0. will do the randomization)

%clear all;close all;
set(0,'DefaultAxesFontSize',14);

startlag = 2; %%%% Why i do this ???? 8.10.2015


%file_cord = 'naz_NC_014248_m0.data';N = 5354700;
%file_cord = 'syg_NC_007776_m0.data';N = 3046682;
%file_cord = 'HIP1cord.NC_019776.m0.Cyanobacterium_aponinum_PCC_10605_uid183340.data.txt'; N = 4114099;
%file_cord = 'result_new_HIP1_cords\cord.m0.Cyanobacterium_aponinum_PCC_10605_uid183340.NC_019776.data.txt';
%file_cord = 'result_new_c2_cords\cord.c2.Cyanobacterium_aponinum_PCC_10605_uid183340.NC_019776.data.txt';
%binCount = 200;
%binSize = 6000;
%motifSize = 8;

N_c = load(file_cord);
N = N_c(1);
c = N_c(2:end);

if (flagCR ~= 0)
	rng(flagCR);
	%size(c,1)
	c = sort(randsample(N,size(c,1))); %%%% RANDOMIZATION of CORDINATES %%%%%%%%%%%%
end

allDist = getCircularDist(c,N);
allDistInRange = allDist(allDist<=(binCount+0)*binSize);
allDist2 = pdist(c);
allDistInRange2 = allDist2(allDist2<=(binCount+0)*binSize);

[histHeight,histCenter] = hist(allDistInRange,binCount);
[histHeight2,histCenter2] = hist(allDistInRange2,binCount);

binCount = binCount - 1; %%%% WHAT IS THIS !!!!????? %%%%%%%%%%%%%
histHeight = histHeight(2:end); %%%% WHAT IS THIS !!!!????? %%%%%%%%%%%%%
histHeight2 = histHeight2(2:end); %%%% WHAT IS THIS !!!!????? %%%%%%%%%%%%%

%%%%%%%% INITIAL RED PLOTS (Inter-motif Distances) %%%%%%%%%%%%%%%%
h1 = figure(1);
set(h1, 'Position', [1500 200 1600 700]);
subplot(4,2,1); % circular distances
plot(1:binCount,histHeight,'r','LineWidth',2); 
ylim([min(histHeight),max(histHeight)]);
xlabel('Bin Number');ylabel('# Observed');
subplot(4,2,3); % circular distances. with smoothing.
%histHeight_smooth = smooth(histHeight,0.03,'loess');
histHeight_smooth = smooth(histHeight,3,'moving');
plot(1:binCount,histHeight_smooth,'r','LineWidth',2); 
ylim([min(histHeight),max(histHeight)]);
xlabel('Bin Number');ylabel('# Observed');
%noSM% histHeight = histHeight_smooth';
subplot(4,2,5); % chromosome assumed to be linear.
%allDistInRange = allDistInRange-motifSize; allDistInRange=allDistInRange(allDistInRange>0);
plot(1:binCount,histHeight2,'r','LineWidth',2); 
ylim([min(histHeight2),max(histHeight2)]);
xlabel('Bin Number');ylabel('# Observed');
xxx = 1:binCount; yyy = histHeight2;ppp = polyfit(xxx,yyy,1);yyy1 = polyval(ppp,xxx);
hold on;plot(xxx,yyy1,'b');hold off;
subplot(4,2,7); % un-skrewe for 4,1,3.
yyy2 = yyy - yyy1;
plot(1:binCount,yyy2,'r','LineWidth',2); 

%subplot(4,2,2); % circular distances
%plot(1:binCount,histHeight,'r','LineWidth',2); 
%ylim([min(histHeight),max(histHeight)]);
%xlabel('Bin Number');ylabel('# Observed');
subplot(4,2,2); % circular distances. with smoothing.
%histHeight_smooth = smooth(histHeight,0.03,'loess');
xxx = 1:binCount;yyy = histHeight_smooth';
ppp = polyfit(xxx,yyy,1);yyy2 = polyval(ppp,xxx);
plot(xxx,yyy,'r','LineWidth',2); 
ylim([min(histHeight),max(histHeight)]);
xlabel('Bin Number');ylabel('# Observed');
hold on;plot(xxx,yyy2,'b');hold off;

subplot(4,2,4); % circular distances. with smoothing.
%histHeight_smooth = smooth(histHeight,0.03,'loess');
xxx = 1:binCount;yyy22 = yyy - yyy2 + mean(yyy);
plot(xxx,yyy22,'r','LineWidth',2); 
ylim([min(histHeight),max(histHeight)]);
xlabel('Bin Number');ylabel('# Observed');


%%%%%%%% Blue PLOTS (Inter-distance) & Fourier Trans %%%%%%%%%%%%%%%%
h2 = figure(2);
subplot(3,1,1);
%L = binCount;
%Fs = binCount;
t = 1:binCount;
y = histHeight;
plot(t,y);
xlabel('Time-domain(t)');ylabel('y(t)');
subplot(3,1,2);
y = y - mean(y);
Y = fft(y);
%f = (0:binCount-1)*1/binCount;
f = (0:binCount-1)*1/binCount;
plot(f,abs(Y));xlim([0,f(floor(length(f)/2))]);
%title('Single-Sided Amplitude Spectrum of y(t)');
xlabel('Frequency-domain(f)');ylabel('|Y(f)|');
subplot(3,1,3);
f_invs = 1./f(2:end);
plot(f_invs,abs(Y(2:end)));%xlim([0,length(x)/2]);
xlim([1,max(f_invs)/4]);
xlabel('Periodicity (1/f)');ylabel('|Y(1/f)|');
[maxY,maxYind] = max(abs(Y(2:end)));
inferred_periodicity_Histo = f_invs(maxYind);
inferred_real_periodicity_Histo = inferred_periodicity_Histo * binSize;
title_str = sprintf('Periodicity (1/f). Peak at %f',inferred_periodicity_Histo);%title(title_str);
xlabel(title_str);


%%%%%%%%% Main PLOTS (ACC and fitting)  %%%%%%%%%%%%%%%%
h3 = figure(3);
set(h3, 'Position', [1700 100 1000 800]);
ylimv = 0.8;
subplot(2,2,1);
%autocorr(y,floor(binCount/2));
[ACF,lags,bounds] = autocorr(histHeight,floor(binCount/2));
plot(lags(startlag:end),ACF(startlag:end));
x = lags(startlag:end);
y = ACF(startlag:end);
p = polyfit(x,y,2);
y1 = polyval(p,x);
hold on;plot(x,y1,'m');hold off;
ylim([-1*ylimv,ylimv]);
subplot(2,2,2);
%autocorr(y,floor(binCount/2));
 [ACF,lags,bounds] = autocorr(histHeight,floor(binCount/2));
 plot(lags(startlag:end),ACF(startlag:end));
 x = lags(startlag:end);
 y = ACF(startlag:end);
 p = polyfit(x,y,2);
 y1 = polyval(p,x);
hold on;plot(x,y1,'m');
%grid;
[sigr resr]=fit_damped_sinewave_my5b_JL(y);
plot(x,sigr,'g','LineWidth',2);hold off;
grid;ylim([-1*ylimv,ylimv]);
subplot(2,2,3);
 y2 = y - y1;
plot(x,y2);
p2 = polyfit(x,y2,2);
y21 = polyval(p2,x);
hold on;plot(x,y21,'m');hold off;
ylim([-1*ylimv,ylimv]);
subplot(2,2,4);
plot(x,y2);hold on;
p2 = polyfit(x,y2,2);
y21 = polyval(p2,x);
hold on;plot(x,y21,'m');%hold off;
   y2 = y2 - mean(y2);
[sigr resr]=fit_damped_sinewave_my5b_JL(y2);
plot(x,sigr,'g','LineWidth',2);hold off;
grid;ylim([-1*ylimv,ylimv]);
SameDMoutput = [resr length(c) binSize*resr(1,1)];

 return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%filename = strcat(parts1{1},'-', parts2{1} ,'(',parts2{2},').overlay2.png');
%set(h1, 'PaperUnits', 'inches');
%x_width=10;y_width=3;%x_width=13;y_width=6.5;
%set(h1, 'PaperPosition', [0 0 x_width y_width]);
%print(h1, '-dpng', filename,'-r100');



%%%%%%%%%%%%%%%%% 
%%%%%%  Blue Plots (Post-acf fourier transformation  %%%%%%%%%%%%%
% h21 = figure(21);
% subplot(3,1,1);
% plot(x,y2);
% xlabel('Time-domain(t)');ylabel('y(t)');
% subplot(3,1,2);
 y2 = y2 - mean(y2);
 Y2 = fft(y2);
% %f = (0:binCount-1)*1/binCount;
 f = (0:length(x)-1)/length(x);
% plot(f,abs(Y2));xlim([0,f(floor(length(x)/2))]);
% %title('Single-Sided Amplitude Spectrum of y(t)');
% xlabel('Frequency-domain(f)');ylabel('|Y(f)|');
% subplot(3,1,3);
 f_invs = 1./f
% plot(f_invs,abs(Y2));%xlim([0,length(x)/2]);
% xlim([1,max(f_invs(2:end))/2]);
% xlabel('Periodicity (1/f)');ylabel('|Y(1/f)|');

 [maxY2,maxY2ind] = max(abs(Y2));
 inferred_periodicity = f_invs(maxY2ind);
 inferred_real_periodicity = inferred_periodicity * binSize


 


