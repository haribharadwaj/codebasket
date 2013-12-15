clear;
clc;
addpath('/home/hari/Models/Zilany2014_new/');
fs = 100e3;
fm = 1000;
m = 0.3;
bw = 200;
dur = 0.4;
rise = 0.025;
phi = 0;
isi = 0.025;

t = 0:(1/fs):(dur -1/fs);
pin = (cos(2*pi*fm*t) > 0.999);
fibertype = 1; % 1- LS, 2- MS, 3-HS

SR = [0.1, 15, 100];
stimdb = 80;
stimrms = sqrt(mean((cos(2*pi*fm*t)).^2));

tuningType = 2;

% 20e-6 RMS Pa is 0 dB SPL..
% pin/stimrms has RMS of 1..
% 20e-6*pi/stimrms has 0 dB SPL
% Hence to get desired dB SPL:
pin = 10^(stimdb/20)*20e-6*pin/stimrms;

load vFreq;
CF_step = 10;
f = vFreq(15:CF_step:end);
for nCF = 1:numel(f)
    CF = f(nCF);
    nrep = 1;
    vihc = model_IHC(pin,CF,nrep,1/fs,dur+isi,1,1,tuningType);
    [synout(nCF,:),synoutvar,psth] = model_Synapse(vihc,CF,nrep,1/fs,fibertype,1,0);
   
    t = (0:(size(synout,2)-1))/fs;
    [Rcn(nCF,:),Ric(nCF,:)] = NelsonCarney2004CNIC(synout(nCF,:),fs);
    fprintf(1,'Done with CF # %d / %d\n',nCF,numel(f));
end
f_step = 25;
f_uniform = min(f):f_step:max(f);

fprintf(1, 'All done! Plotting... Hold on!\n');

% Image plotting requires uniform grid, hence interpolating
[tgrid,fgrid_orig] = meshgrid(t,f);

% Gate out the onset response foor visualization
t_gated = t(t>0.050); 

[tgrid_new, fgrid] = meshgrid(t_gated,f_uniform);
synout_new = interp2(tgrid,fgrid_orig,synout, tgrid_new, fgrid,'spline');
Ric_new = interp2(tgrid,fgrid_orig,Ric, tgrid_new, fgrid,'spline');
Rcn_new = interp2(tgrid,fgrid_orig,Rcn, tgrid_new, fgrid,'spline');
figure;
subplot(2,2,1);
imagesc(t_gated,f_uniform,synout_new, [1.2*SR(fibertype), max(max(synout_new))]);
xlabel('Time (s)','FontSize',20);
ylabel('CF (Hz)','FontSize',20);
title('AN Output','FontSize',20);

subplot(2,2,2);
imagesc(t_gated,f_uniform, Ric_new);
xlabel('Time (s)','FontSize',20);
ylabel('CF (Hz)','FontSize',20);
title('IC MF Cells','FontSize',20);

subplot(2,2,3);
Ric_BR = ((Rcn_new - Ric_new) + abs(Rcn_new - Ric_new))/2;
imagesc(t_gated,f_uniform, Ric_BR);
xlabel('Time (s)','FontSize',20);
ylabel('CF (Hz)','FontSize',20);
title('IC Band Reject Cells','FontSize',20);

subplot(2,2,4);
plot(f_uniform,mean(synout_new,2),'k--','linew',2);
hold on;
plot(f_uniform,mean(Ric_new,2),'b','linew',2);
hold on;
plot(f_uniform,mean(Ric_BR,2),'r','linew',2);
ylabel('Average rate','FontSize',20);
xlabel('CF (Hz)','FontSize',20);

figure;
plot(t,mean(Ric)/max(mean(Ric)),'linewidth',2);
hold on;
plot((0:(numel(pin)-1))/fs, pin/max(pin),'r--','linewidth',2) 