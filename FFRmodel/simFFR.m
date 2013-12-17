clear;
clc;
addpath('./Zilany2014_new/');
fs = 100e3;
fc = 4000;
fm = 100;
m = 0.6;
bw = 200;
dur = 2.0;
rise = 0.025;
phi = 0;
isi = 0.025;
SAM = 1;
notch = 1;
modulated = 0;

plotting = 0;
fiblist = {'Low-SR','Med-SR','High-SR'};
for fibertype = 1:3; % 1- LS, 2- MS, 3-HS
    fib = fiblist{fibertype};
    fprintf(1,'\n ------------Simulating %s fibers :)-----------\n',fib);
    SR = [0.1, 15, 100];
    stimdb = 70;
    
    
    
    synout = 0;
    Ric = 0;
    Rcn = 0;
    
    Ntrials = 8;
    if(SAM)
        pin = makesamtone(fc,fm,m,fs,dur,rise,phi);
    else
        pin = maketranstone(fc,fm,m,bw,fs,dur,rise,phi);
    end
    
    stimrms = sqrt(mean(pin.^2));
    
    for trial = 1:Ntrials
        fprintf (1,'\n########### Doing Trial # %d/%d ############\n',...
            trial, Ntrials);
        if(notch)
            SNR = 30;
            noise = makeNotchNoiseFFT(800,0,fc,dur,fs,rise,0)';
            noiserms = sqrt(mean(noise.^2));
            noise = (noise*stimrms/noiserms)*db2mag(-1*SNR);
            
            if(modulated)
                fmask = 130;
                tmask = (0:(numel(noise)-1))/fs;
                noise = (1+sin(2*pi*fmask*tmask)).*noise;
            end
            pin = pin + noise;
        end
        
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
            nrep = 2;
            vihc = model_IHC(pin,CF,nrep,1/fs,dur+isi,1,1,tuningType);
            [synout_trial(nCF,:),synoutvar,psth] = ...
                model_Synapse(vihc,CF,nrep,1/fs,fibertype,1,0);
            
            t = (0:(size(synout,2)-1))/fs;
            [Rcn_trial(nCF,:),Ric_trial(nCF,:)] = ...
                NelsonCarney2004CNIC(synout_trial(nCF,:),fs);
            
            if(nCF == 1)
                fprintf(1,'Total number of CFs = %d\n',numel(f));
                fprintf(1,'Done with CF #1');
            elseif(mod(nCF,10) == 0)
                fprintf(1,'%d',nCF);
            else
                fprintf(1,'.');
            end
            if(mod(nCF,50) == 0)
                fprintf(1,'\n');
            end
        end
        synout = synout + synout_trial;
        Rcn = Rcn + Rcn_trial;
        Ric = Ric + Ric_trial;
    end
    fprintf(1,'\n');
    
    Ric = Ric/Ntrials;
    Rcn = Rcn/Ntrials;
    synout = synout/Ntrials;
    
    
    % Gate out the onset response foor visualization
    t_gated = t(t>0.050);
    resp = mean(Ric,1);
    [Fresp,fplot] = pmtm(resp(t>0.05),1.25,[],fs);
    
    if(plotting)
        fprintf(1, '\nAll done! Plotting... Hold on!\n');
        
        % Image plotting requires uniform grid, hence interpolating
        [tgrid,fgrid_orig] = meshgrid(t,f);
        
        
        
        [tgrid_new, fgrid] = meshgrid(t_gated,f_uniform);
        synout_new = interp2(tgrid,fgrid_orig,synout, tgrid_new, fgrid,'spline');
        Ric_new = interp2(tgrid,fgrid_orig,Ric, tgrid_new, fgrid,'spline');
        Rcn_new = interp2(tgrid,fgrid_orig,Rcn, tgrid_new, fgrid,'spline');
        
        
        figure;
        subplot(2,2,1);
        imagesc(t_gated,f_uniform,synout_new,...
            [1.2*SR(fibertype), max(max(synout_new))]);
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
        subplot(2,1,1);
        
        plot(t_gated, resp,'linewidth',2);
        xlabel('Time (s)','FontSize',20);
        subplot(2,1,2);
        semilogy(fplot,Fresp,'linewidth',2); xlim([5, 1000]);
        xlabel('Frequency (Hz)','FontSize',20);
    end
    
    pop_resp(fibertype,:) = resp;
    pop_resp_freq(fibertype,:) = Fresp;
end
saveResults = 1;

if(saveResults)
    fname = strcat('RESULTS/SimFFR_results_SNR_',num2str(SNR),...
        '_m',num2str(m*100),'_',num2str(stimdb),'dB_',datestr(clock,'mmm_dd_yyyy_HH_MM'),'.mat');
    save(fname,'t_gated','fplot','pop_resp','m', 'stimdb','SNR','pop_resp_freq','Ntrials');
end

c = 0;
cols = 'rgbyck';
ind = (fplot > 50) & (fplot < 500);
for damage = 0:0.2:1
    w = [(1-damage)*0.15, (1-damage)*0.25, 0.6];
    mixedSR = w*pop_resp;
    c = c+1;
    [pxx,f] = pmtm(mixedSR,2,[],1/mean(diff(t_gated)));
    plot(f(ind),pxx(ind),cols(c),'linew',2);
    xlim([5, 500]);
    hold on;
end