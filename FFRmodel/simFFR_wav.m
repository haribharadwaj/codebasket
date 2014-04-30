clear;
clc;
addpath('./Zilany2014_new/');

wavname = 'four_speakers.wav';
fs = 100e3;

plotting = 1;
fiblist = {'Low-SR','Med-SR','High-SR'};
for fibertype = 1:3; % 1- LS, 2- MS, 3-HS
    fib = fiblist{fibertype};
    fprintf(1,'\n ------------Simulating %s fibers :)-----------\n',fib);
    SR = [0.1, 15, 100];
    stimdb = 75;
    
    
    
    synout = 0;
    Ric = 0;
    Rcn = 0;
    
    Ntrials = 1;
    [pin, fs_wav] = wavread(wavname);
    
    pin = resample(pin,fs,fs_wav)';
    
    dur = size(pin,2)/fs;
    isi = 50e-3; %ms
    
    stimrms = sqrt(mean(pin.^2));
    
    for trial = 1:Ntrials
        fprintf (1,'\n########### Doing Trial # %d/%d ############\n',...
            trial, Ntrials);
        
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
            
            t = (0:(size(synout_trial,2)-1))/fs;
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
    
    if(plotting)
        fprintf(1, '\nAll done! Plotting... Hold on!\n');
        
        % Image plotting requires uniform grid, hence interpolating
        [tgrid,fgrid_orig] = meshgrid(t,f);
        
        f_uniform = linspace(min(f),max(f),100);
        
        [tgrid_new, fgrid] = meshgrid(t_gated,f_uniform);
        synout_new = interp2(tgrid,fgrid_orig,synout, tgrid_new, fgrid,'spline');
        Ric_new = interp2(tgrid,fgrid_orig,Ric, tgrid_new, fgrid,'spline');
        Rcn_new = interp2(tgrid,fgrid_orig,Rcn, tgrid_new, fgrid,'spline');
        
        
        figure;
        imagesc(t_gated,f_uniform,synout_new,...
            [1.2*SR(fibertype), max(max(synout_new))]);
        xlabel('Time (s)','FontSize',20);
        ylabel('CF (Hz)','FontSize',20);
        title('AN Output','FontSize',20);
        
        figure;
        imagesc(t_gated,f_uniform, Ric_new);
        xlabel('Time (s)','FontSize',20);
        ylabel('CF (Hz)','FontSize',20);
        title('IC MF Cells','FontSize',20);
        
%         figure;
%         Ric_BR = ((Rcn_new - Ric_new) + abs(Rcn_new - Ric_new))/2;
%         imagesc(t_gated,f_uniform, Ric_BR);
%         xlabel('Time (s)','FontSize',20);
%         ylabel('CF (Hz)','FontSize',20);
%         title('IC Band Reject Cells','FontSize',20);
        
%         figure;
%         plot(f_uniform,mean(synout_new,2),'k--','linew',2);
%         hold on;
%         plot(f_uniform,mean(Ric_new,2),'b','linew',2);
%         hold on;
%         plot(f_uniform,mean(Ric_BR,2),'r','linew',2);
%         ylabel('Average rate','FontSize',20);
%         xlabel('CF (Hz)','FontSize',20);
%         
    end
   
    IC_all(fibertype,:,:) = Ric_new;
    AN_all(fibertype,:,:) = synout_new;  
    
end

