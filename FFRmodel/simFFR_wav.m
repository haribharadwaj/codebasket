clear;
clc;
addpath('./Zilany2014_new/');

wavname = 'aaa.wav';
wavname_noise = 'restaurant.wav';
%wavname = 'four_speakers.wav';
%wavname = 'tone_4kHz.wav';
%wavname = 'SAM_4kHz_40Hz.wav';
fs = 100e3;

plotting = 0;
fiblist = {'Low-SR','Med-SR','High-SR'};

noisedblist = [50, 60, 70, 80, 90];
Ntrials = 10;

for noisedb = noisedblist
    for trial = 1:Ntrials
        for fibertype = 3:3; % 1- LS, 2- MS, 3-HS
            fib = fiblist{fibertype};
            fprintf(1,'\n ------------Simulating %s fibers :)-----------\n',fib);
            SR = [0.1, 15, 100];
            stimdb = 75;
            
            
            Ntrials = 1;
            [pin_noise, fs_wav] = wavread(wavname_noise);
            tstart = ceil(rand*fs_wav*9); % The file is slightly more than 10 secs long..
            pin_noise = pin_noise(tstart:(tstart+ceil(fs_wav)-1),1);
            
            [pin, fs_wav] = wavread(wavname);
            pin = pin(1:ceil(fs_wav),1);
            
            pin = resample(pin,fs,fs_wav)';
            pin_noise = resample(pin_noise, fs, fs_wav)';
            dur = size(pin,2)/fs;
            isi = 50e-3; %ms
            
            stimrms = rms(pin);
            noiserms = rms(pin_noise);
            
            
            tuningType = 2;
            
            % 20e-6 RMS Pa is 0 dB SPL..
            % pin/stimrms has RMS of 1..
            % 20e-6*pi/stimrms has 0 dB SPL
            % Hence to get desired dB SPL:
            pin = 10^(stimdb/20)*20e-6*pin/stimrms;
            pin_noise = 10^(noisedb/20)*20e-6*pin_noise/noiserms;
            
            load vFreq;
            CF_step = 2;
            f = vFreq(15:CF_step:end);
            noisefac = [0, 1];
            for fac = noisefac
                for nCF = 1:numel(f)
                    CF = f(nCF);
                    nrep = 2;
                    vihc = model_IHC(pin + fac*pin_noise,CF,nrep,1/fs,dur+isi,1,1,tuningType);
                    [synout(nCF,:),synoutvar,psth] = ...
                        model_Synapse(vihc,CF,nrep,1/fs,fibertype,1,0);
                    
                    t = (0:(size(synout,2)-1))/fs;
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
                
                fprintf(1,'\n');
                
                
                % Gate out the onset response foor visualization
                t_gated = t(t>0.050);
                synout = synout(:, t>0.050);
                
                
                
                h_smooth = fspecial('gaussian',[5,5],5);
                synout_new(fac+1,:,:) = imfilter(synout,h_smooth);
                % Rescale to correct rates:
                synout_new(fac+1,:,:) = (synout_new(fac+1,:,:)/2); % Not ideal.. have to fix zilany code
                
                
                if(plotting)
                    fprintf(1, '\nAll done! Plotting... Hold on!\n');
                    
                    subplot(2,1,fac+1);
                    imagesc(t_gated,f,squeeze(synout_new(fac+1,:,:)),...
                        [1.2*SR(fibertype), 400]);
                    xlabel('Time (s)','FontSize',20);
                    ylabel('CF (Hz)','FontSize',20);
                    title('AN Output','FontSize',20);
                    xlim([min(t_gated), min(t_gated)+0.1]);
                    ylim([min(f), 10e3]);
                    
                end
                clear synout
            end
        end
        
        winsamps = round(50e-3*fs);
        nwins = floor(size(synout_new,3)/winsamps);
        tap = dpss(winsamps, 1, 1)';
        
        f_fft = (0:(winsamps-1))*fs/winsamps;
        ind = (f_fft < 2000);
        f_fft = f_fft(ind);
        nconds = 2;
        S = zeros(nconds, numel(f), sum(ind), nwins);
        
        f0 = 140;
        indharms = ((mod(f_fft, f0) < 3) | (mod(f_fft, f0) > (f0-3))) & (f_fft > 20);
        
        for fac = 1:nconds
            for k = 1:nwins
                sig = squeeze(synout_new(fac,:, ((k-1)*winsamps+1):(k*winsamps)));
                wsig = sig.*repmat(tap, numel(f), 1);
                temp = squeeze(abs(fft(wsig, [], 2)));
                temp = temp(:, ind);
                S(fac, :, :, k) = temp;
                SNR(fac, :, k) = 10*log10(sum(temp(:, indharms).^2, 2)./sum(temp(:, ~indharms).^2,2));
            end
        end
        
        
        
        if(plotting)
            cond1 = SNR(1,:,:);
            cond2 = SNR(2,:,:);
            [n1, x1] = hist(cond1(:), 50);
            [n2, x2] = hist(cond2(:), 50);
            figure;
            plot(x1, n1/numel(cond1), 'b-', 'linew', 2);
            hold on;
            plot(x2, n2/numel(cond2), 'r-', 'linew', 2);
            xlabel('Neural SNR in time-frequency atom (dB)', 'FontSize', 20);
            ylabel('Probability of occurrence', 'FontSize', 20);
        end
        fname = strcat('./RESULTS/simNoisy_trial', num2str(trial),'_noise', num2str(noisedb),'dB.mat');
        save(fname, 'S', 'f', 't_gated', 'SNR', 'f_fft', 'f0', 'stimdb', 'noisedb', 'trial');
        clear SNR
    end
end
