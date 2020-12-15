% Super basic DPOAE screening test with 10 dB steps
clear; close all hidden; clc;

fig_num=99;
USB_ch=1;

sID = input('Enter Subject ID (good practice to enter ID_ear (like I13_left):','s');

% Change this to the wherever you put the DPOAE folder
paraDir = 'C:\Experiments\DPOAE_ER10X_SPL\';
if(~exist(strcat(paraDir,'\subjResponses\',sID),'dir'))
    mkdir(strcat(paraDir,'\subjResponses\',sID));
end


IAC = -1;
FS_tag = 3;

Fs = 48828.125;



[f1RZ,RZ,FS]=load_play_circuit(FS_tag,fig_num,USB_ch,0,IAC);

% Noise off
invoke(RZ, 'SoftTrg', 4);

f2 = 4000; %input('Enter frequency of interest(Hz):');
tstart = tic;
L2_list = 0:6:76;
dur = 0.6;
ramp = 0.025;

% Only 1k, 2k, 4k and 8kHz are allowed for now..
% To allow other frequencies, add them to the switch case.. The shift should
% roughly be the difference in SPL between 0 dB HL at that frequency and 0 dB
% HL at 4kHz.. The exact shift is not important.. This is just so that the
% right range on input levels are presented i.e., L2 goes from below
% threshold to well above threshold for most subjects.

switch f2
    case 1000
        shift = 30;
    case 2000
        shift = 25;
    case 4000
        shift = 0;
    case 8000
        shift = 0;
    otherwise
        error('Input frequency not in approved list!');
end

% In-Ear Calbration
calibfreq = 1000;
calibdB = 72;
[gaindB,recorded] = InEarCalib(calibfreq,calibdB,RZ,Fs);
fprintf(1, '########## Gain = %1.2f dB #########\n', gaindB);
if(gaindB > 8)
    error('\nCHECK EARPLUG INSERTION, LEVEL IS OFF BY MORE THAN 8 dB!\n');
end


L2_list = L2_list + shift;

stimTrigger = 1;

figure;
for k = 1:numel(L2_list)
    L2 = L2_list(k);
    [y, f1, L1] = makeDPOAEtrial(f2, L2, Fs, dur, 0.005);
    extrasamples = ceil(Fs*0.025);
    oaelength = size(y, 1) + extrasamples;
    
    
    numtrials = 8;
    
    
    
    jit = rand(numtrials,1)*0.01;
    
    nfloor = 0;
    OAE = [];
    loopcount = 0;
    nfloortarget = -25;
    while (nfloor > nfloortarget) || (loopcount < 2)
        loopcount = loopcount + 1;
        fprintf('Loop count = %d\n', loopcount);
        for j = 1:numtrials
            fprintf(1,'Block %d/%d - Trial Number %d/%d \n',...
                k,numel(L2_list),j,numtrials);
            
            
            nrchannels = size(y,2);
            
            if(nrchannels~=2)
                fprintf(2,'\n HUGE WARNING!! This is a mono wav file. Frequencies may be doubled!!\n');
            end
            
            
            
            rms2 = rms(y(:, 2));
            rms1 = rms(y(:, 1));
            
            digitalAtt = 6; % dB. TO avoid clipping
            %-----------------
            % Why 89?
            %-----------------
            % ER-10X drivers give about 78dB SPL for a 1kHz tone with a 1V-rms drive.
            % Max output is +/-5V peak i.e 3.54V-rms which is 11 dB higher.
            % Thus 89 dB-SPL is the sound level for tones when they occupy full
            % range.
            
            % Full range in MATLAB for a pure tone is +/- 1 which is 0.707 rms and
            % that corresponds to 83 dB SPL at the end. So if we want a signal
            % with rms sigrms to be x dB, then (89 - x) should be
            % db(sigrms/0.707).
            
            drop1 = 89 - L1 + 3 - digitalAtt - gaindB; % The extra for the 0.707 factor
            drop2 = 89 - L2 + 3 - digitalAtt - gaindB;
            x(:, 1) = y(:, 1)*0.95/rms1 * db2mag(-digitalAtt);
            x(:, 2) = y(:, 2)*0.95/rms2 * db2mag(-digitalAtt);
            invoke(RZ, 'SetTagVal', 'attA', drop1);
            invoke(RZ, 'SetTagVal', 'attB', drop2);
            
            %Start dropping from maximum RMS (actual RMS not peak-equivalent)
            
            % In above, we have assumed that both the speakers in the ER-10c give
            % 72 dB SPL for a 1V, 1KHz tone.. At BU, the ER-10Cs do indeed give
            % those.. If yours give you slightly different levels, set these
            % numbers accordingly.
            
            
            
            wavedata = x';
            
            if(any(abs(wavedata(1, :)) > 1))
                error('What did you do!? Sound is clipping!! Cannot Continue!!\n');
            end
            
            %Load the 2ch variable 'x' into the RP2:
            invoke(RZ, 'WriteTagVEX', 'datain', 0, 'I16', (wavedata*2^15));
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Send audio trigger
            invoke(RZ, 'SetTagVal', 'trgname', stimTrigger); % Set Stimulus Trigger
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %Set the delay of the sound
            invoke(RZ, 'SetTagVal', 'onsetdel',100); % onset delay is in ms
            %Start playing from the buffer:
            invoke(RZ, 'SoftTrg', 1); %Playback trigger
            
            currindex = invoke(RZ, 'GetTagVal', 'indexin');
            
            while(currindex < oaelength)
                currindex=invoke(RZ, 'GetTagVal', 'indexin');
            end
            
            invoke(RZ, 'SoftTrg', 2); %Stop stim-play
            
            
            OAE_temp(j,:)=invoke(RZ, 'ReadTagVex', 'dataout', 0, oaelength,'F32','F64',1); %#ok<SAGROW>
            
            if max(max(abs(OAE_temp)))> 2^31/0.96542
                fprintf(2,'\n F*CKING CLIPPING!!! \n')
                fprintf(2, ['\n  Try using a hardware high-pass between'...
                    'mic-out of ER-10X and analog-in of TDT\n']);
            end
            
            pause(jit(j));
            %WaitSecs(jit(j));
            % Get ready for next trial
            invoke(RZ, 'SoftTrg', 8); % Stop and clear OAE buffer
            
            
            %Reset the play index to zero:
            invoke(RZ, 'SoftTrg', 5); %Reset Trigger
            %Clearing I/O memory buffers:
            invoke(RZ,'ZeroTag','datain');
        end
        
        OAE = [OAE; OAE_temp]; %#ok<AGROW>
        % Reject bad trials
        good = [];
        artifactFactor = 1.25;
        for trial = 1:(numtrials*loopcount)
            if rms(OAE(trial, :)) < median(rms(OAE, 2)) * artifactFactor
                good = [good, trial]; %#ok<AGROW>
            end
        end
        if mod(numel(good),2)~=0
            good = good(1:(end-1));
        end
        if (numel(good)>16) && (mod(numel(good),2) == 1)
            good = good(1:(end-1));
        end
        OAEgood = OAE(good, :);
        % Calculate noise level
        OAE_noise = OAEgood;
        OAE_noise(2:2:end, :) = -1 * OAE_noise(2:2:end, :);
        noise = mean(OAE_noise, 1);
        signal = mean(OAEgood, 1);
        w = dpss(numel(noise),1,1)';
        w = w / sum(w);
        f_dp = 2*f1 - f2;
        t_calc = (0:(numel(noise)-1))/Fs;
        wsn = w.*sin(2*pi*f_dp*t_calc);
        wcn = w.*cos(2*pi*f_dp*t_calc);
        
        mic_sens = 50e-3; % V/Pa
        mic_gain = db2mag(40);
        V_to_microV = 1e6;
        ref = 20; % microP
        factor = V_to_microV / (mic_sens * mic_gain * ref);
        % Super crude noise-floor estimate by flipping half the trials
        noisefloor(k) = db(sqrt(sum(wsn.*noise)^2 + sum(wcn.*noise)^2) * factor); %#ok<SAGROW>
        DP(k) = db(sqrt(sum(wsn.*signal)^2 + sum(wcn.*signal)^2) * factor); %#ok<SAGROW>
        nfloor = noisefloor(k);
    end
    hold on;
    plot(L2, noisefloor(k), 'ok', 'linew', 2);
    hold on;
    plot(L2, DP(k), 'or', 'linew', 2);
    xlabel('L2 dB SPL', 'FontSize', 16);
    % Can use probe calibration offline to specify in SPL
    % This is just to see if there is DPOAE above noise floor and if it is
    % growing roughly like it should
    ylabel('DPOAE (dB SPL)', 'FontSize', 16);
    drawnow;
    
    fdir = strcat(paraDir,'/subjResponses/',sID);
    datetag = datestr(clock);
    datetag(strfind(datetag,' ')) = '_';
    datetag(strfind(datetag,'-')) = '_';
    datetag(strfind(datetag,':')) = '_';
    fname = strcat(fdir,'/',sID,'_DPOAE_RAW_',num2str(L2),'_',num2str(round(f2)),'_',datetag,'.mat');
    save(fname, 'L2', 'f1', 'f2', 'sID', 'OAE');
end
set(gca, 'FontSize', 16);
title(['Subject ID: ', strrep(sID, '_', ' - ')], 'FontSize', 16);
datetag = datestr(clock);
datetag(strfind(datetag,' ')) = '_';
datetag(strfind(datetag,'-')) = '_';
datetag(strfind(datetag,':')) = '_';
fdir = strcat(paraDir,'/subjResponses/',sID);
fname = strcat(fdir,'/',sID,'_DPOAE_', num2str(round(f2)),'_',datetag,'.mat');

save(fname,'L2_list', 'f1','f2', 'sID', 'DP', 'noisefloor');
close_play_circuit(f1RZ, RZ);
toc(tstart)



