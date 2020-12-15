function [gaindB,recorded] = InEarCalib(calibfreq,calibdB,RP,Fs)
% Calibrates the sound level in ear at any 1 frequency
%
% USAGE:
%   gainDB = InEarCalib(calibfreq,calibdB,RP,Fs);
%   gainDB = InEarCalib(4000,72,RP,48828.125);
%
%


t = 0:(1/Fs):(0.5 - 1/Fs);

x = rampsound(sin(2*pi*calibfreq*t),Fs,0.05)';

oaelength = numel(x) + ceil(Fs*0.025);
sigrms = sqrt(mean(x.^2));

x = [x zeros(size(x))];

drop = 89 - calibdB + 3;

%Start dropping from maximum RMS (actual RMS not peak-equivalent)
wavedata = x'*db2mag(-1*drop)/sigrms;

recorded = zeros(10,oaelength);

stimTrigger = 1;
numreps = 10;

for j = 1:numreps
    fprintf(1,'\nRunning Calibration Trial #%d/%d\n',j,numreps);
    %Load the 2ch variable 'x' into the RP2:
    invoke(RP, 'WriteTagVEX', 'datain', 0, 'I16', (wavedata*2^15));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Send audio trigger
    invoke(RP, 'SetTagVal', 'trgname', stimTrigger); % Set Stimulus Trigger
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Set the delay of the sound
    invoke(RP, 'SetTagVal', 'onsetdel',100); % onset delay is in ms
    %Start playing from the buffer:
    invoke(RP, 'SoftTrg', 1); %Playback trigger
    
    currindex = invoke(RP, 'GetTagVal', 'indexin');
    
    while(currindex < oaelength)
        currindex=invoke(RP, 'GetTagVal', 'indexin');
    end
    
    invoke(RP, 'SoftTrg', 2); %Stop stim-play
    
    
    recorded(j,:) = invoke(RP, 'ReadTagVex', 'dataout', 0, oaelength,'F32','F64',1);
    if max(max(abs(recorded)))> 2^31/0.96542
        fprintf(2,'F*CKING CLIPPING!!!') 
    end
    pause(0.1);
    %WaitSecs(0.1);
    % Get ready for next trial
    invoke(RP, 'SoftTrg', 8); % Stop and clear OAE buffer
    
    
    %Reset the play index to zero:
    invoke(RP, 'SoftTrg', 5); %Reset Trigger
    %Clearing I/O memory buffers:
    invoke(RP,'ZeroTag','datain');
end

% Microphone characteristics
preAMPgain = 40; % dB

% NOTE:
% 0 dB SPL = 0 dB uV i.e 1 uV rms = 0 dB SPL
% Hence no additional shifts are needed

V_ref = 1e-6;

output = median(recorded,1);

output_V = output;

output_rms_V = rms(output_V);

output_SPL = db(output_rms_V/V_ref) - preAMPgain; 

gaindB = calibdB - output_SPL;





