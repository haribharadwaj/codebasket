% Script to test calibration
clear all;
close all hidden;
clc;

fig_num=99;
USB_ch=1;
IAC = -1;
FS_tag = 3;
Fs = 48828.125;
[f1RP,RP,FS]=load_play_circuit(FS_tag,fig_num,USB_ch,0,IAC);

f2 = 4000; % Hz
t = 0:1/(Fs):(10-1/Fs);
y = sin(2*pi*f2*t);
y = rampsound(y,Fs,0.02);
setpeak = 0.02; 
% This should correspond to the pk-pk of voltage to headphone driver
% View on oscilloscope to confirm

wavedata = setpeak*[y' y']';

Nreps = 10;
%Load the 2ch variable 'x' into the RP2:
invoke(RP, 'WriteTagVEX', 'datain', 0, 'I16', (wavedata*2^15));
WaitSecs(10.0);

for rep = 1:Nreps
    fprintf(1,'\n Doing rep #%d/%d\n',rep,Nreps);
   
    %Set the delay of the sound
    invoke(RP, 'SetTagVal', 'onsetdel',100); % onset delay is in ms
    %Start playing from the buffer:
    invoke(RP, 'SoftTrg', 1); %Playback trigger
    
    WaitSecs(11 + 0.1 + 0.2*rand);
    
    %Reset the play index to zero:
    invoke(RP, 'SoftTrg', 5);
    WaitSecs(0.1);
    
end

%Clearing I/O memory buffers:
invoke(RP,'ZeroTag','datain');
