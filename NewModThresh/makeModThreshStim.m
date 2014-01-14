function y = makeModThreshStim(fc,fm,m,type,SNR,bw,fs,dur,rise,playplot)
% USAGE:
%-------
%  x = makeModThreshStim(fc,fm,m,SNR,bw,fs,dur,rise,playplot);
%
% Parameters
%------------
%   fc - Center Frequency (Hz)
%   fm - Modulation Frequency (Hz)
%   m - Modulation Depth (0 to 1)
%   type - 1 for "transposed", 0 for SAM
%   SNR - SNR in dB (spectrum level)
%   bw - Bandwidth of the noise (same is used for notch of masker) (Hz)
%   fs - Sampling rate
%   dur - Length of the thing (s)
%   rise - Ramp risetime (same is used a falltime at the end) (s)
%   playplot - Whether to play the sounds and plot the signals/spectra
%
% Returns
%-----------
%   y - Sound vector

if(type == 0)
    x = makesamnoise(fc,fm,m,fs,bw,dur,rise,0);
else
    bwnoise = bw;
    bwenv = bw/2;
    x = maketransnoise(fc,fm,m,fs,bwnoise,bwenv,dur,rise,0);
end
n = makeNotchNoiseFFT(bw/2,0,fc,dur,fs,rise,0);

xrms = rms(x)*sqrt(1e3/bw);

% The noise spans 200 Hz to 20 kHz => 19.8 kHz
% It has a notch of width "bw"

nrms = rms(n)*sqrt(1e3/(19.8e3 - bw));

y = scaleSound(db2mag(SNR)*x/xrms + n/nrms);

if(playplot)
    t = (0:(numel(y)-1))/fs;
    subplot(2,1,1);
    plot(t,y,'linew',2);
    xlabel('Time (s)','FontSize',20);
    ylabel('Scaled Signal','FontSize',20);
    
    [pxx, f] = pmtm(y,3,[],fs);
    subplot(2,1,2);
    
    ref = (20e-6)^2;
    plot(f,10*log10(pxx/ref),'linew',2); % pxx is in intensity units
    xlim([0, 16e3]);
    xlabel('Frequency (Hz)','FontSize',20);
    ylabel('Power (dB)','FontSize',20);
    soundsc(y,fs);
end


