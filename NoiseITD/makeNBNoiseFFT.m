function noise = makeNBNoiseFFT(bw,fc,tmax,fs,rampSize,playplot)
% USAGE:
%    noise = makeNBNoiseFFT(bw,fc,tmax,fs,rampSize,playplot);
%  e.g.:
%    noise = makeNBNoiseFFT(50,1000,0.6,48828.125,0.025,1);
% Makes notched noise with different bandwidths. RMS is 0.1 always.
%  bw - Bandwidth of noise in Hz (two-side)
%  tmax - Duration of noise in seconds
%  fs - Sampling rate
%  fc - center frequency in Hz
%  rampSize - seconds
%  playplot - Whether to play the noise (using sound(.)) and plot it's
%                 waveform and spectrum
%-----------------------------------------------------
%% Settings

if(~exist('fs','var'))
    fs = 48828.125; % Sampling Rate
end

if(~exist('tmax','var'))
    tmax = 0.6; % Duration in Seconds
end

if(~exist('rampSize','var'))
    rampSize = 0.025; %In seconds
end

if(~exist('fc','var'))
    fc = 1000;
end


fmin = fc - bw/2;
fmax = fc + bw/2;

if(~exist('playplot','var'))
    playplot = 0;
end

%-----------------------------------------------------
t = 0:(1/fs):(tmax - 1/fs);


%% Making Noise

fstep = 1/tmax; %Frequency bin size

hmin = ceil(fmin/fstep);
hmax = floor(fmax/fstep);


phase = rand(hmax-hmin+1,1)*2*pi;

noiseF = zeros(numel(t),1);
noiseF(hmin:hmax) = exp(1j*phase);
noiseF((end-hmax+1):(end-hmin+1)) = exp(-1*1j*phase);


noise = ifft(noiseF,'symmetric');
noise = rampsound(noise,fs,rampSize);
noiserms = rms(noise);
noise = (noise/noiserms) * 0.1;


if(playplot)
    plot(t,noise);
    [pxx,f] = pwelch(noise,[],[],[],fs);
    figure; plot(f,pxx);
    soundsc(noise,fs);
end

