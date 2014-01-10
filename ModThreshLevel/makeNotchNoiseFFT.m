function noise = makeNotchNoiseFFT(bw,asymm,fc,tmax,fs,rampSize,playplot)
% USAGE:
%    noise = makeNotchNoiseFFT(bw,fc,tmax,fs,rampSize,playplot);
%  e.g.:
%    noise = makeNotchNoiseFFT(50,4000,0.6,48828,0.025,1);
% Makes notched noise with different bandwidths.
%  bw - Bandwidth of notch in Hz (one-side)
%  tmax - Duration of noise in seconds
%  fs - Sampling rate
%  fc - Notch frequency
%  rampSize - seconds
%  playplot - Whether to play the noise (using sound(.)) and plot it's
%                 waveform and spectrum
%-----------------------------------------------------
%% Settings

if(~exist('fs','var'))
    fs = 48828; % Sampling Rate
end

if(~exist('tmax','var'))
    tmax = 0.6; % Duration in Seconds
end

if(~exist('rampSize','var'))
    rampSize = 0.025; %In seconds
end

if(~exist('fc','var'))
    fc = 4000;
end

switch asymm
    case 0
        CutOffLP = fc - bw; % Cut off determining LP, HP
        CutOffHP = fc + bw;
    case 1
        if(bw < 0.2)
            fprintf(2,'WARNING! Filter too-narrow for asymmetric condition!');
        end
        CutOffLP = fc - bw - 0.2*fc; % Cut off determining LP, HP
        CutOffHP = fc + bw;
    case 2
        if(bw < 0.2)
            fprintf(2,'WARNING! Filter too-narrow for asymmetric condition!');
        end
        CutOffLP = fc - bw; % Cut off determining LP, HP
        CutOffHP = fc + bw + 0.2*fc;
end
blobs = 0;
if(blobs)
    blob = 0.25;
    fmax = CutOffHP + blob*fc; % Max frequency in kiloherz
    fmin = CutOffLP - blob*fc;
else
    fmax = 15e3;
    fmin = 200;
end

if(~exist('playplot','var'))
    playplot = 0;
end
%-----------------------------------------------------


t = 0:(1/fs):(tmax - 1/fs);


%% Making Noise

fstep = 1/tmax; %Frequency bin size

%--------------------
% MAKE LP PART
%--------------------

CutOff = CutOffLP;

hmin = ceil(fmin/fstep);
hmax = floor(CutOff/fstep);


phase = rand(hmax-hmin+1,1)*2*pi;

noiseF = zeros(numel(t),1);
noiseF(hmin:hmax) = exp(1j*phase);
noiseF((end-hmax+1):(end-hmin+1)) = exp(-1*1j*phase);


noiseF_LP = noiseF;

%--------------------
% MAKE HP PART
%--------------------
CutOff = CutOffHP;


hmin = ceil(CutOff/fstep);
hmax = floor(fmax/fstep);

phase = rand(hmax-hmin+1,1)*2*pi;

noiseF = zeros(numel(t),1);
noiseF(hmin:hmax) = exp(1j*phase);
noiseF((end-hmax+1):(end-hmin+1)) = exp(-1*1j*phase);

noiseF_HP = noiseF;

noiseF = noiseF_LP + noiseF_HP;
%------------------------------------------------------------------------

noise = ifft(noiseF,'symmetric');
noise = scaleSound(rampsound(noise,fs,rampSize));


if(playplot)
    plot(t,noise);
    [pxx,f] = pwelch(noise,[],[],[],fs);
    figure; plot(f,pxx);
    sound(noise,fs);
end

