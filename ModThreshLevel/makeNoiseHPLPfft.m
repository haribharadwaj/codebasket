
%-----------------------------------------------------
%% Settings
fs = 48828; % Sampling Rate
tmax = 2; % Duration in Seconds
rampSize = 0.025; %In seconds
filter = 1; % 0 = BB, 1 = LP, 2 = HP
CutOff = 2000; % Cut off determining LP, HP
f0 = 100; % Fundamental frequency for clicks
fmax = 10e3; % Max frequency in kiloherz, 10k is a good number
ITD = 200e-6;
playplot = 1;
%-----------------------------------------------------


t = 0:(1/fs):(tmax - 1/fs);


%% Making Noise

fstep = 1/tmax; %Frequency bin size
switch filter
    case 0
        hmin = 1;
        hmax = floor(fmax/fstep);
    case 1
        hmin = 1;
        hmax = floor(CutOff/fstep);
    case 2
        hmin = ceil(CutOff/fstep);
        hmax = floor(fmax/fstep);
    otherwise
        fprintf(2,'\n!WARNING: Filter Setting Unrecognizable!\n');
        fprintf(2,'Defaulting to Broadband!\n');
        hmin = 1;
        hmax = floor(fmax/fstep);
end

phase = rand(hmax-hmin+1,1)*2*pi;

noiseF = zeros(numel(t),1);
noiseF(hmin:hmax) = exp(1j*phase);
noiseF((end-hmax+1):(end-hmin+1)) = exp(-1*1j*phase);

noise = ifft(noiseF,'symmetric');
ramplen = ceil(rampSize*fs);
w = hanning(ramplen)';

noise = scaleSound(rampsound(noise,fs,rampSize));


if(playplot)
    plot(t,noise);
    [pxx,f] = pwelch(noise,[],[],[],fs);
    figure; plot(f,pxx);
    sound(noise,fs);
end

