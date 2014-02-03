function noise = makesamnoise(fc,fm,m,fs,bw,dur,rise,playplot)
% USAGE:
%----------------------------------------------------
%  x = makesamnoise(fc,fm,m,fs,bw,dur,rise,playplot);
%----------------------------------------------------
%    fc - Center frequency of the noise band
%    fm - Modulation frequency
%    m - Modulation Depth
%    fs - Sampling frequency
%    bw - Bandwidth of the noise band
%    dur - Duration of the sound
%    rise - Ramp duration (each side)
%    playplot - Whether to play and plot the sound
%----------------------------------------------------

if(~exist('playplot','var'))
    playplot = 0;
end

fmin = fc - bw/2;
fmax = fc + bw/2;


t = (0:(1/fs):(dur - 1/fs))';


%% Making Noise

fstep = 1/dur; %Frequency bin size

hmin = ceil(fmin/fstep);
hmax = floor(fmax/fstep);


phase = rand(hmax-hmin+1,1)*2*pi;

noiseF = zeros(numel(t),1);
noiseF(hmin:hmax) = exp(1j*phase);
noiseF((end-hmax+1):(end-hmin+1)) = exp(-1*1j*phase);

noise = ifft(noiseF,'symmetric');

% Envelope
env = 1+sin(2*pi*fm*t);
env = env/max(env);

% Imposing depth
if((m < 0) || (m > 1))
    fprintf(2,'WARNING! Making m = 1\n');
    m = 1;
end

env = m*env + (1 - m);

noise = env.*noise;

noise = scaleSound(rampsound(noise,fs,rise));


if(playplot)
    plot(t,noise);
    [pxx,f] = pwelch(noise,[],[],[],fs);
    figure; plot(f,pxx);
    sound(noise,fs);
end

