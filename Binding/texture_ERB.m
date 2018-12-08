function [stim, fs] = texture_ERB(nfreqs, ngroup, rho, seq, fs,...
    dur, playstim)
% USAGE:
%   [stim, fs] = texture_ERB(nfreqs, ngroup, rho, seq, fs,...
%                                 dur, playstim);
%
% Example:
% [stim, fs] = texture_ERB(15, 10, 0.9, [1, 2, 1, 2], 44100,...
%                               1.0, 1.0)
%
%-----------
% nfreqs: Number of tones in mixture (optional, default 25)
% ngroup: Number of tones to be temporally coherent
%         (optional, default ceil(nfreqs * 0.8))
% seq: Sequence of incoherent and ccoherent mixtures.
%      1 is for incoherent, 2 is for coherent. Example: [1, 1, 2] gives
%      two token of incoherent stims follwed by one token of coherent.
%      Note all coherent tokens have the same number of grouped tones.
%      (optional, default [1, 2, 1, 2])
% rho: Correlation between the envelopes of grouped tones
% fs: Sampling rate in Hz (optional, default 44100)
% dur: Duration in seconds of each token in seq (optional, default 1.0)
% playstim: Whether to play the stim
%           (optional: default 1, plays it)
%------------
% Copyright © 2015 Hari Bharadwaj. All rights reserved.
% hari@nmr.mgh.harvard.edu
% Aug 15, 2015
%------------

if ~exist('nfreqs', 'var') || isempty(nfreqs)
    nfreqs = 20;
end

if ~exist('ngroup', 'var') || isempty(ngroup)
    ngroup = ceil(nfreqs * 0.8);
end

if ~exist('seq', 'var') || isempty(seq)
    seq = [1, 2, 1, 2];
end

if ~exist('rho', 'var') || isempty(rho)
    rho = 1.0;
end

if ~exist('fs', 'var') || isempty(fs)
    fs = 48828.125;
end

if ~exist('dur', 'var') || isempty(dur)
    dur = 1.0;
end

if ~exist('playstim', 'var') || isempty(playstim)
    playstim = 1;
end

rise = 0.002;
t = 0: (1/fs): (dur - 1/fs);



f1 = 200;
fmax = 8000;
nERBs = cams(fmax) - cams(f1);
spacing_ERBs = nERBs / (nfreqs - 1);
fprintf(1, 'This stim will have successive tones separated by %2.2f ERBs\n', spacing_ERBs);
if spacing_ERBs < 1.0
    warning('The spacing between tones is LESS THAN 1 ERB!\n');
end

bwfilt = 80;
% Making a filter whose impulse response if purely positive (to avoid phase
% jumps) so that the filtered envelope is purely positive. Using a dpss
% window to minimize sidebands. For a bandwidth of bw, to get the shortest
% filterlength, we need to restrict time-bandwidth product to a minimum.
% Thus we need a length*bw = 2 => length = 2/bw (second). Hence filter
% coefficients are calculated as follows:
b = dpss(floor(2*fs/bwfilt),1,1);  % Using to increase actual bw when rounding
b = b-b(1);
b = b / sum(b);


envrate = 14;
bw = 20;
% The above two together gives 4 to 24 Hz envelop

z = 0;
for k = 1:nfreqs
    env = makeNBNoiseFFT(bw, envrate, dur, fs, rise, 0)';
    env = env.*(env > 0);
    env = filter(b,1,env);
    env = env(1:numel(t));
    x = env .* sin(2*pi*invcams(cams(f1) + spacing_ERBs*(k-1))*t);
    z = z + x;
end
z = rampsound(z, fs, rise);
z = z / rms(z);

y = 0;
group = randperm(nfreqs);
fg = group(1:ngroup);
bg = group((ngroup+1):end);

env1 = makeNBNoiseFFT(bw, envrate, dur, fs, rise, 0)';
env1 = env1.*(env1 > 0);
env1 = filter(b,1,env1);
env1 = env1(1:numel(t));
for k = fg
    env2 = makeNBNoiseFFT(bw, envrate, dur, fs, rise, 0)';
    env2 = env2.*(env2 > 0);
    env2 = filter(b,1,env2);
    env2 = env2(1:numel(t));
    env = sqrt(rho)*env1 + sqrt(1 - rho^2)*env2;
    x = scaleSound(rampsound(env .* sin(2*pi*invcams(cams(f1) + spacing_ERBs*(k-1))*t), fs, rise));
    y = y + x;
end

for k = bg
    env = makeNBNoiseFFT(bw, envrate, dur, fs, rise, 0)';
    env = env.*(env > 0);
    env = filter(b,1,env);
    env = env(1:numel(t));
    x = scaleSound(rampsound(env .* sin(2*pi*invcams(cams(f1) + spacing_ERBs*(k-1))*t), fs, rise));
    y = y + x;
end
y = rampsound(y, fs, rise);

y = y/rms(y);
stim = [];
for k = seq
    if (k == 1)
        stim = [stim, z]; %#ok<*AGROW>
    else
        stim = [stim, y];
    end
end

% Just to start out with a signal amplitude that's less than 1.0 in MATLAB
stim = 0.04 * stim / rms(stim); 

if(playstim)
    sound(stim, fs);
end
