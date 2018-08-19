function [x, itd] = noiseITD(fmin, fmax, itd, fs, nfilts)
% Applies an arbitrary ITD time course to a noise band by vocoding
%
% USAGE:
%   [x, itd] = noiseITD(fmin, fmax, itd, fs, nfilts);
%
% INPUTS:
%   fmin - Low end frequency of noise band (Hz,  fmin > 80)
%   fmax - High end frequency of noise band (Hz, fmax < 8k)
%   itd - The ITD timecourse (milliseconds)
%   fs - Sampling rate of the signals involved
%   nfilts - Number of bands to contruct TF representation with
%
% OUTPUTS:
%   x - two channel signal with specified ITD (same length as 'itd')
%   itd - the itd function input as is or the 5 second-long one generated
%         by default
%
% NOTE:
%  All arguments are optional. The defaults are fmin = 300 Hz, fmax = 3000 Hz,
%  itd is a low-pass noise function with a 40 Hz cutoff and a 0.8 ms max in
%  either direction, nfilts = 128, fs = 48828.125 Hz.
%
%-----------------------
% Copyright Hari Bharadwaj 2018. All rights reserved.
% hbharadwaj@purdue.edu
%-----------------------
if ~exist('fmin', 'var')
    fmin  = 300;
end

if ~exist('fmax', 'var')
    fmax  = 3000;
end

if ~exist('fs', 'var')
    fs = 48828.125;
end

if ~exist('itd', 'var')
    dur = 5;
    itd = makeNBNoiseFFT(40, 20.1, dur, fs, 0.05, 0);
    itd = itd * 0.66/max(abs(itd));
end

if ~exist('nfilts', 'var')
    nfilts  = 128;
end

ms = 1e-3; % ITD is in milliseconds
if any(abs(fmax * itd/2 * ms) > 1)
    warning('Some phase shifts needed exceed a whole cycle!');
end
f_low = 80;
f_high = 8000;

% Equally space CFs on an ERB scale
cfs = invcams(linspace(cams(f_low), cams(f_high), nfilts));


% Make all inputs as column vectors
itd = itd(:);

% Make noise band
dt = 1/fs;
dur = numel(itd)*dt;
bw = fmax - fmin;
fc = (fmin + fmax)/2;
ramp = 0.025;
noise = makeNBNoiseFFT(bw,fc, dur, fs, ramp, 0);

% Extract BM responses for noise carrier
% Phase align the filters so that resynthesizing is trivial

fprintf(1, 'Extracting basilar membrane filter outputs!\n');
bm = gammatoneFast(noise, cfs, fs, true);
bm2 = zeros(size(bm));

% Calculate instantaneous phase and apply CF-dependent phase-shift to get
% the specified ITD

for k = 1:nfilts
    fprintf(1, 'Processing filter # %d / %d \n', k, nfilts);
    analytic = hilbert(bm(:, k));
    amp = abs(analytic);
    phi = unwrap(angle(analytic));
    
    shift = (itd/2) * cfs(k) * ms;  % In cycles
    BW = invcams(cams(cfs(k)) + 0.5) - invcams(cams(cfs(k)) - 0.5);
    if max(abs(diff(shift) * fs)) > BW/2
        msg = ['Phase modulations are too fast for CF = ', num2str(cfs(k))];
        warning(msg);
    end
        
    contra = amp .* cos(phi + shift*2*pi);
    bm2(:, k) = contra;
    ipsi = amp .* cos(phi - shift*2*pi);
    bm2(:, k) = ipsi;
end

x(1, :) = rmsnormalize(sum(bm, 2));
x(2, :) = rmsnormalize(sum(bm2, 2));






