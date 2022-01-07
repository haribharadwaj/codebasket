function [stim, f1, L1] = makeDPOAEtrial(f2,L2,fs,dur,ramp)
%%%%%%%%%%%%%%%%%%%
% USAGE:
%    [stim, f1, L1] = makeDPOAEtrial(f2, L2, fs,dur,ramp)
%    [stim, f1, L1] = makeDPOAEtrial(4000, 50, 48828.125,2,0.025);
%
%  f2: F2 frequency (Hz), F1 is chosen such that F2/F1 = 1.21
%  L2: Level of the primary frequency of interest (i.e) f2 in dB SPL
%  fs: Sampling rate (Hz)
%  dur: Duration of stimulus (s)
%  ramp: ramp duration (i.e each half) (s)
%
% References:
%   Johnson, T. A., Neely, S. T., Garner, C. A., and Gorga, M. P. (2006). In-
% fluence of primary-level and primary-frequency ratios on human distortion
% product otoacoustic emissions,â€? J. Acoust. Soc. Am. 119, 418â€“428.
%
%   Neely, S. T., Johnson, T. A., Kopun, J., Dierking, D. M., & Gorga, M. P. (2009).
% "Distortion-product otoacoustic emission input/output characteristics in 
% normal-hearing and hearing-impaired human ears," J. Acoust. Soc. Am. 126, 728.
%----------------
% Hari Bharadwaj, August 25, 2013
%----------------

t = 0:(1/fs):(dur-1/fs);

Hz2kHz = 1e-3;

f2byf1 = 1.22 + log2(9.6/(f2*Hz2kHz))*(L2/415).^2;
L1 = 80 + 0.137*log2(18/(f2*Hz2kHz))*(L2 - 80);

f1 = f2/f2byf1;


primary2 = scaleSound(rampsound(sin(2*pi*f2*t),fs,ramp));
primary1 = scaleSound(rampsound(sin(2*pi*f1*t),fs,ramp));

stim = [primary1; primary2]';

