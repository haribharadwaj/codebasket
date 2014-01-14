function [x,env] = makesamtone(fc,fm,m,fs,dur,rise,phi)
% Creating a transposed tone
%
% USAGE:
% x = makesamtone(fc,fm,m,fs,dur,rise,phi);
%
% fc: Carrier/center frequency  (Hz)
% fm: Modulation frequency (Hz)
%  m: Modulation depth (0 to 1)
% fs: sampling rate (Hz)
% dur: Duration (rise time will be included within this) in seconds
% rise: Rise time (fall time would be the same) in seconds (dpss ramp)
% phi: Starting phase
%
% -----------------
% Hari Bharadwaj
% hari@nmr.mgh.harvard.edu
%------------------

t = 0:(1/fs):(dur - 1/fs);

env = 1+sin(2*pi*fm*t - phi);
env = env/max(env);

% Imposing depth

if((m < 0) || (m > 1))
    fprintf(2,'WARNING! Making m = 1\n');
    m = 1;
end
env = m*env + (1 - m);

carr = sin(2*pi*fc*t);

x = rampsound(carr.*env, fs, rise);







