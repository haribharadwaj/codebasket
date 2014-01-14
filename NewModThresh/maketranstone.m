function [x,env] = maketranstone(fc,fm,m,bw,fs,dur,rise,phi)
% Creating a transposed tone
%
% USAGE:
% x = maketranstone(fc,fm,m,bw,fs,dur,rise,phi);
%
% fc: Carrier/center frequency  (Hz)
% fm: Modulation frequency (Hz)
%  m: Modulation depth (0 to 1)
% bw: twice the low-pass filter cutoff applied on the envelope (Hz)
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

env = sin(2*pi*fm*t - phi);
env = env.*(env > 0);

% Making a filter whose impulse response if purely positive (to avoid phase
% jumps) so that the filtered envelope is purely positive. Using a dpss
% window to minimize sidebands. For a bandwidth of bw, to get the shortest
% filterlength, we need to restrict time-bandwidth product to a minimum.
% Thus we need a length*bw = 2 => length = 2/bw (second). Hence filter
% coefficients are calculated as follows:
b = dpss(floor(2*fs/bw),1,1);  % Using to increase actual bw when rounding
b = b-b(1);
env = filter(b,1,env);
env = env(1:numel(t));
env = env/max(env);

% Imposing depth

if((m < 0) || (m > 1))
    fprintf(2,'WARNING! Making m = 1\n');
    m = 1;
end
env = m*env + (1 - m);

carr = sin(2*pi*fc*t);

x = rampsound(carr.*env, fs, rise);







