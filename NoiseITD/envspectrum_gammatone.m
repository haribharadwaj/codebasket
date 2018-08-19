function [f_env, cfs, Penvs] = envspectrum_gammatone(x, fs, nfilts)

if ~exist('fs', 'var')
    fs = 48828.125;
end

if ~exist('nfilts', 'var')
    nfilts  = 16;
end

f_low = 80;
f_high = 8000;

% Equally space CFs on an ERB scale
cfs = invcams(linspace(cams(f_low), cams(f_high), nfilts));
bm = gammatoneFast(x, cfs, fs, true);

Penvs = [];
for k = 1:nfilts
    fprintf(1, 'Processing filter # %d / %d \n', k, nfilts);
    env = abs(hilbert(bm(:, k)));
    [Penv, f_env] = pmtm(env, 3.5, [], fs);
    Penvs = [Penvs, pow2db(Penv(:))]; %#ok<AGROW>
end
