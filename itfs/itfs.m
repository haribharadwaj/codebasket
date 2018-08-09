function masked = itfs(noisy, oracle, noise, LC, fs, nfilts)
% Ideal time-frequency segregation using binary masks in TF domain using a
% gammatone filterbank.
%
% USAGE:
%   masked = itfs(noisy, oracle, LC, fs);
%
% INPUTS:
%   noisy - The target-masker mixture
%   oracle - The clean target (same length as noisy)
%   noise - The actual noise added (oracle + noise should be exactly equal
%       to noisy.
%   LC - The local SNR criterion
%   fs - Sampling rate of the signals involved
%   nfilts - Number of bands to contruct TF representation with (optional)
%
% OUTPUTS:
%   masked - The cleaned mixture after ITFS is done
%

if ~exist('nfilts', 'var')
    nfilts  = 128;
end
f_low = 80;
f_high = 8000;

% Equally space CFs on an ERB scale
cfs = invcams(linspace(cams(f_low), cams(f_high), nfilts));


% Make all inputs as column vectors
noisy = noisy(:);
oracle = oracle(:);
noise = noise(:);

% Extract BM responses for mixture and oracle
% Phase align the filters so that resynthesizing is trivial

fprintf(1, 'Extracting basilar membrane filter outputs!\n');
bm_mix = gammatoneFast(noisy, cfs, fs, true);
bm_T = gammatoneFast(oracle, cfs, fs, true);
bm_N = gammatoneFast(noise, cfs, fs, true);


% Window and calculate SNR
win_ms = 20;
win_overlap_ms = 10;

win = blackman(ceil(fs * win_ms * 1e-3));
step = ceil(fs * win_overlap_ms * 1e-3);

masked = zeros(size(bm_mix));
for k = 1:nfilts
    fprintf(1, 'Processing filter # %d / %d \n', k, nfilts);
    [frames_mix, ~] = sig2frames(bm_mix(:, k), win, step);
    [frames_T, ~] = sig2frames(bm_T(:, k), win, step);
    [frames_N, tail] = sig2frames(bm_N(:, k), win, step);
    SNR = pow2db(mean(frames_T.^2, 1)) - ...
        pow2db(mean(frames_N.^2, 1));
    mask = SNR > LC;
    
    frames_masked = frames_mix.* repmat(mask(:)', size(frames_mix, 1), 1);
    masked(:, k) = frames2sig(frames_masked, tail, step);
end

masked = rmsnormalize(sum(masked, 2));





