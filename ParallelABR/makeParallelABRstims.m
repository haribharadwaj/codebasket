function [x, trains] = makeParallelABRstims(fs, dur, n_epochs,...
    burst_rate, fc_list, n_cycles_per_burst)
% Make stimuli for parallel ABR based on Polonenko & Maddox (2019).
%
% USAGE:
%   [x, train] = makeParallelABRstims(fs, dur, n_epochs,...
%    burst_rate, fc_list, n_cycles_per_burst);
%
% INPUTS (all optional):
%   fs: Sampling rate in Hz
%   dur: Duration of each epoch in seconds
%   n_epochs: Number of unique epochs to generate
%   burst_rate: Rate of the poisson process underlying each frequency
%   fc_list: List of tone-burst frequencies
%   n_cycles_per_burst: Number of cycles in each burst (should be odd)
% 
% OUTPUTS:
%   x: The generate stimulus of size (n_epochs x number_of_time_samples)
%   trains: The click trains underlying the burst sequence for each
%       frequency in fc_list. Size is (numel(fc_list) x n_epochs x
%       number_of_time_samples). Note that these trains have sign. To
%       extract responses, the absolute values should be used for circular
%       cross-correlation with EEG measurements.
%
% Reference:
% ---------
% Polonenko, M. J., & Maddox, R. K. (2019). The Parallel Auditory
% Brainstem Response. Trends in Hearing.
% https://doi.org/10.1177/2331216519871395
%
% ---------
% Copyright 2020. Hari Bharadwaj. All rights reserved.
% ---------

if ~exist('fs', 'var')
    fs = 48828.125;
end

if ~exist('dur', 'var')
    dur = 1.0;
end

if ~exist('n_epochs', 'var')
    n_epochs = 30;
end

if ~exist('burst_rate', 'var')
    burst_rate = 40; % For each fc in fc_list
end

if ~exist('fc_list', 'var')
    fc_list = [500, 1000, 2000, 4000, 8000];
end


if ~exist('n_cycles_per_burst', 'var')
    n_cycles_per_burst = 5;
end

%% Allocate space for click trains and stimuli
n_fc = numel(fc_list);
n_samps = floor(fs * dur);
trains = zeros(n_fc, n_epochs, n_samps);
x = zeros(n_fc, n_epochs, n_samps);

%% Draw burst center sample indices for each frequency and (unique) epoch
n_bursts = floor(burst_rate * dur);
click_inds = randi(n_samps, [n_fc, n_epochs, n_bursts]);

%% Populate trains with clicks
for kfc = 1:n_fc
    fc = fc_list(kfc);
    % Make sure odd number of samples in the burst
    N = roundodd(fs * n_cycles_per_burst / fc); 
    t = (-(N-1)/2 : (N-1)/2) / fs; % Center time at 0 exactly
    burst = cos(2 * pi * fc * t) .* blackman(N)';
    
    for kepoch = 1:n_epochs
        % Check if any burst would be cut-off at the ends and remove them
        click_inds_current = squeeze(click_inds(kfc, kepoch, :));
        bad_click_inds = (click_inds_current < (N-1)/2) | ...
            (click_inds_current > (n_samps - (N-1)/2));
        click_inds_current(bad_click_inds) = [];
        
        % Populate trains with 1 or -1 with equal probability
        % Check number of bursts in case some were dropped
        n_bursts_current = numel(click_inds_current);
        trains(kfc, kepoch, click_inds_current) = ...
            sign(rand(n_bursts_current, 1) - 0.5);
        
        % Convolve with burst
        train = squeeze(trains(kfc, kepoch, :));
        % Make group delay zero by using the 'same' argument for conv()
        x(kfc, kepoch, :) = conv(train, burst, 'same');
    end
end

x = scaleSound(squeeze(sum(x, 1)));


