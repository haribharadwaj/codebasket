function [frames, tail] = sig2frames(sig, win, step)

% Extracts windowed frames from continuous signals for overlap add analysis
% and synthesis
%
% USAGE:
%   frames = sig2frames(sig, win, step);
%
% INPUTS:
%   sig - Some signal
%   win - window function (nw x 1)
%   step - step size in samples
%
% OUTPUT:
%   frames - nw x numwins matrix of frames
%   tail - number of samples of unframed signal at the end

nw = numel(win); % Samples per window

% Determine number of windows and extra tail length of unframed signal
numwins = 0;
while ((numwins - 1) * step + nw <= numel(sig))
    numwins = numwins + 1;
end
numwins = numwins - 1;
tail = numel(sig) - ((numwins - 1) * step + nw);


frames = zeros(nw, numwins);

for k = 1:numwins
    frames(:, k) = win.*sig( ((k-1)*step + 1) : ((k-1)*step + nw) );
end

    