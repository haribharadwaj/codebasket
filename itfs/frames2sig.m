function sig = frames2sig(frames, tail, step)

% Resynthesizes signal from frames using overlap add
%
% USAGE:
%   sig = frames2sig(frames, tail, step);
%
% INPUTS:
%   frames - nw x numwins matrix of frames
%   tail - number of samples of unframed signal at the end
%   step - step size in samples
%
% OUTPUT:
%   sig - Resynthesized signal

sig = zeros((size(frames, 2) - 1)*step + tail + size(frames, 1), 1);

nw = size(frames, 1);
numwins = size(frames, 2);


for k = 1:numwins
    sig(((k - 1)*step + 1): ((k-1)*step + nw)) = ...
        sig(((k - 1)*step + 1): ((k-1)*step + nw)) + frames(:, k);
end

    