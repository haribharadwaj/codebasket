function y = scaleSound(x)
% Scales a matrix appropriately to between +/- 0.95

clipbuffer = 0.95;

y = clipbuffer*x./max(abs(x(:)));