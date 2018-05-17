function y = rmsnormalize(x)
% Returns sound vector rmsnormalized to 0.1
r = rms(x(:));
y = x*0.1/r;