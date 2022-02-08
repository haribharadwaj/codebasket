[x, trains] = makeParallelABRstims;
fs = 48828; % Should be same as what's used to create stims
filename_prefix = 'parallelABRtrial';
for k = 1:size(x, 1)
    fname = strcat(filename_prefix, num2str(k), '.wav');
    audiowrite(fname, x(k, :), fs);
end
