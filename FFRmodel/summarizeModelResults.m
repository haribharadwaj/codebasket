% Summarizing the model results
clear all;
close all;
clc;
flist = dir('./RESULTS/*.mat');

Nsims = numel(flist);
for j = 1:Nsims
    fprintf(1,'Reading File # %d / %d\n',j,Nsims);
    fname = strcat('./RESULTS/',flist(j).name);   
    [d(j,:),p(j,:),level(j),m(j)] = getResults(fname);
end

    