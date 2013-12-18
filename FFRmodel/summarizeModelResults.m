% Summarizing the model results

flist = dir('./RESULTS/*.mat');

Nsims = numel(flist);
d = zeros(1,Nsims);
p = zeros(1,Nsims);
level = zeros(1,Nsims);
m = zeros(1,Nsims);
for j = 1:Nsims
    fprintf(1,'Reading File # %d / %d\n',j,Nsims);
    fname = flist(j).name;   
    [d(j),p(j),level(j),m(j)] = getResults(fname);
end

    