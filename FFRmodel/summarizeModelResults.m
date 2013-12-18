% Summarizing the model results
clear all;
close all;
clc;

if(~exist('resultsSummary.mat','file'))
    flist = dir('./RESULTS/*.mat');
    
    Nsims = numel(flist);
    for j = 1:Nsims
        fprintf(1,'Reading File # %d / %d\n',j,Nsims);
        fname = strcat('./RESULTS/',flist(j).name);
        [d(j,:),p(j,:),level(j),m(j)] = getResults(fname);
    end
else
    load('./resultsSummary.mat');
    
    
    figure;
    lev80_ind = (level == 80);
    plot(d(1,:),p(lev80_ind,:),'o-','linew',2)
    legend(num2str(m(lev80_ind)))
    legend(num2str(m(lev80_ind)'))
    set(gca,'FontSize',20);
    figure;
    m60_ind = (m == 60);
    plot(d(1,:),p(m60_ind,:),'o-','linew',2)
    legend(num2str(level(m60_ind)'))
    set(gca,'FontSize',20);
end