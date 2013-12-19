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
    d = d(1,:);
    p_mod = p(lev80_ind,:);
    symblist = 'soxd*<>';
    [sorted_m,sorted_ind] = sort(m(lev80_ind));
    
    for j = 1:numel(sorted_ind)
        plot(d,p_mod(sorted_ind(j),:),strcat('k-',symblist(j)),...
            'linew',2,'MarkerSize',8);
        hold on;
    end
    legend(num2str(sorted_m'))
    set(gca,'FontSize',20);
    figure;
    
    m60_ind = (m == 60);
    p_level = p(m60_ind,:);
    [sorted_level,sorted_ind] = sort(level(m60_ind));
   for j = 1:numel(sorted_ind)
        plot(d,p_level(sorted_ind(j),:),strcat('k-',symblist(j)),...
            'linew',2,'MarkerSize',8);
        hold on;
    end
    legend(num2str(sorted_level'))
    set(gca,'FontSize',20);
end