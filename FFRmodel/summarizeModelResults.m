% Summarizing the model results
clear all;
close all;
clc;

reobtain = 0;
if(~exist('resultsSummaryUnnorm.mat','file') || reobtain)
    flist = dir('./RESULTS/*.mat');
    
    Nsims = numel(flist);
    for j = 1:Nsims
        fprintf(1,'Reading File # %d / %d\n',j,Nsims);
        fname = strcat('./RESULTS/',flist(j).name);
        [d(j,:),p(j,:),level(j),m(j)] = getResults(fname);
    end
else
    load('./resultsSummaryUnnorm.mat');
    
end

% Plot iso-m curves as is
figure;
lev80_ind = (level == 80);
d = d(1,:);
p_mod = p(lev80_ind,:);
symblist = 'soxd*<>';
[sorted_m,sorted_ind] = sort(m(lev80_ind));

for j = 1:numel(sorted_ind)
    plot(d,db(p_mod(sorted_ind(j),:)),strcat('k-',symblist(j)),...
        'linew',2,'MarkerSize',8);
    hold on;
end
legend(num2str(sorted_m'))
xlabel('% Neuropathy','FontSize',20);
ylabel('EFR (arbitrary units)','FontSize',20);
set(gca,'FontSize',20);


% Plot average normalized iso-m curves
figure;
lev80_ind = (level == 80);
d = d(1,:);
p_mod = p(lev80_ind,:);
symblist = 'soxd*<>';
[sorted_m,sorted_ind] = sort(m(lev80_ind));

for j = 1:numel(sorted_ind)
    plot(d,db(p_mod(sorted_ind(j),:)/mean(p_mod(sorted_ind(j),:))),strcat('k-',symblist(j)),...
        'linew',2,'MarkerSize',8);
    hold on;
end
legend(num2str(sorted_m'))
set(gca,'FontSize',20);
xlabel('% Neuropathy','FontSize',20);
ylabel('EFR (dB re: average)','FontSize',20);
set(gca,'FontSize',20);

% Plot self-normalized iso-m curves
figure;
lev80_ind = (level == 80);
d = d(1,:);
p_mod = p(lev80_ind,:);
symblist = 'soxd*<>';
[sorted_m,sorted_ind] = sort(m(lev80_ind));

for j = 1:numel(d)
    plot(sorted_m,db(p_mod(sorted_ind,j)/(p_mod(sorted_ind(end),j))),strcat('k-',symblist(j)),...
        'linew',2,'MarkerSize',8);
    hold on;
end
legend(num2str(d'))
set(gca,'FontSize',20);
xlabel('Modulation Depth','FontSize',20);
ylabel('EFR (dB re: 100% modulation)','FontSize',20);
set(gca,'FontSize',20);

% Get sensitivity versus m
for j = 1:numel(sorted_m)
    sens_m(j) = -mean(diff(db(p_mod(sorted_ind(j),:)))); %#ok<SAGROW>
end


% Plot iso-level curves as is
figure;
m60_ind = (m == 60);
p_level = p(m60_ind,:);
[sorted_level,sorted_ind] = sort(level(m60_ind));
for j = 1:numel(sorted_ind)
    plot(d,db(p_level(sorted_ind(j),:)),strcat('k-',symblist(j)),...
        'linew',2,'MarkerSize',8);
    hold on;
end
legend(num2str(sorted_level'))
set(gca,'FontSize',20);
xlabel('% Neuropathy','FontSize',20);
ylabel('EFR (arbitrary units)','FontSize',20);
set(gca,'FontSize',20);

% Plot normalized iso-level curves as is
figure;
m60_ind = (m == 60);
p_level = p(m60_ind,:);
[sorted_level,sorted_ind] = sort(level(m60_ind));
for j = 1:numel(sorted_ind)
    plot(d,db(p_level(sorted_ind(j),:)/mean(p_level(sorted_ind(j),:))),strcat('k-',symblist(j)),...
        'linew',2,'MarkerSize',8);
    hold on;
end
legend(num2str(sorted_level'))
xlabel('% Neuropathy','FontSize',20);
ylabel('EFR (dB re: average)','FontSize',20);
set(gca,'FontSize',20);

% Plot with level as axis:
figure; 
for j = 1:numel(d)
plot(sorted_level,db(p_level(sorted_ind,j)),strcat('k-',symblist(j)),...
'linew',2,'MarkerSize',8);
hold on;
end
legend(num2str(d'))
set(gca,'FontSize',20);
xlabel('Sound Level (dB SPL)','FontSize',20);
ylabel('EFR (arbirtray units)','FontSize',20);


% Get sensitivity versus level
for j = 1:numel(sorted_level)
    sens_level(j) = -mean(diff(db(p_level(sorted_ind(j),:)))); %#ok<SAGROW>
end


% Plot sensitivity
figure;
subplot(1,2,1);
plot(sorted_m,sens_m,'k--','linew',2);
xlabel('Modulation Depth (%)','FontSize',20);
ylabel('Sensitivity (dB drop per 20% neuropathy)','FontSize',20);
set(gca,'FontSize',20);
subplot(1,2,2);
plot(sorted_level,sens_level,'k--','linew',2);
xlabel('Sound Level (dB SPL)','FontSize',20);
set(gca,'FontSize',20);

