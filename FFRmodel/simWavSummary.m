% clear all;
% close all hidden;
% clc;
% 
% noisedblist = 50:5:90;
% xcommon = -70:10;
% h_all = zeros(numel(noisedblist), numel(xcommon));
% hnull = 0;
% 
% datmat = [];
% 
% for k = 1:numel(noisedblist)
%     fname = strcat('./RESULTS/simNoisy_trial1_noise', num2str(noisedblist(k)),'dB.mat');
%     load(fname);
%     hnull = hnull + interp1(x1, h1, xcommon, 'linear');
%     h_all(k, :) = interp1(x2, h2, xcommon, 'linear');
%     datmat = [datmat; SNR(1,:)'];
%     noise = [noise; noisedblist(k)*zeros(numel(SNR(1,:)),1)];
%     datmat = [datmat; SNR(2,:)'];
%     noise = [noise; noisedblist(k)*ones(numel(SNR(1,:)),1)];
% end
% 
% hnull = hnull/numel(noisedblist);
% 
% figure;
% plot(xcommon, hnull, 'k', 'linew', 2);
% hold on;
% plot(xcommon, h_all, 'linew', 2);
% xlabel('Modulation Domain SNR (dB)', 'FontSize', 20);
% ylabel('Density', 'FontSize', 20);
% 
% 
% M = [datmat, noise];


clear all;
close all hidden;
clc;

noisedblist = 50:5:90;
datmat = [];
noise = [];
for k = 1:numel(noisedblist)
    fname = strcat('./RESULTS/simNoisy_trial1_noise', num2str(noisedblist(k)),'dB.mat');
    load(fname);
    for fac = 1:2
        for l = 1:19
            temp = squeeze(S(fac, :, :, l));
            indsig = 8;
            indnoise = [5,6,10,11];
            SNR_new = 10*log10((temp(:, indsig).^2)./mean(temp(:, indnoise).^2,2));
            datmat = [datmat;  SNR_new(:)];
            noise = [noise; (fac-1)*noisedblist(k)*ones(numel(SNR_new),1)];
        end
    end
end
M = [datmat, noise];