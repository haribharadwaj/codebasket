% Info to pick up
sID = 'Hari_right';
f2 = 4000;

% Typical parameters of recording.. Change if changed during recording
Fs = 48828.125;
L2_list = 0:6:76;
paraDir = 'C:\Experiments\DPOAE_ER10C\';
fdir = strcat(paraDir,'/subjResponses/',sID);

fignum=2;
figure(fignum);
for k = 1: numel(L2_list)
    L2 = L2_list(k);
    fname = strcat(fdir,'/',sID,'_DPOAE_RAW_',num2str(L2),'_',num2str(round(f2)),'_*.mat');
    names = dir(fname);
    if numel(names) ~= 1
        error('Multiple files for the same subject and input level!\n');
    end
    load(strcat(fdir, '/', names.name));
    % Calculate noise level
    good = [];
    for trial = 1:size(OAE,1)
        if ~any(OAE(trial, :) > (median(OAE, 1) + 20*mad(OAE,1)))
            good = [good, trial]; %#ok<AGROW>
        end
    end
    if mod(numel(good),2)~=0
        good = good(1:(end-1));
    end
    OAE = OAE(good, :);
    
    OAE_noise = OAE;
    OAE_noise(2:2:end, :) = -1 * OAE_noise(2:2:end, :);
    noise = mean(OAE_noise, 1);
    signal = mean(OAE, 1);
    w = dpss(numel(noise),1,1)';
    w = w / sum(w);
    f_dp = 2*f1 - f2;
    t_calc = (0:(numel(noise)-1))/Fs;
    wsn = w.*sin(2*pi*f_dp*t_calc);
    wcn = w.*cos(2*pi*f_dp*t_calc);
    
    mic_sens = 50e-3; % V/Pa
    mic_gain = db2mag(40);
    V_to_microV = 1e6;
    ref = 20; % microP
    factor = V_to_microV / (mic_sens * mic_gain * ref);
    % Super crude noise-floor estimate by flipping half the trials
    noisefloor(k) = db(sqrt(sum(wsn.*noise)^2 + sum(wcn.*noise)^2) * factor); %#ok<SAGROW>
    DP(k) = db(sqrt(sum(wsn.*signal)^2 + sum(wcn.*signal)^2) * factor); %#ok<SAGROW>
    
    figure(fignum);
    hold on;
    plot(L2, noisefloor(k), 'ok', 'linew', 2);
    hold on;
    plot(L2, DP(k), 'or', 'linew', 2);
    xlabel('L2 dB SPL', 'FontSize', 16);
    % Can use probe calibration offline to specify in SPL
    % This is just to see if there is DPOAE above noise floor and if it is
    % growing roughly like it should
    ylabel('DPOAE (dB SPL)', 'FontSize', 16);
    drawnow;
end
set(gca, 'FontSize', 16);
title(['Subject ID: ', strrep(sID, '_', ' - ')], 'FontSize', 16);