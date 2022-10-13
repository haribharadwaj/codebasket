% Mouse has high-frequency hearing + we intend to play 32kHz
fs = 97656.25;
dur = 1.0;
n_epochs = 30;
burst_rate  = 40; % 40 burst for each frequency per second
fc_list = [12.14, 30.49]*1e3;
n_cycles_per_burst = 5;
[x, trains] = makeParallelABRstims(fs,  dur, n_epochs, burst_rate,...
    fc_list, n_cycles_per_burst);

y = reshape(x', [1, numel(x)]);
save mouseCentralGainStim y trains fs...
    fc_list n_epochs n_cycles_per_burst dur;