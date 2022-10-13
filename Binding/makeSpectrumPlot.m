[x, fs] = audioread('coh18sample_ERB.wav');
[Pxx, f] = pmtm(x, 1.5, 2^ceil(log2(numel(x))), fs);
semilogx(f, pow2db(Pxx) - max(pow2db(Pxx)), 'k', 'linew', 2);
xlim([100, 10000]);
set(gca, 'XTick', [0.25, 0.5, 1, 2, 4, 8]*1000, 'Xticklabel', {'0.25', '0.5', '1', '2', '4', '8'});
set(gca, 'FontSize', 20);
xlabel('Frequency (kHz)', 'FontSize', 20);
ylabel('Relative Level (dB)', 'FontSize', 20);
ylim([-100, 0]);