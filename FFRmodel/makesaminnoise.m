function x = makesaminnoise(fc,fm,m,fs,dur,rise,phi,bw,asymm,SNR)

x = makesamtone(fc,fm,m,fs,dur,rise,phi);
noise = makeNotchNoiseFFT(bw,asymm,fc,dur,fs,rise,0);

stimrms = sqrt(mean(x.^2));
noiserms = sqrt(mean(noise.^2));
noise = (noise'*stimrms/noiserms)*db2mag(-1*SNR);

x = x + noise;

