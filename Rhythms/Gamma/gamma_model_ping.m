%% Simple PING gamma script

 I0p=20;
 gP=0.2;
 gI=0.5;
 tauI=10;
 T0=100;
 [VP,sP,VI,sI,t] = ping(I0p,gP,gI,tauI,T0);
 
 % Plot Results
 figure(1);
 subplot(2,1,1)
 plot(t,VP);  xlabel('Time [ms]');  ylabel('P-cell [mV]');
 subplot(2,1,2)
 plot(t,VI);  xlabel('Time [ms]');  ylabel('I-cell [mV]');