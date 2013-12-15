function [Rcn,Ric] = NelsonCarney2004CNIC(anrate,fs)


Acn=1.5; % Overall scale factor
Aic=1;
Scn=0.6; % Inhibition to excitation ratio
Sic=1.5;
Dcn=1e-3; % Delay
Dic=2e-3;
Tex=0.5e-3; % Time Constant
Tin=2e-3;
t=0:(1/fs):10*Tin;
dt = 1/fs;



% The 1/T^2 is there to normalize the area under the PSPs
Inhcn = padarray(Scn*(1/Tin^2)*(t).*exp(-(t)/Tin),[0 round(Dcn*fs)],0,'pre');
Inhcn = Inhcn(1:numel(t));
Inhic = padarray(Sic*(1/Tin^2)*(t).*exp(-(t)/Tin),[0 round(Dic*fs)],0,'pre');
Inhic = Inhic(1:numel(t));

Excn = (1/Tex^2)*t.*exp(-t/Tex);
Exic = (1/Tex^2)*t.*exp(-t/Tex);

Rcn = dt*Acn*(conv(Excn,anrate)-conv(Inhcn,anrate));
Rcn = Rcn(1:numel(anrate));
Rcn = Rcn.*(Rcn>0);

Ric = dt*Aic*(conv(Exic,Rcn)-conv(Inhic,Rcn));
Ric = Ric(1:numel(anrate));
Ric = Ric.*(Ric>0);
