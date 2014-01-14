
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [f1,RP,FS]=load_play_circuit(FS_tag,fig_num,USB_ch,Noise_Amp,IAC)
% Loads the TDT circuit and makes actx links necessay
% Legacy code from KC. The TDT matlab syntax has now changed to look more
% like OOPS but this old style is still supported.
%
%------------
% Hari Bharadwaj, September 6, 2010
%------------
warning('off'); %#ok<WNOFF>

CIR_PATH='play_noise_kc_serial_buttonin.rcx'; %The *.rco circuit used to play the files

%Generate the actx control window in a specified figure:
%-------------------------------------------------------
f1=figure(fig_num);
set(f1,'Position',[5 5 30 30],'Visible','off'); %places and hides the ActX Window
RP=actxcontrol('RPco.x',[5 5 30 30],f1); %loads the actx control for using rco circuits
invoke(RP,'ConnectRP2','USB',USB_ch); % opens the RP2 device   *** CHANGED 'USB' to 'GB' TMS 6/2005 ***
% (the USB channel may change for other computer configurations)

% The rco circuit can be run at the following set of discrete sampling
% frequencies (in Hz): 0=6k, 1=12k, 2=25k, 3=50k, 4=100k, 5=200k.
% Use the tags listed above to specify below:
%--------------------------------------------------------------------------
invoke(RP,'LoadCOFsf',CIR_PATH,FS_tag); %loads the circuit using the specified sampling freq.
FS_vec=[6 12 25 50 100 200]*1e3;
FS=FS_vec(FS_tag+1);

invoke(RP,'Run'); %start running the circuit

Status = double(invoke(RP,'GetStatus'));
if bitget(double(Status),1)==0
    error('Error connecting to RP2');
elseif bitget(double(Status),2)==0
    error('Error loading circuit')';
elseif bitget(double(Status),3)==0
    error('Error running circuit')';
end

% Start the background noise
if invoke(RP,'GetStatus')==7  % RP2 circuit status success
    % Set the level of the backgound noise
    invoke(RP, 'SetTagVal', 'noiselev',Noise_Amp);
    invoke(RP, 'SetTagVal', 'phase', IAC);
    invoke(RP, 'SoftTrg', 3); %Playback noise trigger
else
    error('RP2 circuit status failure')
end
end % of load_play_circuit
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%