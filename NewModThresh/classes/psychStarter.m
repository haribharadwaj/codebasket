classdef psychStarter < handle
    % A class to initialize PTB and get screen parameters, disable warning,
    % get white/black index, intialize DAQ device etc. Different paradigms
    % could potentially inherit from this and implement more specialized
    % starters.
    %
    % INITIALIZATION:
    %   PS = psychStarter(useTDT, screenDist, screenWidth, useTrigs);
    %   screenDist and screenWidth should be in same units. useTDT and
    %   useTrigs are 0 or true.
    %
    %--------------------------
    % Hari Bharadwaj, September 6, 2010
    %--------------------------
    properties (SetAccess = private)
        window   % The window on which everything is displayed
        rect  % The dimensions of the window in pixel units
        white
        black
        grey
        red
        blue
        green
        intensityFactor
        centre
        onePercentOfScreen % 1 percent of screen width in pixel units
        oneDeg % 1 degree of visual angle in pixel units
        RP % TDT circuit object
        di % DAQ device port
        % Remembering settings before starter being initialized
        f1 % Figure Number for ActiveX
        oldVisualDebugLevel
        oldSupressAllWarnings
        
        
    end
    
    events
        PTBnotGood
    end
    
    methods
        function PS = psychStarter(useTDT, screenDist, screenWidth, useTrigs,FsampTDT,whichScreen)
            
            % TODO: Set defaults if inputs not given
            
            try
                [versionString, versionStructure] = PsychtoolboxVersion; %#ok<*ASGLU>
                if (versionStructure.major < 3)
                    notify(PS, 'PTBnotGood');
                    fprintf(1,'\n ERROR: Cannot run starter!! \n Psychtoolbox version too old! Please install version 3.xxx or newer!\n') %#ok<*PRTCAL>
                    return;
                end
                
            catch %#ok<CTCH>
                notify(PS, 'PTBnotGood');
                fprintf(1,'\n ERROR: Cannot find Psychtoolbox version number! Probably PTB not installed??\n')
                return;
            end
            
            
            try
                % ---------- Window Setup ----------
                % Screen is able to do a lot of configuration and performance checks on
                % open, and will print out a fair amount of detailed information when
                % it does.  These commands supress that checking behavior and just let
                % the demo go straight into action.  See ScreenTest for an example of
                % how to do detailed checking.
                PS.oldVisualDebugLevel = Screen('Preference', 'VisualDebugLevel', 3);
                PS.oldSupressAllWarnings = Screen('Preference', 'SuppressAllWarnings', 1);
                % Find out how many screens and use largest screen number.
                
                if(~exist('whichScreen','var'))
                    whichScreen = max(Screen('Screens'));
                end
                PS.rect=Screen(whichScreen,'Rect'); % Gets Screen dimensions
                HideCursor;
                
                % Opens a graphics window on the main monitor (screen 0).  If you have
                % multiple monitors connected to your computer, then you can specify a
                % different monitor by supplying a different number in the second
                % argument to OpenWindow, e.g. Screen('OpenWindow', 2).
                PS.window = Screen('OpenWindow', whichScreen);
                
                
                
                % ---------- Color Setup -----------
                % Retrieves color codes for black and white and grey.
                PS.intensityFactor = 1; % CHANGE TO CHANGE INTENSITY
                PS.black = BlackIndex(PS.window);  % Retrieves the CLUT color code for black.
                PS.white = WhiteIndex(PS.window)*PS.intensityFactor;  % Retrieves the CLUT color code for white.
                PS.grey = (PS.black + PS.white) / 2;  % Computes the CLUT color code for grey.
                if round(PS.grey) == PS.white
                    PS.grey = PS.black;
                end
                PS.red = [PS.white 0 0];
                PS.green = [0 PS.white 0];
                PS.blue = [0 0 PS.white];
                
                %--------- Text Setup --------------
                txtSz=40;
                Screen('FillRect',PS.window,PS.black);
                Screen(PS.window,'TextSize',txtSz);
                Screen('DrawText',PS.window,'Starting up ...', ...
                    PS.rect(3)/3,PS.rect(4)/2,PS.white);
                Screen('Flip',PS.window);
                
                %------- Trigger Setup -------------
                if(useTrigs && ~useTDT)
                    PS.di = DaqDeviceIndex; % DaqDevice init only for MEG and MRI
                    port_A = 0;   % Select DIO port A, (A = 0; B = 1)
                    direction_A = 0; % Set port A as output (output = 0; input = 1)
                    DaqDConfigPort(PS.di,port_A,direction_A); % Configure DIO port A0-A7
                    DaqDOut(PS.di,0,0);
                end
                KbName('UnifyKeyNames');
                
                %------- Other properties ----------
                
                PS.centre = [PS.rect(3)/2 PS.rect(4)/2];
                PS.onePercentOfScreen = ceil(PS.rect(3)/100);
                PS.oneDeg = ceil((pi/180)*screenDist*PS.rect(3)/screenWidth);
                
                if(useTDT)
                    % ####### CHECK ########## What are f1, FS, RP?
                    % Perform basic initialization of the TDT
                    Noise_Amp = 10^(-40/20); % This represents 20dB SNR ratio (CHECK!)
                    IAC = -1;
                    if(exist('FsampTDT','var'))
                        FS_tag=FsampTDT;
                    else
                        FS_tag = 3; %Default is 48828
                    end
                    fig_num=99;
                    USB_ch=1;
                    [PS.f1,PS.RP,FS]=load_play_circuit(...
                        FS_tag,fig_num,USB_ch,Noise_Amp,IAC); %#ok<NASGU>
                else
                    % Perform basic initialization of the sound driver
                    InitializePsychSound;
                end
                
            catch %#ok<CTCH>
                
                
                Screen('CloseAll');
                
                % Restores the mouse cursor.
                ShowCursor;
                
                % Restore preferences
                Screen('Preference', 'VisualDebugLevel', PS.oldVisualDebugLevel);
                Screen('Preference', 'SuppressAllWarnings', PS.oldSupressAllWarnings);
                
                % To see error description.
                psychrethrow(psychlasterror);
            end
            
        end
        
    end
end