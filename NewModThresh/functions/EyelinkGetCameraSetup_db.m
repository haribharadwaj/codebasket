function EyelinkGetCameraSetup_db


    
try
    
    %%%%Screen('Preference', 'SkipSyncTests', 1);

    %% OPEN DISPLAY WINDOW, SET BASIC PARAMETERS
    screenNumber=max(Screen('Screens'));
    [w, screenRect]=Screen('OpenWindow', screenNumber, [67, 67, 67], [],32,2);
    [width, height]=Screen('WindowSize', screenNumber);
    
    disp('hello')
    
    % Initialize 'el' eyelink struct with proper defaults for output to
    % window 'w':
    el=EyelinkInitDefaults(w);
    
    % Initialize Eyelink connection (real or dummy). The flag '1' requests
    % use of callback function and eye camera image display:
    if ~EyelinkInit([], 1)
        fprintf('Eyelink Init aborted.\n');
        cleanup;
        return;
    end
    
    %Eyelink('TestSuite');
    
    result = Eyelink('StartSetup',1);
    
    % Perform drift correction: The special flags 1,1,1 request
    % interactive correction with video display:
    % You have to hit esc before return.
    %result = Eyelink('DriftCorrStart',30,30,1,1,1);
    
            
            [keyIsDown,timeSecs,keyCode] = KbCheck(-1);
            char = KbName(keyCode);
                            %if strcmp(char, 'space') || strcmp(char, 'Space')
                            if keyIsDown==1
                                cleanup
                            end

catch
    cleanup
end

cleanup

end

% Cleanup routine:
function cleanup
    Eyelink('Shutdown');
    sca;
    ListenChar(0);
end
