%% Script to test paradigm classes
clear classes
try
    PS = psychStarter(0,0.4,0.3,0);
    tic
    for j = 1:51313
        renderVisFrame(PS,'FIX');
        Screen('Flip',PS.window);
        WaitSecsMex(0.3);
        if(rand > 0.5)
            renderVisFrame(PS,'CUER');
            fprintf(1,'\nRight!\n');
        else
            renderVisFrame(PS,'CUEL');
            fprintf(1,'\nLeft!\n');
        end
        Screen('Flip',PS.window);
        WaitSecsMex(0.15);
        renderVisFrame(PS,'FIX');
        Screen('Flip',PS.window);
        WaitSecsMex(0.5);
%         AlphaTest;
        WaitSecsMex(0.5);
        renderVisFrame(PS,'RESP');
        Screen('Flip',PS.window);
        WaitSecsMex(1.2);
        
    end
    toc
    sca;
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
