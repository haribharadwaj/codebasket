function [respListAll, mlistAll, thresh, bonus] = ...
    getModThresh(sID,bw,soundLevel,nBlocks,useTDT,screenDist,...
    screenWidth,buttonBox)

% ########## DO NOT RUN THIS DIRECTLY: use modscript.m instead ##########
% USAGE (if you insist on calling this):
% [respList, mlist] = getThreshNotch(sID,bw,asymm,soundLevel,nBlocks,useTDT,screenDist,screenWidth,buttonBox)

global globalBlocks
global currBlock
global globalBonus


paraDir = '/home/hari/codebasket/NewModThresh';
% whichScreen = 1;
addpath(genpath(paraDir));
if(~exist(strcat(paraDir,'/subjResponses/',sID),'dir'))
    mkdir(strcat(paraDir,'/subjResponses/',sID));
end
respDir = strcat(paraDir,'/subjResponses/',sID,'/');

FsampTDT = 3; % 48828 Hz
useTrigs = 0;
feedback = 1;

feedbackDuration = 0.2;

PS = psychStarter(useTDT,screenDist,screenWidth,useTrigs,FsampTDT); %,whichScreen);



Nup = 3; % Weighted 1 up-1down with weights of 3:1



if(~exist('bw','var'))
    bw = 3000;
end

NmaxTrials = 200;
NminTrials = 40;
target = (randperm(NmaxTrials) > NmaxTrials/2);


try
    
    if(useTDT)
        % Noise off
        invoke(PS.RP, 'SoftTrg', 4);
        
    end
    bonus = 0;
    
    
    for j = 1:nBlocks
        geometric = 1;
        if(geometric)
            stepDown = 0.8;
            stepUp = (1/stepDown)^Nup;
        else
            stepDown = -0.1;
            stepUp = Nup*(-stepDown);
        end
        m = 1.0; % Starting mod depth, which is adapted
        
        textlocH = PS.rect(3)/4;
        textlocV = PS.rect(4)/3;
        line2line = 50;
        blockNumStr = num2str(currBlock + j);
        totalBlocks = num2str(globalBlocks);
        info = strcat('This is block #',blockNumStr,'/',totalBlocks,'...');
        Screen('DrawText',PS.window,info,textlocH,textlocV,PS.white);
        info = strcat('Press any key twice to begin...');
        Screen('DrawText',PS.window,info,textlocH,textlocV+line2line,PS.white);
        Screen('Flip',PS.window);
        WaitKeyInputs;
        WaitKeyInputs;
        tstart = tic;
        
        
        converged = 0;
        respList = [];
        mlist = [];
        trialCount = 0;
        correctCount = 0;
        
        
        while(~converged)
            renderVisFrame(PS,'FIX');
            Screen('Flip',PS.window);
            
            trialCount = trialCount + 1;
            
            if(trialCount == 1)
                WaitSecs(4);
            end
            
            
            level = soundLevel; % For now
            fc = 8000;
            fs = 48828;
            dur = 0.4;
            rampSize = 0.025;
            fm = 23;
            
            sig = makeModThreshStim(fc,fm,m,1,0,bw,fs,dur,rampSize,0);
            dummy = makeModThreshStim(fc,fm,0,1,0,bw,fs,dur,rampSize,0);
            rmssig = sqrt(mean(sig.^2));
            rmsdummy = sqrt(mean(dummy.^2));
            
            sig = sig/rmssig;
            dummy = dummy/rmsdummy;
            
            
            if(target(trialCount))
                % Correct answer is "1"
                answer = 1;
                y = sig;
                z = dummy;
            else
                % Correct answer is "2"
                answer = 2;
                y = dummy;
                z = sig;
            end
            
            y = [y, y];
            z = [z, z];
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Clear Up buffers for 1st stim
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            if(useTDT)
                invoke(PS.RP, 'SoftTrg', 2); %Stop trigger
                %Reset the play index to zero:
                invoke(PS.RP, 'SoftTrg', 5); %Reset Trigger
                %Clearing I/O memory buffers:
                invoke(PS.RP,'ZeroTag','datain');
            end
            %-----------------
            % Why 122?
            %-----------------
            % ER-1s give about 105dB SPL for a 1kHz tone with a 1V-rms drive.
            % Max output is +/-10V peak i.e 7.07V-rms which is 17 dB higher.
            % Thus 122 dB-SPL is the sound level for tones when they occupy full
            % range.
            
            % Full range in MATLAB for a pure tone is +/- 1 which is 0.707 rms and
            % that corresponds to 122 dB SPL at the end. So if we want a signal
            % with rms sigrms to be x dB, then (122 - x) should be
            % db(sigrms/0.707).  The calculation below is off by 3 dB since it
            % doesn't have the 0.707 factor at the right place => we are actually
            % getting 3dB SPL higher (i.e 78): CHECK! (=> Keep setting on -3dB)
            drop = 122 - level;
            
            %Start dropping from maximum RMS (actual RMS not peak-equivalent)
            wavedata = y'*db2mag(-1*drop); % We want to fix tone SPL
            %-----------------------------------------
            
            
            % The trial flow:
            
            
            if useTDT
                %Load the 2ch variable 'x' into the RP2:
                invoke(PS.RP, 'WriteTagVEX', 'datain', 0, 'I16', (wavedata*2^15));
                %Set the delay of the sound
                invoke(PS.RP, 'SetTagVal', 'onsetdel',100); % onset delay is in ms
                %Start playing from the buffer:
                Screen('DrawText',PS.window,'1',PS.rect(3)/2 - 20,PS.rect(4)/2-20,PS.white);
                Screen('Flip',PS.window);
                invoke(PS.RP, 'SoftTrg', 1); %Playback trigger
            else
                sound(y,fs);
            end
            
            WaitSecs(1.4);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Clear Up buffers for 2nd stim
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            if(useTDT)
                invoke(PS.RP, 'SoftTrg', 2); %Stop trigger
                %Reset the play index to zero:
                invoke(PS.RP, 'SoftTrg', 5); %Reset Trigger
                %Clearing I/O memory buffers:
                invoke(PS.RP,'ZeroTag','datain');
            end
            
            %-----------------
            % Why 122?
            %-----------------
            % ER-1s give about 105dB SPL for a 1kHz tone with a 1V-rms drive.
            % Max output is +/-10V peak i.e 7.07V-rms which is 17 dB higher.
            % Thus 122 dB-SPL is the sound level for tones when they occupy full
            % range.
            
            % Full range in MATLAB for a pure tone is +/- 1 which is 0.707 rms and
            % that corresponds to 122 dB SPL at the end. So if we want a signal
            % with rms sigrms to be x dB, then (122 - x) should be
            % db(sigrms/0.707).  The calculation below is off by 3 dB since it
            % doesn't have the 0.707 factor at the right place => we are actually
            % getting 3dB SPL higher (i.e 78): CHECK! (=> Keep setting on -3dB)
            drop = 122 - level;
            
            %Start dropping from maximum RMS (actual RMS not peak-equivalent)
            wavedata = z'*db2mag(-1*drop); % We want to fix tone SPL
            %-----------------------------------------
            
            if useTDT
                %Load the 2ch variable 'x' into the RP2:
                invoke(PS.RP, 'WriteTagVEX', 'datain', 0, 'I16', (wavedata*2^15));
                %Set the delay of the sound
                invoke(PS.RP, 'SetTagVal', 'onsetdel',100); % onset delay is in ms
                %Start playing from the buffer:
                Screen('DrawText',PS.window,'2',PS.rect(3)/2-20,PS.rect(4)/2-20,PS.white);
                Screen('Flip',PS.window);
                invoke(PS.RP, 'SoftTrg', 1); %Playback trigger
            else
                sound(z,fs);
            end
            
            WaitSecs(1);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %  Response Frame
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            renderVisFrame(PS,'RESP');
            Screen('Flip',PS.window);
            if(buttonBox)
                resp = getResponse(PS.RP);
            else
                resp = getResponseKb;
            end
            fprintf(1,'\n Target = %s, Response = %s',num2str(answer),num2str(resp));
            if((numel(resp)>=1) && ((answer - resp(end)) == 0))
                fprintf(1,'..which is correct!\n');
                respList = [respList, 1]; %#ok<AGROW>
                correct = 1;
                mlist = [mlist, m];
            else
                fprintf(1,'..which is Wrong!\n');
                respList = [respList, 0]; %#ok<AGROW>
                correct = 0;
                mlist = [mlist, m];
            end
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %  Feedback Frame
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            if(feedback)
                if(correct)
                    renderVisFrame(PS,'GO');
                    correctCount = correctCount +1;
                    
                else
                    renderVisFrame(PS,'NOGO');
                    
                end
            end
            if(correct)
                if(geometric)
                    m = m*stepDown;
                else
                    m = m + stepDown;
                end
            else
                if(geometric)
                    m = m*stepUp;
                else
                    m = m + stepUp;
                end
            end
            
            
            if( m > 1)
                m = 1.0;
            end
            
            if (m < 0)
                m = 0.0;
            end
            
            Screen('Flip',PS.window);
            WaitSecs(feedbackDuration + rand*0.1);
            
            % Counting Reversals
            revList = [];
            nReversals = 0;
            for k = 3:numel(mlist)
                if((mlist(k-1) > mlist(k)) && (mlist(k-1) > mlist(k-2)))
                    nReversals = nReversals + 1;  revList = [revList, (k-1)];
                end
                if((mlist(k-1) < mlist(k)) && (mlist(k-1) < mlist(k-2)))
                    nReversals = nReversals + 1;  revList = [revList, (k-1)];
                end
            end
            
            if(nReversals >= 2)
                if(geometric)
                    stepDown = 0.9;
                    stepUp = (1/stepDown)^Nup;
                else
                    stepDown = -0.05;
                    stepUp = Nup*(-stepDown);
                end
                
                
            end
            
            if(nReversals >= 4)
                if(geometric)
                    stepDown = 0.95; %#ok<*UNRCH>
                    stepUp = (1/stepDown)^Nup;
                else
                    stepDown = -0.02;
                    stepUp = Nup*(-stepDown);
                end
                
                
            end
            
            if( (nReversals >= 12) && (trialCount > NminTrials))
                stdev = std(mlist(revList((end-7):end)));
                if( stdev < 4)
                    converged = 1;
                end
            end
            
        end
        
        thresh(j) = mean((mlist(revList((end-7):end)))); %#ok<*AGROW>
        mlistAll{j} = mlist;
        respListAll(j,:) = respList;
        % Save respList
        currentBonus = correctCount*0.01;
        bonus = bonus + currentBonus ;
        bonusDisp = bonus + globalBonus;
        fprintf(2,'\n###### THRESHOLD FOR THIS BLOCK IS %f +/- %f\n',...
            thresh(j),stdev);
        toc(tstart);
        datetag = datestr(clock);
        datetag(strfind(datetag,' ')) = '_';
        datetag(strfind(datetag,'-')) = '_';
        datetag(strfind(datetag,':')) = '_';
        fname_resp = strcat(respDir,sID,'_Block_0',num2str(j),...
            '_mod_',datetag,'.mat');
        save(fname_resp,'mlist','respList',...
            'correctCount','thresh','bonus','bw','soundLevel');
        
        
        % Display Reward
        info = strcat('Done with Block #',blockNumStr,'/',totalBlocks);
        Screen('DrawText',PS.window,info,textlocH,textlocV,PS.white);
        
        info = strcat('This block''s bonus is $',num2str(currentBonus));
        Screen('DrawText',PS.window,info,textlocH,textlocV + line2line,PS.white);
        
        info = strcat('Total bonus is $',num2str(bonusDisp));
        Screen('DrawText',PS.window,info,textlocH,textlocV + 2*line2line,PS.white);
        
        if(j==nBlocks)
            info = strcat('Press any key to exit...');
        else
            info = strcat('Press any key to continue...');
        end
        Screen('DrawText',PS.window,info,textlocH,textlocV + 3*line2line,PS.white);
        Screen('Flip',PS.window);
        WaitKeyInputs;
    end
    sca;
    close_play_circuit(PS.f1,PS.RP);
    
catch %#ok<CTCH>
    
    
    Screen('CloseAll');
    
    % Restores the mouse cursor.
    ShowCursor;
    
    % Save stuff
    
    datetag = datestr(clock);
    datetag(strfind(datetag,' ')) = '_';
    datetag(strfind(datetag,'-')) = '_';
    datetag(strfind(datetag,':')) = '_';
    
    fname_resp = strcat(respDir,sID,'_Block_0',num2str(j),...
        '_crash_mod_',datetag,'.mat');
    save(fname_resp,'mlist','respList',...
        'correctCount','bonus','bw','soundLevel');
    % Restore preferences
    Screen('Preference', 'VisualDebugLevel', PS.oldVisualDebugLevel);
    Screen('Preference', 'SuppressAllWarnings', PS.oldSupressAllWarnings);
    close_play_circuit(PS.f1,PS.RP);
    % To see error description.
    psychrethrow(psychlasterror);
end

