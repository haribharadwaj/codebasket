function resp = getResponseKb(dur)
% USAGE:
%  resp = getResponseKb(dur);
%  OR
% resp = getResponseKb;
%
%  If dur is given (in seconds), it will check for responses for 'dur'
%  seconds...if no argument is given, it will wait for  nresp presses
%
% nresp is hardcoded in the first line of the code

nresp = 1;
if(~exist('dur','var'))
    dur = [];
    waitflag = 1;
else
    waitflag = 0;
end



rate = 0.05;



if(isempty(dur))
    dur = 0.5;
end

resp = [];
carry = 0;
while(numel(resp) < nresp)
    nquery = floor(dur/rate);
    resptemp  = zeros(1,nquery);
    
    for j = 1:nquery
        %     FlushEvents('keyDown'); %This required some stupidass DLL
        [keyIsDown, secs, keyCode, deltaSecs] = KbCheck; %#ok<NASGU,ASGLU>
        if(keyIsDown)
            key = KbName(keyCode);
            
            %Shift options are also returned for the label => need to select
            key = key(1);
            if(iscell(key))
                key = '0';
            end
            
            if((key == '1') || (key == '2'))
                
                resptemp(j) = str2num(key); %#ok<ST2NM>
            else
                fprintf(2,'\n WARNING! Pressed key is not within response set: First label character is %s\n',key);
            end
        end
        WaitSecs(rate);
    end
    resp = [resp resptemp(diff([carry,resptemp])>0)]; %#ok<AGROW>
    carry = resptemp(end);
    if(~waitflag)
        break;
    end
end


