function WaitKeyInputs


if ~IsOSX
    KbWait([], 3);
else
    fprintf(2,'\nWARNING! The OSX case does nothing!\n');
end

