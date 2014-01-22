clear classes;
clc;
global globalBlocks
global currBlock
global globalBonus
globalBlocks = 8;
currBlock = 0;
globalBonus = 0;

sID = input('Enter Subject ID: ','s');
fprintf(1,'\n ### Entered ID is: %s\n', sID);

useButtons = input('UseButtonBox (Enter 1 for yes and 0 for no): ');
if(useButtons)
    fprintf(1,'\n ### Using ButtonBox\n');
else
    fprintf(1,'\n ### NOT Using ButtonBox\n');
end

nBlocksPerLevel = 2;
bw = 3000;

% Level 45
level = 45;
[resp, m ,thresh, bonus] = ...
    getModThresh(sID,bw,level,nBlocksPerLevel,1,0.4,0.4,useButtons);
globalBonus = globalBonus + bonus;
currBlock = currBlock + 2;

while(std(thresh,1) > 0.05)
    globalBlocks = globalBlocks + 2;
    
    [resp, m ,thresh_extra, bonus] = ...
        getModThresh(sID,bw,level,nBlocksPerLevel,1,0.4,0.4,useButtons);
    globalBonus = globalBonus + bonus;
    currBlock = currBlock + 2;
    thresh = [thresh, thresh_extra];
end

% Level 55
level = 55;
[resp, m ,thresh, bonus] = ...
    getModThresh(sID,bw,level,nBlocksPerLevel,1,0.4,0.4,useButtons);
globalBonus = globalBonus + bonus;
currBlock = currBlock + 2;

while(std(thresh,1) > 0.05)
    globalBlocks = globalBlocks + 2;
    
    [resp, m ,thresh_extra, bonus] = ...
        getModThresh(sID,bw,level,nBlocksPerLevel,1,0.4,0.4,useButtons);
    globalBonus = globalBonus + bonus;
    currBlock = currBlock + 2;
    thresh = [thresh, thresh_extra];
end

% Level 65
level = 65;
[resp, m ,thresh, bonus] = ...
    getModThresh(sID,bw,level,nBlocksPerLevel,1,0.4,0.4,useButtons);
globalBonus = globalBonus + bonus;
currBlock = currBlock + 2;

while(std(thresh,1) > 0.05)
    globalBlocks = globalBlocks + 2;
    
    [resp, m ,thresh_extra, bonus] = ...
        getModThresh(sID,bw,level,nBlocksPerLevel,1,0.4,0.4,useButtons);
    globalBonus = globalBonus + bonus;
    currBlock = currBlock + 2;
    thresh = [thresh, thresh_extra];
end

% Level 75
level = 75;
[resp, m ,thresh, bonus] = ...
    getModThresh(sID,bw,level,nBlocksPerLevel,1,0.4,0.4,useButtons);
globalBonus = globalBonus + bonus;
currBlock = currBlock + 2;

while(std(thresh,1) > 0.05)
    globalBlocks = globalBlocks + 2;
    
    [resp, m ,thresh_extra, bonus] = ...
        getModThresh(sID,bw,level,nBlocksPerLevel,1,0.4,0.4,useButtons);
    globalBonus = globalBonus + bonus;
    currBlock = currBlock + 2;
    thresh = [thresh, thresh_extra];
end