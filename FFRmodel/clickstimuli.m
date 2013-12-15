function [stimulus,indcl]=clickstimuli(FS)
stimulus=[0];
for d=1:100
                dur=(randi([90 210]))./1000;
                %dur=500./1000;
                num_samplesforclick=floor(FS*80e-6); %duration of click about 80 microseconds
                num_samplesp=floor(FS*dur);
                stimp=[];
                stimp(1:num_samplesp)=0;
                stim=[];
                stim(1:num_samplesforclick)=0.95;
                %stim(1)=1;
                stimulus=[stimulus,stim,stimp];
end
            
%% find onset of each click
indcl=[];
ind=find(stimulus==max(stimulus));
for i=1:length(ind)
    if (stimulus(ind(i)-1)==0)
        indcl=[indcl;ind(i)];
    end
end 

            %StimGain=sqrt(2)*(10^((80-Lcal)/20))/(0.95); %Gain value of click
            %stimulus(2,:)=0;