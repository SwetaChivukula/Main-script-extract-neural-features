%sample and hold function

function [holdData] = sampleAndHold(rawData, allData, stimTimeIndex, stimIndex, channelSet)
%% Organize the time
fs = 30000;
dt = 1/fs;
%TendSamples = length(NS5.Data(1,:));
%TendSec = TendSamples/fs;
TendSec = rawData.MetaTags.DataDurationSec;
TendSamples = rawData.MetaTags.DataPoints;
TstartSec = 0;
TstartSamples = 1;
timeAll = TstartSec:dt:TendSec-dt;
%% Map the neural data from channel-indexing to electrode-indexing
electrodeNumber = zeros(128,1);
zeroChar = double('0');
for i = 1:128
   elecChar =  rawData.ElectrodesInfo(i).Label(7:9);
   elecDouble = double(elecChar);
   elecDouble2 = elecDouble(elecDouble ~= 0);
   numDigits = size(elecDouble2,2);
   if(numDigits == 1)
       elecNumber = elecDouble2-zeroChar;
   elseif(numDigits == 2)
       elecNumber = sum([10 1] .* (elecDouble2-zeroChar));
   elseif(numDigits == 3)
       elecNumber = sum([100 10 1] .* (elecDouble2-zeroChar));
   else
       disp('Error')
   end
   electrodeNumber(i) = elecNumber;
end
elecChannel = [electrodeNumber,[1:128]'];
elecChannelSorted = sortrows(elecChannel,1);

neuralData = allData; % need to figure out how to map to correct electrode


%% determine the sampling value


chanLength = length(channelSet);

% try to index across electrodes    
for electrodeNum = 1:chanLength
    sample = 0;
    holdData(electrodeNum,:) = neuralData(electrodeNum,:);
    count = 1;
    for counter = 1:length(timeAll)
        if timeAll(1,counter) == stimTimeIndex(1,count)
            if sample == 0
              sample = holdData(electrodeNum,(counter-1)); % determine sample value
            end
        end
    end
    
    % hold for 700 to 800 us to cut out the pulse - try to reduce number of
    % for loops to increase efficiency
    for counting = 1:length(stimTimeIndex) % hold for 700 - 800 us (~ 0.0007)
        clear count;
        count = stimIndex(counting);
        for i = 1:30 % duration of stim from elec 129 - logical indexing, use function diff, repmat
            holdData(electrodeNum,count) = sample;
            count = count + 1;
         end
    end
end
end