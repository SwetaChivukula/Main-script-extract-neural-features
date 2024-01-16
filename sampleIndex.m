% sampling indices function

function [stimIndex, index] = sampleIndex(rawData)
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
allData = (double(rawData.Data)./4).*10^(-6);
neuralData = allData([elecChannelSorted(:,2);129;130],:);

%% Determine duration of stim artifact
thresh = 0.0000133; % determined from plot of electrode 129
for electrode = 129
    stimDuration = timeAll(neuralData(electrode,:) >= thresh);
end

%% Use monitor port to sample and hold each pulse
clear stimIndex
clear index

threshA = 0.0000133; % determined from plot of electrode 129
threshB = 0.0069129;
count = 1;

% detect all pulses
for electrode = 129
    for i = 1:length(timeAll)
        if neuralData(electrode, i) >= threshA && neuralData(electrode, i) <= threshB
            stimIndex(count) = timeAll(i);
            index(count) = i;
            count = count + 1;
        end
    end
end

% extra detection for pulses that are within neural signal
threshC = 0.0000060; % determined from plot of electrode 129
threshD = 0.0000134;
% do not reset count value - append to end of stimIndex and index vectors
stimExtra = [];
stimExtra(1) = stimIndex(1);
stimExtra(2) = stimIndex(end);

for electrode = 129
    for i = 1:length(timeAll)
        if stimExtra(1) <= timeAll(i) && stimExtra(end) >= timeAll(i)
            if neuralData(electrode, i) >= threshC && neuralData(electrode, i) <= threshD
                stimIndex(count) = timeAll(i);
                index(count) = i;
                count = count + 1;
            end
        end
    end
end

end