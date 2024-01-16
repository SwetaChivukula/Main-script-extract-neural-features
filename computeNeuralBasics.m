function neuralDataOutput = computeNeuralBasics(rawData,channelSet,binSize,filterRange,threshold,excludeChannels, stimTimeIndex, stimIndex)
%COMPUTENEURALBASICS is a function designed to extract the threshold
%crossings and the spike power of an NS5 file
%   RawData is the NS5 file, channelSet is a range of channel indices,
%   binSize is the size of the time bins in milliseconds for the binned
%   data, filterRange is the range of frequency the neural data is filtered
%   by in Hz, and threshold is the threshold value as an RMS multiplier.
%   stimTimeIndex is the vector in which the time values are recorded
%   whenever a stim starts according to each electrode 129. stimIndex is
%   the index of the time vector where stims were detected.

    if(nargin<6)
        excludeChannels = [];
    end
    
    if(~isempty(excludeChannels))
        channelSetMask = ones(size(channelSet));
        excludeElecsSubsetA = (excludeChannels<=max(channelSet));
        excludeElecsSubsetB = (excludeChannels>=min(channelSet));
        excludeElecsSubsetIndices = (excludeElecsSubsetA+excludeElecsSubsetB)==2;
        excludeElecsSubset = excludeChannels(excludeElecsSubsetIndices);
        channelSetMask(excludeElecsSubset) = 0;
        includeChannelsIndex = find(channelSetMask==1);
    else
        includeChannelsIndex = 1:length(channelSet);
    end

    voltageConversion = 1/4; % converts to uV instead of 0.25 uV as the NSP records it
    allData = double(rawData.Data(channelSet,:))'.*voltageConversion;

% part two of sample and hold using stim time index and stim index as
% inputs in compute neural basics *adding two inputs into compute neural
% basics

heldData = sampleAndHold(rawData, allData', stimTimeIndex, stimIndex, channelSet)';


    % Decimate the data using frank willet's function
    allData = decimate_frw(heldData,2);
    decSR = 15000;
    
    % Filter the data in the spike band
    [B,A] = butter(4,2*filterRange/decSR);
    allData = filtfilt(B,A,allData);
    
    % Apply CAR to this subset
    allData = allData - mean(allData(:,includeChannelsIndex),2);
    
    % Compute the root mean square of the data
    rmsVal = rms(allData);
    
    % Get spike band power
    sqSignal = allData.^2;

    nSamples = (binSize/1000)*decSR;% samples per bin
    % StartIdx = 0+1;%rawData.MetaTags.Timestamp(end)+1;
    startIdx = rawData.MetaTags.Timestamp(end)+1;
    binIdx = startIdx:(startIdx+nSamples-1);
    nBins = floor((size(sqSignal,1)-startIdx+1)/nSamples);

    spikePow = zeros(nBins, size(sqSignal,2));
    timeAxis = zeros(nBins, 1);
    for n=1:nBins
        binIdx(binIdx>length(sqSignal))=[];
        spikePow(n,:) = mean(sqSignal(binIdx,:));
        timeAxis(n) = (binIdx(1)/decSR);
        binIdx = binIdx + nSamples;
    end
    spikePow = single(spikePow);
    
    % Get TX times (raster data)
    threshVectors = cell(length(threshold),1);
    for t=1:length(threshold) % issue here because threshold only has one value?
        threshVectors{t} = rmsVal*threshold(t);
    end

    timeAxis = (0:(size(allData,1)-1)) * (1/decSR);
    txEvents = cell(size(allData,2), length(threshVectors));
    for t=1:length(threshVectors)
        [txEvents(:,t)] = txEventsFromRawVoltage_parallel( allData, timeAxis, threshVectors{t} );
    end

    nBins = floor((timeAxis(end))*(1000/binSize)); % Why is this nBins ceil vs. the previous nBins is floor? Converted to floor
    binEdges = (0:nBins)*(binSize/1000);
    binTimes = (0:(nBins-1))*binSize;
    
    % Get binned TX times (PSTHs)
    binnedTX = cell(length(threshVectors),1);
    for t=1:length(threshVectors)    % only 1x1 cell array translating to length of 1
        binnedTX{t} = zeros(nBins,size(allData,2));
        for c=1:length(txEvents) 
            binnedTX{t}(:,c) = histcounts(txEvents{c,t}, binEdges); % txEvents full of empty cells leading to histcounts function failing
        end
        binnedTX{t} = uint8(binnedTX{t});
    end
    
    % Map channels to electrodes, assuming map file was loaded correctly
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
    
    electrodeSet = electrodeNumber(channelSet');
    
    neuralDataOutput.binTimes = binTimes;
    neuralDataOutput.spikePower = double(spikePow)';
    neuralDataOutput.binnedTX = double(binnedTX{1,1}')*(1/(binSize/1000));
    neuralDataOutput.TXtimes = txEvents; % not quite sure...
    neuralDataOutput.channels = channelSet';
    neuralDataOutput.electrodes = electrodeSet;
    
end