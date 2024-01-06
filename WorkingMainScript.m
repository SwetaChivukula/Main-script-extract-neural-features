% working main code

% preprocessNeuralDataTactorSingleAuto
% Running computeNeuralBasics on Tactor Single Auto Data

%% Reset the System
close all; clear; clc;

%% Load the Data

%where the NSP neural data is saved
neuralDataDir = 'C:\Users\chivu\OneDrive\Documents\MATLAB\research_project_graczyk\NSP Data';
%where the stimulus data is saved (the .mat files)
stimulationDataDir = 'C:\Users\chivu\OneDrive\Documents\MATLAB\research_project_graczyk\NSP Data';
% where the processed data should be saved
%saveDir = 'C:\Users\bchut\Documents\Research\SensoryFeedbackAnalysis\data\2021-11-19 Tactor Single Preliminary 1\ProcessedData';
neuralDataSaveDir = 'C:\Users\chivu\OneDrive\Documents\MATLAB\research_project_graczyk\NSP Data';
%list of blocks to process (e.g., formatBlockList = [1 2 3 4 5 6 7 8 9 10])
formatBlockList = [1:3]; %1:3 for motor cortex since last file didn't save correctly
excludeElectrodes = cell(3,1);
excludeChannels = cell(3,1);



excludeElectrodes{1} = [98,128]; %update with bad electrodes, 10 and 123 included
excludeChannels{1} = [96,65];%

%% Review the Data
for blockIdx=1:length(formatBlockList)
    % Signal to the command window the status of the reformatting
    disp(['Formatting block ' num2str(formatBlockList(blockIdx))]);
    % blockNumber = formatBlockList(blockIdx);
    %for NSPIdx = 1:3
        % Find the NSP Files
        nspFileSearchString = [neuralDataDir filesep 'RP01_ICMS_2-Contact_Characterization_*_NSP' num2str(blockIdx) '*.ns5']; % note the NSP counter was offset by one
        nspFiles = dir(nspFileSearchString);
        if isempty(nspFiles)
            error(['No NSP files found for block ' num2str(blockIdx) ' that matched the search string: ' nspFileSearchString]);
        end

    for NSPFileIdx = 1:length(nspFiles)
        nspFiles = nspFiles(1:length(nspFiles)); % Take the first 3 assuming the number of blocks goes from 1:99 to avoid confusion between 1 and 11
        % Find the Stimulation File
        blockNumber = NSPFileIdx;
        stimulationFileSearchString = [stimulationDataDir filesep 'RP01_ICMS_2-Contact_Characterization*.mat'];
        stimulationFiles = dir(stimulationFileSearchString);
        % stimulationFile = stimulationFiles(blockNumber);
        if isempty(stimulationFiles)
            error(['No stimulation file found for block ' num2str(blockIdx) ' that matched the search string: ' stimulationFileSearchString]);
            % blockIdx was changed from blockNumber
        end
        
%%
% eventually make a for loop at (1)
        % Load the stimulation .mat file
        stimStruct = load([stimulationDataDir filesep stimulationFiles(NSPFileIdx).name]);
        
        %nspFilename = nspFiles(1).name; % For the sensory cortex
        nspFilename = nspFiles(NSPFileIdx).name; % for any cortex
        % Load the neural data .ns5 file
        %if(blockNumber == 6) % JUST FOR RECEPTIVE FIELD MAPPING DAY
        %NS5 = openNSx([neuralDataDir filesep nspFilename],'read','t:0:7','min'); % Converts the ns5 file to a matlab friendly file
        %else
        NS5 = openNSx([neuralDataDir filesep nspFilename],'read'); % Converts the ns5 file to a matlab friendly file
        %end
        excludeChannelsNSP = excludeChannels{blockIdx}; %blockIdx used to be NSPIdx
        
        %% Analyze neural data - Threshold Crossings (binned & spike times), Spike Power
        threshold = -4.5; % threshold crossing RMS multiplier
        binSize = 10; % milliseconds % can change to 1ms
        %filterRange = [300 5000];
        filterRange = [250 5000];

        % changing electrodeSet, chanSet, electrodeSetFull from 1:128 to
        % 1:130
        electrodeSet = 1:130;
        chanSet = [1:130];
        electrodeSetFull = [1:130];
%         if(~isempty(excludeElecs))
%            electrodeSetMask = ones(size(electrodeSetFull));
%            electrodeSetMask(excludeElecs) = 0;
%            electrodeSet = electrodeSetFull(electrodeSetMask==1);
%         else
%            electrodeSet = electrodeSetFull;
%         end

%blockIdx used to be NSPIdx - checking if this works better
        if(blockIdx == 1)
            % for the sensory cortex
            arrayName = '-Sensory ';
        elseif(blockIdx == 2)
            arrayName = '-AIP-IFG ';
        elseif(blockIdx == 3)
            arrayName = '-Motor ';
        else
            arrayName = '';
        end
        
        if(length(chanSet)<26)
        neuralDataOutputAll = computeNeuralBasics(NS5,chanSet,binSize,filterRange,threshold);
        else
            neuralDataOutputAll = struct('binTimes',[],'spikePower',[],'binnedTX',[],'TXtimes',[],'channels',[],'electrodes',[]);
            numChannelsTotal = length(chanSet);
            numDivisions = 5;
            chanSubsetStart = 1;
            division = 5;
            numChannelsSubset = numChannelsTotal/numDivisions;
            chanSubsetStart = (division-1)*numChannelsSubset+1;
            chanSubsetEnd = division*numChannelsSubset;
            chanSubset = chanSubsetStart:chanSubsetEnd;

            for i = 1:length(chanSubset)
                if chanSubset(i) == 129
                    

            for division = 1:numDivisions
                numChannelsSubset = numChannelsTotal/numDivisions;
                chanSubsetStart = (division-1)*numChannelsSubset+1;
                chanSubsetEnd = division*numChannelsSubset;
                chanSubset = chanSubsetStart:chanSubsetEnd;
                neuralDataOutputChunk = computeNeuralBasics(NS5,chanSubset,binSize,filterRange,threshold,excludeChannelsNSP);
                if(division == 1)
                    neuralDataOutputAll.binTimes = [neuralDataOutputAll.binTimes,neuralDataOutputChunk.binTimes];
                end
                neuralDataOutputAll.spikePower = [neuralDataOutputAll.spikePower; neuralDataOutputChunk.spikePower];
                neuralDataOutputAll.binnedTX = [neuralDataOutputAll.binnedTX; neuralDataOutputChunk.binnedTX];
                neuralDataOutputAll.TXtimes = [neuralDataOutputAll.TXtimes; neuralDataOutputChunk.TXtimes];
                neuralDataOutputAll.channels = [neuralDataOutputAll.channels; neuralDataOutputChunk.channels];
                neuralDataOutputAll.electrodes = [neuralDataOutputAll.electrodes; neuralDataOutputChunk.electrodes];
            end
            
        end
        neuralData = neuralDataOutputAll;
        save([neuralDataSaveDir filesep 'Single Tactor Auto Neural Basics ' num2str(blockIdx) arrayName num2str(blockNumber) '.mat'],"neuralData");
%blockIdx used to be NSPIdx
    end
        end
    end
end