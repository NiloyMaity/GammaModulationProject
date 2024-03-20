% Modified from displayAllElectrodesGRF

% cutoffs: a 4-D array with cutoffs for firingRate, snr, total spikes and absolute change in firing rate.
% Set a particular threshold to zero if you don't want to use it. Default: [5 1.5 0 1]
% added cutoff for difference between baseline and transient firing [5 1.5 0 1 2]

% cutoffs=[1 1.2 0 0.3 2]; For analog recordings,[1 1 0 0.3 2]

% electrodesToUse - an array of electrode numbers. If left empty, the
% program will try to find (1) highRMSElectrodes from [subjectName gridType
% 'RFData.mat'] and (2) badElectrodes from badTrials mat file.

% parameterCombinationVals: Data segments for a particular stimulus
% combination can be selected by providing appropriate values as a 7-D
% array. If left empty, all stimulus repeats are used.

% stimulusPeriod: duration over which firing rate is calculated. Default:[0.25 0.75]
% baselinePeriod: duration over which firing rate is calculated. Default:[-0.5 0]

function [goodSpikeElectrodes,electrodesToUse,firingRate,snr,totalSpikes] = getGoodSpikeElectrodes(subjectName,expDate,protocolName,folderSourceString,cutoffs,badTrialNameStr,electrodesToUse,parameterCombinationVals,stimulusPeriod,baselinePeriod)

if ~exist('folderSourceString','var');  folderSourceString='N:';        end
if ~exist('cutoffs','var');             cutoffs=[];                     end
if ~exist('badTrialNameStr','var');     badTrialNameStr = '_v3';        end %% Please use 'V1' or 'V4' to do area specific analysis [Niloy]
if ~exist('electrodesToUse','var');     electrodesToUse = [];           end
if ~exist('parameterCombinationVals','var'); parameterCombinationVals = []; end
if ~exist('stimulusPeriod','var'); stimulusPeriod = [0.25 0.75];        end
if ~exist('baselinePeriod','var'); baselinePeriod = [-0.5 0];           end
if ~exist('transientPeriod','var'); transientPeriod = [0 0.1];           end

gridType='Microelectrode';

if isempty(cutoffs)
    cutoffs = [5 1.5 0 1]; % cutoffs=[1 1.2 0 0.3]; For analog recordings,[1 1 0 0.3]
end

folderName = fullfile(folderSourceString,'data',subjectName,gridType,expDate,protocolName);

% Get folders
folderExtract = fullfile(folderName,'extractedData');
folderSegment = fullfile(folderName,'segmentedData');
folderSpikes = fullfile(folderSegment,'Spikes');

%%%%%%%%%%%%%%%%% Get information about electrodes %%%%%%%%%%%%%%%%%%%%%%%%
if isempty(electrodesToUse)

    % Find highRMSElectrodes if this data is available
    rfDataFileName = [subjectName gridType 'RFData.mat']; % This file is in DataMap/ReceptiveFieldData/{subjectName} folder and should be in Matlab's path
    if exist(rfDataFileName,'file')
        tmp = load(rfDataFileName);
        SpikeElecs = tmp.highRMSElectrodes;

        if strcmp(badTrialNameStr,'V1')
            electrodesToUse= find(SpikeElecs<49);%%This part added to differentiate between V1 and V4 [Niloy]
        else
            electrodesToUse= find(electrodesToUse>=49);
        end
    else
        electrodesToUse = 1:96; % Take all electrodes
        disp(['Could not find ' rfDataFileName '. Starting with all electrodes.']);
    end
end

%%%%%%%%%%% Find bad electrodes from badTrials file if that exists %%%%%%%%
badTrialsFileName = fullfile(folderSegment,['badTrials' badTrialNameStr '.mat']);
if exist(badTrialsFileName,'file')
    tmp=load(badTrialsFileName);

    if isfield(tmp,'badElecs')
        badElecs = tmp.badElecs;
    else
        badElecs = [];
    end
    badTrials = tmp.badTrials;
else
    badElecs = [];
    badTrials = [];
end

electrodesToUse = setdiff(SpikeElecs(electrodesToUse),badElecs);

[neuralChannelsStored,SourceUnitID] = loadSpikeInfo(folderSpikes);
[electrodesToUse,~,iPos] = intersect(electrodesToUse,neuralChannelsStored);
SourceUnitID = SourceUnitID(iPos);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get Combinations
parameterCombinations = loadParameterCombinations(folderExtract);

if isempty(parameterCombinationVals)
    a = size(parameterCombinations,1);
    e = size(parameterCombinations,2);
    s = size(parameterCombinations,3);
    f = size(parameterCombinations,4);
    o = size(parameterCombinations,5);
    c = size(parameterCombinations,6);
    t = size(parameterCombinations,7);
else
    a = parameterCombinationVals(1);
    e = parameterCombinationVals(2);
    s = parameterCombinationVals(3);
    f = parameterCombinationVals(4);
    o = parameterCombinationVals(5);
    c = parameterCombinationVals(6);
    t = parameterCombinationVals(7);
end

goodPos = setdiff(parameterCombinations{a,e,s,f,o,c,t},badTrials);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N = length(electrodesToUse);

firingRate = zeros(1,N);
snr = zeros(1,N);
totalSpikes = zeros(1,N);
changeInFiringRate = zeros(1,N);

for i=1:N
    disp(i);
    clear spikeData
    x = load(fullfile(folderSpikes,['elec' num2str(SpikeElecs(i)) '_SID' num2str(SourceUnitID(i))]));
    spikeData = x.spikeData;
    H = getSpikeCounts(spikeData(goodPos),stimulusPeriod);
    HBL = getSpikeCounts(spikeData(goodPos),baselinePeriod);
    HT=getSpikeCounts(spikeData(goodPos),transientPeriod);
    
    totalSpikes(i) = sum(H);
    firingRate(i) = mean(H)/diff(stimulusPeriod);

    totalTransientSpike(i)=sum(HT);
    transientfiringRate(i) = mean(HT)/diff(transientPeriod);
    

    changeInFiringRate(i) = abs(firingRate(i) - mean(HBL)/diff(baselinePeriod));
    changeIntransientfiringRate(i) =  transientfiringRate(i)./(mean(HBL)/diff(baselinePeriod));

    clear segmentData
    x = load(fullfile(folderSegment,'Segments',['elec' num2str(electrodesToUse(i))]));
    segmentData = x.segmentData;
    snr(i) = getSNR(segmentData);
end

goodSpikeElectrodes = electrodesToUse((firingRate>=cutoffs(1)) & (snr>=cutoffs(2)) & (totalSpikes>=cutoffs(3)) & (changeInFiringRate>=cutoffs(4))& (changeIntransientfiringRate>=cutoffs(5)));
    goodSpikeElectrodesFilename=fullfile(folderSegment,append('GoodUnits',badTrialNameStr)); %%Added to save the GoodUnits with Area name [Niloy]
    save(goodSpikeElectrodesFilename,'firingRate','snr','goodSpikeElectrodes','totalSpikes','cutoffs')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load Data
function [neuralChannelsStored,SourceUnitID] = loadSpikeInfo(folderSpikes)
fileName = fullfile(folderSpikes,'spikeInfo.mat');
if exist(fileName,'file')
    x=load(fileName);
    neuralChannelsStored=x.neuralChannelsStored;
    SourceUnitID=x.SourceUnitID;
else
    neuralChannelsStored=[];
    SourceUnitID=[];
end
end
function [parameterCombinations,aValsUnique,eValsUnique,sValsUnique,fValsUnique,oValsUnique,cValsUnique,tValsUnique] = loadParameterCombinations(folderExtract)

x=load(fullfile(folderExtract,'parameterCombinations.mat'));
parameterCombinations=x.parameterCombinations;
aValsUnique=x.aValsUnique;
eValsUnique=x.eValsUnique;

if ~isfield(x,'sValsUnique')
    sValsUnique = x.rValsUnique/3;
else
    sValsUnique=x.sValsUnique;
end

fValsUnique=x.fValsUnique;
oValsUnique=x.oValsUnique;

if ~isfield(x,'cValsUnique')
    cValsUnique=[];
else
    cValsUnique=x.cValsUnique;
end

if ~isfield(x,'tValsUnique')
    tValsUnique=[];
else
    tValsUnique=x.tValsUnique;
end
end