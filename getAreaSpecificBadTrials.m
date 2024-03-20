%% This function checks for bad trials for different brain areas and different stimulation type

% Please keep all condition specific protocolLists in Matlab path.

% Accepted range for V1 signal is [-500 300] in general blocks, for Cathodal stim blocks it is [-700 700]
% Accepted range for V4 signal is [-500 400] in general blocks, for Cathodal stim blocks it is [-600 500]
%In case of tACS stim blocks, all will be bad trials, we need to ignore.

%Area=V1 or V4, so AreaFlag= 1 or 2 corresponds to V1 & V4
%StimulationType={'tDCS','tACS'}
%condition={'Stim','Sham'};
%Polarity={'Cathodal','Anodal' or 'SG','FG','Alpha'};
%Session={'single','dual','dual60'};
%SessionID={0,1,2}

function [allBadTrials,badTrials]=getAreaSpecificBadTrials(MonkeyName,folderSourceString,gridType,AreaFlag,StimulationType,condition,Polarity,minLimit,maxLimit,SessionID,saveDataFlag,showElectrodes,checkPeriod,rejectTolerance)

V1=1:48;
V4=49:96;
BrainArea={V1,V4};

if ~exist('checkPeriod','var'); checkPeriod = [-0.7 0.95];end
if ~exist('rejectTolerance','var'); rejectTolerance = 0.5;end
if ~exist('showElectrodes','var'); showElectrodes = [];end
if ~exist('checkPeriod','var'); checkPeriod = [-0.7 0.95];end
if ~exist('saveDataFlag','var'); saveDataFlag = 1;end
if ~exist('threshold','var'); threshold = 6;end

if SessionID==0
    Session='single';
elseif SessionID==1
    Session='dual';
elseif SessionID==2
    Session='dual60';
end

%Getting the right protocols
if SessionID==0
protocollist=append(StimulationType,'_',Polarity,'_',condition);
elseif SessionID==1||2
    protocollist=append(StimulationType,'_',Polarity,'_',condition,'_',Session);
end

[expDates,protocolNames,~] = eval(['allProtocols' protocollist]);

dates=unique(expDates,'stable');

ProtN=size(expDates,2)./size(dates,2); %Gives us the number of Protocols, condition specific

clear protocols
for day=1:length(dates)
    for n=1:ProtN
        protocols{day,n}=protocolNames{1,n+(day-1)*ProtN};%This loop rearranges all the protocols in a nice identifiable structure
    end
end

% Find highRMSElectrodes if this data is available
subjectName=MonkeyName;
rfDataFileName = [subjectName gridType 'RFData.mat']; % This file is in DataMap/ReceptiveFieldData/{subjectName} folder and should be in Matlab's path
if exist(rfDataFileName,'file')
    tmp = load(rfDataFileName);
    electrodesToUse = tmp.highRMSElectrodes;
end

electrodeList = intersect(electrodesToUse,BrainArea{1,AreaFlag});
checkTheseElectrodes = electrodeList;
processAllElectrodes = 0;



for d=1:length(dates)
    for e=1:length(protocols)
        protocolName=protocols{d,e};
        [allBadTrials,badTrials] = findBadTrialsforLFP(MonkeyName,dates{1,d},protocolName,folderSourceString,gridType,checkTheseElectrodes,processAllElectrodes,threshold,maxLimit,minLimit,showElectrodes,saveDataFlag,checkPeriod,rejectTolerance);
    end
end
