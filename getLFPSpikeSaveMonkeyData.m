%% This function saves PSD,Band Powers,TF and PSTH named by SF,Con, Ori and polarity of stimulation

%%Please keep the Stimulation Specific Protocol List in MATLAB path.

%Area=V1 or V4, so AreaFlag= 1 or 2 corresponds to V1 & V4
%StimulationType={'tDCS','tACS'}
%condition={'Stim','Sham'};
%Polarity={'Cathodal','Anodal' or 'SG','FG','Alpha'};
%Session={'single','dual','dual60'};
%SessionID={0,1,2}

function [BandData,ShadeData,PSTHGrid]=get_LFP_Spike_Save_MonkeyData(MonkeyName,folderSourceString,gridType,AreaFlag,StimulationType,condition,Polarity,SessionID,PSDTFFlag,SFVals,ConVals,OriVals,badTrialNameStr)

if ~exist('SFVals','var'); SFVals=[1 2 3 4 5]; end
if ~exist('ConVals','var'); ConVals=[1 2 3 4]; end
if ~exist('OriVals','var'); OriVals=[1 2 3 4 5]; end
if ~exist('badTrialNameStr','var'); badTrialNameStr = 'V1'; end % This is Area specific

V1=1:48;
V4=49:96;
BrainArea={V1,V4};

if strcmp(badTrialNameStr,'V1') % We have to have a grid of different range, to comply with con values

    SGRange={{12 28}, {16,28}, {28,36}, {16,28}};
    FGRange={{32,44}, {32,48}, {48,68}, {36,52}};
elseif strcmp(badTrialNameStr,'V4')
    SGRange={{16 28}, {16,28}, {28,40}, {20,40}};
    FGRange={{32,40}, {32,44}, {60,76}, {52,72}};
end

if SessionID==0
    Session='single';
elseif SessionID==1
    Session='dual';
elseif SessionID==2
    Session='dual60';
end



%% Getting the right protocols
protocollist=append(StimulationType,'_',Polarity,'_',condition);

%Getting the right folders to save data

if SessionID==0
    protocollist=append(StimulationType,'_',Polarity,'_',condition);
    SaveFolder=fullfile(folderSourceString,'single_Stim');
    if ~exist(SaveFolder, 'dir')
        mkdir(SaveFolder)
    end
elseif SessionID==1||2
    protocollist=append(StimulationType,'_',Polarity,'_',condition,'_',Session);
    SaveFolder=fullfile(folderSourceString,'dual_Stim');
    if ~exist(SaveFolder, 'dir')
        mkdir(SaveFolder)
    end
end

[expDates,protocolNames,~] = eval(['allProtocols' protocollist]);

dates=unique(expDates,'stable');

ProtN=size(expDates,2)./size(dates,2); %Gives us the number of Protocols, condition specific

clear protocols
for day=1:length(dates)
    for n=1:ProtN
        protocols{day,n}=protocolNames{1,n+(day-1)*ProtN}; %This loop rearranges all the protocols in a nice identifiable structure
    end
end

% Find highRMSElectrodes if this data is available
rfDataFileName = [MonkeyName gridType 'RFData.mat']; % This file is in DataMap/ReceptiveFieldData/{subjectName} folder and should be in Matlab's path
if exist(rfDataFileName,'file')
    tmp = load(rfDataFileName);
    electrodesToUse = tmp.highRMSElectrodes;
end

electrodeList = intersect(electrodesToUse,BrainArea{1,AreaFlag});

%% Analysis and saving Loop begins
for SF=SFVals(1):SFVals(end)
    for Ori=OriVals(1):OriVals(end)
        for Con=ConVals(1):ConVals(end) %Changes the gamma range

            SGIdx= [(SGRange{1,Con}{1,1}./2)+1 (SGRange{1,Con}{1,2}./2)+1]; %Marking the gamma range according the contrast value
            FGIdx= [(FGRange{1,Con}{1,1}./2)+1 (FGRange{1,Con}{1,2}./2)+1];

            DataSourceString = fullfile(SaveFolder,StimulationType,badTrialNameStr,condition,Polarity);

            if ~exist(DataSourceString, 'dir')
                mkdir(DataSourceString)
            end

            %Finding Good units loop
            for day=1:length(dates)
                respFile = fullfile(folderSourceString,'data',MonkeyName,gridType,dates{1,day},'GRF_001','segmentedData', append('GoodUnits',badTrialNameStr,'.mat'));

                if isfile(respFile)
                    respFileStruct=load(respFile);
                    AllGoodUnit{1,1}=respFileStruct.goodSpikeElectrodes;
                    AreaGoodUnit=intersect(BrainArea{1,AreaFlag},AllGoodUnit{1,1});
                else
                    AreaGoodUnit=[];
                    disp 'No unit found';
                end

                for iProt =1:size(protocols,2)
                    folderName = fullfile(folderSourceString,'data',MonkeyName,gridType,dates{1,day},protocols{day,iProt});
                    folderExtract = fullfile(folderName,'extractedData');
                    folderSegData = fullfile(folderName,'segmentedData');
                    folderLFP = fullfile(folderSegData,'LFP');
                    folderSpikes = fullfile(folderSegData,'Spikes');

                    %%%% load LFP Information
                    [~,timeVals,~,~] = loadlfpInfo(folderLFP);

                    %%%% load Impedance Information
                    badImpedanceCutoff = 2500;
                    impedanceFileName = fullfile(folderSourceString,'data',MonkeyName,gridType,dates{1,day},'impedanceValues.mat');
                    impedanceData = load(impedanceFileName);
                    q = find(impedanceData.impedanceValues>badImpedanceCutoff) ;
                    badImpedanceElecs= intersect(q,BrainArea{1,AreaFlag});

                    %%%% load goodStims file
                    load(fullfile(folderExtract,'goodStimNums.mat'));

                    %%%% load badTrials file
                    load(fullfile(folderSegData,append('badTrials',badTrialNameStr,'.mat')));
                    electrodeList=setdiff(electrodeList,badImpedanceElecs);

                    %%%%Get Combinations
                    [parameterCombinations,aValsUnique,eValsUnique,sValsUnique,fValsUnique,oValsUnique,cValsUnique,tValsUnique] = loadParameterCombinations(folderExtract); %#ok<*ASGLU>
                    goodPos = parameterCombinations{1,1,1,SF,Ori,Con,1};   % Finding parameter combination for gamma; This corresponds to average of all SF and all Ori

                    %%%Steps to check if bad trial file is from tACS stim block or not
                    tf=isequal(length(goodStimNums),length(badTrials));
                    tf=double(tf);

                    if tf==1
                        badTrials=[];
                    end

                    goodPos=setdiff(goodPos,badTrials);

                    %%%% Parameters for MTM
                    % window lengths for MTM; you can vary this to increase frequency resolution in TF plots; window length of 0.2 gives you 5 Hz; 0.25 gives you 4 Hz and accordingly!
                    BWList     = 1;                 % Bandwidths. BW=1 for STFT BW=2 for MTM

                    baselineS  = [-0.5 0];          % Baseline period for computing change
                    stimPeriod = [0.25 0.75];

                    Ts= timeVals(2)-timeVals(1);    % Sampling period
                    Fs = round(1/Ts);               % Sampling Frequency
                    params.pad = -1;                %No padding=-1
                    params.Fs = Fs;
                    params.fpass = [0 100];
                    params.trialave = 1;            %Takes average of trial data if the value is 1, don't if 0
                    params.tapers = [BWList 2*BWList-1];
                    movingWinTF  = [0.250 0.025];

                    range = baselineS;
                    rangePos = round(diff(range)*Fs);
                    blPos = find(timeVals>=baselineS(1),1)+(1:rangePos);
                    stPos = find(round(timeVals,4)>=stimPeriod(1),1)+(1:rangePos);

                    %%%%Spike data saving loop%%%%

                    for chan=1:length(AreaGoodUnit)
                        clear spikeData
                        x=load(fullfile(folderSpikes,['elec' num2str(AreaGoodUnit(1,chan)) '_SID0.mat']));
                        spikeData=x.spikeData;
                        [psthVals,xs] = getPSTH(spikeData(goodPos),10,[timeVals(1) timeVals(end)]);
                        PSTHGrid(chan,:)=(psthVals);
                    end

                    PSTHFileName=fullfile(DataSourceString,append(num2str(SF),'SF','_',num2str(Ori),'Ori','_',num2str(Con),'Con','_',protocols{day,iProt},'_',dates{1,day},'_','PSTH.mat'));
                    save(PSTHFileName,"PSTHGrid","xs");
                    PSTHGrid=[];

                    %%%%loading all electrode data
                 if PSDTFFlag==1
                    for iElec = 1:length(electrodeList)
                        clear x
                        x= load(fullfile(folderLFP,['elec' num2str(electrodeList(1,iElec)) '.mat']));
                        ElecData=x.analogData(goodPos,:);
                        ElecData=ElecData-repmat(mean(ElecData,2),1,size(ElecData,2)); %DC shift correction
                        PSDdataBL = ElecData(:,blPos);
                        PSDdataST = ElecData(:,stPos);

                        %%%%Performing MT for TF
                        [TFSpec,tList,fList] = mtspecgramc(ElecData',movingWinTF,params);
                        tList = tList+ timeVals(1)-1/Fs; % Center the times with respect to the stimulus onset time
                        BlPosList = intersect(find(tList>=baselineS(1)),find(tList<baselineS(2))); %Finds Baseline position

                        SpectraTFBl(iElec,:,:) = mean(TFSpec(BlPosList,:),1);
                        SpectraTFSt(iElec,:,:) = TFSpec; % Stores each electrode's stimulus power spectra


                        %%%%Performing MT for PSD
                        [S1,f1]=mtspectrumc(PSDdataBL',params); %#ok<*SAGROW> %We are averaging across all trials here inside mtspectrumc.
                        [S2,f2]=mtspectrumc(PSDdataST',params); %#ok<*SAGROW>


                        SpectraPSDBl(iElec,:) = S1;% Stores each electrode's baseline power spectra;(total elec x freqpoints, ex: 48x51)
                        SpectraPSDSt(iElec,:) = S2;% Stores each electrode's stimulus power spectra
                    end

                    ShadeData=(10*(log10(SpectraPSDSt))-10*(log10(SpectraPSDBl)));

                    TFBlPower = (squeeze(mean(log10(SpectraTFBl),1))); % Across electrode averaging
                    TFStPower = (squeeze(mean(log10(SpectraTFSt),1)));


                    CollectBlSG(iProt,:,:)=log10(mean(SpectraPSDBl(:,SGIdx(1):SGIdx(2)),2));
                    CollectBlFG(iProt,:,:)=log10(mean(SpectraPSDBl(:,FGIdx(1):FGIdx(2)),2));

                    CollectStSG(iProt,:,:)=log10(mean(SpectraPSDSt(:,SGIdx(1):SGIdx(2)),2));
                    CollectStFG(iProt,:,:)=log10(mean(SpectraPSDSt(:,FGIdx(1):FGIdx(2)),2));

                    TFDeltaPow= 10*(TFStPower-repmat(TFBlPower',length(tList),1));

                    %%%Saving TF file

                    PSDFileName=fullfile(DataSourceString,append(num2str(SF),'SF','_',num2str(Ori),'Ori','_',num2str(Con),'Con','_',protocols{day,iProt},'_',dates{1,day},'_','PSD.mat'));
                    TFFileName=fullfile(DataSourceString,append(num2str(SF),'SF','_',num2str(Ori),'Ori','_',num2str(Con),'Con','_',protocols{day,iProt},'_',dates{1,day},'_','TF.mat'));
                    save(PSDFileName,"ShadeData","f1");
                    save(TFFileName,"TFDeltaPow","TFStPower","fList","tList");
                 end
                end
                
               if PSDTFFlag==1
                CollectSG=10*(CollectStSG-CollectBlSG);
                AvgSG=mean(CollectSG,2);
                CollectFG=10*(CollectStFG-CollectBlFG);
                AvgFG=mean(CollectFG,2);
                BandData=[AvgSG AvgFG];
               
                BandDataFileName=fullfile(DataSourceString,append(num2str(SF),'SF','_',num2str(Ori),'Ori','_',num2str(Con),'Con','_',dates{1,day},'_','BandData.mat'));
                save(BandDataFileName,'BandData','CollectStSG','CollectBlSG','CollectStFG','CollectBlFG','CollectSG','CollectFG');
               end
            end %Date Loop
            clear SpectraPSDSt SpectraPSDBl SpectraTFBl SpectraTFSt CollectStSG CollectBlSG CollectBlFG CollectStFG
        end %Con loop
    end %Ori Loop
end %SF loop
end

%% Accessory Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% load Data
function [analogChannelsStored,timeVals,goodStimPos,analogInputNums] = loadlfpInfo(folderLFP)
load(fullfile(folderLFP,'lfpInfo.mat'));
analogChannelsStored=sort(analogChannelsStored);
if ~exist('analogInputNums','var')
    analogInputNums=[];
end
end

function [neuralChannelsStored,SourceUnitID] = loadspikeInfo(folderSpikes)
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

%Loads Parameter combination
function [parameterCombinations,aValsUnique,eValsUnique,sValsUnique,...
    fValsUnique,oValsUnique,cValsUnique,tValsUnique] = loadParameterCombinations(folderExtract)

load(fullfile(folderExtract,'parameterCombinations.mat')); %#ok<*LOAD>
end

