% In this code, SF/Ori/Con Vals correspond to Parameter values.
% Bad trial is done area specific, that Id will be used to mention the Brain Area

%%Please keep the Stimulation Specific Protocol List in MATLAB path.
%StimulationType={'tDCS','tACS'}
%condition={'Stim','Sham'};
%Polarity={'Cathodal','Anodal' or 'SG','FG','Alpha'};
%SessionID=0,1 or 2 for Single, Dual or Dual with 60 minutes gap

%normalize flag=1 plots normalized PSTH, 0 plots general PSTH
%BadDayFlag=1 takes away the analog days, 0 ignores it.



function plotPSDandSpike(folderSourceString,StimulationType,condition,Polarity,SessionID,normalizeFlag,BadDayFlag,SFVals,ConVals,OriVals,badTrialNameStr)

if ~exist('SFVals','var'); SFVals=5; end
if ~exist('ConVals','var'); ConVals=4; end
if ~exist('OriVals','var'); OriVals=5; end
if ~exist('badTrialNameStr','var'); badTrialNameStr = 'V1'; end % This is Area specific

if SessionID==0
    Session='single';
elseif SessionID==1
    Session='dual';
elseif SessionID==2
    Session='dual60';
end

if strcmp('Sham',condition)
    c=1;
elseif strcmp('Stim',condition)
    c=2;
end


if strcmp(badTrialNameStr,'V1') % We have to have a grid of different range, to comply with con values
    PSDLimsDiff=[-4 15];
elseif strcmp(badTrialNameStr,'V4')
    PSDLimsDiff=[-8 8];
end

timeLimsS= [-0.5 1];
freqLimsHz=[0 100];

baselineS  = [-0.5 0];          % Baseline period for computing change
stimPeriod = [0.25 0.75];
binWidthMS=10;




%% Getting the right protocols

if SessionID==0
    protocollist=append(StimulationType,'_',Polarity,'_',condition);
    DataFolder=fullfile(folderSourceString,'single_Stim');
elseif SessionID==1||2
    protocollist=append(StimulationType,'_',Polarity,'_',condition,'_',Session);
    DataFolder=fullfile(folderSourceString,'dual_Stim');
end

[expDates,protocolNames,~] = eval(['allProtocols' protocollist]);

if BadDayFlag==1
    BadDays=GetBadDays;
    BadDay=intersect(BadDays,expDates);
    IdList = find(contains(expDates,BadDay));
    for Index=1:length(IdList)
        expDates{1,IdList(Index)}=[];
        protocolNames{1,IdList(Index)}=[];
    end
    expDates=expDates(~cellfun('isempty',expDates));
    protocolNames=protocolNames(~cellfun('isempty',protocolNames));
end


dates=unique(expDates,'stable');

ProtN=size(expDates,2)./size(dates,2); %Gives us the number of Protocols, condition specific

clear protocols
for day=1:length(dates)
    for n=1:ProtN
        protocols{day,n}=protocolNames{1,n+(day-1)*ProtN}; %This loop rearranges all the protocols in a nice identifiable structure
    end
end

% Calling the color grid
ColCode=pickIDs(SessionID);

FigSourceString=fullfile(folderSourceString,'Figures');
if ~exist(FigSourceString, 'dir')
    mkdir(FigSourceString)
end

%Finding Stimulation Frequency
if strcmp(Polarity,'SG')
    StimFreq=20;
elseif strcmp(Polarity,'FG')
    StimFreq=40;
elseif strcmp(Polarity,'Alpha')
    StimFreq=10;
end


%Getting Data
DataSourceString = fullfile(DataFolder,StimulationType,badTrialNameStr,condition,Polarity);
dates=unique(expDates,'stable');


%% For Plotting PSD & Spikes
clear RawPSD RawPSTH

%     Collecting the data
for iProt =1:size(protocols,2)
    for j=1:length(SFVals)
        for k=1:length(OriVals)
            for l=1:length(ConVals)
                for day=1:length(dates)

                    DataPSD = dir(fullfile(DataSourceString,append(num2str(SFVals(j)),'SF','_',num2str(OriVals(k)),'Ori','_',num2str(ConVals(l)),'Con','_',protocols{day,iProt},'_',dates{1,day},'*PSD.mat')));
                    DataSpike = dir(fullfile(DataSourceString,append(num2str(SFVals(j)),'SF','_',num2str(OriVals(k)),'Ori','_',num2str(ConVals(l)),'Con','_',protocols{day,iProt},'_',dates{1,day},'*PSTH.mat')));

                    % Loading the spike file
                    if ~isempty(DataSpike)
                        s1 = load(fullfile(DataSourceString,DataSpike.name)); % DataStructure(field).name
                        RawPSTH{day,:,:} = s1.PSTHGrid;
                        SPtList=s1.xs; %Spike timevals
                    else
                        RawPSTH{day,:,:} =[];
                    end

                    % Loading the PSD file
                    d1=load(fullfile(DataSourceString,DataPSD.name)); % DataStructure(field).name
                    RawPSD{day,:,:} = (d1.ShadeData);
                    f1=d1.f1;
                end
                alldayPSTH{j,k,l}=cat(1,RawPSTH{:,1}); %Concatenating across days
                alldayPSD{j,k,l}=cat(1,RawPSD{:,1});
            end
        end
    end

    PSTHData{iProt,:}=mean(cell2mat(permute(alldayPSTH(:),[2 3 1])),3);%Averaging across all condition specific cells
    PSDData{iProt,:}=mean(cell2mat(permute(alldayPSD(:),[2 3 1])),3);%Averaging across all condition specific cells
end

%% Here pulling of all data from required dates has been completed.

for i=1:size(protocols,2)
    for n=1:size(PSTHData{1,1},1)
        maxFR=max(PSTHData{1,1}(n,:)); % takes the max firing rate of the first protocol, per electrode
        normFRData{i,1}(n,:)=PSTHData{i,1}(n,:)./maxFR;
    end
end

% for i=1:size(protocols,2)
%     for n=1:size(PSTHData{1,1},1)
%         maxFR=max(PSTHData{i,1}(n,:));
%         normFRData{i,1}(n,:)=PSTHData{i,1}(n,:)./maxFR;
%     end
% end

%Getting Plot handle
% hPSD=getPlotHandles(1,1,[0.05 0.065 0.4 0.85],0.008,0.02,1);
% hSpikes=getPlotHandles(1,1,[0.57 0.065 0.4 0.85],0.008,0.02,1);

hPSD=getPlotHandles(1,2,[0.09 0.05 0.66 0.3],0.32,0.02,0);
hSpikes=getPlotHandles(1,2,[0.32 0.05 0.66 0.3],0.32,0.02,0);


%% Plot the PSD shade
hold (hPSD(1,c),'on');
for iProt =1:size(protocols,2)
    DeltaShadeData=PSDData{iProt,:};
    options.handle = hPSD(1,c);
    options.color_area = ColCode{1,iProt}{1,1};
    options.alpha      = 0.2;
    options.error      = 'sem';
    options.x_axis= f1;
    plot_areaerrorbar(DeltaShadeData, options);hold on;
end

%% Plot delta PSD plot
hold (hPSD(1,c),'on');
for iProt=1:size(protocols,2)
    h(1)=plot(hPSD(1,c),f1,mean(PSDData{iProt,:},1),LineWidth=1.7,color=ColCode{1,iProt}{1,2});hold (hPSD(1,1),'on');
end

axes(hPSD(1,c))
title(append('\Delta PSD','(',condition,')'),'FontWeight','bold')
yl=yline(0,"--");

if strcmp(StimulationType,'tACS')
    xl=xline(StimFreq,"--",'StimFreq','Color','red',LabelVerticalAlignment='bottom');
end

ylabel('\DeltaPower (dB)','FontSize',12)
text(35,13.4,append('n=',num2str(size(DeltaShadeData,1))),'color','k','FontSize',10,'FontWeight','bold');

axis(hPSD(1,c),[freqLimsHz PSDLimsDiff]);

Xax=gca().XAxis;
Yax=gca().YAxis;
set(Xax,'TickDirection','out');
set(Xax,'TickLength',[0.02 0.025]);
set(Yax,'TickDirection','out');
set(Yax,'TickLength',[0.02 0.025]);
set(gca,'FontWeight','bold');
set(Xax,'FontSize',10);
set(Yax,'FontSize',10);

line([PSDLimsDiff(:,1) freqLimsHz(:,2)], [PSDLimsDiff(:,2) PSDLimsDiff(:,2)],'color','k')
line([freqLimsHz(:,2) freqLimsHz(:,2)], [PSDLimsDiff(:,1) PSDLimsDiff(:,2)],'color','k')

if SessionID==0
    legend(hPSD(1,c),{(''),(''),(''),(''),(''),(''),('Pre-'),append(condition),('Post-'),('Check1'),('Check2'),('Check3')},'Fontsize',6.5);
elseif SessionID==1
    legend(hPSD(1,c),{(''),(''),(''),(''),(''),(''),(''),(''),('Pre-'),append(condition),('Post-'),append(condition),('Post-'),('Check1'),('Check2'),('Check3')},'Fontsize',6.5);
elseif SessionID==2
    legend(hPSD(1,c),{(''),(''),(''),(''),(''),(''),(''),(''),(''),('Pre-'),append(condition),('Post-'),('Post-'),append(condition),('Post-'),('Check1'),('Check2'),('Check3')},'Fontsize',6.5);
end

clear DeltaShadeData
% Plot delta shade PSTH plot
hold (hSpikes(1,c),'on');
for iProt =1:size(protocols,2)
    if normalizeFlag==1
        DeltaShadeData=normFRData{iProt,:};
    else
        DeltaShadeData=PSTHData{iProt,:};
    end
    options.handle = hSpikes(1,c);
    options.color_area = ColCode{1,iProt}{1,1};
    options.alpha      = 0.2;
    options.error      = 'sem';
    options.x_axis= SPtList;
    plot_areaerrorbar(DeltaShadeData, options);hold on;
end

hold (hSpikes(1,c),'on');
for iProt=1:size(protocols,2)
    if ~isempty(PSTHData{iProt,:})
        if normalizeFlag==1
            plot(hSpikes(1,c),SPtList,mean(normFRData{iProt,:},1),LineWidth=1.6, color=ColCode{1,iProt}{1,2});hold (hSpikes(1,1),'on');
        else
            plot(hSpikes(1,c),SPtList,mean(PSTHData{iProt,:},1),LineWidth=1.6, color=ColCode{1,iProt}{1,2});hold (hSpikes(1,1),'on');
        end
    end
end
hold (hSpikes(1,c),'on');

if strcmp(badTrialNameStr,'V1')
    axis(hSpikes(1,c),[timeLimsS [0 1]]);
elseif strcmp(badTrialNameStr,'V4')
    axis(hSpikes(1,c),[timeLimsS [0 8]]);
end

axes(hSpikes(1,c))

Xax=gca().XAxis;
Yax=gca().YAxis;
set(Xax,'TickDirection','out');
set(Xax,'TickLength',[0.02 0.025]);
set(Yax,'TickDirection','out');
set(Yax,'TickLength',[0.02 0.025]);
set(gca,'FontWeight','bold');
set(Xax,'FontSize',10);
set(Yax,'FontSize',10);

line([timeLimsS(:,1) 1], [1 1],'color','k')
line([1 1], [0 1],'color','k')

title(append('PSTH','(',condition,')'),'FontWeight','bold')
text(0.5,0.895,append('n=',num2str(size(DeltaShadeData,1))),'color','k','FontSize',10,'FontWeight','bold');
ylabel('normalized Firing Rate','FontSize',10)

% %% Saving the figure in the fig folder
% folderSave = fullfile(FigSourceString);
% figname1 = fullfile(folderSave,append(Polarity,Session,condition,badTrialNameStr,'PSDSpike1.tiff'));

% set(figure(1), 'Position', get(0, 'Screensize'));
% exportgraphics(figure(1),figname1,"Resolution",600)
% close all
% clear all
