% In this code, SF/Ori/Con Vals correspond to Parameter value.
% Bad trial is done area specific, that Id will be used to mention the Brain Area

%%Please keep the Stimulation Specific Protocol List in MATLAB path.
%StimulationType={'tDCS','tACS'}
%condition={'Stim','Sham'};
%Polarity={'Cathodal','Anodal' or 'SG','FG','Alpha'};
%SessionID=0,1 or 2 for Single, Dual or Dual with 60 minutes gap

%A is area Id, A=1 is V1, A=2 is V4
%BadDayFlag=1 takes away the analog days, 0 ignores it.



function plotTFStimSham(folderSourceString,StimulationType,Polarity,SessionID,A,BadDayFlag,SFVals,ConVals,OriVals)

if ~exist('SFVals','var'); SFVals=5; end
if ~exist('ConVals','var'); ConVals=4; end
if ~exist('OriVals','var'); OriVals=5; end

if SessionID==0
    Session='single';
elseif SessionID==1
    Session='dual';
elseif SessionID==2
    Session='dual60';
end

condition={'Stim','Sham'};

timeLimsS= [-0.5 1];
freqLimsHz=[0 100];

FigSourceString=fullfile(folderSourceString,'Figures');
if ~exist(FigSourceString, 'dir')
    mkdir(FigSourceString)
end

% for A=1:2 %%Area Loop, as this plot will consist both V1 and V4
if A==1
    badTrialNameStr='V1';
elseif A==2
    badTrialNameStr='V4';
end

if strcmp(badTrialNameStr,'V1') % We have to have a grid of different range, to comply with con values
    SGRange={{12 28}, {16,28}, {28,36}, {16,28}};
    FGRange={{32,44}, {32,48}, {48,68}, {36,52}};
    cLimsDiff=[-12 12];
elseif strcmp(badTrialNameStr,'V4')
    SGRange={{16 28}, {16,28}, {28,40}, {20,40}};
    FGRange={{32,40}, {32,44}, {60,76}, {52,72}};
    cLimsDiff=[-8 8];
end

%% Getting the right protocols
for c=1:2

    if A==1
        win=c;
    elseif A==2
        win=c+2;
    end

    if SessionID==0
        protocollist=append(StimulationType,'_',Polarity,'_',condition{1,c});
        DataFolder=fullfile(folderSourceString,'single_Stim');
    elseif SessionID==1||2
        protocollist=append(StimulationType,'_',Polarity,'_',condition{1,c},'_',Session);
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

    %% Getting Plot Handles
    hTF=getPlotHandles(2,size(protocols,2),[0.3 0.71 0.68 0.25],0.008,0.017,0);



    %% Getting Data
    DataSourceString = fullfile(DataFolder,StimulationType,badTrialNameStr,condition{1,c},Polarity);
    dates=unique(expDates,'stable');

    for iProt =1:size(protocols,2)
        for j=1:length(SFVals)
            for k=1:length(OriVals)
                for l=1:length(ConVals)
                    for day=1:length(dates)
                        DataTF = dir(fullfile(DataSourceString,append(num2str(SFVals(j)),'SF','_',num2str(OriVals(k)),'Ori','_',num2str(ConVals(l)),'Con','_',protocols{day,iProt},'_',dates{1,day},'*TF.mat')));
                        % Loading the TF file
                        d2= load(fullfile(DataSourceString,DataTF.name)); % DataStructure(field).name
                        RawTF(day,:,:) = (d2.TFDeltaPow);
                        tList=d2.tList; % Time frequency list
                        fList=d2.fList; % Frequency point list
                    end
                    allday{j,k,l}=squeeze(mean(RawTF,1));%Averaging across days took place here
                end
            end
        end
        TFData(iProt,:,:)=mean(cell2mat(permute(allday(:),[2 3 1])),3);%Averaging across all condition specific cells
    end

    %% Plot TimeFrequency plot
    for iProt=1:size(protocols,2)
        pcolor(tList,fList,squeeze(TFData(iProt,:,:))','Parent',hTF(win,iProt));colormap('jet');
        shading(hTF(win,iProt),'interp');
        clim(hTF(win,iProt),cLimsDiff);
        axes(hTF(win,iProt))
        Xax=gca().XAxis;
        Yax=gca().YAxis;

        % set(Xax,'TickDirection','out');
        % set(Xax,'TickLength',[0.02 0.025]);
        % set(Yax,'TickDirection','out');
        % set(Yax,'TickLength',[0.02 0.025]);
        set(gca,'FontWeight','bold');
        set(Xax,'FontSize',10);
        set(Yax,'FontSize',10);
        yline(cell2mat(SGRange{1,ConVals}),"--");
        yline(cell2mat(FGRange{1,ConVals}),"-");
    end



    if A==1
        axes(hTF(c,1))
    elseif A==2
        axes(hTF(c+2,1))
    end
    ylabel({condition{1,c};'Frequency(Hz)'},'FontSize',9)
end% Condition loop ends
% end% Area loop ends

%% Putting limits
axis(hTF(1:2,:),[timeLimsS freqLimsHz]);
set(hTF(1:2,2:size(protocols,2)),'YTickLabel',[]);
set(hTF(1,1:size(protocols,2)),'XTickLabel',[]);
set(hTF(1:2,1),'YTick',[0,50,100]);
set(hTF(2,1:size(protocols,2)),'XTick',[0 0.5 1]);
set(hTF(1:end,2:end),'ytick',[])
set(hTF(1,:),'xtick',[])


%% Calling the color grid and Protocol titles
[ColCode,titleString]=pickIDs(SessionID);
%% Putting titles
for iProt=1:size(protocols,2)
    t=title(hTF(1,iProt),titleString{1,iProt},'Fontsize',12);
    t.Color=ColCode{1,iProt}{1,2};
end

% %% Saving the figure in the fig folder
% folderSave = fullfile(FigSourceString);
% figname1 = fullfile(folderSave,append(Polarity,Session,'Stim','TF.tiff'));
%
% set(figure(1), 'Position', get(0, 'Screensize'));
% exportgraphics(figure(1),figname1,"Resolution",600)
% close all

