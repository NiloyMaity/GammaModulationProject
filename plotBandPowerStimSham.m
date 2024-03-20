% In this code, SF/Ori/Con Vals correspond to Parameter values.
% Bad trial is done area specific, that Id will be used to mention the Brain Area

%%Please keep the Stimulation Specific Protocol List in MATLAB path.
%StimulationType={'tDCS','tACS'}
%condition={'Stim','Sham'};
%Polarity={'Cathodal','Anodal' or 'SG','FG','Alpha'};
%SessionID=0,1 or 2 for Single, Dual or Dual with 60 minutes gap

%BadDayFlag=1 takes away the analog days, 0 ignores it.


function plotBandPowerStimSham(folderSourceString,StimulationType,Polarity,SessionID,BadDayFlag,SFVals,ConVals,OriVals,badTrialNameStr)

if ~exist('SFVals','var'); SFVals=5; end
if ~exist('ConVals','var'); ConVals=4; end
if ~exist('OriVals','var'); OriVals=5; end
if ~exist('badTrialNameStr','var'); badTrialNameStr = 'V1'; end % This is Area specific

if SessionID==0
    Session='Single';
elseif SessionID==1
    Session='dual';
elseif SessionID==2
    Session='dual60';
end

condition={'Stim','Sham'};

baselineS  = [-0.5 0];          % Baseline period for computing change
stimPeriod = [0.25 0.75];
binWidthMS=10;


% hDelta1=getPlotHandles(1,2,[0.1 0.6 0.3 0.3],0.03,0.02,0);
% set(hDelta1(1,2),'ytick',[])
%
% % hDelta2=getPlotHandles(1,2,[0.46 0.6 0.3 0.3],0.03,0.02,0);
% % % set(hDelta2(1,2),'ytick',[])
%
% hDelta3=getPlotHandles(1,2,[0.1 0.1 0.3 0.3],0.03,0.02,0);
% set(hDelta3(1,2),'ytick',[])
%
% hDelta4=getPlotHandles(1,2,[0.46 0.1 0.3 0.3],0.03,0.02,0);
% set(hDelta4(1,2),'ytick',[])

hDelta1=getPlotHandles(1,1,[0.09 0.71 0.14 0.2],0.015,0.02,0);

hDelta2=getPlotHandles(1,2,[0.09 0.42 0.28 0.2],0.02,0.02,0);
set(hDelta2(1,2),'ytick',[])

hDelta3=getPlotHandles(1,2,[0.4 0.42 0.28 0.2],0.015,0.02,0);
set(hDelta3(1,2),'ytick',[])

hDelta4=getPlotHandles(1,2,[0.7 0.42 0.28 0.2],0.015,0.02,0);
set(hDelta4(1,2),'ytick',[])



% Calling the color grid
[ColCode,titleString,StimblockID]=pickIDs(SessionID);

for ses=1:length(StimblockID)
    ColCode{1,StimblockID(ses)}=[];
    titleString{1,StimblockID(ses)}=[];
end

Col=ColCode(~cellfun('isempty',ColCode));
ColCode=reshape(Col,[1 length(ColCode)-length(StimblockID)]);

%% Getting the right protocols
for c=1:2
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

    for q=1:length(dates)
        for r=1:length(StimblockID)
            protocols{q,StimblockID(r)}=[];
        end
    end

    y=protocols(~cellfun('isempty',protocols));
    protocols=reshape(y,[length(dates) length(protocols)-length(StimblockID)]);

    FigSourceString=fullfile(folderSourceString,'Figures');
    if ~exist(FigSourceString, 'dir')
        mkdir(FigSourceString)
    end

    % Getting Data
    DataSourceString = fullfile(DataFolder,StimulationType,badTrialNameStr,condition{1,c},Polarity);
    dates=unique(expDates,'stable');

    % Collecting Data
    for day=1:length(dates)
        for j=1:length(SFVals)
            for k=1:length(OriVals)
                for l=1:length(ConVals)
                    DataBand = dir(fullfile(DataSourceString,append(num2str(SFVals(j)),'SF','_',num2str(OriVals(k)),'Ori','_',num2str(ConVals(l)),'Con','_',dates{1,day},'*BandData.mat')));

                    clear BandValues
                    BandValues = load(fullfile(DataSourceString,DataBand.name));

                    cSG=BandValues.CollectSG;
                    cFG=BandValues.CollectFG;
                    cBlSG=BandValues.CollectBlSG;
                    cStSG=BandValues.CollectStSG;
                    cBlFG=BandValues.CollectBlFG;
                    cStFG=BandValues.CollectStFG;

                    for r=1:length(StimblockID)
                        cSG(StimblockID(r),:)=0;
                        cFG(StimblockID(r),:)=0;
                        cBlSG(StimblockID(r),:)=0;
                        cStSG(StimblockID(r),:)=0;
                        cBlFG(StimblockID(r),:)=0;
                        cStFG(StimblockID(r),:)=0;
                    end

                    cSG=reshape(nonzeros(cSG),[size(protocols,2) size(cSG,2)]);
                    cFG=reshape(nonzeros(cFG),[size(protocols,2) size(cFG,2)]);
                    cBlSG=reshape(nonzeros(cBlSG),[size(protocols,2) size(cBlSG,2)]);
                    cStSG=reshape(nonzeros(cStSG),[size(protocols,2) size(cStSG,2)]);
                    cBlFG=reshape(nonzeros(cBlFG),[size(protocols,2) size(cBlFG,2)]);
                    cStFG=reshape(nonzeros(cStFG),[size(protocols,2) size(cStFG,2)]);

                    for iProt =1:size(protocols,2)
                        DataSpike = dir(fullfile(DataSourceString,append(num2str(SFVals(j)),'SF','_',num2str(OriVals(k)),'Ori','_',num2str(ConVals(l)),'Con','_',protocols{day,iProt},'_',dates{1,day},'*PSTH.mat')));
                        % Loading the spike file
                        if ~isempty(DataSpike)
                            s1 = load(fullfile(DataSourceString,DataSpike.name)); % DataStructure(field).name
                            RawPSTH{iProt,j,k,l} = s1.PSTHGrid;
                            SPtList=s1.xs; %Spike timevals
                        else
                            RawPSTH{iProt,j,k,l} =[];
                        end

                        blPos = find(SPtList>=baselineS(1),1)+ (1:(diff(baselineS))/(binWidthMS/1000));
                        stPos = find(SPtList>=stimPeriod(1),1)+ (1:(diff(stimPeriod))/(binWidthMS/1000));

                        RawFRBl{iProt,j,k,l}= mean(s1.PSTHGrid(:,blPos),2);
                        RawFRSt{iProt,j,k,l}= mean(s1.PSTHGrid(:,stPos),2);


                        SGAcross{iProt,j,k,l}=cSG(iProt,:)';
                        FGAcross{iProt,j,k,l}=cFG(iProt,:)';

                        SGAcrossBl{iProt,j,k,l}=cBlSG(iProt,:)';
                        SGAcrossSt{iProt,j,k,l}=cStSG(iProt,:)';
                        FGAcrossBl{iProt,j,k,l}=cBlFG(iProt,:)';
                        FGAcrossSt{iProt,j,k,l}=cStFG(iProt,:)';
                    end
                end
            end
        end
        for iProt =1:size(protocols,2)
            RawFRAcrossBl{day,iProt,:}=mean(cell2mat(RawFRBl(iProt,:)),2);
            RawFRAcrossSt{day,iProt,:}=mean(cell2mat(RawFRSt(iProt,:)),2);

            BandSGAcross{day,iProt,:}=mean(cell2mat(SGAcross(iProt,:)),2); %Averaging across conditions
            BandFGAcross{day,iProt,:}=mean(cell2mat(FGAcross(iProt,:)),2);
            BandSGAcrossBl{day,iProt,:}=mean(cell2mat(SGAcrossBl(iProt,:)),2);
            BandFGAcrossBl{day,iProt,:}=mean(cell2mat(FGAcrossBl(iProt,:)),2);
            BandSGAcrossSt{day,iProt,:}=mean(cell2mat(SGAcrossSt(iProt,:)),2);
            BandFGAcrossSt{day,iProt,:}=mean(cell2mat(FGAcrossSt(iProt,:)),2);
        end
    end

    %% Here pulling of all data from required dates has been completed.

    clear SGAcrossAll FGAcrossAll SGAcrossBlAll SGAcrossStAll FGAcrossBlAll FGAcrossStAll FRAcrossBlAll FRAcrossStAll
    for b=1:size(protocols,2)
        FRAcrossBlAll{b,:}=cat(1,RawFRAcrossBl{:,b});
        FRAcrossStAll{b,:}=cat(1,RawFRAcrossSt{:,b});

        SGAcrossAll{b,:}=cat(1,BandSGAcross{:,b});
        FGAcrossAll{b,:}=cat(1,BandFGAcross{:,b});

        SGAcrossBlAll{b,:}=cat(1,BandSGAcrossBl{:,b});
        FGAcrossBlAll{b,:}=cat(1,BandFGAcrossBl{:,b});
        SGAcrossStAll{b,:}=cat(1,BandSGAcrossSt{:,b});
        FGAcrossStAll{b,:}=cat(1,BandFGAcrossSt{:,b});
    end

    clear SGAcross FGAcross SGAcrossBl FGAcrossBl SGAcrossSt FGAcrossSt RawFRAcrossBl RawFRAcrossSt

    % And This brings all electrode in one page
    %% Lets put all electrodes' band specific data in a cell, so that we can later do subtraction between conditions (Stim-Sham)
    CompareCellFRBl{c,:}=FRAcrossBlAll;
    CompareCellFRSt{c,:}=FRAcrossStAll;

    CompareCellSG{c,:}=SGAcrossAll;
    CompareCellFG{c,:}=FGAcrossAll;

    SeparateCellSGBl{c,:}=SGAcrossBlAll;
    SeparateCellSGSt{c,:}=SGAcrossStAll;

    SeparateCellFGBl{c,:}=FGAcrossBlAll;
    SeparateCellFGSt{c,:}=FGAcrossStAll;
end

%% Collecting Stim and Sham data separately for later usage
StimDataFRBl(:,1)=CompareCellFRBl{1,1}; ShamDataFRBl(:,1)=CompareCellFRBl{2,1};
StimDataFRSt(:,1)=CompareCellFRSt{1,1}; ShamDataFRSt(:,1)=CompareCellFRSt{2,1};


StimData(:,1)=CompareCellSG{1,1}; ShamData(:,1)=CompareCellSG{2,1};
StimDataSt(:,1)=SeparateCellSGSt{1,1}; ShamDataSt(:,1)=SeparateCellSGSt{2,1};
StimDataBl(:,1)=SeparateCellSGBl{1,1}; ShamDataBl(:,1)=SeparateCellSGBl{2,1};

StimData(:,2)=CompareCellFG{1,1}; ShamData(:,2)=CompareCellFG{2,1};
StimDataSt(:,2)=SeparateCellFGSt{1,1}; ShamDataSt(:,2)=SeparateCellFGSt{2,1};
StimDataBl(:,2)=SeparateCellFGBl{1,1}; ShamDataBl(:,2)=SeparateCellFGBl{2,1};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for Protnum=1:size(protocols,2)
    AvgFRStim{Protnum}=mean(StimDataFRSt{Protnum},1);
    AvgFRSham{Protnum}=mean(ShamDataFRSt{Protnum},1);
    errorgridFRStim{Protnum}=std((StimDataFRSt{Protnum}))/sqrt(size(StimDataFRSt{Protnum},1));
    errorgridFRSham{Protnum}=std((ShamDataFRSt{Protnum}))/sqrt(size(ShamDataFRSt{Protnum},1));

    AvgdeltaFRStim{Protnum}= mean((StimDataFRSt{Protnum}-StimDataFRBl{Protnum}),1);
    AvgdeltaFRSham{Protnum}=mean((ShamDataFRSt{Protnum}-ShamDataFRBl{Protnum}),1);
    errordeltagridFRStim{Protnum}=std((StimDataFRSt{Protnum}-StimDataFRBl{Protnum}))/sqrt(size(StimDataFRSt{Protnum},1));
    errordeltagridFRSham{Protnum}=std((ShamDataFRSt{Protnum}-ShamDataFRBl{Protnum}))/sqrt(size(ShamDataFRSt{Protnum},1));
end

%plotting the first FR plot
hold on
point=1:size(protocols,2);

hold(hDelta1(1,1),'on')
e1=errorbar(subplot(hDelta1(1,1)),point,[AvgFRSham{:}],[errorgridFRSham{:}],'color','#ED4672','LineWidth',1.5);
for f=1:size(protocols,2)
    w=[AvgFRSham{:}];
    scatter(point(1,f),w(1,f),30,ColCode{1, f}{1, 1},"filled")
end

e2=errorbar(subplot(hDelta1(1,1)),point,[AvgFRStim{:}],[errorgridFRStim{:}],'color',"#CD9A18",'LineWidth',1.5);

for f=1:size(protocols,2)
    w=[AvgFRStim{:}];
    scatter(point(1,f),w(1,f),30,ColCode{1, f}{1, 1},"filled")
end


e1.Marker='.';
e2.Marker='.';
e1.CapSize=3;
e2.CapSize=3;

yline(0,'--')


text(5,8,append('n=',num2str(size(StimDataFRSt{1,1},1))),'color','k','FontWeight','bold');

% title(hDelta1(1,1),'Firing rate','color','#1768ff','Fontsize',10,'FontWeight','bold')
% subtitle(hDelta1(1,1),'Sham(Light) vs Stim(Dark)','Fontsize',10,'FontWeight','bold')


title(hDelta1(1,1),'Firing rate','color','k','Fontsize',10,'FontWeight','bold')

subtitle(hDelta1(1,1),{'\color[rgb]{0.9294, 0.2784, 0.4471}Sham \color[rgb]{0, 0, 0}vs \color[rgb]{0.8039, 0.6039, 0.0941}Stim'},'Fontsize',10,'FontWeight','bold');

if strcmp(badTrialNameStr,'V1')
    axis(hDelta1(1,1),[[0 length(point)+1] [0 10]]);

    % if SessionID==1||2
    % axis(hDelta1(1,1),[[0 length(point)+1] [0 50]]);
    % end
    set(hDelta1(1,1),'XTick',1:length(point));
    Xax=gca().XAxis;
    Yax=gca().YAxis;
    set(Xax,'TickDirection','out');
    set(Xax,'TickLength',[0.02 0.025]);
    set(Yax,'TickDirection','out');
    set(Yax,'TickLength',[0.02 0.025]);
    set(gca,'FontWeight','bold');
    set(Xax,'FontSize',10);
    set(Yax,'FontSize',10);

    line([0 length(point)+1], [10 10],'color','k')
    line([length(point)+1 length(point)+1], [0 10],'color','k')

elseif strcmp(badTrialNameStr,'V4')
    axis(hDelta1(1,1),[[0 length(point)+1] [0 10]]);
    set(hDelta1(1,1),'XTick',1:length(point));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% for Protnum=1:size(protocols,2)
%
%     if size(StimDataFRSt{Protnum},1)~=size(ShamDataFRSt{Protnum},1)%If days of sham and stim does not match
%         elecsize=size(cSG,2);
%         dayStim=size(StimDataFRSt{Protnum},1)./elecsize;
%         daySham=size(ShamDataFRSt{Protnum},1)./elecsize;
%         u{Protnum}=(StimDataFRSt{Protnum});
%         v{Protnum}=(ShamDataFRSt{Protnum});
%
%         %%Averaging Stim days
%         for d=1:dayStim
%             FirstElecs(1,d)=(d*elecsize-elecsize+1);
%             U{1,d}=u{Protnum}(FirstElecs(1,d):FirstElecs(1,d)+elecsize-1);
%         end
%         %%Averaging Sham days
%         for d=1:daySham
%             FirstElecs(1,d)=(d*elecsize-elecsize+1);
%             V{1,d}=v{Protnum}(FirstElecs(1,d):FirstElecs(1,d)+elecsize-1);
%         end
%
%         dFR{Protnum}=mean(cell2mat(U),2)-mean(cell2mat(V),2);
%
%     elseif size(StimDataFRSt{Protnum},1)==size(ShamDataFRSt{Protnum},1)
%         dFR{Protnum}=((StimDataFRSt{Protnum})-(ShamDataFRSt{Protnum}));
%     end
%
%     AvgdFR{Protnum}=mean(dFR{Protnum});
%     deltaPlotDataFR{Protnum}=AvgdFR{Protnum};% not Making the first point zero
%     errordeltadataFR=dFR{Protnum};
%     errordeltagridFR{Protnum}=std(errordeltadataFR)/sqrt(size(errordeltadataFR,1));
% end


%
% %plotting the Second FR plot
% hold on
% point=1:size(protocols,2);
%
% hold(hDelta1(1,2),'on')
% e1=errorbar(subplot(hDelta1(1,2)),point,[deltaPlotDataFR{:}],[errordeltagridFR{:}],'color',"#E35B00",'LineWidth',1.5);
% for f=1:size(protocols,2)
%     w=[AvgdFR{:}];
%     scatter(point(1,f),w(1,f),30,ColCode{1, f}{1, 1},"filled")
% end
%
%
% e1.Marker='.';
% e1.CapSize=3;
%
% yline(0,'--')
%
% text(5,0.8,append('n=',num2str(size(StimDataFRSt{1,1},1))),'color','k','FontSize',10,'FontWeight','bold');
%
% title(hDelta1(1,2),'Firing rate','color','#1768ff','Fontsize',10,'FontWeight','bold')
% subtitle(hDelta1(1,2),'Stim-Sham','Fontsize',9,'FontWeight','bold')
%
%
%
% if strcmp(badTrialNameStr,'V1')
%     axis(hDelta1(1,2),[[0 length(point)+1] [-5 5]]);
%     set(hDelta1(1,2),'XTick',1:length(point));
% elseif strcmp(badTrialNameStr,'V4')
%     axis(hDelta1(1,2),[[0 length(point)+1] [-5 5]]);
%     set(hDelta1(1,2),'XTick',1:length(point));
%
% end
%
% for ProtN=2:size(protocols,2)
%     x1=dFR{1,1};% Post-Pre
%     y1=dFR{ProtN};
%     [p1,h1]=ranksum(y1,x1); %This checks if gamma power in a protocol is significantly different than the pre-stim/sham protocol
%     SGdeltaStat(1,ProtN)=double(h1);
%     SGdeltaStat(2,ProtN)=p1;
% end
%
%
%
% % Plotting the significance star
%
%
% v=[deltaPlotDataFR{:}];
% Asterixdata=(v);
% Xcentres=point;
% hold on
%
%
% stat=SGdeltaStat(1,:);
%
%
% for d=1:size(protocols,2)
%     if stat(1,d)==1
%         hold on
%         subplot(hDelta1(1,2))
%         ypoint=(Asterixdata(:,d));
%         if ypoint<0
%             text(Xcentres(:,d)-0.1,ypoint-1.3,'\ast','fontWeight','bold',HandleVisibility='off'); hold on
%         elseif ypoint>0
%             text(Xcentres(:,d)-0.2,ypoint+1.25,'\ast','fontWeight','bold',HandleVisibility='off'); hold on
%         end
%     end
% end
%




%plotting the Second FR plot
hold on
point=1:size(protocols,2);

hold(hDelta2(1,1),'on')
e1=errorbar(subplot(hDelta2(1,1)),point,[AvgdeltaFRSham{:}]-[AvgdeltaFRSham{1}],[errordeltagridFRSham{:}],'color','#ED4672','LineWidth',1.5);
for f=1:size(protocols,2)
    w=[AvgdeltaFRSham{:}]-[AvgdeltaFRSham{1}];
    scatter(point(1,f),w(1,f),30,ColCode{1, f}{1, 1},"filled")
end

e2=errorbar(subplot(hDelta2(1,1)),point,[AvgdeltaFRStim{:}]-[AvgdeltaFRStim{1}],[errordeltagridFRStim{:}],'color',"#CD9A18",'LineWidth',1.5);

for f=1:size(protocols,2)
    w=[AvgdeltaFRStim{:}]-[AvgdeltaFRStim{1}];
    scatter(point(1,f),w(1,f),30,ColCode{1, f}{1, 1},"filled")
end


e1.Marker='.';
e2.Marker='.';
e1.CapSize=3;
e2.CapSize=3;

yline(0,'--')


text(5,1.8,append('n=',num2str(size(StimDataFRSt{1,1},1))),'color','k','FontSize',9,'FontWeight','bold');

title(hDelta2(1,1),'\Delta Firing rate','color','k','Fontsize',9,'FontWeight','bold')
% subtitle(hDelta2(1,1),'Sham(Light) vs Stim(Dark)','Fontsize',9,'FontWeight','bold')
subtitle(hDelta2(1,1),{'\color[rgb]{0.9294, 0.2784, 0.4471}Sham \color[rgb]{0, 0, 0}vs \color[rgb]{0.8039, 0.6039, 0.0941}Stim'},'Fontsize',10,'FontWeight','bold');



if strcmp(badTrialNameStr,'V1')
    axis(hDelta2(1,1),[[0 length(point)+1] [-2 2]]);
    axis(hDelta2(1,2),[[0 length(point)+1] [-1 1]]);
    set(hDelta2(1,1),'XTick',1:length(point));
    set(hDelta2(1,2),'YTick',[-1 0 1]);
    Xax=gca().XAxis;
    Yax=gca().YAxis;
    set(Xax,'TickDirection','out');
    set(Xax,'TickLength',[0.02 0.025]);
    set(Yax,'TickDirection','out');
    set(Yax,'TickLength',[0.02 0.025]);
    set(gca,'FontWeight','bold');
    set(Xax,'FontSize',10);
    set(Yax,'FontSize',10);

    line([0 length(point)+1], [2 2],'color','k')
    line([length(point)+1 length(point)+1], [-2 2],'color','k')

elseif strcmp(badTrialNameStr,'V4')
    axis(hDelta2(1,1),[[0 length(point)+1] [-0.2 0.2]]);
    axis(hDelta2(1,2),[[0 length(point)+1] [-0.2 0.2]]);
    % set(hDelta1(1,:),'XTick',[]);
    set(hDelta2(1,1),'XTick',1:length(point));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%

if SessionID==0
    for Protnum=1:size(protocols,2)

        if size(StimDataFRSt{Protnum},1)~=size(ShamDataFRSt{Protnum},1)%If days of sham and stim does not match
            elecsize=size(cSG,2);
            dayStim=size(StimDataFRSt{Protnum},1)./elecsize;
            daySham=size(ShamDataFRSt{Protnum},1)./elecsize;
            u{Protnum}=(StimDataFRSt{Protnum}-StimDataFRBl{Protnum});
            v{Protnum}=(ShamDataFRSt{Protnum}-ShamDataFRBl{Protnum});

            %%Averaging Stim days
            for d=1:dayStim
                FirstElecs(1,d)=(d*elecsize-elecsize+1);
                U{1,d}=u{Protnum}(FirstElecs(1,d):FirstElecs(1,d)+elecsize-1);
            end
            %%Averaging Sham days
            for d=1:daySham
                FirstElecs(1,d)=(d*elecsize-elecsize+1);
                V{1,d}=v{Protnum}(FirstElecs(1,d):FirstElecs(1,d)+elecsize-1);
            end

            dFR{Protnum}=mean(cell2mat(U),2)-mean(cell2mat(V),2);

        elseif size(StimDataFRSt{Protnum},1)==size(ShamDataFRSt{Protnum},1)
            dFR{Protnum}=((StimDataFRSt{Protnum}-StimDataFRBl{Protnum})-(ShamDataFRSt{Protnum}-ShamDataFRBl{Protnum}));
        end

        AvgdFR{Protnum}=mean(dFR{Protnum});
        deltaPlotDataFR{Protnum}=AvgdFR{Protnum}-AvgdFR{1};%Making the first point zero
        errordeltadataFR=dFR{Protnum};
        errordeltagridFR{Protnum}=std(errordeltadataFR)/sqrt(size(errordeltadataFR,1));
    end



    %plotting the fourth FR plot
    hold on
    point=1:size(protocols,2);
    hold(hDelta2(1,2),'on')
    e=errorbar(subplot(hDelta2(1,2)),point,[deltaPlotDataFR{:}],[errordeltagridFR{:}],'color',"#CD9A18",'LineWidth',1.5);

    for f=1:size(protocols,2)
        w=[deltaPlotDataFR{:}];
        scatter(point(1,f),w(1,f),30,ColCode{1, f}{1, 1},"filled")
    end

    Xax=gca().XAxis;
    Yax=gca().YAxis;

    set(Xax,'TickDirection','out');
    set(Xax,'TickLength',[0.02 0.025]);
    set(Yax,'TickDirection','out');
    set(Yax,'TickLength',[0.02 0.025]);
    set(gca,'FontWeight','bold');
    set(Xax,'FontSize',10);
    set(Yax,'FontSize',10);

    line([0 length(point)+1], [1 1],'color','k')
    line([length(point)+1 length(point)+1], [-1 1],'color','k')

    e.Marker='.';
    e.CapSize=3;
    yline(deltaPlotDataFR{1,1},'--')
    text(3.4,0.8,append('n=',num2str(size(dFR{1,1},1)),', p=0.05'),'color','k','FontSize',9,'FontWeight','bold');

    title(hDelta2(1,2),'\Delta Firing rate','color','k','Fontsize',9,'FontWeight','bold')

    % subtitle(hDelta2(1,2),'(Stim-Sham)','Fontsize',9,'FontWeight','bold')
    subtitle(hDelta2(1,2),{'\color[rgb]{0.8039, 0.6039, 0.0941}Stim \color[rgb]{0, 0, 0}- \color[rgb]{0.9294, 0.2784, 0.4471}Sham' },'Fontsize',10,'FontWeight','bold');

    if strcmp(badTrialNameStr,'V1')
        axis(hDelta2(1,2),[[0 length(point)+1] [-1 1]]);
        set(hDelta2(1,2),'XTick',1:length(point));
    elseif strcmp(badTrialNameStr,'V4')
        axis(hDelta2(1,2),[[0 length(point)+1] [-1 1]]);
        axis(hDelta2(1,2),[[0 length(point)+1] [-1 1]]);
        % set(hDelta1(1,:),'XTick',[]);
        set(hDelta2(1,2),'XTick',1:length(point));
    end

    %In case of SG
    for ProtN=2:size(protocols,2)
        x1=dFR{1,1};% Post-Pre
        y1=dFR{ProtN};
        [p1,h1]=ranksum(y1,x1); %This checks if gamma power in a protocol is significantly different than the pre-stim/sham protocol
        SGdeltaStat(1,ProtN)=double(h1);
        SGdeltaStat(2,ProtN)=p1;
    end



    % Plotting the significance star


    v=[deltaPlotDataFR{:}];
    Asterixdata=(v);
    Xcentres=point;
    hold on


    stat=SGdeltaStat(1,:);


    for d=1:size(protocols,2)
        if stat(1,d)==1
            hold on
            subplot(hDelta2(1,2))
            ypoint=(Asterixdata(:,d));
            if ypoint<0
                text(Xcentres(:,d)-0.1,ypoint-0.3,'\ast','fontWeight','bold',HandleVisibility='off'); hold on
            elseif ypoint>0
                text(Xcentres(:,d)-0.2,ypoint+0.25,'\ast','fontWeight','bold',HandleVisibility='off'); hold on
            end
        end
    end
end




%% For band plot

for Band=1:2
    for Protnum=1:size(protocols,2)
        AvgdeltaStim{Protnum,Band}= mean((StimDataSt{Protnum,Band}-StimDataBl{Protnum,Band}),1);
        AvgdeltaSham{Protnum,Band}=mean((ShamDataSt{Protnum,Band}-ShamDataBl{Protnum,Band}),1);
        errordeltagridStim{Protnum,Band}=std((StimDataSt{Protnum,Band}-StimDataBl{Protnum,Band}))/sqrt(size(StimDataSt{Protnum,Band},1));
        errordeltagridSham{Protnum,Band}=std((ShamDataSt{Protnum,Band}-ShamDataBl{Protnum,Band}))/sqrt(size(ShamDataSt{Protnum,Band},1));
    end
end

%plotting First two band plot
hold on
point=1:size(protocols,2);
for Band=1:2 %%% For slow and fast gamma band

    hold(hDelta3(1,Band),'on')

    e1=errorbar(subplot(hDelta3(1,Band)),point,[AvgdeltaSham{:,Band}]-[AvgdeltaSham{1,Band}],[errordeltagridSham{:,Band}],'color','#ED4672','LineWidth',1.5);
    for f=1:size(protocols,2)
        w=[AvgdeltaSham{:,Band}]-[AvgdeltaSham{1,Band}];
        scatter(point(1,f),w(1,f),30,ColCode{1, f}{1, 1},"filled")
    end

    Xax=gca().XAxis;
    Yax=gca().YAxis;

    set(Xax,'TickDirection','out');
    set(Xax,'TickLength',[0.02 0.025]);
    set(Yax,'TickDirection','out');
    set(Yax,'TickLength',[0.02 0.025]);
    set(gca,'FontWeight','bold');
    set(Xax,'FontSize',10);
    set(Yax,'FontSize',10);

    line([0 length(point)+1], [0.2 0.2],'color','k')
    line([length(point)+1 length(point)+1], [-0.2 0.2],'color','k')


    e2=errorbar(subplot(hDelta3(1,Band)),point,[AvgdeltaStim{:,Band}]-[AvgdeltaStim{1,Band}],[errordeltagridStim{:,Band}],'color',"#CD9A18",'LineWidth',1.5);

    for f=1:size(protocols,2)
        w=[AvgdeltaStim{:,Band}]-[AvgdeltaStim{1,Band}];
        scatter(point(1,f),w(1,f),30,ColCode{1, f}{1, 1},"filled")
    end



    e1.Marker='.';
    e2.Marker='.';
    e1.CapSize=3;
    e2.CapSize=3;

    yline(0,'--');




    text(4.7,0.18,append('n=',num2str(size(StimData{1,1},1))),'color','k','FontSize',9,'FontWeight','bold');
    if Band==1
        title(hDelta3(1,Band),append('SG  ','\Delta(log Power)'),'color','k','Fontsize',9,'FontWeight','bold')

        subtitle(hDelta3(1,Band),{'\color[rgb]{0.9294, 0.2784, 0.4471}Sham \color[rgb]{0, 0, 0}vs \color[rgb]{0.8039, 0.6039, 0.0941}Stim'},'Fontsize',10,'FontWeight','bold');
    elseif Band==2
        title(hDelta3(1,Band),append('FG  ','\Delta(log Power)'),'color','k','Fontsize',9,'FontWeight','bold')
        subtitle(hDelta3(1,Band),{'\color[rgb]{0.9294, 0.2784, 0.4471}Sham \color[rgb]{0, 0, 0}vs \color[rgb]{0.8039, 0.6039, 0.0941}Stim'},'Fontsize',10,'FontWeight','bold');
    end
end

if strcmp(badTrialNameStr,'V1')
    axis(hDelta3(1,1),[[0 length(point)+1] [-0.2 0.2]]);
    axis(hDelta3(1,2),[[0 length(point)+1] [-0.2 0.2]]);
    set(hDelta3(1,:),'XTick',1:length(point));
elseif strcmp(badTrialNameStr,'V4')
    axis(hDelta3(1,1),[[0 length(point)+1] [-0.2 0.2]]);
    axis(hDelta3(1,2),[[0 length(point)+1] [-0.2 0.2]]);
    set(hDelta3(1,:),'XTick',[]);
    set(hDelta3(1,1),'XTick',1:length(point));
end

%% For Third and fourth Plot

for Band=1:2
    for Protnum=1:size(protocols,2)

        if size(StimDataSt{Protnum,Band},1)~=size(ShamDataSt{Protnum,Band},1)%If days of sham and stim does not match
            elecsize=size(cSG,2);
            dayStim=size(StimDataSt{Protnum,Band},1)./elecsize;
            daySham=size(ShamDataSt{Protnum,Band},1)./elecsize;
            u{Protnum,Band}=(StimDataSt{Protnum,Band}-StimDataBl{Protnum,Band});
            v{Protnum,Band}=(ShamDataSt{Protnum,Band}-ShamDataBl{Protnum,Band});

            %%Averaging Stim days
            for d=1:dayStim
                FirstElecs(1,d)=(d*elecsize-elecsize+1);
                U{1,d}=u{Protnum,Band}(FirstElecs(1,d):FirstElecs(1,d)+elecsize-1);
            end
            %%Averaging Sham days
            for d=1:daySham
                FirstElecs(1,d)=(d*elecsize-elecsize+1);
                V{1,d}=v{Protnum,Band}(FirstElecs(1,d):FirstElecs(1,d)+elecsize-1);
            end

            dPower{Protnum,Band}=mean(cell2mat(U),2)-mean(cell2mat(V),2);

        elseif size(StimDataSt{Protnum,Band},1)==size(ShamDataSt{Protnum,Band},1)
            dPower{Protnum,Band}=((StimDataSt{Protnum,Band}-StimDataBl{Protnum,Band})-(ShamDataSt{Protnum,Band}-ShamDataBl{Protnum,Band}));
        end

        AvgdPower{Protnum,Band}=mean(dPower{Protnum,Band});
        deltaPlotDataF{Protnum,Band}=AvgdPower{Protnum,Band}-AvgdPower{1,Band};%Making the first point zero
        errordeltadataF=10*dPower{Protnum,Band}; %Log power converted to dB
        errordeltagridF{Protnum,Band}=std(errordeltadataF)/sqrt(size(errordeltadataF,1));
    end
end


%Plotting last two band plot
hold on
point=1:size(protocols,2);
for Band=1:2 %%% For slow and fast gamma band
    hold(hDelta4(1,Band),'on')
    e=errorbar(subplot(hDelta4(1,Band)),point,10*[deltaPlotDataF{:,Band}],[errordeltagridF{:,Band}],'color',"#CD9A18",'LineWidth',1.5);

    for f=1:size(protocols,2)
        w=10*[deltaPlotDataF{:,Band}];
        scatter(point(1,f),w(1,f),30,ColCode{1, f}{1, 1},"filled")
    end

    Xax=gca().XAxis;
    Yax=gca().YAxis;

    set(Xax,'TickDirection','out');
    set(Xax,'TickLength',[0.02 0.025]);
    set(Yax,'TickDirection','out');
    set(Yax,'TickLength',[0.02 0.025]);
    set(gca,'FontWeight','bold');
    set(Xax,'FontSize',10);
    set(Yax,'FontSize',10);


    line([0 length(point)+1], [2 2],'color','k')
    line([length(point)+1 length(point)+1], [-2 2],'color','k')

    e.Marker='.';
    e.CapSize=3;
    yline(deltaPlotDataF{1,1},'--')
    text(3.4,1.8,append('n=',num2str(size(dPower{1,1},1)),', p=0.05'),'color','k','FontSize',9,'FontWeight','bold');
    if Band==1
        title(hDelta4(1,Band),append('SG  ','\Delta power(dB)'),'color','k','Fontsize',9,'FontWeight','bold')
        subtitle(hDelta4(1,Band),{'\color[rgb]{0.8039, 0.6039, 0.0941}Stim \color[rgb]{0, 0, 0}- \color[rgb]{0.9294, 0.2784, 0.4471}Sham' },'Fontsize',10,'FontWeight','bold');
    elseif Band==2
        title(hDelta4(1,Band),append('FG  ','\Delta power(dB)'),'color','k','Fontsize',9,'FontWeight','bold')
        subtitle(hDelta4(1,Band),{'\color[rgb]{0.8039, 0.6039, 0.0941}Stim \color[rgb]{0, 0, 0}- \color[rgb]{0.9294, 0.2784, 0.4471}Sham' },'Fontsize',10,'FontWeight','bold');
    end

    if strcmp(badTrialNameStr,'V1')
        axis(hDelta4(1,1),[[0 length(point)+1] [-2 2]]);
        axis(hDelta4(1,2),[[0 length(point)+1] [-2 2]]);
        set(hDelta4(1,:),'XTick',[]);
        set(hDelta4(1,1),'XTick',1:length(point));
        set(hDelta4(1,2),'XTick',1:length(point));

    elseif strcmp(badTrialNameStr,'V4')
        axis(hDelta4(1,1),[[0 length(point)+1] [-2 2]]);
        axis(hDelta4(1,2),[[0 length(point)+1] [-2 2]]);
        set(hDelta4(1,:),'XTick',[]);
        set(hDelta4(1,:),'XTick',1:length(point));
    end
end

%In case of SG
for Band=1:2
    for ProtN=2:size(protocols,2)
        x1=dPower{1,Band};% Post-Pre
        y1=dPower{ProtN,Band};
        [p1,h1]=ranksum(y1,x1); %This checks if gamma power in a protocol is significantly different than the pre-stim/sham protocol

        if Band==1
            SGdeltaStat(1,ProtN)=double(h1);
            SGdeltaStat(2,ProtN)=p1;
            clear x1 y1 h1 p1
        elseif Band==2
            FGdeltaStat(1,ProtN)=double(h1);
            FGdeltaStat(2,ProtN)=p1;
            clear x1 y1 h1 p1
        end
    end
end


% Plotting the significance star

for freqstat=1:2 %Band specific identity
    v=10*[deltaPlotDataF{:,freqstat}];
    Asterixdata=(v);
    Xcentres=point;
    hold on

    if freqstat==1
        stat=SGdeltaStat(1,:);
    elseif freqstat==2
        stat=FGdeltaStat(1,:);
    end


    for d=1:size(protocols,2)
        if stat(1,d)==1
            hold on
            subplot(hDelta4(1,freqstat))
            ypoint=(Asterixdata(:,d));
            if ypoint<0
                text(Xcentres(:,d)-0.1,ypoint-0.3,'\ast','fontWeight','bold',HandleVisibility='off'); hold on
            elseif ypoint>0
                text(Xcentres(:,d)-0.2,ypoint+0.3,'\ast','fontWeight','bold',HandleVisibility='off'); hold on
            end
        end
    end
end

%% Saving the figure in the fig folder
% folderSave = fullfile(FigSourceString);
% figname1 = fullfile(folderSave,append(Polarity,Session,'Stim',badTrialNameStr,'Band.tiff'));

% set(figure(1), 'Position', get(0, 'Screensize'));
% exportgraphics(figure(1),figname1,"Resolution",600)
% close all
