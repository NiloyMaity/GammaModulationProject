% Bad trial is done area specific, that Id will be used to mention the Brain Area

%%Please keep the Stimulation Specific Protocol List in MATLAB path.
%StimulationType={'tDCS','tACS'}
%condition={'Stim','Sham'};
%Polarity={'Cathodal','Anodal' or 'SG','FG','Alpha'};
%SessionID=0,1 or 2 for Single, Dual or Dual with 60 minutes gap

%normalize flag=1 plots normalized PSTH, 0 plots general PSTH
%BadDayFlag=1 takes away the analog days, 0 ignores it.

%A is area Id, A=1 is V1, A=2 is V4



function plotAll(folderSourceString,StimulationType,Polarity,SessionID,normalizeFlag,BadDayFlag,A)

FigSourceString=fullfile(folderSourceString,'Figures');
if ~exist(FigSourceString, 'dir')
    mkdir(FigSourceString)
end

if SessionID==0
    Session='single';
elseif SessionID==1
    Session='dual';
elseif SessionID==2
    Session='dual60';
end

Allcondition={'Sham','Stim'};

%% Plots The PSD and PSTH
for c=1:2
    condition=Allcondition{1,c};
plotPSDandSpike(folderSourceString,StimulationType,condition,Polarity,SessionID,normalizeFlag,BadDayFlag);
end

%% PLots the band power and firing rate
plotBandPowerStimSham(folderSourceString,StimulationType,Polarity,SessionID,BadDayFlag)

%% Plots the time frequency plot
plotTFStimSham(folderSourceString,StimulationType,Polarity,SessionID,A,BadDayFlag)

%% Saving the figure in the fig folder

folderSave = fullfile(FigSourceString);
figname1 = fullfile(folderSave,append(Polarity,Session,'.tiff'));

set(figure(1), 'Position', get(0, 'Screensize'));
exportgraphics(figure(1),figname1,"Resolution",600)
close all
end
