%% This small code provides right color ID for different protocol, associated protocol type/title, and stimulation block ids

function [ColCode,titleString,StimblockID]=pickIDs(SessionId)
   
    % %% RGB codes
    % preStimCol{1,1} = [62 199 11]./255; %Green
    % StimCol1{1,1} = [255, 165, 0]./255;%orange-ish
    % postStimCol1{1,1} = [210 19 18]./255;%Red
    % postStimCol2{1,1} = [176 74 90]./255;%Carmine Pink
    % StimCol2{1,1}= [241, 196, 15]./255;%ButterCup
    % postStimCol3{1,1} = [0,0,0]./255;%Black
    % check1Col{1,1} = [255, 0, 204]./255;%Hot magenta
    % check2Col{1,1} = [88, 0, 255]./255;%Electric Indigo
    % check3Col{1,1} = [22 160 133]./255;%Mountain meadow
    % 
    % %% HEX codes
    % preStimCol{1,2}="#3EC70B";%Green 62, 199, 11
    % StimCol1{1,2}="#FFA500";%Orange-ish 255, 165, 0
    % postStimCol1{1,2}="#D21312";%Red 210, 19, 18
    % postStimCol2{1,2}="#B04A5A";%Carmine pink 176, 74, 90
    % StimCol2{1,2}="#F1C40F";%Buttercup 241, 196, 15
    % postStimCol3{1,2}="#000000";%Black 0, 0, 0
    % check1Col{1,2}="#FF00CC"; %Hot Magenta 255, 0, 204
    % check2Col{1,2}="#5800FF";%Electric Indigo 88, 0, 255
    % check3Col{1,2}="#16A085";%Mountain Medow 22, 160, 133




    %% RGB codes
    preStimCol{1,1} = [151 28 81]./255; %Wine Grape

    StimCol1{1,1} = [255, 165, 0]./255;%orange-ish

    postStimCol1{1,1} = [210 19 18]./255;%Red
    postStimCol2{1,1} = [176 74 90]./255;%Carmine Pink

    StimCol2{1,1}= [241, 196, 15]./255;%ButterCup

    postStimCol3{1,1} = [8 48 35]./255;%Cyprus
    check1Col{1,1} = [45 93 78]./255;%brunswick green
    check2Col{1,1} = [94 141 37]./255;%Vida loca
    check3Col{1,1} = [156 176 23]./255;%Stipule Green

    %% HEX codes
    preStimCol{1,2}="#971C51";% Wine Grape 151,28, 81

    StimCol1{1,2}="#FFA500";%Orange-ish 255, 165, 0

    postStimCol1{1,2}="#D21312";%Red 210, 19, 18
    postStimCol2{1,2}="#B04A5A";%Carmine pink 176, 74, 90

    StimCol2{1,2}="#F1C40F";%Buttercup 241, 196, 15

    postStimCol3{1,2}="#083023";%Cyprus 8, 48, 35
    check1Col{1,2}="#2D5D4E"; %brunswick green 45, 93, 78
    check2Col{1,2}="#5E8D25";%Vida loca 94, 141, 37
    check3Col{1,2}="#9CB017";%Stipule Green 156, 176, 23

    if SessionId==0   
    ColCode={preStimCol StimCol1 postStimCol3 check1Col check2Col check3Col};
    titleString={'pre' 'Stim/Sham' 'post' 'check1' 'check2' 'check3'};
    StimblockID=2;
    elseif SessionId==1
    ColCode={preStimCol StimCol1 postStimCol1 StimCol2 postStimCol3 check1Col check2Col check3Col};
    titleString={'pre' 'Stim/Sham' 'post' 'Stim/Sham' 'postStim' 'check1' 'check2' 'check3'};
    StimblockID=[2 4];
    elseif SessionId==2
    ColCode={preStimCol StimCol1 postStimCol1 postStimCol2 StimCol2 postStimCol3 check1Col check2Col check3Col};
    titleString={'pre' 'Stim/Sham' 'post1' 'post2' 'Stim/Sham' 'post1' 'check1' 'check2' 'check3'};
    StimblockID=[2 5];
    end


