MainFolder= 'Working Codes';
folderSourceString='F:';
for DisplayType=2

    if DisplayType==1
        DataFolder= 'New code Revert TF';
        df = 4;
        NumberOfPoints = 26;
        freqLimsHz = [0 100];   

     elseif DisplayType==2
        DataFolder= 'New code baseline';
        df = 2;
        NumberOfPoints = 51;
        freqLimsHz = [0 100];
    end

    Path = fullfile(folderSourceString,'Programs',MainFolder,'June',DataFolder);
    % colorNames = [{'b'} {'g'} {'r'}];
    % Path= 'F:\Programs\Working Codes\TF_Data';

    hFig = figure(1);
    set(hFig,'PaperPositionMode','auto')
    tickLen = 1.3*[0.02 0.05];
    
    fontSize = 9;
    lineWidth = 1.2;
    tickLabel1 = [0 0.5];

    for iprot=1:5
        if iprot==1
            DataStructure = dir(fullfile(Path,'ShamPre*.mat')); % Struct directory where the .mat files are saved
        elseif iprot==2
            DataStructure = dir(fullfile(Path,'ShamPost*.mat')); % Struct directory where the .mat files are saved
        elseif iprot==3
            DataStructure = dir(fullfile(Path,'StimPre*.mat')); % Struct directory where the .mat files are saved
        elseif iprot==4
            DataStructure = dir(fullfile(Path,'StimPost*.mat')); % Struct directory where the .mat files are saved
        elseif iprot==5
            DataStructure = dir(fullfile(Path,'Check*.mat')); % Struct directory where the .mat files are saved
        end

        SubjN =  numel(DataStructure); %Gets number of subjects

        ProtocolCell = cell(1,SubjN); % preallocate cell array for Pre/Post protocol



        %%%%%%%%%%%%%%%%%%%%% Collects Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
            clear RawProtocolCell
        

        for SubjId = 1:SubjN
            SubjectData = load(fullfile(Path,DataStructure(SubjId).name));
            RawProtocolCell(SubjId) = struct2cell(SubjectData); % Copies the struct data to usable preallocated cell in each loop. RawProtocol has stpower/blpower.
        end

         if DisplayType==2
            ProtocolCell{1,1}= 10*(log10(RawProtocolCell{1,1}));
            ProtocolCell{1,2}= 10*(log10(RawProtocolCell{1,2}));
            ProtocolCell{1,3}= 10*(log10(RawProtocolCell{1,3}));
            

            AllMatrix = cat(3,ProtocolCell{:}); % join all of the original matrices into a 1 x fList x NUMBER OF SUBJECT array
            MeanData = mean(AllMatrix,3); % Average plottable data

            ProtocolCell{1,SubjN+1} = MeanData; % Add average data in the data grid.


        end  
        

        

        RawAllMatrix = cat(3,RawProtocolCell{:}); % join all of the original matrices into a 1 x fList x NUMBER OF SUBJECT array
        RawMeanData = mean(RawAllMatrix,3); % Average plottable data
        SubjN = SubjN+1;
        RawProtocolCell{1,SubjN} = RawMeanData; % Add average data in the data grid.


         %%%%%%%%%%%%%%% Sets up timepoints and Frequencypoints %%%%%%%%%%%%%%
        fStart = 0;
        frequencyPoints = fStart + (0:NumberOfPoints-1)*df; %Gets Frequency list

          tStart = -0.7220; % Depending on multitaper method
            n = 72;
            timeLimsS = [0 0.7];

        dt = 0.025;       % Depending on windowstep
        timePoints = tStart + (0:n-1)*dt;  % Gets timepoint list


        % time interval to be displayed in seconds. Stimulus onset is at 0.
        cLimsDiff = [-2 3];  % colormap limits for change in power
        deltaPowerRange = [-8 12];

       if DisplayType==2 
        data_shade(1,:) = ProtocolCell{1,1};
        data_shade(2,:) = ProtocolCell{1,2};
        data_shade(3,:) = ProtocolCell{1,3};
        
        options.handle     = figure(1);
        subplot(1,2,1)
         if iprot==1
        options.color_area = [153 255 255]./255;    % Blue theme
        options.color_line = [0.3010 0.7450 0.9330];
          elseif iprot==2 
        options.color_area = [153 153 255]./255;    % Orange theme
        options.color_line = [0 0.4470 0.7410];
          elseif iprot==3
         options.color_area = [153 255 153]./255;    % Orange theme
         options.color_line = [0.4660 0.6740 0.1880];
         elseif iprot==4
         options.color_area = [255 153 153]./255;    % Orange theme
         options.color_line = [0.6350 0.0780 0.1840];
         elseif iprot==5
         options.color_area = [224 224 224]./255;    % Orange theme
         options.color_line = [0.9290 0.6940 0.1250]; 
         end
%         legend('PreSham','PostSham','PreStim','PostStim','(+10)');
        

        options.alpha      = 0.5;
        options.line_width = 2;
        options.error      = 'sem';
        plot_areaerrorbar(data_shade, options)
        
       end
       hold on
    end
    title('Delta PSD')
    xlabel('Frequency(Hz)')
    ylabel('DeltaPower(dB)')
end
%        legend('SEM1','SEM2','SEM3','SEM4','SEM5',Location='south');
% legend('PreSham','PostSham','PreStim','PostStim','(+10)',Location='north');

p= categorical({'Alpha','Slow Gamma','Fast Gamma'});
p= reordercats(p,{'Alpha','Slow Gamma','Fast Gamma'});
            load Tdata.mat
            load CollectAvgT.mat
            
%     SEM_data = std(data, [], 2)./ sqrt(size(data,2));
            subplot(1,2,2)
            B=bar(p,data,(1),'grouped');

            B(1).FaceColor = [0.3010 0.7450 0.9330];
            B(2).FaceColor = [0 0.4470 0.7410];
            B(3).FaceColor = [0.4660 0.6740 0.1880];
            B(4).FaceColor = [0.6350 0.0780 0.1840];
            B(5).FaceColor = [0.9290 0.6940 0.1250]; 

            ylim([-4 4]);

% 
           
            hold on


Calcul_grid = CollectAvg(1:3,:);
error_grid = std(Calcul_grid)/sqrt(size(Calcul_grid,1));

error_grid = [error_grid(1,1),error_grid(1,2),error_grid(1,3),error_grid(1,4),error_grid(1,5)
              error_grid(1,6),error_grid(1,7),error_grid(1,8),error_grid(1,9),error_grid(1,10)
              error_grid(1,11),error_grid(1,12),error_grid(1,13),error_grid(1,14),error_grid(1,15)];

Xcentres=vertcat(B.XEndPoints);
data=data';
error_grid=error_grid';
errorbar(Xcentres(:),data(:),error_grid(:),'k.')
title('Change in band specific Power')
    xlabel('Frequency(Hz)')
    ylabel('DeltaPower(dB)')
hold off

