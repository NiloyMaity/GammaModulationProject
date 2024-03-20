
%% 0. Keep all stimulation specific protocol list in the matlab path
% To pull the protocols accordingly

%% 1. After extracting the data, we should do the bad trial anlysis per brain area.
% getAreaSpecificBadtrials do this. We get bad trial file with 'V1' or 'V4' string.

%% 2. After Bad trial analysis is done, we find Good units.
% using GetGoodSpikeElectrodes.

%% 3. With Number of good units and good trials, we proceed to save PSD,TF, Band Powers, and PSTH
% get_LFP_Spike_Save_MonkeyData

%% 4. Plotting Delta TF plots
% plotTFStimSham

%% 5. Plotting PSD and spikes in a plot
% PlotPSDandSpike 

%% 6. Plotting BandPower across
% plotBandPowerStimSham

