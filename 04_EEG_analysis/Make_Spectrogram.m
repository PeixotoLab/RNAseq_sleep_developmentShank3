function fig_handle = Make_Spectrogram(options)
%%%% plot spectrogram and hypnogram of selected animals, this script is
%%%% specific for the AD experiemnts from SK, March 2023, CM
% Modified from CM_Spectrogam_Individuals.m by MJR July/August 2023


arguments
    options.BL_file
    options.SD_file 
    options.path_to_BL_file
    options.path_to_SD_file 
    options.epoch_dur
    options.window_length_min
    options.hypnogram_pooling_window_min  
end 

BL_file             = options.BL_file;
SD_file             = options.SD_file;
path_to_BL_file     = options.path_to_BL_file;
path_to_SD_file     = options.path_to_SD_file;
epoch_duration_secs = options.epoch_dur;
window_length_min   = options.window_length_min;
hypnogram_pooling_window_min = options.hypnogram_pooling_window_min;



%%% create a directory for the output/figures
%path2='/Users/christinemuheim/Desktop/Sabrina/Figures/';
%path2='/Users/michael/code/Peixoto_lab/Lizzy/EEG_HeatMap_Figures/';

%%% get the filenames liste
%file={info(contains({info.name},'.mat')).name};  
%numfile=(size(file,2));

%startepoch=23;
%startrow=4;  % MJR: I think you mean column here??

% for the spectrogram, define how long the window for the running mean
% should be
%k=300; %running mean window, 150 = 10min, 300=20min
k = window_length_min*(60/epoch_duration_secs);    


% for the hypnogram, pool segments (to make less noisy looking?)
int_hyp = hypnogram_pooling_window_min*60/epoch_duration_secs;




% for mac
BLdata = load([path_to_BL_file,BL_file]);  % NOTE:  BLdata and SDdata are Tables!!
SDdata = load([path_to_SD_file,SD_file]);

BLdata = BLdata.EEG;
SDdata = SDdata.EEG;


% %     % for win
% %     data=importdata([path+"\"+file{anim}]);
    

%get the frequency range
%Frequency=data(startepoch-1,startrow:end-1);


%%% extract only the FFT data
namesBL    = BLdata.Properties.VariableNames;
namesSD    = SDdata.Properties.VariableNames;
startcolBL = min(find(contains(namesBL,'Hz')));
startcolSD = min(find(contains(namesSD,'Hz')));
if startcolBL ~= startcolSD
    error('In Make_Spectrogram: the BL file and SD do not have EEG data in the same columns.')
end 

EEG_BL_table = BLdata(1:end-1,startcolBL:end);
EEG_SD_table = SDdata(1:end-1,startcolSD:end);  

% convert from cell to mat
% EEG_BL = cell2mat(table2array(EEG_BL));
% EEG_SD = cell2mat(table2array(EEG_SD));

EEG_BL = table2array(EEG_BL_table);
if iscell(EEG_BL_table{1,4})
    EEG_BL = cell2mat(EEG_BL);
    %EEG_BL = cat(1,EEG_BL); % faster?
end
EEG_SD = table2array(EEG_SD_table);
if iscell(EEG_SD_table{1,4})
    EEG_SD = cell2mat(EEG_SD);
end 

% get reference values, average accross all times and frequencies
RefEEG_BL = mean(EEG_BL,'all','omitnan');
RefEEG_SD = mean(EEG_SD,'all','omitnan');


%normalize the entire EEGs by those values
EEGrel_BL = EEG_BL./RefEEG_BL;
EEGrel_SD = EEG_SD./RefEEG_SD;

%% Spectrograms
%use a running mean in the time domain, for 10min= 150epoch
EEG_smoothed_BL = movmean(EEGrel_BL,k,1);
EEG_smoothed_SD = movmean(EEGrel_SD,k,1);


%% Hypnogram
% extract the scoring values
Scores_BL = BLdata.Stage(1:end-1);
Scores_SD = SDdata.Stage(1:end-1);

%find all the epochs per state, do not differentiate artifacts
IndxW_BL = find(strcmp(Scores_BL,'W*') | strcmp(Scores_BL,'W') | strcmp(Scores_BL,'WA'));
IndxN_BL = find(strcmp(Scores_BL,'N*') | strcmp(Scores_BL,'N') | strcmp(Scores_BL,'NR'));
IndxR_BL = find(strcmp(Scores_BL,'R*') | strcmp(Scores_BL,'R'));

IndxW_SD = find(strcmp(Scores_SD,'W*') | strcmp(Scores_SD,'W') | strcmp(Scores_SD,'WA'));
IndxN_SD = find(strcmp(Scores_SD,'N*') | strcmp(Scores_SD,'N') | strcmp(Scores_SD,'NR'));
IndxR_SD = find(strcmp(Scores_SD,'R*') | strcmp(Scores_SD,'R'));

Hypno1_BL = NaN(length(Scores_BL),1);
Hypno1_BL(IndxW_BL) = 1;
Hypno1_BL(IndxN_BL) = 0;
Hypno1_BL(IndxR_BL) = -1;

Hypno1_SD = NaN(length(Scores_SD),1);
Hypno1_SD(IndxW_SD) = 1;
Hypno1_SD(IndxN_SD) = 0;
Hypno1_SD(IndxR_SD) = -1;

% use mode to find the stage that is most often occuring in a 2min
% window
clear i
% int=30; %how many epochs are pooled
segments_BL = size(Hypno1_BL,1)/int_hyp;
segments_SD = size(Hypno1_SD,1)/int_hyp;
for i=1:segments_BL
     Hypno_BL(i)=mode(Hypno1_BL(i*int_hyp-int_hyp+1:i*int_hyp));
end
for i=1:segments_SD
     Hypno_SD(i)=mode(Hypno1_SD(i*int_hyp-int_hyp+1:i*int_hyp));
end



%AnimalID=extractBefore(file{anim},'_Calcium.mat');
%%% plotting, take two baseline files and plot first the spectrogram, then
%%% the hypnogram underneath it
% xlabelvector_heatmap = [1800 3600 5400 7200 9000 10800 12600 14400 16200 18000 19800 21600];
% xlabelvector_hypno   = [60 120 180 240 300 360 420 480 540 600 660 720];
labels    = namesBL;
labels    = labels(4:end-1);
labels2   = strrep(labels,'Hz','');
Hz_labels = round(cellfun(@str2num, labels2),2);


%[d,three_Hz_in_columns]=min(abs(Hz_labels-3));
%ylabelvector_heatmap = three_Hz_in_columns:three_Hz_in_columns
ylabelvector_heatmap = 3:3:max(Hz_labels);  %[3 6 9 12 15 18 21 24];

two_hours_in_epochs  = (2*60*60)/epoch_duration_secs;
%xlabelvector_heatmap = two_hours_in_epochs:two_hours_in_epochs:12*two_hours_in_epochs;
xlabelvector_heatmap = 0:two_hours_in_epochs:12*two_hours_in_epochs;

%xlabelvector_hypno   = int_hyp*2:int_hyp*2:int_hyp*24;
xlabelvector_hypno = 0:(2*length(Hypno_BL))/24:length(Hypno_BL);

Freq_as_numbers = Hz_labels(ylabelvector_heatmap);  % AHA: this will only work if using 1-Hz freq bins.  
Freq = arrayfun(@num2str, Freq_as_numbers, 'UniformOutput', 0);
%Freq=({'1.6' '3.9' '6.3' '8.6' '10.9' '13.3' '15.6' '17.9'});

Time=({'0','2','4','6','8','10','12','14','16','18','20','22','24'});
clims=[0 4];
%to get the titles for the plot:
BL_filename_for_title = extractBefore(BL_file,'_c_FFT');
SD_filename_for_title = extractBefore(SD_file,'_c_FFT');
% figurename1=fullfile([path2,'Spectrogram_BL_SD_',animal1])
% figurename2=fullfile([path2,'Spectrogram_BL_SD_',animal2])


% --- Make the figure ---------------------------------
F1=figure; 
F1.Name = 'Spectrogram';
F1.Renderer = 'painters';
F1.Position = [81 303 1277 535];
myfontsize = 12;

% BL spectrogram
subplot(3,6,[1:3 7:9])
imagesc((EEG_smoothed_BL'),clims)
ax=F1.Children;
%ax.YTick=ylabelvector_heatmap;
%ax.YTickLabel=Freq;
ax.YTickLabel = Hz_labels(ax.YTick);
ax.XTick=xlabelvector_heatmap;
ax.XTickLabel=Time;
ax.YDir='normal';
ax.YLabel.String = 'Frequency (Hz)';
ax.FontSize = myfontsize;
title(BL_filename_for_title,'Spectrogram for Baseline','interpreter','none')

% BL hypnogram
subplot(3,6,[13:15])
plot(Hypno_BL,'LineWidth',1,'Color','k')
ax=F1.Children(1);
ax.YTick=[-1 0 1];
ax.YTickLabel={'REM' 'NREM' 'Wake'};
% ax.XTick=xlabelvector_hypno;
ax.XTick=xlabelvector_hypno;
ax.XTickLabel=Time;
ax.XLabel.String='ZT';
ax.YLim=[-1.5 1.5];
ax.XLim = [0 length(Hypno_BL)];
ax.FontSize = myfontsize;

% SD spectrogram
subplot(3,6,[4:6 10:12])
imagesc((EEG_smoothed_SD'),clims)
ax=F1.Children(1);
%ax.YTick=ylabelvector_heatmap;
ax.YTickLabel=[];
ax.XTick=xlabelvector_heatmap;
ax.XTickLabel=Time;
ax.YDir='normal';
ax.FontSize = myfontsize;
title(SD_filename_for_title,'Spectrogram for SD','interpreter','none')
a=colorbar;
ylabel(a,'Rel EEG power','FontSize',10,'Rotation',270);
a.Label.Position(1) = 2.2;
a.Position=[0.92 0.42 0.03 0.49];
a.FontSize = myfontsize;

% SD hypnogram
subplot(3,6,[16:18])
plot(Hypno_SD,'LineWidth',1,'Color','k')
ax=F1.Children(1);
ax.YTick=[-1 0 1];
ax.YTickLabel={'REM' 'NREM' 'Wake'};
ax.XTick=xlabelvector_hypno;
ax.XTickLabel=Time;
ax.XLabel.String='ZT';
ax.YLim=[-1.5 1.5];
ax.XLim = [0 length(Hypno_SD)];
ax.FontSize = myfontsize;
% saveas(gcf,figurename1,'fig')
% saveas(gcf,figurename1,'tif')
fig_handle = gcf;

