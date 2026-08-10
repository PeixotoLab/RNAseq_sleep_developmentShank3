%%%% plot spectrogram and hypnogram of selected animals, this script is
%%%% specific for the AD experiemnts from SK, March 2023, CM

%% 
clear;
close all hidden;
clc;
%% file location
path = uigetdir(); 
info = dir(path); %contents of selected directory

%%% create a directory for the output/figures
%path2='/Users/christinemuheim/Desktop/Sabrina/Figures/';
path2='/Users/michael/code/Peixoto_lab/Lizzy/EEG_HeatMap_Figures/';

%%% get the filenames liste
file={info(contains({info.name},'.mat')).name};  
numfile=(size(file,2));
startepoch=23;
startrow=4;
% for the spectrogram, define how long the window for the running mean
% should be
k=300; %running mean window, 10min
%define a matrix for the new hypnogram
Hypno1=NaN(21600,1);
%%% for each file, plot the spectrogram, normalized to total EEG power,
%%% make 10min bins, also plot the hypnogram underneath each spectrogram.
%%% safe as .fig file and the load later again to plot in the same window

for anim=1:numfile
    % for mac
    data=importdata([path+"/"+file{anim}]);
% %     % for win
% %     data=importdata([path+"\"+file{anim}]);
    %get the frequency range
    Frequency=data(startepoch-1,startrow:end-1);
    %%% extract only the FFT data
    EEG=data(startepoch:end-1,startrow:end-1);%%%% where the epoch data starts, 21600 epochs, 26 Hz bins
    % convert from cell to mat
    EEG=cell2mat(EEG);
    % get reference value, average accross all times and frequencies
    RefEEG=mean(EEG,'all','omitnan');
    %normalize the entire EEG by that value
    EEGrel=EEG./RefEEG;

    %% Spectrograms
    %use a running mean in the time domain, for 10min= 150epoch
    EEGave(:,:,anim)=movmean(EEGrel,k,1);

    %% Hypnogram
    % extract the scoring values
    Score=data(startepoch:end-1,2);
    %find all the epochs per state, do not differentiate artifacts
    IndxW=find(strcmp(Score,'W*')| strcmp(Score,'W'));
    IndxN=find(strcmp(Score,'N*')| strcmp(Score,'NR'));
    IndxR=find(strcmp(Score,'R*')| strcmp(Score,'R'));

    Hypno1(IndxW)=1;
    Hypno1(IndxN)=0;
    Hypno1(IndxR)=-1;

    % use mode to find the stage that is most often occuring in a 2min
    % window
    clear i
    int=30; %how many epochs are pooled
    segments=size(Hypno1,1)/int;
    for i=1:segments
         Hypno(i,anim)=mode(Hypno1(i*int-int+1:i*int));
    end

end

AnimalID=extractBefore(file{anim},'_Calcium.mat');
%%% plotting, take two baseline files and plot first the spectrogram, then
%%% the hypnogram underneath it
ylabelvector=[3 6 9 12 15 18 21 24];
xlabelvector=[1800 3600 5400 7200 9000 10800 12600 14400 16200 18000 19800 21600];
xlabelvector2=[60 120 180 240 300 360 420 480 540 600 660 720];

Freq=({'1.6' '3.9' '6.3' '8.6' '10.9' '13.3' '15.6' '17.9'});  % you can get these by looking at the frequency ranges by column.  ie column 3 has 1.6Hz, etc.
Time=({'2','4','6','8','10','12','14','16','18','20','22','24'});
clims=[0 4]
%to get the titles for the plot:
animal1=extractBefore(file{1,1},'_BL_cFFT.mat');
animal2=extractBefore(file{1,3},'_BL_cFFT.mat');
figurename1=fullfile([path2,'Spectrogram_BL_SD_',animal1])
figurename2=fullfile([path2,'Spectrogram_BL_SD_',animal2])

F1=figure; %first animal
subplot(3,6,[1:3 7:9])
imagesc((EEGave(:,:,1)'),clims)
set(gca,'YTick',ylabelvector,'YTickLabel',Freq);
set(gca,'XTick',xlabelvector,'XTickLabel',Time);
set(gca,'YDir','normal')
ylabel('Frequency')
title(animal1,'Spectrogram for Baseline')
subplot(3,6,[13:15])
plot(Hypno(:,1),'LineWidth',1,'Color','k')
set(gca,'YTick',[-1 0 1],'YTickLabel',{'REM' 'NREM' 'Wake'})
set(gca,'XTick',xlabelvector2,'XTickLabel',Time);
xlabel('ZT')
ylim([-1.5 1.5])
xlim([0 720])
subplot(3,6,[4:6 10:12])
imagesc((EEGave(:,:,2)'),clims)
set(gca,'YTick',ylabelvector,'YTickLabel',[]);
set(gca,'XTick',xlabelvector,'XTickLabel',Time);
set(gca,'YDir','normal')
title(animal1,'Spectrogram for SD')
a=colorbar;
ylabel(a,'Rel EEG power','FontSize',10,'Rotation',270);
a.Label.Position(1) = 2.7;
a.Position=[0.92 0.42 0.03 0.49];
subplot(3,6,[16:18])
plot(Hypno(:,2),'LineWidth',1,'Color','k')
set(gca,'XTick',xlabelvector2,'XTickLabel',Time);
xlabel('ZT')
ylim([-1.5 1.5])
xlim([0 720])
saveas(gcf,figurename1,'fig')
saveas(gcf,figurename1,'tif')

F2=figure;
subplot(3,6,[1:3 7:9])
imagesc((EEGave(:,:,3)'),clims)
set(gca,'YTick',ylabelvector,'YTickLabel',Freq);
set(gca,'XTick',xlabelvector,'XTickLabel',Time);
set(gca,'YDir','normal')
ylabel('Frequency')
title(animal2,'Spectrogram for Baseline')
subplot(3,6,[13:15])
plot(Hypno(:,3),'LineWidth',1,'Color','k')
set(gca,'YTick',[-1 0 1],'YTickLabel',{'REM' 'NREM' 'Wake'})
set(gca,'XTick',xlabelvector2,'XTickLabel',Time);
xlabel('ZT')
ylim([-1.5 1.5])
xlim([0 720])

subplot(3,6,[4:6 10:12])
imagesc((EEGave(:,:,4)'),clims)
set(gca,'YTick',ylabelvector,'YTickLabel',[]);
set(gca,'XTick',xlabelvector,'XTickLabel',Time);
set(gca,'YDir','normal')
title(animal2,'Spectrogram for SD')
a=colorbar;
ylabel(a,'Rel EEG power','FontSize',10,'Rotation',270);
a.Label.Position(1) = 2.7;
a.Position=[0.92 0.42 0.03 0.49];
subplot(3,6,[16:18])
plot(Hypno(:,4),'LineWidth',1,'Color','k')
set(gca,'XTick',xlabelvector2,'XTickLabel',Time);
xlabel('ZT')
ylim([-1.5 1.5])
xlim([0 720])
saveas(gcf,figurename2,'fig')
saveas(gcf,figurename2,'tif')

%%% plotting, take two baseline files and plot first
