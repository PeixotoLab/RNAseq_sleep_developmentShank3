% November 2022- update to LM_SleepAnalysis_100622.m 
% Michael Rempe  
 
close all hidden;
clear all;


%  --- Load all parameters from InputParameters.m (Change this line to your input file name saved in the same directory. ie if your file is called MyInputParams.m you will have S=MyInputParams;)
S = InputParamsCAITLIN_all_forspectra_p30;
%S = InputParamsCAITLIN_all_forspectra_P45;
%S = InputParamsCAITLIN_all_notforspectra_p24;



% ----- You shouldn't need to change anything below this line ------------------------------------------
% ------------------------------------------------------------------------------------------------------
disp(['I am using the input parameters file ',S.FileName,' which contains the function ',S.FuncName])
disp('If this is not what you intended, please cancel and start over.')
epochs = S.epochs;               
numhrs = S.numhrs;
%Light  = S.Light;       
%Dark   = S.Dark;        
swa    = S.swa;         
% WakeEpocsDWTBL      = S.WakeEpocsDWTBL;         
% NREMEpocsDWTBL      = S.NREMEpocsDWTBL;
% REMEpocsDWTBL       = S.REMEpocsDWTBL;
SD_length_hrs       = S.SD_length_hrs;
epoch_duration_secs = S.epoch_duration_secs;
firstNREM_episode_duration_epochs = S.firstNREM_episode_duration_epochs;  
NREM_char = S.NREM_char; 
EEGLowerLimit_Hz = S.EEGLowerLimit_Hz; 

Bout_Minimums.W = S.Bout_Minimums.W;        
Bout_Minimums.N = S.Bout_Minimums.N;
Bout_Minimums.R = S.Bout_Minimums.R;
LegendLabels    = S.LegendLabels; 


% If Analyze_TIS_DP_6hr_segments is not there, assume a reasonable default
if isfield(S,'Analyze_TIS_DP_6hr_segments')
    Analyze_TIS_DP_6hr_segments = S.Analyze_TIS_DP_6hr_segments;
else
    Analyze_TIS_DP_6hr_segments = false;
end 

if isfield(S,'Analyze_Bouts_DP_6hr_segments')
    Analyze_Bouts_DP_6hr_segments = S.Analyze_Bouts_DP_6hr_segments;
else
    Analyze_Bouts_DP_6hr_segments = false;
end 

if isfield(S,'EEGLowerLimit_Hz')
    EEGLowerLimit_Hz = S.EEGLowerLimit_Hz;
else
    EEGLowerLimit_Hz = 0.25;
end 

if isfield(S,'Normalization')
    Normalization = S.Normalization;
else
    Normalization = 'MeanPowerAllStates'; %'AreaUnderCurve'; %3.22.24: Marcos said to make MeanPower the default
end 

if isfield(S,'PlotHourlyBaselineNREM_Delta')
    PlotHourlyBaselineNREM_Delta = S.PlotHourlyBaselineNREM_Delta;
else
    PlotHourlyBaselineNREM_Delta = false;
end

if isfield(S,'SeparateSpectralIntoLPDP')
    SeparateSpectralIntoLPDP = S.SeparateSpectralIntoLPDP;
else
    SeparateSpectralIntoLPDP = false;
end


if isfield(S,'Analyze_Recovery_2hr_bins')
    Analyze_Recovery_2hr_bins = S.Analyze_Recovery_2hr_bins;
else
    Analyze_Recovery_2hr_bins = false;
end

if isfield(S,'Sexes')
    Sexes = S.Sexes;
else 
    Sexes = {'Male','Female'};
end 

if isfield(S,'SexVarName')
    SexVarName = S.SexVarName;
else 
    SexVarName = 'Sex';
end 


% Paths to data
path_to_files_WT_BL.(Sexes{1})  = S.path_to_files_WT_BL.(Sexes{1});
path_to_files_WT_SD.(Sexes{1})  = S.path_to_files_WT_SD.(Sexes{1});
path_to_files_Mut_BL.(Sexes{1}) = S.path_to_files_Mut_BL.(Sexes{1});
path_to_files_Mut_SD.(Sexes{1}) = S.path_to_files_Mut_SD.(Sexes{1});
path_to_files_WT_BL.(Sexes{2})  = S.path_to_files_WT_BL.(Sexes{2});
path_to_files_WT_SD.(Sexes{2})  = S.path_to_files_WT_SD.(Sexes{2});
path_to_files_Mut_BL.(Sexes{2}) = S.path_to_files_Mut_BL.(Sexes{2});
path_to_files_Mut_SD.(Sexes{2}) = S.path_to_files_Mut_SD.(Sexes{2});

load_data_from_mat_file_instead = S.load_data_from_mat_file_instead; 
MatFileContainingData = S.MatFileContainingData; 



if ~isempty(path_to_files_WT_BL.(Sexes{1}))
    ffl_list_WT_BL.(Sexes{1})  = {dir(fullfile(path_to_files_WT_BL.(Sexes{1}), '*.mat')).name};  else ffl_list_WT_BL.(Sexes{1})  = {}; end
if ~isempty(path_to_files_WT_SD.(Sexes{1}))
    ffl_list_WT_SD.(Sexes{1})  = {dir(fullfile(path_to_files_WT_SD.(Sexes{1}), '*.mat')).name};  else ffl_list_WT_SD.(Sexes{1})  = {}; end
if ~isempty(path_to_files_Mut_BL.(Sexes{1}))
    ffl_list_Mut_BL.(Sexes{1}) = {dir(fullfile(path_to_files_Mut_BL.(Sexes{1}),'*.mat')).name};  else ffl_list_Mut_BL.(Sexes{1}) = {}; end
if ~isempty(path_to_files_Mut_SD.(Sexes{1}))
    ffl_list_Mut_SD.(Sexes{1}) = {dir(fullfile(path_to_files_Mut_SD.(Sexes{1}),'*.mat')).name};  else ffl_list_Mut_SD.(Sexes{1}) = {}; end
if ~isempty(path_to_files_WT_BL.(Sexes{2}))
    ffl_list_WT_BL.(Sexes{2})  = {dir(fullfile(path_to_files_WT_BL.(Sexes{2}), '*.mat')).name};  else ffl_list_WT_BL.(Sexes{2})  = {}; end
if ~isempty(path_to_files_WT_SD.(Sexes{2}))
    ffl_list_WT_SD.(Sexes{2})  = {dir(fullfile(path_to_files_WT_SD.(Sexes{2}), '*.mat')).name};  else ffl_list_WT_SD.(Sexes{2})  = {}; end 
if ~isempty(path_to_files_Mut_BL.(Sexes{2}))
    ffl_list_Mut_BL.(Sexes{2}) =  {dir(fullfile(path_to_files_Mut_BL.(Sexes{2}),'*.mat')).name}; else ffl_list_Mut_BL.(Sexes{2}) = {}; end
if ~isempty(path_to_files_Mut_SD.(Sexes{2}))
    ffl_list_Mut_SD.(Sexes{2}) =  {dir(fullfile(path_to_files_Mut_SD.(Sexes{2}),'*.mat')).name}; else ffl_list_Mut_SD.(Sexes{2}) = {}; end

total_recordings = length(ffl_list_WT_BL.(Sexes{1}))  + length(ffl_list_WT_SD.(Sexes{1}))  + ...
                   length(ffl_list_Mut_BL.(Sexes{1})) + length(ffl_list_Mut_SD.(Sexes{1})) + ...
                   length(ffl_list_WT_BL.(Sexes{2}))  + length(ffl_list_WT_SD.(Sexes{2}))  + ...
                   length(ffl_list_Mut_BL.(Sexes{2})) + length(ffl_list_Mut_SD.(Sexes{2}));                
disp(['You are running the code on ',num2str(total_recordings),' recordings total.'])
% -------------------------------------------------------------------------------------------------------------------------------------
% -------------------------------------------------------------------------------------------------------------------------------------

disp('Extracting EEG data ...')


% If you can just read it from a .mat file, do that and save the time of doing all those 'load' commands over OneDrive
if load_data_from_mat_file_instead
    disp(['WARNING: Loading data from a .mat file (',S.MatFileContainingData,' rather than from individual *_c_FFT files. If this is not what you want, set load_data_from_mat_file to 0.'])
    load(S.MatFileContainingData)
    %load('Abrar_Hemi_data_Jan2024.mat')         % Abrar's data
    %load('Adult_Mut_WT_BL_SD_M_F_Aug2023.mat') % Lizzy data from Aug 2023

else 


    % Males
    [swaNREMMutSDhourly,swaREMMutSDhourly,swaWakeMutSDhourly,EEGMutSDWake,EEGMutSDNREM,EEGMutSDREM,Wake24hMutSD,NREM24hMutSD,REM24hMutSD,Scores.(Sexes{1}).Mut.SD,EEG_bin_edgesMutSDM] = Extract_EEG_data(path_to_files_Mut_SD.(Sexes{1}),ffl_list_Mut_SD.(Sexes{1}),epochs,numhrs,swa,EEGLowerLimit_Hz); % Male Mut SleepDep
    [swaNREMMutBLhourly,swaREMMutBLhourly,swaWakeMutBLhourly,EEGMutBLWake,EEGMutBLNREM,EEGMutBLREM,Wake24hMutBL,NREM24hMutBL,REM24hMutBL,Scores.(Sexes{1}).Mut.BL,EEG_bin_edgesMutBLM] = Extract_EEG_data(path_to_files_Mut_BL.(Sexes{1}),ffl_list_Mut_BL.(Sexes{1}),epochs,numhrs,swa,EEGLowerLimit_Hz); % Male Mut BL
    [swaNREMWTSDhourly, swaREMWTSDhourly, swaWakeWTSDhourly, EEGWTSDWake, EEGWTSDNREM, EEGWTSDREM, Wake24hWTSD, NREM24hWTSD, REM24hWTSD,Scores.(Sexes{1}).WT.SD,EEG_bin_edgesWTSDM]    = Extract_EEG_data(path_to_files_WT_SD.(Sexes{1}),ffl_list_WT_SD.(Sexes{1}),epochs,numhrs,swa,EEGLowerLimit_Hz);   % Male WT SD
    [swaNREMWTBLhourly, swaREMWTBLhourly, swaWakeWTBLhourly, EEGWTBLWake, EEGWTBLNREM, EEGWTBLREM, Wake24hWTBL, NREM24hWTBL, REM24hWTBL,Scores.(Sexes{1}).WT.BL,EEG_bin_edgesWTBLM]    = Extract_EEG_data(path_to_files_WT_BL.(Sexes{1}),ffl_list_WT_BL.(Sexes{1}),epochs,numhrs,swa,EEGLowerLimit_Hz);   % Male WT BL

    % Females
    [swaNREMMutSDhourly_F,swaREMMutSDhourly_F,swaWakeMutSDhourly_F,EEGMutSDWake_F,EEGMutSDNREM_F,EEGMutSDREM_F,Wake24hMutSD_F,NREM24hMutSD_F,REM24hMutSD_F,Scores.(Sexes{2}).Mut.SD,EEG_bin_edgesMutSDF] = Extract_EEG_data(path_to_files_Mut_SD.(Sexes{2}),ffl_list_Mut_SD.(Sexes{2}),epochs,numhrs,swa,EEGLowerLimit_Hz); % Female Mut SleepDep
    [swaNREMMutBLhourly_F,swaREMMutBLhourly_F,swaWakeMutBLhourly_F,EEGMutBLWake_F,EEGMutBLNREM_F,EEGMutBLREM_F,Wake24hMutBL_F,NREM24hMutBL_F,REM24hMutBL_F,Scores.(Sexes{2}).Mut.BL,EEG_bin_edgesMutBLF] = Extract_EEG_data(path_to_files_Mut_BL.(Sexes{2}),ffl_list_Mut_BL.(Sexes{2}),epochs,numhrs,swa,EEGLowerLimit_Hz); % Female Mut BL
    [swaNREMWTSDhourly_F, swaREMWTSDhourly_F, swaWakeWTSDhourly_F, EEGWTSDWake_F, EEGWTSDNREM_F, EEGWTSDREM_F, Wake24hWTSD_F, NREM24hWTSD_F, REM24hWTSD_F,Scores.(Sexes{2}).WT.SD,EEG_bin_edgesWTSDF]    = Extract_EEG_data(path_to_files_WT_SD.(Sexes{2}),ffl_list_WT_SD.(Sexes{2}),epochs,numhrs,swa,EEGLowerLimit_Hz);  % Female WT SD
    [swaNREMWTBLhourly_F, swaREMWTBLhourly_F, swaWakeWTBLhourly_F, EEGWTBLWake_F, EEGWTBLNREM_F, EEGWTBLREM_F, Wake24hWTBL_F, NREM24hWTBL_F, REM24hWTBL_F,Scores.(Sexes{2}).WT.BL,EEG_bin_edgesWTBLF]    = Extract_EEG_data(path_to_files_WT_BL.(Sexes{2}),ffl_list_WT_BL.(Sexes{2}),epochs,numhrs,swa,EEGLowerLimit_Hz);  % Female WT BL
end

if ~load_data_from_mat_file_instead & isempty(ffl_list_WT_BL.(Sexes{1}))
    error('You did not request that I load data from a .mat file, and I can''t find any files in the directories. Do you need to connect to FS1?')
end 
disp('Done extracting EEG data.')
% ------------------------------------------------------------------------------------------------------
% ------------------------------------------------------------------------------------------------------
% ------------------------------------------------------------------------------------------------------
% --- END OF EXTRACTING DATA ---------------------------------------------------------------------------


% --- Sanity check #1:  calc sleep dep efficiency, to make sure no files were mixed up between BL and SD --
efficiency_struct = calc_sleep_dep_efficiency(Scores,Sexes,SD_length_hrs,epoch_duration_secs);
SD_efficiency_Table = table;
SD_efficiency_Table.Sex           = [repmat('M',length(efficiency_struct.M.WT) + length(efficiency_struct.M.Mut),1); repmat('F',length(efficiency_struct.F.WT) + length(efficiency_struct.F.Mut),1)]; 
SD_efficiency_Table.Genotype      = [repmat('WT ',length(efficiency_struct.M.WT),1); repmat('Mut',length(efficiency_struct.M.Mut),1); repmat('WT ',length(efficiency_struct.F.WT),1); repmat('Mut',length(efficiency_struct.F.Mut),1)]; 
SD_efficiency_Table.Filename      = [ffl_list_WT_SD.(Sexes{1})'; ffl_list_Mut_SD.(Sexes{1})'; ffl_list_WT_SD.(Sexes{2})'; ffl_list_Mut_SD.(Sexes{2})'];
SD_efficiency_Table.SD_Efficiency = [efficiency_struct.M.WT'; efficiency_struct.M.Mut'; efficiency_struct.F.WT'; efficiency_struct.F.Mut']; 
if sum(SD_efficiency_Table.SD_Efficiency<0.90)>0  % If you have any recordings with SD efficiency under 90% tell the user 
    answer = questdlg('At least one of the recordings had a sleep dep efficiency less than 90 percent. Would you like me to stop and tell you which file(s) have low sleep dep efficiency?','Sleep Dep Efficiency Issue','STOP. And tell me.','Keep Going.  It''s fine.','Keep Going.  It''s fine.');
    drawnow; pause(0.05);  % this prevents questdlg from hanging
    if strcmp(answer,'STOP. And tell me.')
        disp('Here are the problematic file(s):')
        disp(SD_efficiency_Table(SD_efficiency_Table.SD_Efficiency<0.90,:))
        error('Stopped because of low sleep dep efficiency.')
    end 
end 


% --- Sanity check #2:  What percentage of epochs are scored as artifact? 
[artefact_percentages_Table,artefact_perc_av,ex_code] = calc_artifact_percentages(Scores,Sexes);
if ex_code
    artefact_percentages_Table.Filename      = [ffl_list_WT_BL_Males'; ffl_list_WT_SD_Males'; ffl_list_Mut_BL_Males'; ffl_list_Mut_SD_Males'; ffl_list_WT_BL_Females'; ffl_list_WT_SD_Females'; ffl_list_Mut_BL_Females'; ffl_list_Mut_SD_Females'];
    display(artefact_percentages_Table)
    error('Stopped because of high levels of artifacts.')
end

%profile on
% ------------------------------------------------------------------------------------------------------
% --- Set up structures to pass into functions
% - Time in State - 
% Each row is an hour (24) and each column is for an animal.  24x5.  Values represent minutes per hour in that stage
TimeInState = SetUpTimeInStateStruct(Wake24hWTBL,   NREM24hWTBL,   REM24hWTBL,   Wake24hWTSD,   NREM24hWTSD,   REM24hWTSD, ...
                                     Wake24hMutBL,  NREM24hMutBL,  REM24hMutBL,  Wake24hMutSD,  NREM24hMutSD,  REM24hMutSD, ...
                                     Wake24hWTBL_F, NREM24hWTBL_F, REM24hWTBL_F, Wake24hWTSD_F, NREM24hWTSD_F, REM24hWTSD_F, ...
                                     Wake24hMutBL_F,NREM24hMutBL_F,REM24hMutBL_F,Wake24hMutSD_F,NREM24hMutSD_F,REM24hMutSD_F,SD_length_hrs,Sexes);




% Set up the Spectral Data struct EEG_struct
EEG_struct = SetUpEEGstruct(EEGWTBLWake,   EEGWTBLNREM,   EEGWTBLREM,   EEGWTSDWake,   EEGWTSDNREM,   EEGWTSDREM,...
                            EEGMutBLWake,  EEGMutBLNREM,  EEGMutBLREM,  EEGMutSDWake,  EEGMutSDNREM,  EEGMutSDREM, ...
                            EEGWTBLWake_F, EEGWTBLNREM_F, EEGWTBLREM_F, EEGWTSDWake_F, EEGWTSDNREM_F, EEGWTSDREM_F,...
                            EEGMutBLWake_F,EEGMutBLNREM_F,EEGMutBLREM_F,EEGMutSDWake_F,EEGMutSDNREM_F,EEGMutSDREM_F,Sexes);



% --- ANALYSIS -----------------------------------------------------------------------------------------
% ------------------------------------------------------------------------------------------------------
% ------------------------------------------------------------------------------------------------------


% ------------------------------------------------------------------------------------------------------
% ------- Time In State Analysis -----------------------------------------------------------------------
disp('Beginning Time in State Analysis...')
% NOTE: Rep Measures for the SD day first 12 hours ignores the actual sleep dep time
%       First12hrSD is actually the first 12 hours (LP) excluding the actual sleep dep 

% Set up struct
Twelve_hour_avg = struct();
Twelve_hour_avg.percentages.(Sexes{1})    = [];
Twelve_hour_avg.percentages.(Sexes{2})    = [];
Twelve_hour_avg.bout_counts.(Sexes{1})    = [];
Twelve_hour_avg.bout_counts.(Sexes{2})    = [];
Twelve_hour_avg.bout_durations.(Sexes{1}) = [];
Twelve_hour_avg.bout_durations.(Sexes{2}) = [];

% Males
[Twelve_hour_avg.percentages.(Sexes{1}),HourlyStatePercentages.(Sexes{1}),RepMeas_ANOVA_tables.(Sexes{1}).hourly_percentage,~]     = Time_In_State_Analysis(TimeInState=TimeInState,Sex=Sexes{1},LegendLabels=LegendLabels,SD_length_hrs=SD_length_hrs);
[~,                               ~,                        RepMeas_ANOVA_tables.(Sexes{1}).hourly_percentage.First12hrBL, rm_model_BL_LP] = Time_In_State_Analysis(TimeInState=TimeInState,TimeFrameRMANOVA=[1 12], inBLorSD='BL',Sex=Sexes{1},SD_length_hrs=SD_length_hrs);
[~,                               ~,                        RepMeas_ANOVA_tables.(Sexes{1}).hourly_percentage.Last12hrBL,  rm_model_BL_DP] = Time_In_State_Analysis(TimeInState=TimeInState,TimeFrameRMANOVA=[13 24],inBLorSD='BL',Sex=Sexes{1},SD_length_hrs=SD_length_hrs);
[~,                               ~,                        RepMeas_ANOVA_tables.(Sexes{1}).hourly_percentage.First12hrSD, rm_model_SD_LP] = Time_In_State_Analysis(TimeInState=TimeInState,TimeFrameRMANOVA=[SD_length_hrs+1 12], inBLorSD='SD',Sex=Sexes{1},SD_length_hrs=SD_length_hrs);
[~,                               ~,                        RepMeas_ANOVA_tables.(Sexes{1}).hourly_percentage.Last12hrSD,  rm_model_SD_DP] = Time_In_State_Analysis(TimeInState=TimeInState,TimeFrameRMANOVA=[13 24], inBLorSD='SD',Sex=Sexes{1},SD_length_hrs=SD_length_hrs);
[~,                               ~,                        RepMeas_ANOVA_tables.(Sexes{1}).hourly_percentage.First6hrDPBL,rm_mod_BL_DPf6] = Time_In_State_Analysis(TimeInState=TimeInState,TimeFrameRMANOVA=[13 18], inBLorSD='BL',Sex=Sexes{1},SD_length_hrs=SD_length_hrs); % these 4 lines are for 6-hour segments during DP
[~,                               ~,                        RepMeas_ANOVA_tables.(Sexes{1}).hourly_percentage.Last6hrDPBL, rm_mod_BL_DPl6] = Time_In_State_Analysis(TimeInState=TimeInState,TimeFrameRMANOVA=[19 24], inBLorSD='BL',Sex=Sexes{1},SD_length_hrs=SD_length_hrs);
[~,                               ~,                        RepMeas_ANOVA_tables.(Sexes{1}).hourly_percentage.First6hrDPSD,rm_mod_SD_DPf6] = Time_In_State_Analysis(TimeInState=TimeInState,TimeFrameRMANOVA=[13 18], inBLorSD='SD',Sex=Sexes{1},SD_length_hrs=SD_length_hrs);
[~,                               ~,                        RepMeas_ANOVA_tables.(Sexes{1}).hourly_percentage.Last6hrDPSD, rm_mod_SD_DPl6] = Time_In_State_Analysis(TimeInState=TimeInState,TimeFrameRMANOVA=[19 24], inBLorSD='SD',Sex=Sexes{1},SD_length_hrs=SD_length_hrs);


% % Females
[Twelve_hour_avg.percentages.(Sexes{2}),HourlyStatePercentages.(Sexes{2}),RepMeas_ANOVA_tables.(Sexes{2}).hourly_percentage,~]       = Time_In_State_Analysis(TimeInState=TimeInState,Sex=Sexes{2},MakeFigures=0,LegendLabels=LegendLabels,SD_length_hrs=SD_length_hrs);
[~,                               ~,                        RepMeas_ANOVA_tables.(Sexes{2}).hourly_percentage.First12hrBL,rm_model_F_BL_LP] = Time_In_State_Analysis(TimeInState=TimeInState,TimeFrameRMANOVA=[1 12], inBLorSD='BL',Sex=Sexes{2},SD_length_hrs=SD_length_hrs);
[~,                               ~,                        RepMeas_ANOVA_tables.(Sexes{2}).hourly_percentage.Last12hrBL, rm_model_F_BL_DP] = Time_In_State_Analysis(TimeInState=TimeInState,TimeFrameRMANOVA=[13 24],inBLorSD='BL',Sex=Sexes{2},SD_length_hrs=SD_length_hrs);
[~,                               ~,                        RepMeas_ANOVA_tables.(Sexes{2}).hourly_percentage.First12hrSD,rm_model_F_SD_LP] = Time_In_State_Analysis(TimeInState=TimeInState,TimeFrameRMANOVA=[SD_length_hrs+1 12], inBLorSD='SD',Sex=Sexes{2},SD_length_hrs=SD_length_hrs);
[~,                               ~,                        RepMeas_ANOVA_tables.(Sexes{2}).hourly_percentage.Last12hrSD, rm_model_F_SD_DP] = Time_In_State_Analysis(TimeInState=TimeInState,TimeFrameRMANOVA=[13 24], inBLorSD='SD',Sex=Sexes{2},SD_length_hrs=SD_length_hrs);


RepMeas_ANOVA_tables.(Sexes{1}).TIS.BL.First12hr  = RepMeas_ANOVA_tables.(Sexes{1}).hourly_percentage.First12hrBL.BL_custom;
RepMeas_ANOVA_tables.(Sexes{1}).TIS.BL.Last12hr   = RepMeas_ANOVA_tables.(Sexes{1}).hourly_percentage.Last12hrBL.BL_custom;
RepMeas_ANOVA_tables.(Sexes{1}).TIS.SD.First12hr  = RepMeas_ANOVA_tables.(Sexes{1}).hourly_percentage.First12hrSD.SD_custom;
RepMeas_ANOVA_tables.(Sexes{1}).TIS.SD.Last12hr   = RepMeas_ANOVA_tables.(Sexes{1}).hourly_percentage.Last12hrSD.SD_custom;
% if you want to separate DP into 6-hour segments
if Analyze_TIS_DP_6hr_segments
    % BL 
    RepMeas_ANOVA_tables.(Sexes{1}).TIS.BL.First6hrDP = RepMeas_ANOVA_tables.(Sexes{1}).hourly_percentage.First6hrDPBL.BL_custom;
    RepMeas_ANOVA_tables.(Sexes{1}).TIS.BL.Last6hrDP  = RepMeas_ANOVA_tables.(Sexes{1}).hourly_percentage.Last6hrDPBL.BL_custom;

    % SD
    RepMeas_ANOVA_tables.(Sexes{1}).TIS.SD.First6hrDP = RepMeas_ANOVA_tables.(Sexes{1}).hourly_percentage.First6hrDPSD.SD_custom;
    RepMeas_ANOVA_tables.(Sexes{1}).TIS.SD.Last6hrDP  = RepMeas_ANOVA_tables.(Sexes{1}).hourly_percentage.Last6hrDPSD.SD_custom;

end 

RM_Models.(Sexes{1}).BL.LP = rm_model_BL_LP;
RM_Models.(Sexes{1}).BL.DP = rm_model_BL_DP;
RM_Models.(Sexes{1}).SD.LP = rm_model_SD_LP;
RM_Models.(Sexes{1}).SD.DP = rm_model_SD_DP;
% if you want to separate DP into 6-hr segments
if Analyze_TIS_DP_6hr_segments
    % BL
    RM_Models.(Sexes{1}).BL.First6hrDP = rm_mod_BL_DPf6;
    RM_Models.(Sexes{1}).BL.Last6hrDP  = rm_mod_BL_DPl6;

    % SD
    RM_Models.(Sexes{1}).SD.First6hrDP = rm_mod_SD_DPf6;
    RM_Models.(Sexes{1}).SD.Last6hrDP  = rm_mod_SD_DPl6;
end  

% if ~isempty(ffl_list_WT_BL_Females)
%     RepMeas_ANOVA_tables_Females_TIS.BL.First12hr = RepMeas_ANOVA_tables_Females_hourly_percentageFirst12hrBL.BL_custom;
%     RepMeas_ANOVA_tables_Females_TIS.BL.Last12hr  = RepMeas_ANOVA_tables_Females_hourly_percentageLast12hrBL.BL_custom;
%     RepMeas_ANOVA_tables_Females_TIS.SD.First12hr = RepMeas_ANOVA_tables_Females_hourly_percentageFirst12hrSD.SD_custom;
%     RepMeas_ANOVA_tables_Females_TIS.SD.Last12hr  = RepMeas_ANOVA_tables_Females_hourly_percentageLast12hrSD.SD_custom;

%     RM_Models_F.BL.LP = rm_model_F_BL_LP;
%     RM_Models_F.BL.DP = rm_model_F_BL_DP;
%     RM_Models_F.SD.LP = rm_model_F_SD_LP;
%     RM_Models_F.SD.DP = rm_model_F_SD_DP;
% end 

if isempty(ffl_list_WT_BL.(Sexes{2}))  % If no female data, make the 3-panel figures for male data (using RepMeas_ANOVA_tables_Males_TIS)
    [RepMeas_ANOVA_TIS_p_values,PostHocs_TIS.(Sexes{1})] = plot_TimeInState(TimeInState=TimeInState,Sex='Male',  LegendLabels=LegendLabels,SD_length_hrs=SD_length_hrs,RepMeas_ANOVA_tables=RepMeas_ANOVA_tables.(Sexes{1}).TIS,RM_Models=RM_Models.(Sexes{1}));
end 

% if ~isempty(ffl_list_WT_BL_Females) 
%     PostHocs_TIS = plot_TimeInState(TimeInState=TimeInState,Sex='Female',LegendLabels=LegendLabels,SD_length_hrs=SD_length_hrs,RepMeas_ANOVA_tables=RepMeas_ANOVA_tables_Males_TIS,RM_Models=RM_Models_F);
% end 


% make a 9-panel figure: time in state panels, all 3 states, WT BL vs SD, Mut BL vs SD, WT SD vs Mut SD
[RepMeas_ANOVA_TimeInState_Recovery.(Sexes{1}), PostHocs_9_Panel_TimeInState.(Sexes{1})] = Make_9_Panel_TimeInState_Figure(HourlyStatePercentages.(Sexes{1}),TimeInState,SD_length_hrs,LegendLabels,Sexes{1});    % Male 
[RepMeas_ANOVA_TimeInState_Recovery.(Sexes{2}), PostHocs_9_Panel_TimeInState.(Sexes{2})] = Make_9_Panel_TimeInState_Figure(HourlyStatePercentages.(Sexes{2}),TimeInState,SD_length_hrs,LegendLabels,Sexes{2});  % Female 


TimeInStateResults = Build_Time_In_State_Results_Struct(RepMeas_ANOVA_tables.(Sexes{1}).hourly_percentage.First12hrBL,RepMeas_ANOVA_tables.(Sexes{1}).hourly_percentage.Last12hrBL, ...
                                                        RepMeas_ANOVA_tables.(Sexes{1}).hourly_percentage.First12hrSD,RepMeas_ANOVA_tables.(Sexes{1}).hourly_percentage.Last12hrSD, ...
                                                        RepMeas_ANOVA_tables.(Sexes{2}).hourly_percentage.First12hrBL,RepMeas_ANOVA_tables.(Sexes{2}).hourly_percentage.Last12hrBL, ...
                                                        RepMeas_ANOVA_tables.(Sexes{2}).hourly_percentage.First12hrSD,RepMeas_ANOVA_tables.(Sexes{2}).hourly_percentage.Last12hrSD);

% Interaction of Sex and Genotype
if ~all(structfun(@isempty,TimeInState.WT.(Sexes{2}).BL))  % If you have female data, do a RM anova with Sex and Genotype as factors
    %[RepMeasures_ANOVA_tables_BothSexes_hourly_percentage,~]  = Time_In_State_Analysis_Both_Genders(TimeInState=TimeInState);
    
    [RepMeas_ANOVA_tables_BothSexes_hourly_percentageFirst12hrBL, RM_Models_BothSexes.BL.LP] = Time_In_State_Analysis_Both_Genders(TimeInState=TimeInState,Sexes=Sexes,SexVarName=SexVarName,TimeFrameRMANOVA=[1 12],inBLorSD='BL');
    [RepMeas_ANOVA_tables_BothSexes_hourly_percentageLast12hrBL,  RM_Models_BothSexes.BL.DP] = Time_In_State_Analysis_Both_Genders(TimeInState=TimeInState,Sexes=Sexes,SexVarName=SexVarName,TimeFrameRMANOVA=[13 24],inBLorSD='BL');
    [RepMeas_ANOVA_tables_BothSexes_hourly_percentageFirst12hrSD, RM_Models_BothSexes.SD.LP] = Time_In_State_Analysis_Both_Genders(TimeInState=TimeInState,Sexes=Sexes,SexVarName=SexVarName,TimeFrameRMANOVA=[SD_length_hrs+1 12],inBLorSD='SD');
    [RepMeas_ANOVA_tables_BothSexes_hourly_percentageLast12hrSD,  RM_Models_BothSexes.SD.DP] = Time_In_State_Analysis_Both_Genders(TimeInState=TimeInState,Sexes=Sexes,SexVarName=SexVarName,TimeFrameRMANOVA=[13 24],inBLorSD='SD');


    RepMeas_ANOVA_tables_BothSexes_TIS.BL.First12hr = RepMeas_ANOVA_tables_BothSexes_hourly_percentageFirst12hrBL.BL_custom_interaction;
    RepMeas_ANOVA_tables_BothSexes_TIS.BL.Last12hr  = RepMeas_ANOVA_tables_BothSexes_hourly_percentageLast12hrBL.BL_custom_interaction;
    RepMeas_ANOVA_tables_BothSexes_TIS.SD.First12hr = RepMeas_ANOVA_tables_BothSexes_hourly_percentageFirst12hrSD.SD_custom_interaction;
    RepMeas_ANOVA_tables_BothSexes_TIS.SD.Last12hr  = RepMeas_ANOVA_tables_BothSexes_hourly_percentageLast12hrSD.SD_custom_interaction;


    [RepMeas_ANOVA_TIS_BothSexes_p_values,PostHocs_TIS_MandF] = plot_TimeInState_Both_Genders(TimeInState=TimeInState,Sexes=Sexes,SexVarName=SexVarName,LegendLabels=LegendLabels,SD_length_hrs=SD_length_hrs,RepMeas_ANOVA_tables=RepMeas_ANOVA_tables_BothSexes_TIS,RM_Models=RM_Models_BothSexes);


end 

% Perform 2-way ANOVA for time in each behavioral state (Sex and Genotype) using 12-hour averages
[TwoWayANOVAs_TimeInState_p_vals,TwoWayANOVA_TimeInState_sig_differences,ANOVA_tables_TimeInState] = perform_All_2way_ANOVAS_TISorBouts(Twelve_hour_avg.percentages.(Sexes{1}),Twelve_hour_avg.percentages.(Sexes{2}),12,Sexes,SexVarName);

if ~all(structfun(@isempty,TimeInState.WT.(Sexes{2}).BL))  % if female data present
    ANOVA_TimeInState_p_vals = RepMeas_ANOVA_TIS_BothSexes_p_values;
else 
    ANOVA_TimeInState_p_vals = RepMeas_ANOVA_TIS_p_values;
end

% Make Time in State plots grouped by sex and genotype and separated by light and dark
% BL
%W_TIS_BL_Light_Sex_Gen_fig_handle    = plot_time_in_state_M_vs_F_WT_vs_Mut(Twelve_hour_avg.percentages.Male,Twelve_hour_avg.percentages.Female,TwoWayANOVAs_TimeInState_p_vals,'BL','First12hrs','Wake');
W_TIS_BL_Light_Sex_Gen_fig_handle    = plot_output_var_M_vs_F_WT_vs_Mut(Twelve_hour_avg.percentages.(Sexes{1}),Twelve_hour_avg.percentages.(Sexes{2}),'Percentages',ANOVA_TimeInState_p_vals,'BL','First12hrs','Wake',Sexes,LegendLabels);
NREM_TIS_BL_Light_Sex_Gen_fig_handle = plot_output_var_M_vs_F_WT_vs_Mut(Twelve_hour_avg.percentages.(Sexes{1}),Twelve_hour_avg.percentages.(Sexes{2}),'Percentages',ANOVA_TimeInState_p_vals,'BL','First12hrs','NREM',Sexes,LegendLabels);
REM_TIS_BL_Light_Sex_Gen_fig_handle  = plot_output_var_M_vs_F_WT_vs_Mut(Twelve_hour_avg.percentages.(Sexes{1}),Twelve_hour_avg.percentages.(Sexes{2}),'Percentages',ANOVA_TimeInState_p_vals,'BL','First12hrs','REM',Sexes,LegendLabels);

W_TIS_BL_Dark_Sex_Gen_fig_handle    = plot_output_var_M_vs_F_WT_vs_Mut(Twelve_hour_avg.percentages.(Sexes{1}),Twelve_hour_avg.percentages.(Sexes{2}),'Percentages',ANOVA_TimeInState_p_vals,'BL','Last12hrs','Wake',Sexes,LegendLabels);
NREM_TIS_BL_Dark_Sex_Gen_fig_handle = plot_output_var_M_vs_F_WT_vs_Mut(Twelve_hour_avg.percentages.(Sexes{1}),Twelve_hour_avg.percentages.(Sexes{2}),'Percentages',ANOVA_TimeInState_p_vals,'BL','Last12hrs','NREM',Sexes,LegendLabels);
REM_TIS_BL_Dark_Sex_Gen_fig_handle  = plot_output_var_M_vs_F_WT_vs_Mut(Twelve_hour_avg.percentages.(Sexes{1}),Twelve_hour_avg.percentages.(Sexes{2}),'Percentages',ANOVA_TimeInState_p_vals,'BL','Last12hrs','REM',Sexes,LegendLabels);

% SD
W_TIS_SD_Light_Sex_Gen_fig_handle    = plot_output_var_M_vs_F_WT_vs_Mut(Twelve_hour_avg.percentages.(Sexes{1}),Twelve_hour_avg.percentages.(Sexes{2}),'Percentages',ANOVA_TimeInState_p_vals,'SDexcSD','First12hrs','Wake',Sexes,LegendLabels);
NREM_TIS_SD_Light_Sex_Gen_fig_handle = plot_output_var_M_vs_F_WT_vs_Mut(Twelve_hour_avg.percentages.(Sexes{1}),Twelve_hour_avg.percentages.(Sexes{2}),'Percentages',ANOVA_TimeInState_p_vals,'SDexcSD','First12hrs','NREM',Sexes,LegendLabels);
REM_TIS_SD_Light_Sex_Gen_fig_handle  = plot_output_var_M_vs_F_WT_vs_Mut(Twelve_hour_avg.percentages.(Sexes{1}),Twelve_hour_avg.percentages.(Sexes{2}),'Percentages',ANOVA_TimeInState_p_vals,'SDexcSD','First12hrs','REM',Sexes,LegendLabels);

W_TIS_SD_Dark_Sex_Gen_fig_handle    = plot_output_var_M_vs_F_WT_vs_Mut(Twelve_hour_avg.percentages.(Sexes{1}),Twelve_hour_avg.percentages.(Sexes{2}),'Percentages',ANOVA_TimeInState_p_vals,'SD','Last12hrs','Wake',Sexes,LegendLabels);
NREM_TIS_SD_Dark_Sex_Gen_fig_handle = plot_output_var_M_vs_F_WT_vs_Mut(Twelve_hour_avg.percentages.(Sexes{1}),Twelve_hour_avg.percentages.(Sexes{2}),'Percentages',ANOVA_TimeInState_p_vals,'SD','Last12hrs','NREM',Sexes,LegendLabels);
REM_TIS_SD_Dark_Sex_Gen_fig_handle  = plot_output_var_M_vs_F_WT_vs_Mut(Twelve_hour_avg.percentages.(Sexes{1}),Twelve_hour_avg.percentages.(Sexes{2}),'Percentages',ANOVA_TimeInState_p_vals,'SD','Last12hrs','REM',Sexes,LegendLabels);



%W_TIS_BL_Light_Sex_Gen_fig_handle    = plot_output_var_M_vs_F_WT_vs_Mut(Twelve_hour_avg.percentages.Male,Twelve_hour_avg.percentages.Female,[],TwoWayANOVAs_TimeInState_p_vals,'BL','First12hrs','Wake');

%W_TIS_SDexcSD_Light_Sex_Gen_fig_handle = plot_time_in_state_M_vs_F_WT_vs_Mut(Twelve_hour_avg_percentages_Male,Twelve_hour_avg_percentages_Female,'SDexcSD','First12hrs','Wake');


% Comment this line if you want to keep the individual panels as separate figures (rather then combined into one figure with 6 panels)
combine_time_in_state_panels_into_one_big_figure('Time In State')


disp('Done with Time in State Analysis.')
% ----------------------------------------------------------------------------------------------------------------------------------
% --- End of Time in State Analysis ------------------------------------------------------------------------------------------------

% ----------------------------------------------------------------------------------------------------------------------------------
% --- Bout Number Analysis   ---
% Males
[Twelve_hour_avg.bout_counts.(Sexes{1}),Twelve_hour_avg.bout_durations.(Sexes{1})] = Bout_Analysis(Scores=Scores,epoch_duration=epoch_duration_secs,sleep_dep_length=SD_length_hrs,Sex=Sexes{1},BoutMinimums=Bout_Minimums,WindowLengthHrs=12);


% Females
[Twelve_hour_avg.bout_counts.(Sexes{2}),Twelve_hour_avg.bout_durations.(Sexes{2})] = Bout_Analysis(Scores=Scores,epoch_duration=epoch_duration_secs,sleep_dep_length=SD_length_hrs,Sex=Sexes{2},BoutMinimums=Bout_Minimums,WindowLengthHrs=12);

if Analyze_Bouts_DP_6hr_segments
    % Males
    [Six_hour_avg.bout_counts.(Sexes{1}),Six_hour_avg.bout_durations.(Sexes{1})] = Bout_Analysis(Scores=Scores,epoch_duration=epoch_duration_secs,sleep_dep_length=SD_length_hrs,Sex=Sexes{1},BoutMinimums=Bout_Minimums,WindowLengthHrs=6);

    % Females
    [Six_hour_avg.bout_counts.(Sexes{2}),Six_hour_avg.bout_durations.(Sexes{2})] = Bout_Analysis(Scores=Scores,epoch_duration=epoch_duration_secs,sleep_dep_length=SD_length_hrs,Sex=Sexes{2},BoutMinimums=Bout_Minimums,WindowLengthHrs=6);

    % Perform 2-way ANOVA for bout counts of each behavioral state (Sex and Genotype)
    [TwoWayANOVAs_BoutCounts6hr_p_vals,TwoWayANOVA_BoutCounts6hr_sig_differences,ANOVA_tables_BoutNumber_6hr] = perform_All_2way_ANOVAS_TISorBouts(Six_hour_avg.bout_counts.(Sexes{1}),Six_hour_avg.bout_counts.(Sexes{2}),6,Sexes,SexVarName);

end  


% Perform 2-way ANOVA for bout counts of each behavioral state (Sex and Genotype)
[TwoWayANOVAs_BoutCounts_p_vals,TwoWayANOVA_BoutCounts_sig_differences,ANOVA_tables_BoutNumber] = perform_All_2way_ANOVAS_TISorBouts(Twelve_hour_avg.bout_counts.(Sexes{1}),Twelve_hour_avg.bout_counts.(Sexes{2}),12,Sexes,SexVarName);


% Make plots 
% W_boutnum_BL_Light_Sex_Gen_fig_handle    = plot_output_var_M_vs_F_WT_vs_Mut(Twelve_hour_avg.bout_counts.Male,Twelve_hour_avg.bout_counts.Female,'Bout_Counts',TwoWayANOVAs_BoutCounts_p_vals,'BL','First12hrs','Wake');
% BL Light period
W_boutnum_BL_Light_Sex_Gen_fig_handle    = plot_output_var_M_vs_F_WT_vs_Mut(Twelve_hour_avg.bout_counts.(Sexes{1}),Twelve_hour_avg.bout_counts.(Sexes{2}),'Bout_Counts',TwoWayANOVAs_BoutCounts_p_vals,'BL','First12hrs','Wake',Sexes,LegendLabels);
NREM_boutnum_BL_Light_Sex_Gen_fig_handle = plot_output_var_M_vs_F_WT_vs_Mut(Twelve_hour_avg.bout_counts.(Sexes{1}),Twelve_hour_avg.bout_counts.(Sexes{2}),'Bout_Counts',TwoWayANOVAs_BoutCounts_p_vals,'BL','First12hrs','NREM',Sexes,LegendLabels);
REM_boutnum_BL_Light_Sex_Gen_fig_handle  = plot_output_var_M_vs_F_WT_vs_Mut(Twelve_hour_avg.bout_counts.(Sexes{1}),Twelve_hour_avg.bout_counts.(Sexes{2}),'Bout_Counts',TwoWayANOVAs_BoutCounts_p_vals,'BL','First12hrs','REM',Sexes,LegendLabels);

% BL dark period
W_boutnum_BL_Dark_Sex_Gen_fig_handle    = plot_output_var_M_vs_F_WT_vs_Mut(Twelve_hour_avg.bout_counts.(Sexes{1}),Twelve_hour_avg.bout_counts.(Sexes{2}),'Bout_Counts',TwoWayANOVAs_BoutCounts_p_vals,'BL','Last12hrs','Wake',Sexes,LegendLabels);
NREM_boutnum_BL_Dark_Sex_Gen_fig_handle = plot_output_var_M_vs_F_WT_vs_Mut(Twelve_hour_avg.bout_counts.(Sexes{1}),Twelve_hour_avg.bout_counts.(Sexes{2}),'Bout_Counts',TwoWayANOVAs_BoutCounts_p_vals,'BL','Last12hrs','NREM',Sexes,LegendLabels);
REM_boutnum_BL_Dark_Sex_Gen_fig_handle  = plot_output_var_M_vs_F_WT_vs_Mut(Twelve_hour_avg.bout_counts.(Sexes{1}),Twelve_hour_avg.bout_counts.(Sexes{2}),'Bout_Counts',TwoWayANOVAs_BoutCounts_p_vals,'BL','Last12hrs','REM',Sexes,LegendLabels);

% SD (first 12 hours)
% W_boutnum_SD_Light_Sex_Gen_fig_handle    = plot_output_var_M_vs_F_WT_vs_Mut(Twelve_hour_avg.bout_counts.Male,Twelve_hour_avg.bout_counts.Female,'Bout_Counts',TwoWayANOVAs_BoutCounts_p_vals,'SD','First12hrs','Wake');
% NREM_boutnum_SD_Light_Sex_Gen_fig_handle = plot_output_var_M_vs_F_WT_vs_Mut(Twelve_hour_avg.bout_counts.Male,Twelve_hour_avg.bout_counts.Female,'Bout_Counts',TwoWayANOVAs_BoutCounts_p_vals,'SD','First12hrs','NREM');
% REM_boutnum_SD_Light_Sex_Gen_fig_handle  = plot_output_var_M_vs_F_WT_vs_Mut(Twelve_hour_avg.bout_counts.Male,Twelve_hour_avg.bout_counts.Female,'Bout_Counts',TwoWayANOVAs_BoutCounts_p_vals,'SD','First12hrs','REM');

% SDexcSD (first 12 hours)
W_boutnum_SD_Light_Sex_Gen_fig_handle    = plot_output_var_M_vs_F_WT_vs_Mut(Twelve_hour_avg.bout_counts.(Sexes{1}),Twelve_hour_avg.bout_counts.(Sexes{2}),'Bout_Counts',TwoWayANOVAs_BoutCounts_p_vals,'SDexcSD','First12hrs','Wake',Sexes,LegendLabels);
NREM_boutnum_SD_Light_Sex_Gen_fig_handle = plot_output_var_M_vs_F_WT_vs_Mut(Twelve_hour_avg.bout_counts.(Sexes{1}),Twelve_hour_avg.bout_counts.(Sexes{2}),'Bout_Counts',TwoWayANOVAs_BoutCounts_p_vals,'SDexcSD','First12hrs','NREM',Sexes,LegendLabels);
REM_boutnum_SD_Light_Sex_Gen_fig_handle  = plot_output_var_M_vs_F_WT_vs_Mut(Twelve_hour_avg.bout_counts.(Sexes{1}),Twelve_hour_avg.bout_counts.(Sexes{2}),'Bout_Counts',TwoWayANOVAs_BoutCounts_p_vals,'SDexcSD','First12hrs','REM',Sexes,LegendLabels);

% SD (last 12 hours)
W_boutnum_SD_Dark_Sex_Gen_fig_handle    = plot_output_var_M_vs_F_WT_vs_Mut(Twelve_hour_avg.bout_counts.(Sexes{1}),Twelve_hour_avg.bout_counts.(Sexes{2}),'Bout_Counts',TwoWayANOVAs_BoutCounts_p_vals,'SD','Last12hrs','Wake',Sexes,LegendLabels);
NREM_boutnum_SD_Dark_Sex_Gen_fig_handle = plot_output_var_M_vs_F_WT_vs_Mut(Twelve_hour_avg.bout_counts.(Sexes{1}),Twelve_hour_avg.bout_counts.(Sexes{2}),'Bout_Counts',TwoWayANOVAs_BoutCounts_p_vals,'SD','Last12hrs','NREM',Sexes,LegendLabels);
REM_boutnum_SD_Dark_Sex_Gen_fig_handle  = plot_output_var_M_vs_F_WT_vs_Mut(Twelve_hour_avg.bout_counts.(Sexes{1}),Twelve_hour_avg.bout_counts.(Sexes{2}),'Bout_Counts',TwoWayANOVAs_BoutCounts_p_vals,'SD','Last12hrs','REM',Sexes,LegendLabels);

combine_time_in_state_panels_into_one_big_figure('Bout Number')
                                                                                                      
if Analyze_Bouts_DP_6hr_segments
    
    % BL Dark period (fist 6 hours)
    W_boutnum_BL_Darkfirst6_Sex_Gen_fig_handle    = plot_output_var_M_vs_F_WT_vs_Mut(Six_hour_avg.bout_counts.(Sexes{1}),Six_hour_avg.bout_counts.(Sexes{2}),'Bout_Counts',TwoWayANOVAs_BoutCounts6hr_p_vals,'BL','DPfirst6','Wake',Sexes,LegendLabels);
    NREM_boutnum_BL_Darkfirst6_Sex_Gen_fig_handle = plot_output_var_M_vs_F_WT_vs_Mut(Six_hour_avg.bout_counts.(Sexes{1}),Six_hour_avg.bout_counts.(Sexes{2}),'Bout_Counts',TwoWayANOVAs_BoutCounts6hr_p_vals,'BL','DPfirst6','NREM',Sexes,LegendLabels);
    REM_boutnum_BL_Darkfirst6_Sex_Gen_fig_handle  = plot_output_var_M_vs_F_WT_vs_Mut(Six_hour_avg.bout_counts.(Sexes{1}),Six_hour_avg.bout_counts.(Sexes{2}),'Bout_Counts',TwoWayANOVAs_BoutCounts6hr_p_vals,'BL','DPfirst6','REM',Sexes,LegendLabels);

    % BL Dark period (last 6 hours)
    W_boutnum_BL_Dark_Sexlast6_Gen_fig_handle    = plot_output_var_M_vs_F_WT_vs_Mut(Six_hour_avg.bout_counts.(Sexes{1}),Six_hour_avg.bout_counts.(Sexes{2}),'Bout_Counts',TwoWayANOVAs_BoutCounts6hr_p_vals,'BL','DPlast6','Wake',Sexes,LegendLabels);
    NREM_boutnum_BL_Darklast6_Sex_Gen_fig_handle = plot_output_var_M_vs_F_WT_vs_Mut(Six_hour_avg.bout_counts.(Sexes{1}),Six_hour_avg.bout_counts.(Sexes{2}),'Bout_Counts',TwoWayANOVAs_BoutCounts6hr_p_vals,'BL','DPlast6','NREM',Sexes,LegendLabels);
    REM_boutnum_BL_Darklast6_Sex_Gen_fig_handle  = plot_output_var_M_vs_F_WT_vs_Mut(Six_hour_avg.bout_counts.(Sexes{1}),Six_hour_avg.bout_counts.(Sexes{2}),'Bout_Counts',TwoWayANOVAs_BoutCounts6hr_p_vals,'BL','DPlast6','REM',Sexes,LegendLabels);

    % give the figures better names
    W_boutnum_BL_Darkfirst6_Sex_Gen_fig_handle.Name    = 'Bout Number Wake BL DP First 6 hrs';
    NREM_boutnum_BL_Darkfirst6_Sex_Gen_fig_handle.Name = 'Bout Number NREM BL DP First 6 hrs';
    REM_boutnum_BL_Darkfirst6_Sex_Gen_fig_handle.Name  = 'Bout Number REM BL DP First 6 hrs';
    W_boutnum_BL_Dark_Sexlast6_Gen_fig_handle.Name     = 'Bout Number Wake BL DP Last 6 hrs';
    NREM_boutnum_BL_Darklast6_Sex_Gen_fig_handle.Name  = 'Bout Number NREM BL DP Last 6 hrs';
    REM_boutnum_BL_Darklast6_Sex_Gen_fig_handle.Name   = 'Bout Number REM BL DP Last 6 hrs';

    % SD Dark period (fist 6 hours)
    W_boutnum_SD_Darkfirst6_Sex_Gen_fig_handle    = plot_output_var_M_vs_F_WT_vs_Mut(Six_hour_avg.bout_counts.(Sexes{1}),Six_hour_avg.bout_counts.(Sexes{2}),'Bout_Counts',TwoWayANOVAs_BoutCounts6hr_p_vals,'SD','DPfirst6','Wake',Sexes,LegendLabels);
    NREM_boutnum_SD_Darkfirst6_Sex_Gen_fig_handle = plot_output_var_M_vs_F_WT_vs_Mut(Six_hour_avg.bout_counts.(Sexes{1}),Six_hour_avg.bout_counts.(Sexes{2}),'Bout_Counts',TwoWayANOVAs_BoutCounts6hr_p_vals,'SD','DPfirst6','NREM',Sexes,LegendLabels);
    REM_boutnum_SD_Darkfirst6_Sex_Gen_fig_handle  = plot_output_var_M_vs_F_WT_vs_Mut(Six_hour_avg.bout_counts.(Sexes{1}),Six_hour_avg.bout_counts.(Sexes{2}),'Bout_Counts',TwoWayANOVAs_BoutCounts6hr_p_vals,'SD','DPfirst6','REM',Sexes,LegendLabels);

    % SD Dark period (last 6 hours)
    W_boutnum_SD_Dark_Sexlast6_Gen_fig_handle    = plot_output_var_M_vs_F_WT_vs_Mut(Six_hour_avg.bout_counts.(Sexes{1}),Six_hour_avg.bout_counts.(Sexes{2}),'Bout_Counts',TwoWayANOVAs_BoutCounts6hr_p_vals,'SD','DPlast6','Wake',Sexes,LegendLabels);
    NREM_boutnum_SD_Darklast6_Sex_Gen_fig_handle = plot_output_var_M_vs_F_WT_vs_Mut(Six_hour_avg.bout_counts.(Sexes{1}),Six_hour_avg.bout_counts.(Sexes{2}),'Bout_Counts',TwoWayANOVAs_BoutCounts6hr_p_vals,'SD','DPlast6','NREM',Sexes,LegendLabels);
    REM_boutnum_SD_Darklast6_Sex_Gen_fig_handle  = plot_output_var_M_vs_F_WT_vs_Mut(Six_hour_avg.bout_counts.(Sexes{1}),Six_hour_avg.bout_counts.(Sexes{2}),'Bout_Counts',TwoWayANOVAs_BoutCounts6hr_p_vals,'SD','DPlast6','REM',Sexes,LegendLabels);

    % give the figures better names
    W_boutnum_SD_Darkfirst6_Sex_Gen_fig_handle.Name    = 'Bout Number Wake SD DP First 6 hrs';
    NREM_boutnum_SD_Darkfirst6_Sex_Gen_fig_handle.Name = 'Bout Number NREM SD DP First 6 hrs';
    REM_boutnum_SD_Darkfirst6_Sex_Gen_fig_handle.Name  = 'Bout Number REM SD DP First 6 hrs';
    W_boutnum_SD_Dark_Sexlast6_Gen_fig_handle.Name     = 'Bout Number Wake SD DP Last 6 hrs';
    NREM_boutnum_SD_Darklast6_Sex_Gen_fig_handle.Name  = 'Bout Number NREM SD DP Last 6 hrs';
    REM_boutnum_SD_Darklast6_Sex_Gen_fig_handle.Name   = 'Bout Number REM SD DP Last 6 hrs';

end 
% ----------------------------------------------------------------------------------------------------------------------------------
% --- End of Bout Number Analysis --------------------------------------------------------------------------------------------------



% --- Bout Duration Analysis ---------------------------------------------------------------------
% bout durations have already been computed above in Bout_Analysis.m

if Analyze_Bouts_DP_6hr_segments
    
    % Perform 2-way ANOVA for bout durations of each behavioral state (Sex and Genotype)
    [TwoWayANOVAs_BoutDurations6hr_p_vals,TwoWayANOVA_BoutDurations6hr_sig_differences,ANOVA_tables_BoutDuration_6hr] = perform_All_2way_ANOVAS_TISorBouts(Six_hour_avg.bout_durations.(Sexes{1}),Six_hour_avg.bout_durations.(Sexes{2}),6,Sexes,SexVarName);
end  




% Perform 2-way ANOVA for bout durations of each behavioral state (Sex and Genotype)
[TwoWayANOVAs_BoutDurations_p_vals,TwoWayANOVA_BoutDurations_sig_differences,ANOVA_tables_BoutDuration] = perform_All_2way_ANOVAS_TISorBouts(Twelve_hour_avg.bout_durations.(Sexes{1}),Twelve_hour_avg.bout_durations.(Sexes{2}),12,Sexes,SexVarName);



% Make plots 
% BL Light period
W_boutdur_BL_Light_Sex_Gen_fig_handle    = plot_output_var_M_vs_F_WT_vs_Mut(Twelve_hour_avg.bout_durations.(Sexes{1}),Twelve_hour_avg.bout_durations.(Sexes{2}),'Bout_Durations',TwoWayANOVAs_BoutDurations_p_vals,'BL','First12hrs','Wake',Sexes,LegendLabels);
NREM_boutdur_BL_Light_Sex_Gen_fig_handle = plot_output_var_M_vs_F_WT_vs_Mut(Twelve_hour_avg.bout_durations.(Sexes{1}),Twelve_hour_avg.bout_durations.(Sexes{2}),'Bout_Durations',TwoWayANOVAs_BoutDurations_p_vals,'BL','First12hrs','NREM',Sexes,LegendLabels);
REM_boutdur_BL_Light_Sex_Gen_fig_handle  = plot_output_var_M_vs_F_WT_vs_Mut(Twelve_hour_avg.bout_durations.(Sexes{1}),Twelve_hour_avg.bout_durations.(Sexes{2}),'Bout_Durations',TwoWayANOVAs_BoutDurations_p_vals,'BL','First12hrs','REM',Sexes,LegendLabels);

% BL dark period
W_boutdur_BL_Dark_Sex_Gen_fig_handle    = plot_output_var_M_vs_F_WT_vs_Mut(Twelve_hour_avg.bout_durations.(Sexes{1}),Twelve_hour_avg.bout_durations.(Sexes{2}),'Bout_Durations',TwoWayANOVAs_BoutDurations_p_vals,'BL','Last12hrs','Wake',Sexes,LegendLabels);
NREM_boutdur_BL_Dark_Sex_Gen_fig_handle = plot_output_var_M_vs_F_WT_vs_Mut(Twelve_hour_avg.bout_durations.(Sexes{1}),Twelve_hour_avg.bout_durations.(Sexes{2}),'Bout_Durations',TwoWayANOVAs_BoutDurations_p_vals,'BL','Last12hrs','NREM',Sexes,LegendLabels);
REM_boutdur_BL_Dark_Sex_Gen_fig_handle  = plot_output_var_M_vs_F_WT_vs_Mut(Twelve_hour_avg.bout_durations.(Sexes{1}),Twelve_hour_avg.bout_durations.(Sexes{2}),'Bout_Durations',TwoWayANOVAs_BoutDurations_p_vals,'BL','Last12hrs','REM',Sexes,LegendLabels);

% SDexcSD (first 12 hours)
W_boutdur_SD_Light_Sex_Gen_fig_handle    = plot_output_var_M_vs_F_WT_vs_Mut(Twelve_hour_avg.bout_durations.(Sexes{1}),Twelve_hour_avg.bout_durations.(Sexes{2}),'Bout_Durations',TwoWayANOVAs_BoutDurations_p_vals,'SDexcSD','First12hrs','Wake',Sexes,LegendLabels);
NREM_boutdur_SD_Light_Sex_Gen_fig_handle = plot_output_var_M_vs_F_WT_vs_Mut(Twelve_hour_avg.bout_durations.(Sexes{1}),Twelve_hour_avg.bout_durations.(Sexes{2}),'Bout_Durations',TwoWayANOVAs_BoutDurations_p_vals,'SDexcSD','First12hrs','NREM',Sexes,LegendLabels);
REM_boutdur_SD_Light_Sex_Gen_fig_handle  = plot_output_var_M_vs_F_WT_vs_Mut(Twelve_hour_avg.bout_durations.(Sexes{1}),Twelve_hour_avg.bout_durations.(Sexes{2}),'Bout_Durations',TwoWayANOVAs_BoutDurations_p_vals,'SDexcSD','First12hrs','REM',Sexes,LegendLabels);

% SDexcSD (last 12 hours)
W_boutdur_SD_Dark_Sex_Gen_fig_handle    = plot_output_var_M_vs_F_WT_vs_Mut(Twelve_hour_avg.bout_durations.(Sexes{1}),Twelve_hour_avg.bout_durations.(Sexes{2}),'Bout_Durations',TwoWayANOVAs_BoutDurations_p_vals,'SD','Last12hrs','Wake',Sexes,LegendLabels);
NREM_boutdur_SD_Dark_Sex_Gen_fig_handle = plot_output_var_M_vs_F_WT_vs_Mut(Twelve_hour_avg.bout_durations.(Sexes{1}),Twelve_hour_avg.bout_durations.(Sexes{2}),'Bout_Durations',TwoWayANOVAs_BoutDurations_p_vals,'SD','Last12hrs','NREM',Sexes,LegendLabels);
REM_boutdur_SD_Dark_Sex_Gen_fig_handle  = plot_output_var_M_vs_F_WT_vs_Mut(Twelve_hour_avg.bout_durations.(Sexes{1}),Twelve_hour_avg.bout_durations.(Sexes{2}),'Bout_Durations',TwoWayANOVAs_BoutDurations_p_vals,'SD','Last12hrs','REM',Sexes,LegendLabels);

combine_time_in_state_panels_into_one_big_figure('Bout Duration')

if Analyze_Bouts_DP_6hr_segments
    
    % BL Dark period (fist 6 hours)
    W_boutdur_BL_Darkfirst6_Sex_Gen_fig_handle    = plot_output_var_M_vs_F_WT_vs_Mut(Six_hour_avg.bout_durations.(Sexes{1}),Six_hour_avg.bout_durations.(Sexes{2}),'Bout_Durations',TwoWayANOVAs_BoutDurations6hr_p_vals,'BL','DPfirst6','Wake',Sexes,LegendLabels);
    NREM_boutdur_BL_Darkfirst6_Sex_Gen_fig_handle = plot_output_var_M_vs_F_WT_vs_Mut(Six_hour_avg.bout_durations.(Sexes{1}),Six_hour_avg.bout_durations.(Sexes{2}),'Bout_Durations',TwoWayANOVAs_BoutDurations6hr_p_vals,'BL','DPfirst6','NREM',Sexes,LegendLabels);
    REM_boutdur_BL_Darkfirst6_Sex_Gen_fig_handle  = plot_output_var_M_vs_F_WT_vs_Mut(Six_hour_avg.bout_durations.(Sexes{1}),Six_hour_avg.bout_durations.(Sexes{2}),'Bout_Durations',TwoWayANOVAs_BoutDurations6hr_p_vals,'BL','DPfirst6','REM',Sexes,LegendLabels);

    % BL Dark period (last 6 hours)
    W_boutdur_BL_Dark_Sexlast6_Gen_fig_handle    = plot_output_var_M_vs_F_WT_vs_Mut(Six_hour_avg.bout_durations.(Sexes{1}),Six_hour_avg.bout_durations.(Sexes{2}),'Bout_Durations',TwoWayANOVAs_BoutDurations6hr_p_vals,'BL','DPlast6','Wake',Sexes,LegendLabels);
    NREM_boutdur_BL_Darklast6_Sex_Gen_fig_handle = plot_output_var_M_vs_F_WT_vs_Mut(Six_hour_avg.bout_durations.(Sexes{1}),Six_hour_avg.bout_durations.(Sexes{2}),'Bout_Durations',TwoWayANOVAs_BoutDurations6hr_p_vals,'BL','DPlast6','NREM',Sexes,LegendLabels);
    REM_boutdur_BL_Darklast6_Sex_Gen_fig_handle  = plot_output_var_M_vs_F_WT_vs_Mut(Six_hour_avg.bout_durations.(Sexes{1}),Six_hour_avg.bout_durations.(Sexes{2}),'Bout_Durations',TwoWayANOVAs_BoutDurations6hr_p_vals,'BL','DPlast6','REM',Sexes,LegendLabels);

    % give the figures better names
    W_boutdur_BL_Darkfirst6_Sex_Gen_fig_handle.Name    = 'Bout Duration Wake BL DP First 6 hrs';
    NREM_boutdur_BL_Darkfirst6_Sex_Gen_fig_handle.Name = 'Bout Duration NREM BL DP First 6 hrs';
    REM_boutdur_BL_Darkfirst6_Sex_Gen_fig_handle.Name  = 'Bout Duration REM BL DP First 6 hrs';
    W_boutdur_BL_Dark_Sexlast6_Gen_fig_handle.Name     = 'Bout Duration Wake BL DP Last 6 hrs';
    NREM_boutdur_BL_Darklast6_Sex_Gen_fig_handle.Name  = 'Bout Duration NREM BL DP Last 6 hrs';
    REM_boutdur_BL_Darklast6_Sex_Gen_fig_handle.Name   = 'Bout Duration REM BL DP Last 6 hrs';

    % SD Dark period (fist 6 hours)
    W_boutdur_SD_Darkfirst6_Sex_Gen_fig_handle    = plot_output_var_M_vs_F_WT_vs_Mut(Six_hour_avg.bout_durations.(Sexes{1}),Six_hour_avg.bout_durations.(Sexes{2}),'Bout_Durations',TwoWayANOVAs_BoutDurations6hr_p_vals,'SD','DPfirst6','Wake',Sexes,LegendLabels);
    NREM_boutdur_SD_Darkfirst6_Sex_Gen_fig_handle = plot_output_var_M_vs_F_WT_vs_Mut(Six_hour_avg.bout_durations.(Sexes{1}),Six_hour_avg.bout_durations.(Sexes{2}),'Bout_Durations',TwoWayANOVAs_BoutDurations6hr_p_vals,'SD','DPfirst6','NREM',Sexes,LegendLabels);
    REM_boutdur_SD_Darkfirst6_Sex_Gen_fig_handle  = plot_output_var_M_vs_F_WT_vs_Mut(Six_hour_avg.bout_durations.(Sexes{1}),Six_hour_avg.bout_durations.(Sexes{2}),'Bout_Durations',TwoWayANOVAs_BoutDurations6hr_p_vals,'SD','DPfirst6','REM',Sexes,LegendLabels);

    % SD Dark period (last 6 hours)
    W_boutdur_SD_Dark_Sexlast6_Gen_fig_handle    = plot_output_var_M_vs_F_WT_vs_Mut(Six_hour_avg.bout_durations.(Sexes{1}),Six_hour_avg.bout_durations.(Sexes{2}),'Bout_Durations',TwoWayANOVAs_BoutDurations6hr_p_vals,'SD','DPlast6','Wake',Sexes,LegendLabels);
    NREM_boutdur_SD_Darklast6_Sex_Gen_fig_handle = plot_output_var_M_vs_F_WT_vs_Mut(Six_hour_avg.bout_durations.(Sexes{1}),Six_hour_avg.bout_durations.(Sexes{2}),'Bout_Durations',TwoWayANOVAs_BoutDurations6hr_p_vals,'SD','DPlast6','NREM',Sexes,LegendLabels);
    REM_boutdur_SD_Darklast6_Sex_Gen_fig_handle  = plot_output_var_M_vs_F_WT_vs_Mut(Six_hour_avg.bout_durations.(Sexes{1}),Six_hour_avg.bout_durations.(Sexes{2}),'Bout_Durations',TwoWayANOVAs_BoutDurations6hr_p_vals,'SD','DPlast6','REM',Sexes,LegendLabels);

    % give the figures better names
    W_boutdur_SD_Darkfirst6_Sex_Gen_fig_handle.Name    = 'Bout Duration Wake SD DP First 6 hrs';
    NREM_boutdur_SD_Darkfirst6_Sex_Gen_fig_handle.Name = 'Bout Duration NREM SD DP First 6 hrs';
    REM_boutdur_SD_Darkfirst6_Sex_Gen_fig_handle.Name  = 'Bout Duration REM SD DP First 6 hrs';
    W_boutdur_SD_Dark_Sexlast6_Gen_fig_handle.Name     = 'Bout Duration Wake SD DP Last 6 hrs';
    NREM_boutdur_SD_Darklast6_Sex_Gen_fig_handle.Name  = 'Bout Duration NREM SD DP Last 6 hrs';
    REM_boutdur_SD_Darklast6_Sex_Gen_fig_handle.Name   = 'Bout Duration REM SD DP Last 6 hrs';

end                   

% ----------------------------------------------------------------------------------------------------------------------------------
% --- End of Bout Duration Analysis ------------------------------------------------------------------------------------------------


% ----------------------------------------------------------------------------------------------------------------------------------
% --- Spectral Analysis        -----------------------------------------------------------------------------------------------------
disp('Beginning spectral analysis...')
allEEG_M_bins_equal = isequal(length(EEG_bin_edgesWTBLM),length(EEG_bin_edgesWTSDM),length(EEG_bin_edgesMutBLM),length(EEG_bin_edgesMutSDM));
allEEG_F_bins_equal = isequal(length(EEG_bin_edgesWTBLF),length(EEG_bin_edgesWTSDF),length(EEG_bin_edgesMutBLF),length(EEG_bin_edgesMutSDF));
if ~allEEG_M_bins_equal | ~allEEG_F_bins_equal
    error('Just before call to spectral_analysis: The EEG bins are not the same for some combo of WT/Mut and BL/SD')
else 
    EEG_bin_edges = EEG_bin_edgesWTBLM;
end 

[spectral_pvals] = spectral_analysis(EEG_struct=EEG_struct,SD_length_hrs=SD_length_hrs,LegendLabels=LegendLabels, ...
                                     EEG_bin_edges=EEG_bin_edges,EEGLowerLimit_Hz=EEGLowerLimit_Hz,Normalization=Normalization,...
                                     SeparateSpectralIntoLPDP=SeparateSpectralIntoLPDP,Sexes=Sexes);

spectrogram = Make_Spectrogram(BL_file=ffl_list_Mut_BL.(Sexes{1}){2},SD_file=ffl_list_Mut_SD.(Sexes{1}){2},path_to_BL_file=path_to_files_Mut_BL.(Sexes{1}),...
                               path_to_SD_file = path_to_files_Mut_SD.(Sexes{1}),epoch_dur=epoch_duration_secs,window_length_min=10,hypnogram_pooling_window_min=2);

disp('Done with spectral analysis.')
% ----------------------------------------------------------------------------------------------------------------------------------
% --- End of Spectral Analysis -----------------------------------------------------------------------------------------------------


% ----------------------------------------------------------------------------------------------------------------------------------
% --- Sleep latency analysis ------------------------------------------------------------------------------------------------------- 

all_scores_table = build_scores_table(Scores,Sexes); 
[latency_struct,latency_fig_handles,latency_p_values,latency_sig_difference_bool,ANOVA_table_SleepLatency] = sleep_latency_analysis(all_scores_table,SD_length_hrs,epoch_duration_secs,firstNREM_episode_duration_epochs,NREM_char,LegendLabels,Sexes,SexVarName);
% ----------------------------------------------------------------------------------------------------------------------------------
% --- End of Sleep latency   -------------------------------------------------------------------------------------------------------


% --- NREM Delta power after SD (compare between genotypes for both Male and Female) ----
% ---------------------------------------------------------------------------------------
% NOTE: delta is defined as 0.9-3.9 Hz in Extract_EEG_data.m.  This is different than what Lizzy does in poster (0.8-3.9 Hz)
Hourly_Delta_inNREM.WT.(Sexes{1}).BL  = swaNREMWTBLhourly;
Hourly_Delta_inNREM.Mut.(Sexes{1}).BL = swaNREMMutBLhourly;
Hourly_Delta_inNREM.WT.(Sexes{2}).BL  = swaNREMWTBLhourly_F;
Hourly_Delta_inNREM.Mut.(Sexes{2}).BL = swaNREMMutBLhourly_F;

Hourly_Delta_inNREM.WT.(Sexes{1}).SD  = swaNREMWTSDhourly;
Hourly_Delta_inNREM.Mut.(Sexes{1}).SD = swaNREMMutSDhourly;
Hourly_Delta_inNREM.WT.(Sexes{2}).SD  = swaNREMWTSDhourly_F;
Hourly_Delta_inNREM.Mut.(Sexes{2}).SD = swaNREMMutSDhourly_F;

% Hourly_Sigma_inWAKE.WT.Male.BL = ...  % Fill these in then call the function with Hourly_Sigma_inWake

% -- 1-hr bins ---
[RepMeas_ANOVA_tables_SleepPressure,PostHocs_SleepPressure] = hourly_NREM_Delta_or_Sigma_afterSD(Hourly_power_inSTATE=Hourly_Delta_inNREM,DeltaOrTheta='Delta',WorNREM='NREM', BinSizeHrs = 1, Sexes=Sexes, Plot_BL_too = PlotHourlyBaselineNREM_Delta, SD_length_hrs= SD_length_hrs, ...
                                                                ffl_list_WT_BL_Males   = ffl_list_WT_BL.(Sexes{1}), ffl_list_Mut_BL_Males   = ffl_list_Mut_BL.(Sexes{1}), ...
                                                                ffl_list_WT_BL_Females = ffl_list_WT_BL.(Sexes{2}), ffl_list_Mut_BL_Females = ffl_list_Mut_BL.(Sexes{2}), ...
                                                                ffl_list_WT_SD_Males   = ffl_list_WT_SD.(Sexes{1}), ffl_list_Mut_SD_Males   = ffl_list_Mut_SD.(Sexes{1}), ...
                                                                ffl_list_WT_SD_Females = ffl_list_WT_SD.(Sexes{2}), ffl_list_Mut_SD_Females = ffl_list_Mut_SD.(Sexes{2}));

% -- 2-hr bins ---
if Analyze_Recovery_2hr_bins
    [RepMeas_ANOVA_tbls_SleepPress_2hrBins,PostHocs_SleepPress_2hrBins] = hourly_NREM_Delta_or_Sigma_afterSD(Hourly_power_inSTATE=Hourly_Delta_inNREM,DeltaOrTheta='Delta',WorNREM='NREM', BinSizeHrs = 2, Sexes=Sexes,Plot_BL_too = PlotHourlyBaselineNREM_Delta, SD_length_hrs=SD_length_hrs, ...
                                                                            ffl_list_WT_BL_Males   = ffl_list_WT_BL.(Sexes{1}), ffl_list_Mut_BL_Males   = ffl_list_Mut_BL.(Sexes{1}), ...
                                                                            ffl_list_WT_BL_Females = ffl_list_WT_BL.(Sexes{2}), ffl_list_Mut_BL_Females = ffl_list_Mut_BL.(Sexes{2}), ...
                                                                            ffl_list_WT_SD_Males   = ffl_list_WT_SD.(Sexes{1}), ffl_list_Mut_SD_Males   = ffl_list_Mut_SD.(Sexes{1}), ...
                                                                            ffl_list_WT_SD_Females = ffl_list_WT_SD.(Sexes{2}), ffl_list_Mut_SD_Females = ffl_list_Mut_SD.(Sexes{2}));
end 
%profile viewer
