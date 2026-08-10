function S = InputParams_p24_p30_for_spectra


% -- These are the parameters that used to be set in the main script -----------------
% -- These are only the files that were used based on the excel document
% for eeg analysis in lizzy's paper ----
% ------------------------------------------------------------------------------------
S.epochs = 900;               % number of epochs for the spectra bins, 900=1h, 1800=2h
S.numhrs = 24;
S.Light  = ([1:10800]);       % first 12h
S.Dark   = ([10801:21600]);   % second12h ****
S.swa    = ([900,24]);          
S.WakeEpocsDWTBL = [];        
S.NREMEpocsDWTBL = [];
S.REMEpocsDWTBL  = [];

S.SD_length_hrs  = 3;						% The length of the sleep dep in hours
S.epoch_duration_secs = 4;
S.firstNREM_episode_duration_epochs = 7;  	% the number of NREM epochs to count as a NREM bout for calculating sleep latency
S.EEGLowerLimit_Hz = 0; 				% the lower limit of frequencies read in from .mat files
S.NREM_char = {'N','NR','N*'};            	% the character used to label NREM sleep in sleep latency analysis  
S.Bout_Minimums.W = 2;        				% You need at least 2 epochs in a row to be counted as a bout.
S.Bout_Minimums.N = 2;
S.Bout_Minimums.R = 2;
S.LegendLabels = {'WildType','Shank3^{\DeltaC}'}; % {'WildType','Shank3^{\DeltaC}'} or {'WT','Mut'} or whatever
S.Sexes    	   = {'P24','P30'};    % {'Male','Female'} or {'P24','P30'} This only changes labels.  For if you want to analyze two ages instead of two sexes, for instance.
								   % Please list both 'Male' and 'Female' even if you only have Male data.   
S.SexVarName = 'Age';  			   % Make this 'Sex' if you have male and female data or only male data

S.Analyze_TIS_DP_6hr_segments = true; 	% true if you want to add an additional analysis of the time in state data: 
										% the first 6 hours of the dark period and the last 6 hours of the dark period.   

S.Normalization = 'MeanPowerAllStates';	% How to normalize the spectral curves.  Options are 'MeanPowerAllStates' or 'AreaUnderCurve'
S.PlotHourlyBaselineNREM_Delta = false;	% By default the code plots hourly NREM Delta after SD.  Do you also want to see that during baseline? 

S.SeparateSpectralIntoLPDP = true;		% Do you want to separate the spectral analysis into light period and dark period?   	

S.Analyze_Recovery_2hr_bins = false; 	% In addition to analyzing recovery sleep in 1-hr bins, do you also want to do it using 2-hr bins?  

% Paths to data
%p24 male p30 female
S.path_to_files_WT_BL.(S.Sexes{1})    = '/Users/caitlinottaway/Library/CloudStorage/OneDrive-WashingtonStateUniversity(email.wsu.edu)/Pexioto/EEG_Analysis/P24_P30_P45_P60/P24_MAT_forspectra/WT_BL/';
S.path_to_files_WT_SD.(S.Sexes{1})    = '/Users/caitlinottaway/Library/CloudStorage/OneDrive-WashingtonStateUniversity(email.wsu.edu)/Pexioto/EEG_Analysis/P24_P30_P45_P60/P24_MAT_forspectra/WT_SD/';
S.path_to_files_Mut_BL.(S.Sexes{1})   = '/Users/caitlinottaway/Library/CloudStorage/OneDrive-WashingtonStateUniversity(email.wsu.edu)/Pexioto/EEG_Analysis/P24_P30_P45_P60/P24_MAT_forspectra/S3_BL/';
S.path_to_files_Mut_SD.(S.Sexes{1})   = '/Users/caitlinottaway/Library/CloudStorage/OneDrive-WashingtonStateUniversity(email.wsu.edu)/Pexioto/EEG_Analysis/P24_P30_P45_P60/P24_MAT_forspectra/S3_SD/';
S.path_to_files_WT_BL.(S.Sexes{2})    = '/Users/caitlinottaway/Library/CloudStorage/OneDrive-WashingtonStateUniversity(email.wsu.edu)/Pexioto/EEG_Analysis/P24_P30_P45_P60/P30_MAT_forspectra/WT_BL/';
S.path_to_files_WT_SD.(S.Sexes{2})    = '/Users/caitlinottaway/Library/CloudStorage/OneDrive-WashingtonStateUniversity(email.wsu.edu)/Pexioto/EEG_Analysis/P24_P30_P45_P60/P30_MAT_forspectra/WT_SD/';
S.path_to_files_Mut_BL.(S.Sexes{2})   = '/Users/caitlinottaway/Library/CloudStorage/OneDrive-WashingtonStateUniversity(email.wsu.edu)/Pexioto/EEG_Analysis/P24_P30_P45_P60/P30_MAT_forspectra/S3_BL/';
S.path_to_files_Mut_SD.(S.Sexes{2})   = '/Users/caitlinottaway/Library/CloudStorage/OneDrive-WashingtonStateUniversity(email.wsu.edu)/Pexioto/EEG_Analysis/P24_P30_P45_P60/P30_MAT_forspectra/S3_SD/';
% -- IMPORTANT:  Do you want to read data from a .mat file instead of the directories above? (to save time)
S.load_data_from_mat_file_instead = 0; % 0 for no, 1 for yes
S.MatFileContainingData = 'P24_P30.mat';



% You don't need to mess with this
S.FileName = mfilename;
st = dbstack;
S.FuncName = st.name; 
