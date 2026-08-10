function [p_vals,PostHoc_tables] = plot_TimeInState(options)
%
% USAGE: plot_TimeInState(TimeInState=TimeInState,Sex='Male',LegendLabels=LegendLabels,SD_length_hrs=SD_length_hrs,RepMeas_ANOVA_tables=RepMeas_ANOVA_tables_Males_TIS
% 
% This function simply plots the time in state plots (3 panels, one for each state) and adds asterisks if the repeated measures 
% ANOVAs are significant.  Post-hocs?  




arguments
    options.TimeInState
    %options.SD_length_hrs    = 5
    options.SD_length_hrs    = S.SD_length_hrs
    options.Sex              = 'Male'
    options.LegendLabels     = []
    options.RepMeas_ANOVA_tables = []
    options.RM_Models            = []
    options.PlotIndividualFigs   = 0
end 

TimeInState_struct 	 = options.TimeInState;
SD_length_hrs      	 = options.SD_length_hrs;
Sex                	 = options.Sex;
LegendLabels       	 = options.LegendLabels;
RepMeas_ANOVA_tables = options.RepMeas_ANOVA_tables;
RM_Models            = options.RM_Models;
PlotIndividualFigs   = options.PlotIndividualFigs;

% If Female data was requested, and you don't have female data, just stop
if strcmp(Sex,'Female') & isempty(TimeInState_struct.WT.Female.BL.Wake)
	PostHoc_tables = [];
	return
end 


% -- Unpack the TimeInState struct so we can use the plotting code below in a less messy way

% These are MxN matrices.  M rows for M hours, N columns for N animals.  The values are the number of minutes per hour spent 
% in that state by that animal.   
Wake24hWTBL      = TimeInState_struct.WT.(Sex).BL.Wake;
Wake24hWTSD      = TimeInState_struct.WT.(Sex).SD.Wake;
Wake24hWTSDexcSD = TimeInState_struct.WT.(Sex).SDexcSD.Wake;
NREM24hWTBL      = TimeInState_struct.WT.(Sex).BL.NREM;
NREM24hWTSD      = TimeInState_struct.WT.(Sex).SD.NREM;
NREM24hWTSDexcSD = TimeInState_struct.WT.(Sex).SDexcSD.NREM;
REM24hWTBL       = TimeInState_struct.WT.(Sex).BL.REM;
REM24hWTSD       = TimeInState_struct.WT.(Sex).SD.REM;
REM24hWTSDexcSD  = TimeInState_struct.WT.(Sex).SDexcSD.REM;


Wake24hMutBL      = TimeInState_struct.Mut.(Sex).BL.Wake;
Wake24hMutSD      = TimeInState_struct.Mut.(Sex).SD.Wake;
Wake24hMutSDexcSD = TimeInState_struct.Mut.(Sex).SDexcSD.Wake;
NREM24hMutBL      = TimeInState_struct.Mut.(Sex).BL.NREM;
NREM24hMutSD      = TimeInState_struct.Mut.(Sex).SD.NREM;
NREM24hMutSDexcSD = TimeInState_struct.Mut.(Sex).SDexcSD.NREM; 
REM24hMutBL       = TimeInState_struct.Mut.(Sex).BL.REM;
REM24hMutSD       = TimeInState_struct.Mut.(Sex).SD.REM;
REM24hMutSDexcSD  = TimeInState_struct.Mut.(Sex).SDexcSD.REM;

%% time in stage per hour normalized to recording time (aka, relative to 60min
WakeWTBLRelmean = mean(Wake24hWTBL./60*100,2,'omitnan')';
WakeWTBLRelstd  = std((Wake24hWTBL./60*100)','omitnan')/sqrt(size(Wake24hWTBL,2));
NREMWTBLRelmean = mean(NREM24hWTBL./60*100,2,'omitnan')';
NREMWTBLRelstd  = std((NREM24hWTBL./60*100)','omitnan')/sqrt(size(NREM24hWTBL,2));
REMWTBLRelmean  = mean(REM24hWTBL./60*100,2,'omitnan')';
REMWTBLRelstd   = std((REM24hWTBL./60*100)','omitnan')/sqrt(size(REM24hWTBL,2));

% --- Mutant Baseline time in stage per hour normalized
WakeMutBLRelmean = mean(Wake24hMutBL./60*100,2,'omitnan')';
WakeMutBLRelstd  = std((Wake24hMutBL./60*100)','omitnan')/sqrt(size(Wake24hMutBL,2));
NREMMutBLRelmean = mean(NREM24hMutBL./60*100,2,'omitnan')';
NREMMutBLRelstd  = std((NREM24hMutBL./60*100)','omitnan')/sqrt(size(NREM24hMutBL,2));
REMMutBLRelmean  = mean(REM24hMutBL./60*100,2,'omitnan')';
REMMutBLRelstd   = std((REM24hMutBL./60*100)','omitnan')/sqrt(size(REM24hMutBL,2));

% --- WT Sleep Dep time in stage per hour
WakeWTSDRelmean = mean(Wake24hWTSD./60*100,2,'omitnan')';
WakeWTSDRelstd  = std((Wake24hWTSD./60*100)','omitnan')/sqrt(size(Wake24hWTSD,2));
NREMWTSDRelmean = mean(NREM24hWTSD./60*100,2,'omitnan')';
NREMWTSDRelstd  = std((NREM24hWTSD./60*100)','omitnan')/sqrt(size(NREM24hWTSD,2));
REMWTSDRelmean  = mean(REM24hWTSD./60*100,2,'omitnan')';
REMWTSDRelstd   = std((REM24hWTSD./60*100)','omitnan')/sqrt(size(REM24hWTSD,2));

% --- Mutant Sleep Dep time in stage per hour 
WakeMutSDRelmean = mean(Wake24hMutSD./60*100,2,'omitnan')';
WakeMutSDRelstd  = std((Wake24hMutSD./60*100)','omitnan')/sqrt(size(Wake24hMutSD,2));
NREMMutSDRelmean = mean(NREM24hMutSD./60*100,2,'omitnan')';
NREMMutSDRelstd  = std((NREM24hMutSD./60*100)','omitnan')/sqrt(size(NREM24hMutSD,2));
REMMutSDRelmean  = mean(REM24hMutSD./60*100,2,'omitnan')';
REMMutSDRelstd   = std((REM24hMutSD./60*100)','omitnan')/sqrt(size(REM24hMutSD,2));
%%
% ---------------------------------------------------------------------------------




% Figure out which ANOVAs need a post-hoc done, so we know where 
% to add asterisks to the plot. 
% ------ BL first 12 hours ------
% - Wake - 
if RepMeas_ANOVA_tables.BL.First12hr.Wake.pValue('Genotype') < 0.05 | RepMeas_ANOVA_tables.BL.First12hr.Wake.pValue('Genotype:Time') < 0.05
	PostHoc_bool.BLFirst12hr.W = 1;
	mc_table.BL.First12hr.Wake 		  = multcompare(RM_Models.BL.LP.Wake,'Genotype','By','Time','ComparisonType','lsd');
	
	mc_table.BL.First12hr.Wake = renamevars(mc_table.BL.First12hr.Wake,'pValue','pValueRAW');
	mc_table.BL.First12hr.Wake.pValueAdj = mafdr(mc_table.BL.First12hr.Wake.pValueRAW,'BHFDR',true);
	[~, ia] = unique(mc_table.BL.First12hr.Wake.Time(:), 'stable');  % remove duplicate rows
	mc_table.BL.First12hr.Wake = mc_table.BL.First12hr.Wake(ia,:);
	mc_table.BL.First12hr.Wake = movevars(mc_table.BL.First12hr.Wake, 'pValueAdj', 'After', 'pValueRAW');
	

	xlocs_sig_effects.BL.WT_vs_Mut_LP_W = double(string(mc_table.BL.First12hr.Wake.Time(mc_table.BL.First12hr.Wake.pValueAdj<0.05))); % these are times (h), not indices

else 
	PostHoc_bool.BLFirst12hr.W = 0;
	mc_table.BL.First12hr.Wake = [];
end 

% - NREM -
if RepMeas_ANOVA_tables.BL.First12hr.NREM.pValue('Genotype') < 0.05 | RepMeas_ANOVA_tables.BL.First12hr.NREM.pValue('Genotype:Time') < 0.05
	PostHoc_bool.BLFirst12hr.N = 1; 
	mc_table.BL.First12hr.NREM 		  = multcompare(RM_Models.BL.LP.NREM,'Genotype','By','Time','ComparisonType','lsd');
	
	mc_table.BL.First12hr.NREM = renamevars(mc_table.BL.First12hr.NREM,'pValue','pValueRAW');
	mc_table.BL.First12hr.NREM.pValueAdj = mafdr(mc_table.BL.First12hr.NREM.pValueRAW,'BHFDR',true); 
	[~, ia] = unique(mc_table.BL.First12hr.NREM.Time(:), 'stable');  % remove duplicate rows
	mc_table.BL.First12hr.NREM = mc_table.BL.First12hr.NREM(ia,:);
	mc_table.BL.First12hr.NREM = movevars(mc_table.BL.First12hr.NREM,'pValueAdj', 'After', 'pValueRAW');

	xlocs_sig_effects.BL.WT_vs_Mut_LP_N = double(string(mc_table.BL.First12hr.NREM.Time(mc_table.BL.First12hr.NREM.pValueAdj<0.05))); % these are times (h), not indices
else 
	PostHoc_bool.BLFirst12hr.N = 0; 
	mc_table.BL.First12hr.NREM = [];
end 

% - REM - 
if RepMeas_ANOVA_tables.BL.First12hr.REM.pValue('Genotype') < 0.05 | RepMeas_ANOVA_tables.BL.First12hr.REM.pValue('Genotype:Time') < 0.05
	PostHoc_bool.BLFirst12hr.R = 1;
	mc_table.BL.First12hr.REM 		 = multcompare(RM_Models.BL.LP.REM,'Genotype','By','Time','ComparisonType','lsd');
	
	mc_table.BL.First12hr.REM = renamevars(mc_table.BL.First12hr.REM,'pValue','pValueRAW');
	mc_table.BL.First12hr.REM.pValueAdj = mafdr(mc_table.BL.First12hr.REM.pValueRAW,'BHFDR',true);
	[~, ia] = unique(mc_table.BL.First12hr.REM.Time(:), 'stable');  % remove duplicate rows
	mc_table.BL.First12hr.REM = mc_table.BL.First12hr.REM(ia,:);
	mc_table.BL.First12hr.REM = movevars(mc_table.BL.First12hr.REM,'pValueAdj', 'After', 'pValueRAW');

	xlocs_sig_effects.BL.WT_vs_Mut_LP_R = double(string(mc_table.BL.First12hr.REM.Time(mc_table.BL.First12hr.REM.pValueAdj<0.05))); % these are times (h), not indices
else 
	PostHoc_bool.BLFirst12hr.R = 0;
	mc_table.BL.First12hr.REM = [];
end 

% ------ BL last 12 hours ------
% - Wake - 
if RepMeas_ANOVA_tables.BL.Last12hr.Wake.pValue('Genotype') < 0.05 | RepMeas_ANOVA_tables.BL.Last12hr.Wake.pValue('Genotype:Time') < 0.05
	PostHoc_bool.BLLast12hr.W = 1;
	mc_table.BL.Last12hr.Wake 		 = multcompare(RM_Models.BL.DP.Wake,'Genotype','By','Time','ComparisonType','lsd');
	
	mc_table.BL.Last12hr.Wake = renamevars(mc_table.BL.Last12hr.Wake,'pValue','pValueRAW');
	mc_table.BL.Last12hr.Wake.pValueAdj = mafdr(mc_table.BL.Last12hr.Wake.pValueRAW,'BHFDR',true); 
	[~, ia] = unique(mc_table.BL.Last12hr.Wake.Time(:), 'stable');  % remove duplicate rows
	mc_table.BL.Last12hr.Wake = mc_table.BL.Last12hr.Wake(ia,:);
	mc_table.BL.Last12hr.Wake = movevars(mc_table.BL.Last12hr.Wake,'pValueAdj', 'After', 'pValueRAW');

	xlocs_sig_effects.BL.WT_vs_Mut_DP_W = double(string(mc_table.BL.Last12hr.Wake.Time(mc_table.BL.Last12hr.Wake.pValueAdj<0.05))); % these are times (h), not indices
else 
	PostHoc_bool.BLLast12hr.W = 0;
	mc_table.BL.Last12hr.Wake = [];
end 

% - NREM - 
if RepMeas_ANOVA_tables.BL.Last12hr.NREM.pValue('Genotype') < 0.05 | RepMeas_ANOVA_tables.BL.Last12hr.NREM.pValue('Genotype:Time') < 0.05
	PostHoc_bool.BLLast12hr.N = 1;
	mc_table.BL.Last12hr.NREM 		 = multcompare(RM_Models.BL.DP.NREM,'Genotype','By','Time','ComparisonType','lsd');
	
	mc_table.BL.Last12hr.NREM = renamevars(mc_table.BL.Last12hr.NREM,'pValue','pValueRAW');
	mc_table.BL.Last12hr.NREM.pValueAdj = mafdr(mc_table.BL.Last12hr.NREM.pValueRAW,'BHFDR',true); 
	[~, ia] = unique(mc_table.BL.Last12hr.NREM.Time(:), 'stable');  % remove duplicate rows
	mc_table.BL.Last12hr.NREM = mc_table.BL.Last12hr.NREM(ia,:);
	mc_table.BL.Last12hr.NREM = movevars(mc_table.BL.Last12hr.NREM,'pValueAdj', 'After', 'pValueRAW');

	xlocs_sig_effects.BL.WT_vs_Mut_DP_N = double(string(mc_table.BL.Last12hr.NREM.Time(mc_table.BL.Last12hr.NREM.pValueAdj<0.05))); % these are times (h), not indices
else 
	PostHoc_bool.BLLast12hr.N = 0;
	mc_table.BL.Last12hr.NREM = [];
end 

% - REM - 
if RepMeas_ANOVA_tables.BL.Last12hr.REM.pValue('Genotype') < 0.05 | RepMeas_ANOVA_tables.BL.Last12hr.REM.pValue('Genotype:Time') < 0.05
	PostHoc_bool.BLLast12hr.R = 1;
	mc_table.BL.Last12hr.REM 		= multcompare(RM_Models.BL.DP.REM,'Genotype','By','Time','ComparisonType','lsd');
	
	mc_table.BL.Last12hr.REM = renamevars(mc_table.BL.Last12hr.REM,'pValue','pValueRAW');
	mc_table.BL.Last12hr.REM.pValueAdj = mafdr(mc_table.BL.Last12hr.REM.pValueRAW,'BHFDR',true); 
	[~, ia] = unique(mc_table.BL.Last12hr.REM.Time(:), 'stable');  % remove duplicate rows
	mc_table.BL.Last12hr.REM = mc_table.BL.Last12hr.REM(ia,:);
	mc_table.BL.Last12hr.REM = movevars(mc_table.BL.Last12hr.REM,'pValueAdj', 'After', 'pValueRAW');

	xlocs_sig_effects.BL.WT_vs_Mut_DP_R = double(string(mc_table.BL.Last12hr.REM.Time(mc_table.BL.Last12hr.REM.pValueAdj<0.05))); % these are times (h), not indices
else 
	PostHoc_bool.BLLast12hr.R = 0;
	mc_table.BL.Last12hr.REM  = [];
end 

% ------ SD first 12 hours ------
% - Wake - 
if RepMeas_ANOVA_tables.SD.First12hr.Wake.pValue('Genotype') < 0.05 | RepMeas_ANOVA_tables.SD.First12hr.Wake.pValue('Genotype:Time') < 0.05
	PostHoc_bool.SDFirst12hr.W = 1;
	mc_table.SD.First12hr.Wake 		  = multcompare(RM_Models.SD.LP.Wake,'Genotype','By','Time','ComparisonType','lsd');
	
	mc_table.SD.First12hr.Wake = renamevars(mc_table.SD.First12hr.Wake,'pValue','pValueRAW');
	mc_table.SD.First12hr.Wake.pValueAdj = mafdr(mc_table.SD.First12hr.Wake.pValueRAW,'BHFDR',true);
	[~, ia] = unique(mc_table.SD.First12hr.Wake.Time(:), 'stable');  % remove duplicate rows
	mc_table.SD.First12hr.Wake = mc_table.SD.First12hr.Wake(ia,:);
	mc_table.SD.First12hr.Wake = movevars(mc_table.SD.First12hr.Wake,'pValueAdj', 'After', 'pValueRAW');

	xlocs_sig_effects.SD.WT_vs_Mut_LP_W = double(string(mc_table.SD.First12hr.Wake.Time(mc_table.SD.First12hr.Wake.pValueAdj<0.05))); % these are times (h), not indices
else 
	PostHoc_bool.SDFirst12hr.W = 0;
	mc_table.SD.First12hr.Wake = [];
end 

% - NREM -
if RepMeas_ANOVA_tables.SD.First12hr.NREM.pValue('Genotype') < 0.05 | RepMeas_ANOVA_tables.SD.First12hr.NREM.pValue('Genotype:Time') < 0.05
	PostHoc_bool.SDFirst12hr.N = 1; 
	mc_table.SD.First12hr.NREM 		  = multcompare(RM_Models.SD.LP.NREM,'Genotype','By','Time','ComparisonType','lsd');
	
	mc_table.SD.First12hr.NREM = renamevars(mc_table.SD.First12hr.NREM,'pValue','pValueRAW');
	mc_table.SD.First12hr.NREM.pValueAdj = mafdr(mc_table.SD.First12hr.NREM.pValueRAW,'BHFDR',true);
	[~, ia] = unique(mc_table.SD.First12hr.NREM.Time(:), 'stable');  % remove duplicate rows
	mc_table.SD.First12hr.NREM = mc_table.SD.First12hr.NREM(ia,:);
	mc_table.SD.First12hr.NREM = movevars(mc_table.SD.First12hr.NREM,'pValueAdj', 'After', 'pValueRAW');

	xlocs_sig_effects.SD.WT_vs_Mut_LP_N = double(string(mc_table.SD.First12hr.NREM.Time(mc_table.SD.First12hr.NREM.pValueAdj<0.05))); % these are times (h), not indices
else 
	PostHoc_bool.SDFirst12hr.N = 0;
	mc_table.SD.First12hr.NREM = []; 
end 

% - REM - 
if RepMeas_ANOVA_tables.SD.First12hr.REM.pValue('Genotype') < 0.05 | RepMeas_ANOVA_tables.SD.First12hr.REM.pValue('Genotype:Time') < 0.05
	PostHoc_bool.SDFirst12hr.R = 1;
	mc_table.SD.First12hr.REM 		 = multcompare(RM_Models.SD.LP.REM,'Genotype','By','Time','ComparisonType','lsd');
	
	mc_table.SD.First12hr.REM = renamevars(mc_table.SD.First12hr.REM,'pValue','pValueRAW');
	mc_table.SD.First12hr.REM.pValueAdj = mafdr(mc_table.SD.First12hr.REM.pValueRAW,'BHFDR',true);
	[~, ia] = unique(mc_table.SD.First12hr.REM.Time(:), 'stable');  % remove duplicate rows
	mc_table.SD.First12hr.REM = mc_table.SD.First12hr.REM(ia,:);
	mc_table.SD.First12hr.REM = movevars(mc_table.SD.First12hr.REM,'pValueAdj', 'After', 'pValueRAW');

	xlocs_sig_effects.SD.WT_vs_Mut_LP_R = double(string(mc_table.SD.First12hr.REM.Time(mc_table.SD.First12hr.REM.pValueAdj<0.05))); % these are times (h), not indices
else 
	PostHoc_bool.SDFirst12hr.R = 0;
	mc_table.SD.First12hr.REM  = [];
end 

% ------ SD last 12 hours ------
% - Wake - 
if RepMeas_ANOVA_tables.SD.Last12hr.Wake.pValue('Genotype') < 0.05 | RepMeas_ANOVA_tables.SD.Last12hr.Wake.pValue('Genotype:Time') < 0.05
	PostHoc_bool.SDLast12hr.W = 1;
	mc_table.SD.Last12hr.Wake 		 = multcompare(RM_Models.SD.DP.Wake,'Genotype','By','Time','ComparisonType','lsd');
	
	mc_table.SD.Last12hr.Wake = renamevars(mc_table.SD.Last12hr.Wake,'pValue','pValueRAW');
	mc_table.SD.Last12hr.Wake.pValueAdj = mafdr(mc_table.SD.Last12hr.Wake.pValueRAW,'BHFDR',true);
	[~, ia] = unique(mc_table.SD.Last12hr.Wake.Time(:), 'stable');  % remove duplicate rows
	mc_table.SD.Last12hr.Wake = mc_table.SD.Last12hr.Wake(ia,:);
	mc_table.SD.Last12hr.Wake = movevars(mc_table.SD.Last12hr.Wake,'pValueAdj', 'After', 'pValueRAW');

	xlocs_sig_effects.SD.WT_vs_Mut_DP_W = double(string(mc_table.SD.Last12hr.Wake.Time(mc_table.SD.Last12hr.Wake.pValueAdj<0.05))); % these are times (h), not indices
else 
	PostHoc_bool.SDLast12hr.W = 0;
	mc_table.SD.Last12hr.Wake = [];
end 

% - NREM - 
if RepMeas_ANOVA_tables.SD.Last12hr.NREM.pValue('Genotype') < 0.05 | RepMeas_ANOVA_tables.SD.Last12hr.NREM.pValue('Genotype:Time') < 0.05
	PostHoc_bool.SDLast12hr.N = 1;
	mc_table.SD.Last12hr.NREM 		 = multcompare(RM_Models.SD.DP.NREM,'Genotype','By','Time','ComparisonType','lsd');
	
	mc_table.SD.Last12hr.NREM = renamevars(mc_table.SD.Last12hr.NREM,'pValue','pValueRAW');
	mc_table.SD.Last12hr.NREM.pValueAdj = mafdr(mc_table.SD.Last12hr.NREM.pValueRAW,'BHFDR',true);
	[~, ia] = unique(mc_table.SD.Last12hr.NREM.Time(:), 'stable');  % remove duplicate rows
	mc_table.SD.Last12hr.NREM = mc_table.SD.Last12hr.NREM(ia,:);
	mc_table.SD.Last12hr.NREM = movevars(mc_table.SD.Last12hr.NREM,'pValueAdj', 'After', 'pValueRAW');

	xlocs_sig_effects.SD.WT_vs_Mut_DP_N = double(string(mc_table.SD.Last12hr.NREM.Time(mc_table.SD.Last12hr.NREM.pValueAdj<0.05))); % these are times (h), not indices
else 
	PostHoc_bool.SDLast12hr.N = 0;
	mc_table.SD.Last12hr.NREM = [];
end 

% - REM - 
if RepMeas_ANOVA_tables.SD.Last12hr.REM.pValue('Genotype') < 0.05 | RepMeas_ANOVA_tables.SD.Last12hr.REM.pValue('Genotype:Time') < 0.05
	PostHoc_bool.SDLast12hr.R = 1;
	mc_table.SD.Last12hr.REM 		= multcompare(RM_Models.SD.DP.REM,'Genotype','By','Time','ComparisonType','lsd');
	
	mc_table.SD.Last12hr.REM = renamevars(mc_table.SD.Last12hr.REM,'pValue','pValueRAW');
	mc_table.SD.Last12hr.REM.pValueAdj = mafdr(mc_table.SD.Last12hr.REM.pValueRAW,'BHFDR',true);
	[~, ia] = unique(mc_table.SD.Last12hr.REM.Time(:), 'stable');  % remove duplicate rows
	mc_table.SD.Last12hr.REM = mc_table.SD.Last12hr.REM(ia,:);
	mc_table.SD.Last12hr.REM = movevars(mc_table.SD.Last12hr.REM,'pValueAdj', 'After', 'pValueRAW');

	xlocs_sig_effects.SD.WT_vs_Mut_DP_R = double(string(mc_table.SD.Last12hr.REM.Time(mc_table.SD.Last12hr.REM.pValueAdj<0.05))); % these are times (h), not indices
else 
	PostHoc_bool.SDLast12hr.R = 0;
	mc_table.SD.Last12hr.REM  = [];
end 
% -- end of SD last 12 hours -----------------


% ------- BL First 6hr Dark Period -----------
if isfield(RepMeas_ANOVA_tables.BL,'First6hrDP')
	% - Wake - 
	if RepMeas_ANOVA_tables.BL.First6hrDP.Wake.pValue('Genotype') < 0.05 | RepMeas_ANOVA_tables.BL.First6hrDP.Wake.pValue('Genotype:Time') < 0.05
		PostHoc_bool.BLFirst6hrDP.W = 1;
		mc_table.BL.First6hrDP.Wake = multcompare(RM_Models.BL.First6hrDP.Wake,'Genotype','By','Time','ComparisonType','lsd');
		
		mc_table.BL.First6hrDP.Wake = renamevars(mc_table.BL.First6hrDP.Wake,'pValue','pValueRAW');
		mc_table.BL.First6hrDP.Wake.pValueAdj = mafdr(mc_table.BL.First6hrDP.Wake.pValueRAW,'BHFDR',true);
		[~, ia] = unique(mc_table.BL.First6hrDP.Wake.Time(:), 'stable');  % remove duplicate rows
		mc_table.BL.First6hrDP.Wake = mc_table.BL.First6hrDP.Wake(ia,:);
		mc_table.BL.First6hrDP.Wake = movevars(mc_table.BL.First6hrDP.Wake,'pValueAdj', 'After', 'pValueRAW');

		xlocs_sig_effects.BL.WT_vs_Mut_First6hrDP_W = double(string(mc_table.BL.First6hrDP.Wake.Time(mc_table.BL.First6hrDP.Wake.pValueAdj<0.05))); % these are times (h), not indices
	else 
		PostHoc_bool.BLFirst6hrDP.W = 0;
		mc_table.BL.First6hrDP.Wake = [];
	end 

	% - NREM - 
	if RepMeas_ANOVA_tables.BL.First6hrDP.NREM.pValue('Genotype') < 0.05 | RepMeas_ANOVA_tables.BL.First6hrDP.NREM.pValue('Genotype:Time') < 0.05
		PostHoc_bool.BLFirst6hrDP.N = 1;
		mc_table.BL.First6hrDP.NREM 	   = multcompare(RM_Models.BL.First6hrDP.NREM,'Genotype','By','Time','ComparisonType','lsd');
		
		mc_table.BL.First6hrDP.NREM = renamevars(mc_table.BL.First6hrDP.NREM,'pValue','pValueRAW');
		mc_table.BL.First6hrDP.NREM.pValueAdj = mafdr(mc_table.BL.First6hrDP.NREM.pValueRAW,'BHFDR',true);
		[~, ia] = unique(mc_table.BL.First6hrDP.NREM.Time(:), 'stable');  % remove duplicate rows
		mc_table.BL.First6hrDP.NREM = mc_table.BL.First6hrDP.NREM(ia,:);
		mc_table.BL.First6hrDP.NREM = movevars(mc_table.BL.First6hrDP.NREM,'pValueAdj', 'After', 'pValueRAW');

		xlocs_sig_effects.BL.WT_vs_Mut_First6hrDP_N = double(string(mc_table.BL.First6hrDP.NREM.Time(mc_table.BL.First6hrDP.NREM.pValueAdj<0.05))); % these are times (h), not indices
	else 
		PostHoc_bool.BLFirst6hrDP.N = 0;
		mc_table.BL.First6hrDP.NREM = [];
	end 

	% - REM - 
	if RepMeas_ANOVA_tables.BL.First6hrDP.REM.pValue('Genotype') < 0.05 | RepMeas_ANOVA_tables.BL.First6hrDP.REM.pValue('Genotype:Time') < 0.05
		PostHoc_bool.BLFirst6hrDP.R = 1;
		mc_table.BL.First6hrDP.REM 		  = multcompare(RM_Models.BL.First6hrDP.REM,'Genotype','By','Time','ComparisonType','lsd');
		
		mc_table.BL.First6hrDP.REM = renamevars(mc_table.BL.First6hrDP.REM,'pValue','pValueRAW');
		mc_table.BL.First6hrDP.REM.pValueAdj = mafdr(mc_table.BL.First6hrDP.REM.pValueRAW,'BHFDR',true);
		[~, ia] = unique(mc_table.BL.First6hrDP.REM.Time(:), 'stable');  % remove duplicate rows
		mc_table.BL.First6hrDP.REM = mc_table.BL.First6hrDP.REM(ia,:);
		mc_table.BL.First6hrDP.REM = movevars(mc_table.BL.First6hrDP.REM,'pValueAdj', 'After', 'pValueRAW');

		xlocs_sig_effects.BL.WT_vs_Mut_First6hrDP_R = double(string(mc_table.BL.First6hrDP.REM.Time(mc_table.BL.First6hrDP.REM.pValueAdj<0.05))); % these are times (h), not indices
	else 
		PostHoc_bool.BLFirst6hrDP.R = 0;
		mc_table.BL.First6hrDP.REM = [];
	end 
end % end of BL first 6hr Dark Period

% ------- BL Last 6hr Dark Period -----------
if isfield(RepMeas_ANOVA_tables.BL,'Last6hrDP')
	% - Wake - 
	if RepMeas_ANOVA_tables.BL.Last6hrDP.Wake.pValue('Genotype') < 0.05 | RepMeas_ANOVA_tables.BL.Last6hrDP.Wake.pValue('Genotype:Time') < 0.05
		PostHoc_bool.BLLast6hrDP.W = 1;
		mc_table.BL.Last6hrDP.Wake 		  = multcompare(RM_Models.BL.Last6hrDP.Wake,'Genotype','By','Time','ComparisonType','lsd');
		
		mc_table.BL.Last6hrDP.Wake = renamevars(mc_table.BL.Last6hrDP.Wake,'pValue','pValueRAW');
		mc_table.BL.Last6hrDP.Wake.pValueAdj = mafdr(mc_table.BL.Last6hrDP.Wake.pValueRAW,'BHFDR',true);
		[~, ia] = unique(mc_table.BL.Last6hrDP.Wake.Time(:), 'stable');  % remove duplicate rows
		mc_table.BL.Last6hrDP.Wake = mc_table.BL.Last6hrDP.Wake(ia,:);
		mc_table.BL.Last6hrDP.Wake = movevars(mc_table.BL.Last6hrDP.Wake,'pValueAdj', 'After', 'pValueRAW');

		xlocs_sig_effects.BL.WT_vs_Mut_Last6hrDP_W = double(string(mc_table.BL.Last6hrDP.Wake.Time(mc_table.BL.Last6hrDP.Wake.pValueAdj<0.05))); % these are times (h), not indices
	else 
		PostHoc_bool.BLLast6hrDP.W = 0;
		mc_table.BL.Last6hrDP.Wake = [];
	end 

	% - NREM - 
	if RepMeas_ANOVA_tables.BL.Last6hrDP.NREM.pValue('Genotype') < 0.05 | RepMeas_ANOVA_tables.BL.Last6hrDP.NREM.pValue('Genotype:Time') < 0.05
		PostHoc_bool.BLLast6hrDP.N = 1;
		mc_table.BL.Last6hrDP.NREM 		  = multcompare(RM_Models.BL.Last6hrDP.NREM,'Genotype','By','Time','ComparisonType','lsd');
		
		mc_table.BL.Last6hrDP.NREM = renamevars(mc_table.BL.Last6hrDP.NREM,'pValue','pValueRAW');
		mc_table.BL.Last6hrDP.NREM.pValueAdj = mafdr(mc_table.BL.Last6hrDP.NREM.pValueRAW,'BHFDR',true);
		[~, ia] = unique(mc_table.BL.Last6hrDP.NREM.Time(:), 'stable');  % remove duplicate rows
		mc_table.BL.Last6hrDP.NREM = mc_table.BL.Last6hrDP.NREM(ia,:);
		mc_table.BL.Last6hrDP.NREM = movevars(mc_table.BL.Last6hrDP.NREM,'pValueAdj', 'After', 'pValueRAW');

		xlocs_sig_effects.BL.WT_vs_Mut_Last6hrDP_N = double(string(mc_table.BL.Last6hrDP.NREM.Time(mc_table.BL.Last6hrDP.NREM.pValueAdj<0.05))); % these are times (h), not indices
	else 
		PostHoc_bool.BLLast6hrDP.N = 0;
		mc_table.BL.Last6hrDP.NREM = [];
	end 

	% - REM - 
	if RepMeas_ANOVA_tables.BL.Last6hrDP.REM.pValue('Genotype') < 0.05 | RepMeas_ANOVA_tables.BL.Last6hrDP.REM.pValue('Genotype:Time') < 0.05
		PostHoc_bool.BLLast6hrDP.R = 1;
		mc_table.BL.Last6hrDP.REM 		 = multcompare(RM_Models.BL.Last6hrDP.REM,'Genotype','By','Time','ComparisonType','lsd');
		
		mc_table.BL.Last6hrDP.REM = renamevars(mc_table.BL.Last6hrDP.REM,'pValue','pValueRAW');
		mc_table.BL.Last6hrDP.REM.pValueAdj = mafdr(mc_table.BL.Last6hrDP.REM.pValueRAW,'BHFDR',true);
		[~, ia] = unique(mc_table.BL.Last6hrDP.REM.Time(:), 'stable');  % remove duplicate rows
		mc_table.BL.Last6hrDP.REM = mc_table.BL.Last6hrDP.REM(ia,:);
		mc_table.BL.Last6hrDP.REM = movevars(mc_table.BL.Last6hrDP.REM,'pValueAdj', 'After', 'pValueRAW');

		xlocs_sig_effects.BL.WT_vs_Mut_Last6hrDP_R = double(string(mc_table.BL.Last6hrDP.REM.Time(mc_table.BL.Last6hrDP.REM.pValueAdj<0.05))); % these are times (h), not indices
	else 
		PostHoc_bool.BLLast6hrDP.R = 0;
		mc_table.BL.Last6hrDP.REM = [];
	end 
end % end of BL last 6hr Dark Period



% ------- SD First 6hr Dark Period -----------
if isfield(RepMeas_ANOVA_tables.SD,'First6hrDP')
	% - Wake - 
	if RepMeas_ANOVA_tables.SD.First6hrDP.Wake.pValue('Genotype') < 0.05 | RepMeas_ANOVA_tables.SD.First6hrDP.Wake.pValue('Genotype:Time') < 0.05
		PostHoc_bool.SDFirst6hrDP.W = 1;
		mc_table.SD.First6hrDP.Wake 	   = multcompare(RM_Models.SD.First6hrDP.Wake,'Genotype','By','Time','ComparisonType','lsd');
		
		mc_table.SD.First6hrDP.Wake = renamevars(mc_table.SD.First6hrDP.Wake,'pValue','pValueRAW');
		mc_table.SD.First6hrDP.Wake.pValueAdj = mafdr(mc_table.SD.First6hrDP.Wake.pValueRAW,'BHFDR',true);
		[~, ia] = unique(mc_table.SD.First6hrDP.Wake.Time(:), 'stable');  % remove duplicate rows
		mc_table.SD.First6hrDP.Wake = mc_table.SD.First6hrDP.Wake(ia,:);
		mc_table.SD.First6hrDP.Wake = movevars(mc_table.SD.First6hrDP.Wake,'pValueAdj', 'After', 'pValueRAW');

		xlocs_sig_effects.SD.WT_vs_Mut_First6hrDP_W = double(string(mc_table.SD.First6hrDP.Wake.Time(mc_table.SD.First6hrDP.Wake.pValueAdj<0.05))); % these are times (h), not indices
	else 
		PostHoc_bool.SDFirst6hrDP.W = 0;
		mc_table.SD.First6hrDP.Wake = [];
	end 

	% - NREM - 
	if RepMeas_ANOVA_tables.SD.First6hrDP.NREM.pValue('Genotype') < 0.05 | RepMeas_ANOVA_tables.SD.First6hrDP.NREM.pValue('Genotype:Time') < 0.05
		PostHoc_bool.SDFirst6hrDP.N = 1;
		mc_table.SD.First6hrDP.NREM 	   = multcompare(RM_Models.SD.First6hrDP.NREM,'Genotype','By','Time','ComparisonType','lsd');
		
		mc_table.SD.First6hrDP.NREM = renamevars(mc_table.SD.First6hrDP.NREM,'pValue','pValueRAW');
		mc_table.SD.First6hrDP.NREM.pValueAdj = mafdr(mc_table.SD.First6hrDP.NREM.pValueRAW,'BHFDR',true);
		[~, ia] = unique(mc_table.SD.First6hrDP.NREM.Time(:), 'stable');  % remove duplicate rows
		mc_table.SD.First6hrDP.NREM = mc_table.SD.First6hrDP.NREM(ia,:);
		mc_table.SD.First6hrDP.NREM = movevars(mc_table.SD.First6hrDP.NREM,'pValueAdj', 'After', 'pValueRAW');

		xlocs_sig_effects.SD.WT_vs_Mut_First6hrDP_N = double(string(mc_table.SD.First6hrDP.NREM.Time(mc_table.SD.First6hrDP.NREM.pValueAdj<0.05))); % these are times (h), not indices
	else 
		PostHoc_bool.SDFirst6hrDP.N = 0;
		mc_table.SD.First6hrDP.NREM = [];
	end 

	% - REM - 
	if RepMeas_ANOVA_tables.SD.First6hrDP.REM.pValue('Genotype') < 0.05 | RepMeas_ANOVA_tables.SD.First6hrDP.REM.pValue('Genotype:Time') < 0.05
		PostHoc_bool.SDFirst6hrDP.R = 1;
		mc_table.SD.First6hrDP.REM 		  = multcompare(RM_Models.SD.First6hrDP.REM,'Genotype','By','Time','ComparisonType','lsd');
		
		mc_table.SD.First6hrDP.REM = renamevars(mc_table.SD.First6hrDP.REM,'pValue','pValueRAW');
		mc_table.SD.First6hrDP.REM.pValueAdj = mafdr(mc_table.SD.First6hrDP.REM.pValueRAW,'BHFDR',true);
		[~, ia] = unique(mc_table.SD.First6hrDP.REM.Time(:), 'stable');  % remove duplicate rows
		mc_table.SD.First6hrDP.REM = mc_table.SD.First6hrDP.REM(ia,:);
		mc_table.SD.First6hrDP.REM = movevars(mc_table.SD.First6hrDP.REM,'pValueAdj', 'After', 'pValueRAW');

		xlocs_sig_effects.SD.WT_vs_Mut_First6hrDP_R = double(string(mc_table.SD.First6hrDP.REM.Time(mc_table.SD.First6hrDP.REM.pValueAdj<0.05))); % these are times (h), not indices
	else 
		PostHoc_bool.SDFirst6hrDP.R = 0;
		mc_table.SD.First6hrDP.REM = [];
	end 
end % end of SD first 6hr Dark Period


% ------- SD Last 6hr Dark Period -----------
if isfield(RepMeas_ANOVA_tables.SD,'Last6hrDP')
	% - Wake - 
	if RepMeas_ANOVA_tables.SD.Last6hrDP.Wake.pValue('Genotype') < 0.05 | RepMeas_ANOVA_tables.SD.Last6hrDP.Wake.pValue('Genotype:Time') < 0.05
		PostHoc_bool.SDLast6hrDP.W = 1;
		mc_table.SD.Last6hrDP.Wake 		  = multcompare(RM_Models.SD.Last6hrDP.Wake,'Genotype','By','Time','ComparisonType','lsd');
		
		mc_table.SD.Last6hrDP.Wake = renamevars(mc_table.SD.Last6hrDP.Wake,'pValue','pValueRAW');
		mc_table.SD.Last6hrDP.Wake.pValueAdj = mafdr(mc_table.SD.Last6hrDP.Wake.pValueRAW,'BHFDR',true);
		[~, ia] = unique(mc_table.SD.Last6hrDP.Wake.Time(:), 'stable');  % remove duplicate rows
		mc_table.SD.Last6hrDP.Wake = mc_table.SD.Last6hrDP.Wake(ia,:);
		mc_table.SD.Last6hrDP.Wake = movevars(mc_table.SD.Last6hrDP.Wake,'pValueAdj', 'After', 'pValueRAW');

		xlocs_sig_effects.SD.WT_vs_Mut_Last6hrDP_W = double(string(mc_table.SD.Last6hrDP.Wake.Time(mc_table.SD.Last6hrDP.Wake.pValueAdj<0.05))); % these are times (h), not indices
	else 
		PostHoc_bool.SDLast6hrDP.W = 0;
		mc_table.SD.Last6hrDP.Wake = [];
	end 

	% - NREM - 
	if RepMeas_ANOVA_tables.SD.Last6hrDP.NREM.pValue('Genotype') < 0.05 | RepMeas_ANOVA_tables.SD.Last6hrDP.NREM.pValue('Genotype:Time') < 0.05
		PostHoc_bool.SDLast6hrDP.N = 1;
		mc_table.SD.Last6hrDP.NREM 		  = multcompare(RM_Models.SD.Last6hrDP.NREM,'Genotype','By','Time','ComparisonType','lsd');
		
		mc_table.SD.Last6hrDP.NREM = renamevars(mc_table.SD.Last6hrDP.NREM,'pValue','pValueRAW');
		mc_table.SD.Last6hrDP.NREM.pValueAdj = mafdr(mc_table.SD.Last6hrDP.NREM.pValueRAW,'BHFDR',true);
		[~, ia] = unique(mc_table.SD.Last6hrDP.NREM.Time(:), 'stable');  % remove duplicate rows
		mc_table.SD.Last6hrDP.NREM = mc_table.SD.Last6hrDP.NREM(ia,:);
		mc_table.SD.Last6hrDP.NREM = movevars(mc_table.SD.Last6hrDP.NREM,'pValueAdj', 'After', 'pValueRAW');

		xlocs_sig_effects.SD.WT_vs_Mut_Last6hrDP_N = double(string(mc_table.SD.Last6hrDP.NREM.Time(mc_table.SD.Last6hrDP.NREM.pValueAdj<0.05))); % these are times (h), not indices
	else 
		PostHoc_bool.SDLast6hrDP.N = 0;
		mc_table.SD.Last6hrDP.NREM = [];
	end 

	% - REM - 
	if RepMeas_ANOVA_tables.SD.Last6hrDP.REM.pValue('Genotype') < 0.05 | RepMeas_ANOVA_tables.SD.Last6hrDP.REM.pValue('Genotype:Time') < 0.05
		PostHoc_bool.SDLast6hrDP.R = 1;
		mc_table.SD.Last6hrDP.REM 		 = multcompare(RM_Models.SD.Last6hrDP.REM,'Genotype','By','Time','ComparisonType','lsd');
		
		mc_table.SD.Last6hrDP.REM = renamevars(mc_table.SD.Last6hrDP.REM,'pValue','pValueRAW');
		mc_table.SD.Last6hrDP.REM.pValueAdj = mafdr(mc_table.SD.Last6hrDP.REM.pValueRAW,'BHFDR',true);
		[~, ia] = unique(mc_table.SD.Last6hrDP.REM.Time(:), 'stable');  % remove duplicate rows
		mc_table.SD.Last6hrDP.REM = mc_table.SD.Last6hrDP.REM(ia,:);
		mc_table.SD.Last6hrDP.REM = movevars(mc_table.SD.Last6hrDP.REM,'pValueAdj', 'After', 'pValueRAW');

		xlocs_sig_effects.SD.WT_vs_Mut_Last6hrDP_R = double(string(mc_table.SD.Last6hrDP.REM.Time(mc_table.SD.Last6hrDP.REM.pValueAdj<0.05))); % these are times (h), not indices
	else 
		PostHoc_bool.SDLast6hrDP.R = 0;
		mc_table.SD.Last6hrDP.REM = [];
	end 
end % end of SD last 6hr Dark Period



% ---------------------------------------------------------------------------------
% --- MAKE THE FIGURES ------------------------------------------------------------
posA=([12 0 12 100]);
posB=([0 0 6 100]);

myfontsize   = 16;
mylinewidth  = 2;
mymarkersize = 5; 
astr_size    = 24;

LPboxpos = ([0 -5 12 5]);  %[x y w h] Light period box along bottom of plot
DPboxpos = ([12 -5 12 5]); %[x y w h] Dark period box along bottom of plot
SDboxpos = ([0 -1 SD_length_hrs 100]);  % [x y w h] Sleep Dep hatching

LPboxposREM = ([0 -1 12 1]);  %[x y w h] Light period box along bottom of plot REM VERSION: SHORTER
DPboxposREM = ([12 -1 12 1]); %[x y w h] Dark period box along bottom of plot  REM VERSION: SHORTER

Sig12hrColor = [0.8 0.8 0.8];  % RGB for color of bar that indicates significant RM ANova over 12 hours


if PlotIndividualFigs
	FTIS_BL_RelativetoTRTa=figure;
	FTIS_BL_RelativetoTRTa.Name = 'BL Wake relative to Total Recording Time';

	
	% --- BL Wake ---
	rectangle('Position',LPboxpos,'FaceColor','w','EdgeColor','k','LineWidth',2);     % LP rectangle
	hold on
	rectangle('Position',DPboxpos,'FaceColor','k','EdgeColor','k','LineWidth',2);    % DP rectangle
	e1=errorbar(WakeWTBLRelmean,WakeWTBLRelstd,'Color','k', 'LineWidth',mylinewidth,'Marker','o','MarkerFaceColor','k','MarkerSize',mymarkersize,'CapSize',0);
	e2=errorbar(WakeMutBLRelmean,WakeMutBLRelstd,'Color','r','LineWidth',mylinewidth,'Marker','o','MarkerFaceColor','r','MarkerSize',mymarkersize,'CapSize',0);
	hold off 
	set(gca, 'Layer','top')
	%title('Wake, Baseline');
	ylabel({'WAKE','Time in Wake (% TRT)'});
	xlabel('Time (hours)');
	ax=gca;
	l = legend([e1 e2],LegendLabels);
	l.Box = 'off';
	l.Location = 'None';
	l.Position = [0.7326 0.1649 0.1634 0.0798];
	l.FontSize = 12;
	set(gca, 'XTick',[0 6 12 18 24]); %[1 4 7 10 13 16 19 22]
	set(gca, 'XTickLabels',{'0' '6' '12' '18' '24'}); %{'1' '4' '7' '10' '13' '16' '19' '22'}
	xlim([0 24])
	ylim([-5 100])
	ax.FontSize = myfontsize;
	ax.LineWidth = mylinewidth;
	% -- Add stars if sig diffs (LP)
	if PostHoc_bool.BLFirst12hr.W
		hold on 
		
		rectangle('Position',[0 ax.YLim(2)-5 12 5],'FaceColor',Sig12hrColor,'EdgeColor','none'); 
		bool_xlocsBL_LP_sig_values = ismember(e1.XData,xlocs_sig_effects.BL.WT_vs_Mut_LP_W);
		largest_y1_vals = e1.YData(bool_xlocsBL_LP_sig_values) + e1.YPositiveDelta(bool_xlocsBL_LP_sig_values);
		largest_y2_vals = e2.YData(bool_xlocsBL_LP_sig_values) + e2.YPositiveDelta(bool_xlocsBL_LP_sig_values);
		ylocs_sig_effects = max([largest_y1_vals; largest_y2_vals]) + 10*ones(1,length(largest_y1_vals));
		ylocs_sig_effects(ylocs_sig_effects>ax.YLim(2))=ax.YLim(2);
		text(xlocs_sig_effects.BL.WT_vs_Mut_LP_W, ylocs_sig_effects, '*', 'FontSize',astr_size,'HorizontalAlignment','center','VerticalAlignment','middle');
		hold off 
	end 
	% -- Add stars if sig diffs (DP)
	if PostHoc_bool.BLLast12hr.W
		hold on 
		
		rectangle('Position',[12 ax.YLim(2)-5 12 5],'FaceColor',Sig12hrColor,'EdgeColor','none'); 
		bool_xlocsBL_DP_sig_values = ismember(e1.XData,xlocs_sig_effects.BL.WT_vs_Mut_DP_W);
		largest_y1_vals = e1.YData(bool_xlocsBL_DP_sig_values) + e1.YPositiveDelta(bool_xlocsBL_DP_sig_values);
		largest_y2_vals = e2.YData(bool_xlocsBL_DP_sig_values) + e2.YPositiveDelta(bool_xlocsBL_DP_sig_values);
		ylocs_sig_effects = max([largest_y1_vals; largest_y2_vals]) + 10*ones(1,length(largest_y1_vals));
		ylocs_sig_effects(ylocs_sig_effects>ax.YLim(2))=ax.YLim(2);
		text(xlocs_sig_effects.BL.WT_vs_Mut_DP_W, ylocs_sig_effects, '*', 'FontSize',astr_size,'HorizontalAlignment','center','VerticalAlignment','middle');
		hold off 
	end 

	% --- BL NREM ---
	h=figure;
	h.Name = 'BL NREM relative to Total Recording Time';
	
	rectangle('Position',LPboxpos,'FaceColor','w','EdgeColor','k','LineWidth',2);     % LP rectangle
	hold on
	rectangle('Position',DPboxpos,'FaceColor','k','EdgeColor','k','LineWidth',2);    % DP rectangle

	e1=errorbar(NREMMutBLRelmean,NREMMutBLRelstd,'Color','r','LineWidth',mylinewidth,'Marker','o','MarkerFaceColor','r','MarkerSize',mymarkersize,'CapSize',0);
	e2=errorbar(NREMWTBLRelmean,NREMWTBLRelstd,'Color','k','LineWidth',mylinewidth,'Marker','o','MarkerFaceColor','k','MarkerSize',mymarkersize,'CapSize',0);
	hold off
	set(gca, 'Layer','top')
	
	set(gca, 'XTick',[0 6 12 18 24]);  %[1 4 7 10 13 16 19 22]
	set(gca, 'XTickLabels',{'0' '6' '12' '18' '24'}); %{'1' '4' '7' '10' '13' '16' '19' '22'}
	xlim([0 24])
	ylim([-5 100])
	ylabel({'NREM','Time in NREM (% TRT)'});
	xlabel('Time (hours)')
	ax = gca;
	ax.FontSize = myfontsize;
	ax.LineWidth = mylinewidth;
	% -- Add stars if sig diffs (LP)
	if PostHoc_bool.BLFirst12hr.N
		hold on 
		
		rectangle('Position',[0 ax.YLim(2)-5 12 5],'FaceColor',Sig12hrColor,'EdgeColor','none'); 
		bool_xlocsWTBL_LP_sig_values = ismember(e1.XData,xlocs_sig_effects.BL.WT_vs_Mut_LP_N);
		largest_y1_vals = e1.YData(bool_xlocsWTBL_LP_sig_values) + e1.YPositiveDelta(bool_xlocsWTBL_LP_sig_values);
		largest_y2_vals = e2.YData(bool_xlocsWTBL_LP_sig_values) + e2.YPositiveDelta(bool_xlocsWTBL_LP_sig_values);
		ylocs_sig_effects = max([largest_y1_vals; largest_y2_vals]) + 10*ones(1,length(largest_y1_vals));
		ylocs_sig_effects(ylocs_sig_effects>ax.YLim(2))=ax.YLim(2);
		text(xlocs_sig_effects.BL.WT_vs_Mut_LP_N, ylocs_sig_effects, '*', 'FontSize',astr_size,'HorizontalAlignment','center','VerticalAlignment','middle');
		hold off 
	end 
	% -- Add stars if sig diffs (DP)
	if PostHoc_bool.BLLast12hr.N
		hold on 
		
		rectangle('Position',[12 ax.YLim(2)-5 12 5],'FaceColor',Sig12hrColor,'EdgeColor','none'); 
		bool_xlocsWTBL_DP_sig_values = ismember(e1.XData,xlocs_sig_effects.BL.WT_vs_Mut_DP_N);
		largest_y1_vals = e1.YData(bool_xlocsWTBL_DP_sig_values) + e1.YPositiveDelta(bool_xlocsWTBL_DP_sig_values);
		largest_y2_vals = e2.YData(bool_xlocsWTBL_DP_sig_values) + e2.YPositiveDelta(bool_xlocsWTBL_DP_sig_values);
		ylocs_sig_effects = max([largest_y1_vals; largest_y2_vals]) + 10*ones(1,length(largest_y1_vals));
		ylocs_sig_effects(ylocs_sig_effects>ax.YLim(2))=ax.YLim(2);
		text(xlocs_sig_effects.BL.WT_vs_Mut_DP_N, ylocs_sig_effects, '*', 'FontSize',astr_size,'HorizontalAlignment','center','VerticalAlignment','middle');
		hold off 
	end 

	% --- BL REM ---
	H=figure;
	H.Name = 'BL REM relative to Total Recording Time';
	 
	rectangle('Position',LPboxposREM,'FaceColor','w','EdgeColor','k','LineWidth',2);     % LP rectangle
	hold on
	rectangle('Position',DPboxposREM,'FaceColor','k','EdgeColor','k','LineWidth',2);    % DP rectangle

	e1=errorbar(REMMutBLRelmean,REMMutBLRelstd,'Color','r','LineWidth',mylinewidth,'Marker','o','MarkerFaceColor','r','MarkerSize',mymarkersize,'CapSize',0);
	e2=errorbar(REMWTBLRelmean,REMWTBLRelstd,'Color','k','LineWidth',mylinewidth,'Marker','o','MarkerFaceColor','k','MarkerSize',mymarkersize,'CapSize',0);
	hold off
	set(gca, 'Layer','top')
	
	set(gca, 'XTick',[0 6 12 18 24]); %[1 4 7 10 13 16 19 22]
	set(gca, 'XTickLabels',{'0' '6' '12' '18' '24'}); %{'1' '4' '7' '10' '13' '16' '19' '22'}
	xlim([0 24])
	ylim ([-1 30])
	xlabel('Time (hours)')
	ylabel({'REM','Time in REM (% TRT)'});
	ax = gca;
	ax.FontSize = myfontsize;
	ax.LineWidth = mylinewidth;
	% -- Add stars if sig diffs (LP)
	if PostHoc_bool.BLFirst12hr.R
		hold on 
		
		rectangle('Position',[0 ax.YLim(2)-1 12 1],'FaceColor',Sig12hrColor,'EdgeColor','none'); 
		bool_xlocsWTBL_LP_sig_values = ismember(e1.XData,xlocs_sig_effects.BL.WT_vs_Mut_LP_R);
		largest_y1_vals = e1.YData(bool_xlocsWTBL_LP_sig_values) + e1.YPositiveDelta(bool_xlocsWTBL_LP_sig_values);
		largest_y2_vals = e2.YData(bool_xlocsWTBL_LP_sig_values) + e2.YPositiveDelta(bool_xlocsWTBL_LP_sig_values);
		ylocs_sig_effects = max([largest_y1_vals; largest_y2_vals]) + 10*ones(1,length(largest_y1_vals));
		ylocs_sig_effects(ylocs_sig_effects>ax.YLim(2))=ax.YLim(2);
		text(xlocs_sig_effects.BL.WT_vs_Mut_LP_R, ylocs_sig_effects, '*', 'FontSize',astr_size,'HorizontalAlignment','center','VerticalAlignment','middle');
		hold off 
	end 
	% -- Add stars if sig diffs (DP)
	if PostHoc_bool.BLLast12hr.R
		hold on 
		
		rectangle('Position',[12 ax.YLim(2)-1 12 1],'FaceColor',Sig12hrColor,'EdgeColor','none'); 
		bool_xlocsWTBL_DP_sig_values = ismember(e1.XData,xlocs_sig_effects.BL.WT_vs_Mut_DP_R);
		largest_y1_vals = e1.YData(bool_xlocsWTBL_DP_sig_values) + e1.YPositiveDelta(bool_xlocsWTBL_DP_sig_values);
		largest_y2_vals = e2.YData(bool_xlocsWTBL_DP_sig_values) + e2.YPositiveDelta(bool_xlocsWTBL_DP_sig_values);
		ylocs_sig_effects = max([largest_y1_vals; largest_y2_vals]) + 10*ones(1,length(largest_y1_vals));
		ylocs_sig_effects(ylocs_sig_effects>ax.YLim(2))=ax.YLim(2);
		text(xlocs_sig_effects.BL.WT_vs_Mut_DP_R, ylocs_sig_effects, '*', 'FontSize',astr_size,'HorizontalAlignment','center','VerticalAlignment','middle');
		hold off 
	end 
end % end of if PlotIndividualFigs

% ------ 3-panel figure for all of BL ---------
H=figure;
H.Name = ['Three-Panel BL Time in State Figure ',Sex];
H.Renderer = 'painters';
H.Position = [429 118 678 698];

t= tiledlayout(3,1,'TileSpacing','compact','Padding','compact','TileIndexing', 'columnmajor');
nexttile

% -- BL W -- 
rectangle('Position',LPboxpos,'FaceColor','w','EdgeColor','k','LineWidth',2);     % LP rectangle
hold on
rectangle('Position',DPboxpos,'FaceColor','k','EdgeColor','k','LineWidth',2);    % DP rectangle
e1=errorbar(WakeWTBLRelmean,WakeWTBLRelstd,'Color','k', 'LineWidth',mylinewidth,'Marker','o','MarkerFaceColor','k','MarkerSize',mymarkersize,'CapSize',0);
e2=errorbar(WakeMutBLRelmean,WakeMutBLRelstd,'Color','r','LineWidth',mylinewidth,'Marker','o','MarkerFaceColor','r','MarkerSize',mymarkersize,'CapSize',0);
hold off 
set(gca, 'Layer','top')
ylabel({'WAKE','Time in Wake (% TRT)'});
l = legend([e1 e2],LegendLabels);
ax=gca;
l.Box = 'off';
l.Location = 'None';
l.Location = 'SouthEast';
l.FontSize = 12;
set(gca,'XTick',[])
xlim([0 24])
ylim([-5 100])
ax.FontSize = myfontsize;
ax.LineWidth = mylinewidth;
% -- Add stars if sig diffs (LP)
if PostHoc_bool.BLFirst12hr.W
	hold on 
	
	rectangle('Position',[0 ax.YLim(2)-5 12 5],'FaceColor',Sig12hrColor,'EdgeColor','none'); 
	bool_xlocsBL_LP_sig_values = ismember(e1.XData,xlocs_sig_effects.BL.WT_vs_Mut_LP_W);
	largest_y1_vals = e1.YData(bool_xlocsBL_LP_sig_values) + e1.YPositiveDelta(bool_xlocsBL_LP_sig_values);
	largest_y2_vals = e2.YData(bool_xlocsBL_LP_sig_values) + e2.YPositiveDelta(bool_xlocsBL_LP_sig_values);
	ylocs_sig_effects = max([largest_y1_vals; largest_y2_vals]) + 10*ones(1,length(largest_y1_vals));
	ylocs_sig_effects(ylocs_sig_effects>ax.YLim(2))=ax.YLim(2);
	text(xlocs_sig_effects.BL.WT_vs_Mut_LP_W, ylocs_sig_effects, '*', 'FontSize',astr_size,'HorizontalAlignment','center','VerticalAlignment','middle');
	hold off 
end 
% -- Add stars if sig diffs (DP)
if PostHoc_bool.BLLast12hr.W
	hold on 
	
	rectangle('Position',[12 ax.YLim(2)-5 12 5],'FaceColor',Sig12hrColor,'EdgeColor','none'); 
	bool_xlocsBL_DP_sig_values = ismember(e1.XData,xlocs_sig_effects.BL.WT_vs_Mut_DP_W);
	largest_y1_vals = e1.YData(bool_xlocsBL_DP_sig_values) + e1.YPositiveDelta(bool_xlocsBL_DP_sig_values);
	largest_y2_vals = e2.YData(bool_xlocsBL_DP_sig_values) + e2.YPositiveDelta(bool_xlocsBL_DP_sig_values);
	ylocs_sig_effects = max([largest_y1_vals; largest_y2_vals]) + 10*ones(1,length(largest_y1_vals));
	ylocs_sig_effects(ylocs_sig_effects>ax.YLim(2))=ax.YLim(2);
	text(xlocs_sig_effects.BL.WT_vs_Mut_DP_W, ylocs_sig_effects, '*', 'FontSize',astr_size,'HorizontalAlignment','center','VerticalAlignment','middle');
	hold off 
end 



% -- BL NREM --
nexttile
rectangle('Position',LPboxpos,'FaceColor','w','EdgeColor','k','LineWidth',2);     % LP rectangle
hold on
rectangle('Position',DPboxpos,'FaceColor','k','EdgeColor','k','LineWidth',2);    % DP rectangle

e1=errorbar(NREMMutBLRelmean,NREMMutBLRelstd,'Color','r','LineWidth',mylinewidth,'Marker','o','MarkerFaceColor','r','MarkerSize',mymarkersize,'CapSize',0);
e2=errorbar(NREMWTBLRelmean,NREMWTBLRelstd,'Color','k','LineWidth',mylinewidth,'Marker','o','MarkerFaceColor','k','MarkerSize',mymarkersize,'CapSize',0);
hold off
set(gca, 'Layer','top')

set(gca,'XTick',[])
xlim([0 24])
ylim([-5 100])
ylabel({'NREM','Time in NREM (% TRT)'});
ax = gca;
ax.FontSize = myfontsize;
ax.LineWidth = mylinewidth;
% -- Add stars if sig diffs (LP)
if PostHoc_bool.BLFirst12hr.N
	hold on 
	
	rectangle('Position',[0 ax.YLim(2)-5 12 5],'FaceColor',Sig12hrColor,'EdgeColor','none'); 
	bool_xlocsWTBL_LP_sig_values = ismember(e1.XData,xlocs_sig_effects.BL.WT_vs_Mut_LP_N);
	largest_y1_vals = e1.YData(bool_xlocsWTBL_LP_sig_values) + e1.YPositiveDelta(bool_xlocsWTBL_LP_sig_values);
	largest_y2_vals = e2.YData(bool_xlocsWTBL_LP_sig_values) + e2.YPositiveDelta(bool_xlocsWTBL_LP_sig_values);
	ylocs_sig_effects = max([largest_y1_vals; largest_y2_vals]) + 10*ones(1,length(largest_y1_vals));
	ylocs_sig_effects(ylocs_sig_effects>ax.YLim(2))=ax.YLim(2);
	text(xlocs_sig_effects.BL.WT_vs_Mut_LP_N, ylocs_sig_effects, '*', 'FontSize',astr_size,'HorizontalAlignment','center','VerticalAlignment','middle');
	hold off 
end 
% -- Add stars if sig diffs (DP)
if PostHoc_bool.BLLast12hr.N
	hold on 
	
	rectangle('Position',[12 ax.YLim(2)-5 12 5],'FaceColor',Sig12hrColor,'EdgeColor','none'); 
	bool_xlocsWTBL_DP_sig_values = ismember(e1.XData,xlocs_sig_effects.BL.WT_vs_Mut_DP_N);
	largest_y1_vals = e1.YData(bool_xlocsWTBL_DP_sig_values) + e1.YPositiveDelta(bool_xlocsWTBL_DP_sig_values);
	largest_y2_vals = e2.YData(bool_xlocsWTBL_DP_sig_values) + e2.YPositiveDelta(bool_xlocsWTBL_DP_sig_values);
	ylocs_sig_effects = max([largest_y1_vals; largest_y2_vals]) + 10*ones(1,length(largest_y1_vals));
	ylocs_sig_effects(ylocs_sig_effects>ax.YLim(2))=ax.YLim(2);
	text(xlocs_sig_effects.BL.WT_vs_Mut_DP_N, ylocs_sig_effects, '*', 'FontSize',astr_size,'HorizontalAlignment','center','VerticalAlignment','middle');
	hold off 
end 

% -- BL REM --
nexttile
rectangle('Position',LPboxposREM,'FaceColor','w','EdgeColor','k','LineWidth',2);     % LP rectangle
hold on
rectangle('Position',DPboxposREM,'FaceColor','k','EdgeColor','k','LineWidth',2);    % DP rectangle

e1=errorbar(REMMutBLRelmean,REMMutBLRelstd,'Color','r','LineWidth',mylinewidth,'Marker','o','MarkerFaceColor','r','MarkerSize',mymarkersize,'CapSize',0);
e2=errorbar(REMWTBLRelmean,REMWTBLRelstd,'Color','k','LineWidth',mylinewidth,'Marker','o','MarkerFaceColor','k','MarkerSize',mymarkersize,'CapSize',0);
hold off
set(gca, 'Layer','top')

set(gca, 'XTick',[0 6 12 18 24]); %[1 4 7 10 13 16 19 22]
set(gca, 'XTickLabels',{'0' '6' '12' '18' '24'}); %{'1' '4' '7' '10' '13' '16' '19' '22'}
xlim([0 24])
ylim ([-1 30])
xlabel('Time (hours)')
ylabel({'REM','Time in REM (% TRT)'});
ax = gca;
ax.FontSize = myfontsize;
ax.LineWidth = mylinewidth;
% -- Add stars if sig diffs (LP)
if PostHoc_bool.BLFirst12hr.R
	hold on 
	
	rectangle('Position',[0 ax.YLim(2)-1 12 1],'FaceColor',Sig12hrColor,'EdgeColor','none'); 
	bool_xlocsWTBL_LP_sig_values = ismember(e1.XData,xlocs_sig_effects.BL.WT_vs_Mut_LP_R);
	largest_y1_vals = e1.YData(bool_xlocsWTBL_LP_sig_values) + e1.YPositiveDelta(bool_xlocsWTBL_LP_sig_values);
	largest_y2_vals = e2.YData(bool_xlocsWTBL_LP_sig_values) + e2.YPositiveDelta(bool_xlocsWTBL_LP_sig_values);
	ylocs_sig_effects = max([largest_y1_vals; largest_y2_vals]) + 10*ones(1,length(largest_y1_vals));
	ylocs_sig_effects(ylocs_sig_effects>ax.YLim(2))=ax.YLim(2);
	text(xlocs_sig_effects.BL.WT_vs_Mut_LP_R, ylocs_sig_effects, '*', 'FontSize',astr_size,'HorizontalAlignment','center','VerticalAlignment','middle');
	hold off 
end 
% -- Add stars if sig diffs (DP)
if PostHoc_bool.BLLast12hr.R
	hold on 
	
	rectangle('Position',[12 ax.YLim(2)-1 12 1],'FaceColor',Sig12hrColor,'EdgeColor','none'); 
	bool_xlocsWTBL_DP_sig_values = ismember(e1.XData,xlocs_sig_effects.BL.WT_vs_Mut_DP_R);
	largest_y1_vals = e1.YData(bool_xlocsWTBL_DP_sig_values) + e1.YPositiveDelta(bool_xlocsWTBL_DP_sig_values);
	largest_y2_vals = e2.YData(bool_xlocsWTBL_DP_sig_values) + e2.YPositiveDelta(bool_xlocsWTBL_DP_sig_values);
	ylocs_sig_effects = max([largest_y1_vals; largest_y2_vals]) + 10*ones(1,length(largest_y1_vals));
	ylocs_sig_effects(ylocs_sig_effects>ax.YLim(2))=ax.YLim(2);
	text(xlocs_sig_effects.BL.WT_vs_Mut_DP_R, ylocs_sig_effects, '*', 'FontSize',astr_size,'HorizontalAlignment','center','VerticalAlignment','middle');
	hold off 
end 


% ------ SD Figures ----------------------------------------------------------------
if PlotIndividualFigs
	FRS_TIS_relativetoTRT=figure;
	FRS_TIS_relativetoTRT.Name = 'SD Wake relative to Total Recording Time';



	% --- SD Wake ---
	rectangle('Position',LPboxpos,'FaceColor','w','EdgeColor','k','LineWidth',2);     % LP rectangle
	hold on 
	rectangle('Position',DPboxpos,'FaceColor','k','EdgeColor','k','LineWidth',2);    % DP rectangle
	rectangle('Position',[0 -5 SD_length_hrs 5],'FaceColor',[0.5 0.5 0.5],'EdgeColor','k','LineWidth',2);  % SD rectangle
	ax=gca;
	e1=errorbar(WakeWTSDRelmean,WakeWTSDRelstd,'Color','k','LineWidth',mylinewidth,'Marker','o','MarkerFaceColor','k','MarkerSize',mymarkersize,'CapSize',0);
	e2=errorbar(WakeMutSDRelmean,WakeMutSDRelstd,'Color','r','LineWidth',mylinewidth,'Marker','o','MarkerFaceColor','r','MarkerSize',mymarkersize,'CapSize',0);
	hold off 
	set(gca, 'Layer','top')
	ylabel({'WAKE','Time in WAKE (% TRT)'});
	xlabel('Time (hours)')
	l = legend([e1 e2],LegendLabels);

	l.Box = 'off';
	l.Location = 'None';
	l.Position = [0.7326 0.1649 0.1634 0.0798];
	l.FontSize = 12;
	set(gca, 'XTick',[0 6 12 18 24]); %[1 4 7 10 13 16 19 22]
	set(gca, 'XTickLabels',{'0' '6' '12' '18' '24'}); %{'1' '4' '7' '10' '13' '16' '19' '22'}
	xlim([0 24])
	ylim([-5 100])
	ax.FontSize = myfontsize;
	ax.LineWidth = mylinewidth;
	% -- Add stars if sig diffs (LP)
	if PostHoc_bool.SDFirst12hr.W
		hold on 
		
		rectangle('Position',[0 ax.YLim(2)-5 12 5],'FaceColor',Sig12hrColor,'EdgeColor','none'); 
		bool_xlocsSD_LP_sig_values = ismember(e1.XData,xlocs_sig_effects.SD.WT_vs_Mut_LP_W);
		largest_y1_vals = e1.YData(bool_xlocsSD_LP_sig_values) + e1.YPositiveDelta(bool_xlocsSD_LP_sig_values);
		largest_y2_vals = e2.YData(bool_xlocsSD_LP_sig_values) + e2.YPositiveDelta(bool_xlocsSD_LP_sig_values);
		ylocs_sig_effects = max([largest_y1_vals; largest_y2_vals]) + 10*ones(1,length(largest_y1_vals));
		ylocs_sig_effects(ylocs_sig_effects>ax.YLim(2))=ax.YLim(2);
		text(xlocs_sig_effects.SD.WT_vs_Mut_LP_W, ylocs_sig_effects, '*', 'FontSize',astr_size,'HorizontalAlignment','center','VerticalAlignment','middle');
		hold off 
	end 
	% -- Add stars if sig diffs (DP)
	if PostHoc_bool.SDLast12hr.W
		hold on 
		
		rectangle('Position',[12 ax.YLim(2)-5 12 5],'FaceColor',Sig12hrColor,'EdgeColor','none'); 
		bool_xlocsSD_DP_sig_values = ismember(e1.XData,xlocs_sig_effects.SD.WT_vs_Mut_DP_W);
		largest_y1_vals = e1.YData(bool_xlocsSD_DP_sig_values) + e1.YPositiveDelta(bool_xlocsSD_DP_sig_values);
		largest_y2_vals = e2.YData(bool_xlocsSD_DP_sig_values) + e2.YPositiveDelta(bool_xlocsSD_DP_sig_values);
		ylocs_sig_effects = max([largest_y1_vals; largest_y2_vals]) + 10*ones(1,length(largest_y1_vals));
		ylocs_sig_effects(ylocs_sig_effects>ax.YLim(2))=ax.YLim(2);
		text(xlocs_sig_effects.SD.WT_vs_Mut_DP_W, ylocs_sig_effects, '*', 'FontSize',astr_size,'HorizontalAlignment','center','VerticalAlignment','middle');
		hold off 
	end 

	% -- Add stars if sig diffs (First 6 hours of DP)
	if isfield(PostHoc_bool,'SDFirst6hrDP') & PostHoc_bool.SDFirst6hrDP.W
		hold on 
		
		rectangle('Position',[12 ax.YLim(2)-5 6 5],'FaceColor',[0.9 0.9 0.9],'EdgeColor','none'); 
		bool_xlocsSD_First6hrDP_sig_values = ismember(e1.XData,xlocs_sig_effects.SD.WT_vs_Mut_First6hrDP_W);
		largest_y1_vals = e1.YData(bool_xlocsSD_First6hrDP_sig_values) + e1.YPositiveDelta(bool_xlocsSD_First6hrDP_sig_values);
		largest_y2_vals = e2.YData(bool_xlocsSD_First6hrDP_sig_values) + e2.YPositiveDelta(bool_xlocsSD_First6hrDP_sig_values);
		ylocs_sig_effects = max([largest_y1_vals; largest_y2_vals]) + 10*ones(1,length(largest_y1_vals));
		ylocs_sig_effects(ylocs_sig_effects>ax.YLim(2))=ax.YLim(2);
		text(xlocs_sig_effects.SD.WT_vs_Mut_First6hrDP_W, ylocs_sig_effects, '*', 'FontSize',astr_size,'HorizontalAlignment','center','VerticalAlignment','middle','Color',[0.5 0.5 0.5]);
		hold off 
	end 
	% -- Add stars if sig diffs (Last 6 hours of DP)
	if isfield(PostHoc_bool,'SDLast6hrDP') & PostHoc_bool.SDLast6hrDP.W
		hold on 
		
		rectangle('Position',[18 ax.YLim(2)-5 6 5],'FaceColor',[0.9 0.9 0.9],'EdgeColor','none'); 
		bool_xlocsSD_Last6hrDP_sig_values = ismember(e1.XData,xlocs_sig_effects.SD.WT_vs_Mut_Last6hrDP_W);
		largest_y1_vals = e1.YData(bool_xlocsSD_Last6hrDP_sig_values) + e1.YPositiveDelta(bool_xlocsSD_Last6hrDP_sig_values);
		largest_y2_vals = e2.YData(bool_xlocsSD_Last6hrDP_sig_values) + e2.YPositiveDelta(bool_xlocsSD_Last6hrDP_sig_values);
		ylocs_sig_effects = max([largest_y1_vals; largest_y2_vals]) + 10*ones(1,length(largest_y1_vals));
		ylocs_sig_effects(ylocs_sig_effects>ax.YLim(2))=ax.YLim(2);
		text(xlocs_sig_effects.SD.WT_vs_Mut_Last6hrDP_W, ylocs_sig_effects, '*', 'FontSize',astr_size,'HorizontalAlignment','center','VerticalAlignment','middle','Color',[0.5 0.5 0.5]);
		hold off 
	end 

	% -- SD NREM ---

	H=figure;
	H.Name = 'SD NREM relative to Total Recording Time';
	
	rectangle('Position',LPboxpos,'FaceColor','w','EdgeColor','k','LineWidth',2);     % LP rectangle
	hold on
	rectangle('Position',DPboxpos,'FaceColor','k','EdgeColor','k','LineWidth',2);    % DP rectangle
	rectangle('Position',[0 -5 SD_length_hrs 5],'FaceColor',[0.5 0.5 0.5],'EdgeColor','k','LineWidth',2);  % SD rectangle

	e1=errorbar(NREMMutSDRelmean,NREMMutSDRelstd,'Color','r','LineWidth',mylinewidth,'Marker','o','MarkerFaceColor','r','MarkerSize',mymarkersize,'CapSize',0);
	e2=errorbar(NREMWTSDRelmean,NREMWTSDRelstd,'Color','k','LineWidth',mylinewidth,'Marker','o','MarkerFaceColor','k','MarkerSize',mymarkersize,'CapSize',0);
	hold off
	set(gca, 'Layer','top')
	
	set(gca, 'XTick',[0 6 12 18 24]);  %[1 4 7 10 13 16 19 22]
	set(gca, 'XTickLabels',{'0' '6' '12' '18' '24'}); %{'1' '4' '7' '10' '13' '16' '19' '22'}
	xlim([0 24])
	ylim([-5 100])
	xlabel('Time (hours)')
	ylabel({'NREM','Time in NREM (% TRT)'});
	ax = gca;
	ax.FontSize = myfontsize;
	ax.LineWidth = mylinewidth;
	% -- Add stars if sig diffs (LP)
	if PostHoc_bool.SDFirst12hr.N
		hold on 
		
		rectangle('Position',[0 ax.YLim(2)-5 12 5],'FaceColor',Sig12hrColor,'EdgeColor','none');
		bool_xlocsSD_LP_sig_values = ismember(e1.XData,xlocs_sig_effects.SD.WT_vs_Mut_LP_N);
		largest_y1_vals = e1.YData(bool_xlocsSD_LP_sig_values) + e1.YPositiveDelta(bool_xlocsSD_LP_sig_values);
		largest_y2_vals = e2.YData(bool_xlocsSD_LP_sig_values) + e2.YPositiveDelta(bool_xlocsSD_LP_sig_values);
		ylocs_sig_effects = max([largest_y1_vals; largest_y2_vals]) + 10*ones(1,length(largest_y1_vals));
		ylocs_sig_effects(ylocs_sig_effects>ax.YLim(2))=ax.YLim(2);
		text(xlocs_sig_effects.SD.WT_vs_Mut_LP_N, ylocs_sig_effects, '*', 'FontSize',astr_size,'HorizontalAlignment','center','VerticalAlignment','middle');
		hold off 
	end 
	% -- Add stars if sig diffs (DP)
	if PostHoc_bool.SDLast12hr.N
		hold on 
		
		rectangle('Position',[12 ax.YLim(2)-5 12 5],'FaceColor',Sig12hrColor,'EdgeColor','none');
		bool_xlocsSD_DP_sig_values = ismember(e1.XData,xlocs_sig_effects.SD.WT_vs_Mut_DP_N);
		largest_y1_vals = e1.YData(bool_xlocsSD_DP_sig_values) + e1.YPositiveDelta(bool_xlocsSD_DP_sig_values);
		largest_y2_vals = e2.YData(bool_xlocsSD_DP_sig_values) + e2.YPositiveDelta(bool_xlocsSD_DP_sig_values);
		ylocs_sig_effects = max([largest_y1_vals; largest_y2_vals]) + 10*ones(1,length(largest_y1_vals));
		ylocs_sig_effects(ylocs_sig_effects>ax.YLim(2))=ax.YLim(2);
		text(xlocs_sig_effects.SD.WT_vs_Mut_DP_N, ylocs_sig_effects, '*', 'FontSize',astr_size,'HorizontalAlignment','center','VerticalAlignment','middle');
		hold off 
	end 
	% -- Add stars if sig diffs (First 6 hours of DP)
	if isfield(PostHoc_bool,'SDFirst6hrDP') & PostHoc_bool.SDFirst6hrDP.N
		hold on 
		rectangle('Position',[12 ax.YLim(2)-5 6 5],'FaceColor',[0.9 0.9 0.9],'EdgeColor','none');
		bool_xlocsSD_First6hrDP_sig_values = ismember(e1.XData,xlocs_sig_effects.SD.WT_vs_Mut_First6hrDP_N);
		largest_y1_vals = e1.YData(bool_xlocsSD_First6hrDP_sig_values) + e1.YPositiveDelta(bool_xlocsSD_First6hrDP_sig_values);
		largest_y2_vals = e2.YData(bool_xlocsSD_First6hrDP_sig_values) + e2.YPositiveDelta(bool_xlocsSD_First6hrDP_sig_values);
		ylocs_sig_effects = max([largest_y1_vals; largest_y2_vals]) + 10*ones(1,length(largest_y1_vals));
		ylocs_sig_effects(ylocs_sig_effects>ax.YLim(2))=ax.YLim(2);
		text(xlocs_sig_effects.SD.WT_vs_Mut_First6hrDP_N, ylocs_sig_effects, '*', 'FontSize',astr_size,'HorizontalAlignment','center','VerticalAlignment','middle','Color',[0.5 0.5 0.5]);
		hold off 
	end 
	% -- Add stars if sig diffs (Last 6 hours of DP)
	if isfield(PostHoc_bool,'SDLast6hrDP') & PostHoc_bool.SDLast6hrDP.N
		hold on 
		
		rectangle('Position',[18 ax.YLim(2)-5 6 5],'FaceColor',[0.9 0.9 0.9],'EdgeColor','none');
		bool_xlocsSD_Last6hrDP_sig_values = ismember(e1.XData,xlocs_sig_effects.SD.WT_vs_Mut_Last6hrDP_N);
		largest_y1_vals = e1.YData(bool_xlocsSD_Last6hrDP_sig_values) + e1.YPositiveDelta(bool_xlocsSD_Last6hrDP_sig_values);
		largest_y2_vals = e2.YData(bool_xlocsSD_Last6hrDP_sig_values) + e2.YPositiveDelta(bool_xlocsSD_Last6hrDP_sig_values);
		ylocs_sig_effects = max([largest_y1_vals; largest_y2_vals]) + 10*ones(1,length(largest_y1_vals));
		ylocs_sig_effects(ylocs_sig_effects>ax.YLim(2))=ax.YLim(2);
		text(xlocs_sig_effects.SD.WT_vs_Mut_Last6hrDP_N, ylocs_sig_effects, '*', 'FontSize',astr_size,'HorizontalAlignment','center','VerticalAlignment','middle','Color',[0.5 0.5 0.5]);
		hold off 
	end 

	% -- SD REM ---
	H=figure;
	H.Name = 'SD REM relative to Total Recording Time';
	
	rectangle('Position',LPboxposREM,'FaceColor','w','EdgeColor','k','LineWidth',2);     % LP rectangle
	hold on
	rectangle('Position',DPboxposREM,'FaceColor','k','EdgeColor','k','LineWidth',2);    % DP rectangle
	rectangle('Position',[0 -1 SD_length_hrs 1],'FaceColor',[0.5 0.5 0.5],'EdgeColor','k','LineWidth',2);  % SD rectangle

	e1=errorbar(REMMutSDRelmean,REMMutSDRelstd,'Color','r','LineWidth',mylinewidth,'Marker','o','MarkerFaceColor','r','MarkerSize',mymarkersize,'CapSize',0);
	e2=errorbar(REMWTSDRelmean,REMWTSDRelstd,'Color','k','LineWidth',mylinewidth,'Marker','o','MarkerFaceColor','k','MarkerSize',mymarkersize,'CapSize',0);
	hold off 

	set(gca, 'Layer','top')
	
	set(gca, 'XTick',[0 6 12 18 24]);   %[1 4 7 10 13 16 19 22]
	set(gca, 'XTickLabels',{'0' '6' '12' '18' '24'}); %{'1' '4' '7' '10' '13' '16' '19' '22'}
	xlim([0 24])
	ylim([-1 30])
	xlabel('Time (hours)')
	ylabel({'REM','Time in REM (% TRT)'});
	ax = gca;
	ax.FontSize = myfontsize;
	ax.LineWidth = mylinewidth;
	% -- Add stars if sig diffs (LP)
	if PostHoc_bool.SDFirst12hr.R
		hold on 
		
		rectangle('Position',[0 ax.YLim(2)-1 12 1],'FaceColor',Sig12hrColor,'EdgeColor','none');
		bool_xlocsSD_LP_sig_values = ismember(e1.XData,xlocs_sig_effects.SD.WT_vs_Mut_LP_R);
		largest_y1_vals = e1.YData(bool_xlocsSD_LP_sig_values) + e1.YPositiveDelta(bool_xlocsSD_LP_sig_values);
		largest_y2_vals = e2.YData(bool_xlocsSD_LP_sig_values) + e2.YPositiveDelta(bool_xlocsSD_LP_sig_values);
		ylocs_sig_effects = max([largest_y1_vals; largest_y2_vals]) + 10*ones(1,length(largest_y1_vals));
		ylocs_sig_effects(ylocs_sig_effects>ax.YLim(2))=ax.YLim(2);
		text(xlocs_sig_effects.SD.WT_vs_Mut_LP_R, ylocs_sig_effects, '*', 'FontSize',astr_size,'HorizontalAlignment','center','VerticalAlignment','middle');
		hold off 
	end 
	% -- Add stars if sig diffs (DP)
	if PostHoc_bool.SDLast12hr.R
		hold on 
		
		rectangle('Position',[12 ax.YLim(2)-1 12 1],'FaceColor',Sig12hrColor,'EdgeColor','none');
		bool_xlocsSD_DP_sig_values = ismember(e1.XData,xlocs_sig_effects.SD.WT_vs_Mut_DP_R);
		largest_y1_vals = e1.YData(bool_xlocsSD_DP_sig_values) + e1.YPositiveDelta(bool_xlocsSD_DP_sig_values);
		largest_y2_vals = e2.YData(bool_xlocsSD_DP_sig_values) + e2.YPositiveDelta(bool_xlocsSD_DP_sig_values);
		ylocs_sig_effects = max([largest_y1_vals; largest_y2_vals]) + 10*ones(1,length(largest_y1_vals));
		ylocs_sig_effects(ylocs_sig_effects>ax.YLim(2))=ax.YLim(2);
		text(xlocs_sig_effects.SD.WT_vs_Mut_DP_R, ylocs_sig_effects, '*', 'FontSize',astr_size,'HorizontalAlignment','center','VerticalAlignment','middle');
		hold off 
	end 
	% -- Add stars if sig diffs (First 6 hours of DP)
	if isfield(PostHoc_bool,'SDFirst6hrDP') & PostHoc_bool.SDFirst6hrDP.R
		hold on 
		
		rectangle('Position',[12 ax.YLim(2)-1 6 1],'FaceColor',[0.9 0.9 0.9],'EdgeColor','none');
		bool_xlocsSD_First6hrDP_sig_values = ismember(e1.XData,xlocs_sig_effects.SD.WT_vs_Mut_First6hrDP_R);
		largest_y1_vals = e1.YData(bool_xlocsSD_First6hrDP_sig_values) + e1.YPositiveDelta(bool_xlocsSD_First6hrDP_sig_values);
		largest_y2_vals = e2.YData(bool_xlocsSD_First6hrDP_sig_values) + e2.YPositiveDelta(bool_xlocsSD_First6hrDP_sig_values);
		ylocs_sig_effects = max([largest_y1_vals; largest_y2_vals]) + 10*ones(1,length(largest_y1_vals));
		ylocs_sig_effects(ylocs_sig_effects>ax.YLim(2))=ax.YLim(2);
		text(xlocs_sig_effects.SD.WT_vs_Mut_First6hrDP_R, ylocs_sig_effects, '*', 'FontSize',astr_size,'HorizontalAlignment','center','VerticalAlignment','middle','Color',[0.5 0.5 0.5]);
		hold off 
	end 
	% -- Add stars if sig diffs (Last 6 hours of DP)
	if isfield(PostHoc_bool,'SDFirst6hrDP') & PostHoc_bool.SDLast6hrDP.R
		hold on 
		
		rectangle('Position',[18 ax.YLim(2)-1 6 1],'FaceColor',[0.9 0.9 0.9],'EdgeColor','none');
		bool_xlocsSD_Last6hrDP_sig_values = ismember(e1.XData,xlocs_sig_effects.SD.WT_vs_Mut_Last6hrDP_R);
		largest_y1_vals = e1.YData(bool_xlocsSD_Last6hrDP_sig_values) + e1.YPositiveDelta(bool_xlocsSD_Last6hrDP_sig_values);
		largest_y2_vals = e2.YData(bool_xlocsSD_Last6hrDP_sig_values) + e2.YPositiveDelta(bool_xlocsSD_Last6hrDP_sig_values);
		ylocs_sig_effects = max([largest_y1_vals; largest_y2_vals]) + 10*ones(1,length(largest_y1_vals));
		ylocs_sig_effects(ylocs_sig_effects>ax.YLim(2))=ax.YLim(2);
		text(xlocs_sig_effects.SD.WT_vs_Mut_Last6hrDP_R, ylocs_sig_effects, '*', 'FontSize',astr_size,'HorizontalAlignment','center','VerticalAlignment','middle','Color',[0.5 0.5 0.5]);
		hold off 
	end 
 
end

% ------ 3-panel figure for all of SD ---------

H=figure;
H.Name = ['Three-Panel SD Time in State Figure ',Sex];
H.Renderer = 'painters';
H.Position = [429 118 678 698];

t= tiledlayout(3,1,'TileSpacing','compact','Padding','compact','TileIndexing', 'columnmajor');
nexttile
%-- SD W as part of 3-panel --
rectangle('Position',LPboxpos,'FaceColor','w','EdgeColor','k','LineWidth',2);     % LP rectangle
hold on
rectangle('Position',DPboxpos,'FaceColor','k','EdgeColor','k','LineWidth',2);    % DP rectangle
rectangle('Position',[0 -5 SD_length_hrs 5],'FaceColor',[0.5 0.5 0.5],'EdgeColor','k','LineWidth',2);  % SD rectangle

e1=errorbar(WakeWTSDRelmean,WakeWTSDRelstd,'Color','k','LineWidth',mylinewidth,'Marker','o','MarkerFaceColor','k','MarkerSize',mymarkersize,'CapSize',0);
e2=errorbar(WakeMutSDRelmean,WakeMutSDRelstd,'Color','r','LineWidth',mylinewidth,'Marker','o','MarkerFaceColor','r','MarkerSize',mymarkersize,'CapSize',0);
hold off 
set(gca, 'Layer','top')
ylabel({'WAKE','Time in WAKE (% TRT)'});
xlabel('Time (hours)')
xlabel('')
l = legend([e1 e2],LegendLabels);
l.Box = 'off';
l.Location = 'SouthEast';
l.FontSize = 12;
set(gca,'XTick',[])
set(gca, 'XTickLabels',{'' '' '' '' '' '' '' ''})

xlim([0 24])
ylim([-5 100])
ax = gca;
ax.FontSize = myfontsize;
ax.LineWidth = mylinewidth;
% -- Add stars if sig diffs (LP)
if PostHoc_bool.SDFirst12hr.W
	hold on 
	rectangle('Position',[0 ax.YLim(2)-5 12 5],'FaceColor',Sig12hrColor,'EdgeColor','none');
	bool_xlocsSD_LP_sig_values = ismember(e1.XData,xlocs_sig_effects.SD.WT_vs_Mut_LP_W);
	largest_y1_vals = e1.YData(bool_xlocsSD_LP_sig_values) + e1.YPositiveDelta(bool_xlocsSD_LP_sig_values);
	largest_y2_vals = e2.YData(bool_xlocsSD_LP_sig_values) + e2.YPositiveDelta(bool_xlocsSD_LP_sig_values);
	ylocs_sig_effects = max([largest_y1_vals; largest_y2_vals]) + 10*ones(1,length(largest_y1_vals));
	ylocs_sig_effects(ylocs_sig_effects>ax.YLim(2))=ax.YLim(2);
	text(xlocs_sig_effects.SD.WT_vs_Mut_LP_W, ylocs_sig_effects, '*', 'FontSize',astr_size,'HorizontalAlignment','center','VerticalAlignment','middle');
	hold off 
end 
% -- Add stars if sig diffs (DP)
if PostHoc_bool.SDLast12hr.W
	hold on 
	rectangle('Position',[12 ax.YLim(2)-5 12 5],'FaceColor',Sig12hrColor,'EdgeColor','none');
	bool_xlocsSD_DP_sig_values = ismember(e1.XData,xlocs_sig_effects.SD.WT_vs_Mut_DP_W);
	largest_y1_vals = e1.YData(bool_xlocsSD_DP_sig_values) + e1.YPositiveDelta(bool_xlocsSD_DP_sig_values);
	largest_y2_vals = e2.YData(bool_xlocsSD_DP_sig_values) + e2.YPositiveDelta(bool_xlocsSD_DP_sig_values);
	ylocs_sig_effects = max([largest_y1_vals; largest_y2_vals]) + 10*ones(1,length(largest_y1_vals));
	ylocs_sig_effects(ylocs_sig_effects>ax.YLim(2))=ax.YLim(2);
	text(xlocs_sig_effects.SD.WT_vs_Mut_DP_W, ylocs_sig_effects, '*', 'FontSize',astr_size,'HorizontalAlignment','center','VerticalAlignment','middle');
	hold off 
end 

% -- Add stars if sig diffs (First 6 hours of DP)
if isfield(PostHoc_bool,'SDFirst6hrDP') & PostHoc_bool.SDFirst6hrDP.W
	hold on 
	rectangle('Position',[12 ax.YLim(2)-5 6 5],'FaceColor',[0.9 0.9 0.9],'EdgeColor','none');
	bool_xlocsSD_First6hrDP_sig_values = ismember(e1.XData,xlocs_sig_effects.SD.WT_vs_Mut_First6hrDP_W);
	largest_y1_vals = e1.YData(bool_xlocsSD_First6hrDP_sig_values) + e1.YPositiveDelta(bool_xlocsSD_First6hrDP_sig_values);
	largest_y2_vals = e2.YData(bool_xlocsSD_First6hrDP_sig_values) + e2.YPositiveDelta(bool_xlocsSD_First6hrDP_sig_values);
	ylocs_sig_effects = max([largest_y1_vals; largest_y2_vals]) + 10*ones(1,length(largest_y1_vals));
	ylocs_sig_effects(ylocs_sig_effects>ax.YLim(2))=ax.YLim(2);
	text(xlocs_sig_effects.SD.WT_vs_Mut_First6hrDP_W, ylocs_sig_effects, '*', 'FontSize',astr_size,'HorizontalAlignment','center','VerticalAlignment','middle','Color',[0.5 0.5 0.5]);
	hold off 
end 
% -- Add stars if sig diffs (Last 6 hours of DP)
if isfield(PostHoc_bool,'SDLast6hrDP') & PostHoc_bool.SDLast6hrDP.W
	hold on 
	rectangle('Position',[18 ax.YLim(2)-5 6 5],'FaceColor',[0.9 0.9 0.9],'EdgeColor','none');
	bool_xlocsSD_Last6hrDP_sig_values = ismember(e1.XData,xlocs_sig_effects.SD.WT_vs_Mut_Last6hrDP_W);
	largest_y1_vals = e1.YData(bool_xlocsSD_Last6hrDP_sig_values) + e1.YPositiveDelta(bool_xlocsSD_Last6hrDP_sig_values);
	largest_y2_vals = e2.YData(bool_xlocsSD_Last6hrDP_sig_values) + e2.YPositiveDelta(bool_xlocsSD_Last6hrDP_sig_values);
	ylocs_sig_effects = max([largest_y1_vals; largest_y2_vals]) + 10*ones(1,length(largest_y1_vals));
	ylocs_sig_effects(ylocs_sig_effects>ax.YLim(2))=ax.YLim(2);
	text(xlocs_sig_effects.SD.WT_vs_Mut_Last6hrDP_W, ylocs_sig_effects, '*', 'FontSize',astr_size,'HorizontalAlignment','center','VerticalAlignment','middle','Color',[0.5 0.5 0.5]);
	hold off 
end 

%-- SD NREM as part of 3-panel --
nexttile 
rectangle('Position',LPboxpos,'FaceColor','w','EdgeColor','k','LineWidth',2);     % LP rectangle
hold on
rectangle('Position',DPboxpos,'FaceColor','k','EdgeColor','k','LineWidth',2);    % DP rectangle
rectangle('Position',[0 -5 SD_length_hrs 5],'FaceColor',[0.5 0.5 0.5],'EdgeColor','k','LineWidth',2);  % SD rectangle

e1=errorbar(NREMMutSDRelmean,NREMMutSDRelstd,'Color','r','LineWidth',mylinewidth,'Marker','o','MarkerFaceColor','r','MarkerSize',mymarkersize,'CapSize',0);
e2=errorbar(NREMWTSDRelmean,NREMWTSDRelstd,'Color','k','LineWidth',mylinewidth,'Marker','o','MarkerFaceColor','k','MarkerSize',mymarkersize,'CapSize',0);
hold off
set(gca, 'Layer','top')

set(gca,'XTick',[])

xlim([0 24])
ylim([-5 100])
xlabel('')
ylabel({'NREM','Time in NREM (% TRT)'});
ax = gca;
ax.FontSize = myfontsize;
ax.LineWidth = mylinewidth;

% -- Add stars if sig diffs (LP)
if PostHoc_bool.SDFirst12hr.N
	hold on 
	rectangle('Position',[0 ax.YLim(2)-5 12 5],'FaceColor',Sig12hrColor,'EdgeColor','none');
	bool_xlocsSD_LP_sig_values = ismember(e1.XData,xlocs_sig_effects.SD.WT_vs_Mut_LP_N);
	largest_y1_vals = e1.YData(bool_xlocsSD_LP_sig_values) + e1.YPositiveDelta(bool_xlocsSD_LP_sig_values);
	largest_y2_vals = e2.YData(bool_xlocsSD_LP_sig_values) + e2.YPositiveDelta(bool_xlocsSD_LP_sig_values);
	ylocs_sig_effects = max([largest_y1_vals; largest_y2_vals]) + 10*ones(1,length(largest_y1_vals));
	ylocs_sig_effects(ylocs_sig_effects>ax.YLim(2))=ax.YLim(2);
	text(xlocs_sig_effects.SD.WT_vs_Mut_LP_N, ylocs_sig_effects, '*', 'FontSize',astr_size,'HorizontalAlignment','center','VerticalAlignment','middle');
	hold off 
end 
% -- Add stars if sig diffs (DP)
if PostHoc_bool.SDLast12hr.N
	hold on 
	rectangle('Position',[12 ax.YLim(2)-5 12 5],'FaceColor',Sig12hrColor,'EdgeColor','none');
	bool_xlocsSD_DP_sig_values = ismember(e1.XData,xlocs_sig_effects.SD.WT_vs_Mut_DP_N);
	largest_y1_vals = e1.YData(bool_xlocsSD_DP_sig_values) + e1.YPositiveDelta(bool_xlocsSD_DP_sig_values);
	largest_y2_vals = e2.YData(bool_xlocsSD_DP_sig_values) + e2.YPositiveDelta(bool_xlocsSD_DP_sig_values);
	ylocs_sig_effects = max([largest_y1_vals; largest_y2_vals]) + 10*ones(1,length(largest_y1_vals));
	ylocs_sig_effects(ylocs_sig_effects>ax.YLim(2))=ax.YLim(2);
	text(xlocs_sig_effects.SD.WT_vs_Mut_DP_N, ylocs_sig_effects, '*', 'FontSize',astr_size,'HorizontalAlignment','center','VerticalAlignment','middle');
	hold off 
end 
% -- Add stars if sig diffs (First 6 hours of DP)
if isfield(PostHoc_bool,'SDFirst6hrDP') & PostHoc_bool.SDFirst6hrDP.N
	hold on 
	rectangle('Position',[12 ax.YLim(2)-5 6 5],'FaceColor',[0.9 0.9 0.9],'EdgeColor','none');
	bool_xlocsSD_First6hrDP_sig_values = ismember(e1.XData,xlocs_sig_effects.SD.WT_vs_Mut_First6hrDP_N);
	largest_y1_vals = e1.YData(bool_xlocsSD_First6hrDP_sig_values) + e1.YPositiveDelta(bool_xlocsSD_First6hrDP_sig_values);
	largest_y2_vals = e2.YData(bool_xlocsSD_First6hrDP_sig_values) + e2.YPositiveDelta(bool_xlocsSD_First6hrDP_sig_values);
	ylocs_sig_effects = max([largest_y1_vals; largest_y2_vals]) + 10*ones(1,length(largest_y1_vals));
	ylocs_sig_effects(ylocs_sig_effects>ax.YLim(2))=ax.YLim(2);
	text(xlocs_sig_effects.SD.WT_vs_Mut_First6hrDP_N, ylocs_sig_effects, '*', 'FontSize',astr_size,'HorizontalAlignment','center','VerticalAlignment','middle','Color',[0.5 0.5 0.5]);
	hold off 
end 
% -- Add stars if sig diffs (Last 6 hours of DP)
if isfield(PostHoc_bool,'SDLast6hrDP') & PostHoc_bool.SDLast6hrDP.N
	hold on 
	rectangle('Position',[18 ax.YLim(2)-5 6 5],'FaceColor',[0.9 0.9 0.9],'EdgeColor','none');
	bool_xlocsSD_Last6hrDP_sig_values = ismember(e1.XData,xlocs_sig_effects.SD.WT_vs_Mut_Last6hrDP_N);
	largest_y1_vals = e1.YData(bool_xlocsSD_Last6hrDP_sig_values) + e1.YPositiveDelta(bool_xlocsSD_Last6hrDP_sig_values);
	largest_y2_vals = e2.YData(bool_xlocsSD_Last6hrDP_sig_values) + e2.YPositiveDelta(bool_xlocsSD_Last6hrDP_sig_values);
	ylocs_sig_effects = max([largest_y1_vals; largest_y2_vals]) + 10*ones(1,length(largest_y1_vals));
	ylocs_sig_effects(ylocs_sig_effects>ax.YLim(2))=ax.YLim(2);
	text(xlocs_sig_effects.SD.WT_vs_Mut_Last6hrDP_N, ylocs_sig_effects, '*', 'FontSize',astr_size,'HorizontalAlignment','center','VerticalAlignment','middle','Color',[0.5 0.5 0.5]);
	hold off 
end 

% -- SD REM as part of 3-panel --
nexttile
rectangle('Position',LPboxposREM,'FaceColor','w','EdgeColor','k','LineWidth',2);     % LP rectangle
hold on
rectangle('Position',DPboxposREM,'FaceColor','k','EdgeColor','k','LineWidth',2);    % DP rectangle
rectangle('Position',[0 -1 SD_length_hrs 1],'FaceColor',[0.5 0.5 0.5],'EdgeColor','k','LineWidth',2);  % SD rectangle

e1=errorbar(REMMutSDRelmean,REMMutSDRelstd,'Color','r','LineWidth',mylinewidth,'Marker','o','MarkerFaceColor','r','MarkerSize',mymarkersize,'CapSize',0);
e2=errorbar(REMWTSDRelmean,REMWTSDRelstd,'Color','k','LineWidth',mylinewidth,'Marker','o','MarkerFaceColor','k','MarkerSize',mymarkersize,'CapSize',0);
hold off 

set(gca, 'Layer','top')

set(gca, 'XTick',[0 6 12 18 24]); %[1 4 7 10 13 16 19 22]
set(gca, 'XTickLabels',{'0' '6' '12' '18' '24'});  %{'1' '4' '7' '10' '13' '16' '19' '22'}
xlim([0 24])
ylim([-1 30])
xlabel('Time (hours)')
ylabel({'REM','Time in REM (% TRT)'});
ax = gca;
ax.FontSize = myfontsize;
ax.LineWidth = mylinewidth;

% -- Add stars if sig diffs (LP)
if PostHoc_bool.SDFirst12hr.R
	hold on 
	rectangle('Position',[0 ax.YLim(2)-1 12 1],'FaceColor',Sig12hrColor,'EdgeColor','none');
	bool_xlocsSD_LP_sig_values = ismember(e1.XData,xlocs_sig_effects.SD.WT_vs_Mut_LP_R);
	largest_y1_vals = e1.YData(bool_xlocsSD_LP_sig_values) + e1.YPositiveDelta(bool_xlocsSD_LP_sig_values);
	largest_y2_vals = e2.YData(bool_xlocsSD_LP_sig_values) + e2.YPositiveDelta(bool_xlocsSD_LP_sig_values);
	ylocs_sig_effects = max([largest_y1_vals; largest_y2_vals]) + 10*ones(1,length(largest_y1_vals));
	ylocs_sig_effects(ylocs_sig_effects>ax.YLim(2))=ax.YLim(2);
	text(xlocs_sig_effects.SD.WT_vs_Mut_LP_R, ylocs_sig_effects, '*', 'FontSize',astr_size,'HorizontalAlignment','center','VerticalAlignment','middle');
	hold off 
end 
% -- Add stars if sig diffs (DP)
if PostHoc_bool.SDLast12hr.R
	hold on 
	rectangle('Position',[12 ax.YLim(2)-1 12 1],'FaceColor',Sig12hrColor,'EdgeColor','none');
	bool_xlocsSD_DP_sig_values = ismember(e1.XData,xlocs_sig_effects.SD.WT_vs_Mut_DP_R);
	largest_y1_vals = e1.YData(bool_xlocsSD_DP_sig_values) + e1.YPositiveDelta(bool_xlocsSD_DP_sig_values);
	largest_y2_vals = e2.YData(bool_xlocsSD_DP_sig_values) + e2.YPositiveDelta(bool_xlocsSD_DP_sig_values);
	ylocs_sig_effects = max([largest_y1_vals; largest_y2_vals]) + 10*ones(1,length(largest_y1_vals));
	ylocs_sig_effects(ylocs_sig_effects>ax.YLim(2))=ax.YLim(2);
	text(xlocs_sig_effects.SD.WT_vs_Mut_DP_R, ylocs_sig_effects, '*', 'FontSize',astr_size,'HorizontalAlignment','center','VerticalAlignment','middle');
	hold off 
end 
% -- Add stars if sig diffs (First 6 hours of DP)
if isfield(PostHoc_bool,'SDFirst6hrDP') & PostHoc_bool.SDFirst6hrDP.R
	hold on 
	rectangle('Position',[12 ax.YLim(2)-1 6 1],'FaceColor',[0.9 0.9 0.9],'EdgeColor','none');
	bool_xlocsSD_First6hrDP_sig_values = ismember(e1.XData,xlocs_sig_effects.SD.WT_vs_Mut_First6hrDP_R);
	largest_y1_vals = e1.YData(bool_xlocsSD_First6hrDP_sig_values) + e1.YPositiveDelta(bool_xlocsSD_First6hrDP_sig_values);
	largest_y2_vals = e2.YData(bool_xlocsSD_First6hrDP_sig_values) + e2.YPositiveDelta(bool_xlocsSD_First6hrDP_sig_values);
	ylocs_sig_effects = max([largest_y1_vals; largest_y2_vals]) + 10*ones(1,length(largest_y1_vals));
	ylocs_sig_effects(ylocs_sig_effects>ax.YLim(2))=ax.YLim(2);
	text(xlocs_sig_effects.SD.WT_vs_Mut_First6hrDP_R, ylocs_sig_effects, '*', 'FontSize',astr_size,'HorizontalAlignment','center','VerticalAlignment','middle','Color',[0.5 0.5 0.5]);
	hold off 
end 
% -- Add stars if sig diffs (Last 6 hours of DP)
if isfield(PostHoc_bool,'SDLast6hrDP') & PostHoc_bool.SDLast6hrDP.R
	hold on 
	rectangle('Position',[18 ax.YLim(2)-1 6 1],'FaceColor',[0.9 0.9 0.9],'EdgeColor','none');
	bool_xlocsSD_Last6hrDP_sig_values = ismember(e1.XData,xlocs_sig_effects.SD.WT_vs_Mut_Last6hrDP_R);
	largest_y1_vals = e1.YData(bool_xlocsSD_Last6hrDP_sig_values) + e1.YPositiveDelta(bool_xlocsSD_Last6hrDP_sig_values);
	largest_y2_vals = e2.YData(bool_xlocsSD_Last6hrDP_sig_values) + e2.YPositiveDelta(bool_xlocsSD_Last6hrDP_sig_values);
	ylocs_sig_effects = max([largest_y1_vals; largest_y2_vals]) + 10*ones(1,length(largest_y1_vals));
	ylocs_sig_effects(ylocs_sig_effects>ax.YLim(2))=ax.YLim(2);
	text(xlocs_sig_effects.SD.WT_vs_Mut_Last6hrDP_R, ylocs_sig_effects, '*', 'FontSize',astr_size,'HorizontalAlignment','center','VerticalAlignment','middle','Color',[0.5 0.5 0.5]);
	hold off 
end 
% ---------------------------------------------------------------------------------------
% ---------------------------------------------------------------------------------------


% ---- Set up the Post hoc tables to return ----
if ~isempty(mc_table.BL.First12hr.Wake)  PostHoc_tables.BL.First12hr.Wake  = mc_table.BL.First12hr.Wake;  else PostHoc_tables.BL.First12hr.Wake  = []; end  
if ~isempty(mc_table.BL.First12hr.NREM)  PostHoc_tables.BL.First12hr.NREM  = mc_table.BL.First12hr.NREM;  else PostHoc_tables.BL.First12hr.NREM  = []; end  
if ~isempty(mc_table.BL.First12hr.REM)   PostHoc_tables.BL.First12hr.REM   = mc_table.BL.First12hr.REM;   else PostHoc_tables.BL.First12hr.REM   = []; end
if ~isempty(mc_table.BL.Last12hr.Wake)   PostHoc_tables.BL.Last12hr.Wake   = mc_table.BL.Last12hr.Wake;   else PostHoc_tables.BL.Last12hr.Wake   = []; end
if ~isempty(mc_table.BL.Last12hr.NREM)   PostHoc_tables.BL.Last12hr.NREM   = mc_table.BL.Last12hr.NREM;   else PostHoc_tables.BL.Last12hr.NREM   = []; end
if ~isempty(mc_table.BL.Last12hr.REM)    PostHoc_tables.BL.Last12hr.REM    = mc_table.BL.Last12hr.REM;    else PostHoc_tables.BL.Last12hr.REM    = []; end 
if ~isempty(mc_table.SD.First12hr.Wake)  PostHoc_tables.SD.First12hr.Wake  = mc_table.SD.First12hr.Wake;  else PostHoc_tables.SD.First12hr.Wake  = []; end
if ~isempty(mc_table.SD.First12hr.NREM)  PostHoc_tables.SD.First12hr.NREM  = mc_table.SD.First12hr.NREM;  else PostHoc_tables.SD.First12hr.NREM  = []; end
if ~isempty(mc_table.SD.First12hr.REM)   PostHoc_tables.SD.First12hr.REM   = mc_table.SD.First12hr.REM;   else PostHoc_tables.SD.First12hr.REM   = []; end
if ~isempty(mc_table.SD.Last12hr.Wake)   PostHoc_tables.SD.Last12hr.Wake   = mc_table.SD.Last12hr.Wake;   else PostHoc_tables.SD.Last12hr.Wake   = []; end
if ~isempty(mc_table.SD.Last12hr.NREM)   PostHoc_tables.SD.Last12hr.NREM   = mc_table.SD.Last12hr.NREM;   else PostHoc_tables.SD.Last12hr.NREM   = []; end
if ~isempty(mc_table.SD.Last12hr.REM)    PostHoc_tables.SD.Last12hr.REM    = mc_table.SD.Last12hr.REM;    else PostHoc_tables.SD.Last12hr.REM    = []; end 

% If you are doing 6-hour segments of the DP
% -- BL --
if isfield(mc_table.BL,'First6hrDP')
	if ~isempty(mc_table.BL.First6hrDP.Wake) PostHoc_tables.BL.First6hrDP.Wake = mc_table.BL.First6hrDP.Wake; else PostHoc_tables.BL.First6hrDP.Wake = []; end
	if ~isempty(mc_table.BL.First6hrDP.NREM) PostHoc_tables.BL.First6hrDP.NREM = mc_table.BL.First6hrDP.NREM; else PostHoc_tables.BL.First6hrDP.NREM = []; end
	if ~isempty(mc_table.BL.First6hrDP.REM)  PostHoc_tables.BL.First6hrDP.REM  = mc_table.BL.First6hrDP.REM;  else PostHoc_tables.BL.First6hrDP.REM  = []; end 
end
if isfield(mc_table.BL,'Last6hrDP')
	if ~isempty(mc_table.BL.Last6hrDP.Wake)  PostHoc_tables.BL.Last6hrDP.Wake  = mc_table.BL.Last6hrDP.Wake;  else PostHoc_tables.BL.Last6hrDP.Wake  = []; end
	if ~isempty(mc_table.BL.Last6hrDP.NREM)  PostHoc_tables.BL.Last6hrDP.NREM  = mc_table.BL.Last6hrDP.NREM;  else PostHoc_tables.BL.Last6hrDP.NREM  = []; end
	if ~isempty(mc_table.BL.Last6hrDP.REM)   PostHoc_tables.BL.Last6hrDP.REM   = mc_table.BL.Last6hrDP.REM;   else PostHoc_tables.BL.Last6hrDP.REM   = []; end 
end 

% -- SD --
if isfield(mc_table.SD,'First6hrDP')
	if ~isempty(mc_table.SD.First6hrDP.Wake) PostHoc_tables.SD.First6hrDP.Wake = mc_table.SD.First6hrDP.Wake; else PostHoc_tables.SD.First6hrDP.Wake = []; end
	if ~isempty(mc_table.SD.First6hrDP.NREM) PostHoc_tables.SD.First6hrDP.NREM = mc_table.SD.First6hrDP.NREM; else PostHoc_tables.SD.First6hrDP.NREM = []; end
	if ~isempty(mc_table.SD.First6hrDP.REM)  PostHoc_tables.SD.First6hrDP.REM  = mc_table.SD.First6hrDP.REM;  else PostHoc_tables.SD.First6hrDP.REM  = []; end 
end
if isfield(mc_table.SD,'Last6hrDP')
	if ~isempty(mc_table.SD.Last6hrDP.Wake)  PostHoc_tables.SD.Last6hrDP.Wake  = mc_table.SD.Last6hrDP.Wake;  else PostHoc_tables.SD.Last6hrDP.Wake  = []; end
	if ~isempty(mc_table.SD.Last6hrDP.NREM)  PostHoc_tables.SD.Last6hrDP.NREM  = mc_table.SD.Last6hrDP.NREM;  else PostHoc_tables.SD.Last6hrDP.NREM  = []; end
	if ~isempty(mc_table.SD.Last6hrDP.REM)   PostHoc_tables.SD.Last6hrDP.REM   = mc_table.SD.Last6hrDP.REM;   else PostHoc_tables.SD.Last6hrDP.REM   = []; end 
end 

% Set up p_vals struct to be used in plot_output_var_M_vs_F_WT_vs_Mut
% --- BL ---
% -- LP --
% - W - 
p_vals.Male.BL.LP.Wake.Posthocs.WTvsMut       = RepMeas_ANOVA_tables.BL.First12hr.Wake.pValue('Genotype');
p_vals.Male.BL.LP.Wake.Posthocs.GenotypeXTime = RepMeas_ANOVA_tables.BL.First12hr.Wake.pValue('Genotype:Time');


% - NREM -
p_vals.Male.BL.LP.NREM.Posthocs.WTvsMut        = RepMeas_ANOVA_tables.BL.First12hr.NREM.pValue('Genotype');
p_vals.Male.BL.LP.NREM.Posthocs.GenotypeXTime  = RepMeas_ANOVA_tables.BL.First12hr.NREM.pValue('Genotype:Time');

% - REM - 
p_vals.Male.BL.LP.REM.Posthocs.WTvsMut         = RepMeas_ANOVA_tables.BL.First12hr.REM.pValue('Genotype');
p_vals.Male.BL.LP.REM.Posthocs.GenotypeXTime   = RepMeas_ANOVA_tables.BL.First12hr.REM.pValue('Genotype:Time');

% -- DP --
% - W - 
p_vals.Male.BL.DP.Wake.Posthocs.WTvsMut         = RepMeas_ANOVA_tables.BL.Last12hr.Wake.pValue('Genotype');
p_vals.Male.BL.DP.Wake.Posthocs.GenotypeXTime   = RepMeas_ANOVA_tables.BL.Last12hr.Wake.pValue('Genotype:Time');

% - NREM -
p_vals.Male.BL.DP.NREM.Posthocs.WTvsMut         = RepMeas_ANOVA_tables.BL.Last12hr.NREM.pValue('Genotype');
p_vals.Male.BL.DP.NREM.Posthocs.GenotypeXTime   = RepMeas_ANOVA_tables.BL.Last12hr.NREM.pValue('Genotype:Time');

% - REM - 
p_vals.Male.BL.DP.REM.Posthocs.WTvsMut         = RepMeas_ANOVA_tables.BL.Last12hr.REM.pValue('Genotype');
p_vals.Male.BL.DP.REM.Posthocs.GenotypeXTime   = RepMeas_ANOVA_tables.BL.Last12hr.REM.pValue('Genotype:Time');

% --- SDexcSD ---
% -- LP --
% - W - 
p_vals.Male.SDexcSD.LP.Wake.Posthocs.WTvsMut         = RepMeas_ANOVA_tables.SD.First12hr.Wake.pValue('Genotype');
p_vals.Male.SDexcSD.LP.Wake.Posthocs.GenotypeXTime   = RepMeas_ANOVA_tables.SD.First12hr.Wake.pValue('Genotype:Time');

% - NREM -
p_vals.Male.SDexcSD.LP.NREM.Posthocs.WTvsMut         = RepMeas_ANOVA_tables.SD.First12hr.NREM.pValue('Genotype');
p_vals.Male.SDexcSD.LP.NREM.Posthocs.GenotypeXTime   = RepMeas_ANOVA_tables.SD.First12hr.NREM.pValue('Genotype:Time');

% - REM - 
p_vals.Male.SDexcSD.LP.REM.Posthocs.WTvsMut          = RepMeas_ANOVA_tables.SD.First12hr.REM.pValue('Genotype');
p_vals.Male.SDexcSD.LP.REM.Posthocs.GenotypeXTime    = RepMeas_ANOVA_tables.SD.First12hr.REM.pValue('Genotype:Time');


% --- SD ---
% -- DP --
% - W - 
p_vals.Male.SD.DP.Wake.Posthocs.WTvsMut         = RepMeas_ANOVA_tables.SD.Last12hr.Wake.pValue('Genotype');
p_vals.Male.SD.DP.Wake.Posthocs.GenotypeXTime   = RepMeas_ANOVA_tables.SD.Last12hr.Wake.pValue('Genotype:Time');

% - NREM -
p_vals.Male.SD.DP.NREM.Posthocs.WTvsMut         = RepMeas_ANOVA_tables.SD.Last12hr.NREM.pValue('Genotype');
p_vals.Male.SD.DP.NREM.Posthocs.GenotypeXTime   = RepMeas_ANOVA_tables.SD.Last12hr.NREM.pValue('Genotype:Time');

% - REM - 
p_vals.Male.SD.DP.REM.Posthocs.WTvsMut         = RepMeas_ANOVA_tables.SD.Last12hr.REM.pValue('Genotype');
p_vals.Male.SD.DP.REM.Posthocs.GenotypeXTime   = RepMeas_ANOVA_tables.SD.Last12hr.REM.pValue('Genotype:Time');






































