function [p_vals,PostHoc_tables] = plot_TimeInState_Both_Genders(options)
%
% USAGE: PostHoc_tables = plot_TimeInState_Both_Genders(TimeInState=TimeInState,LegendLabels=LegendLabels,SD_length_hrs=SD_length_hrs,RepMeas_ANOVA_tables=RepMeas_ANOVA_tables_BothMandF_TIS)
% 
% This function simply plots the time in state plots (3 panels, one for each state) and adds asterisks if the repeated measures 
% ANOVAs are significant.  Post-hocs?  




arguments
    options.TimeInState
    %options.SD_length_hrs    = 5
    options.SD_length_hrs    = S.SD_length_hrs
    options.LegendLabels     = []
    options.RepMeas_ANOVA_tables = []
    options.RM_Models            = []
    options.PlotIndividualFigs   = 0
    options.Plot3panelFigs       = 0
    options.Sexes                = []
    options.SexVarName           = []
end 

TimeInState_struct 	 = options.TimeInState;
SD_length_hrs      	 = options.SD_length_hrs;
LegendLabels       	 = options.LegendLabels;
RepMeas_ANOVA_tables = options.RepMeas_ANOVA_tables;
RM_Models            = options.RM_Models;
PlotIndividualFigs   = options.PlotIndividualFigs;
Plot3panelFigs       = options.Plot3panelFigs;
Sexes                = options.Sexes;
SexVarName     	     = options.SexVarName;


% -- Unpack the TimeInState struct so we can use the plotting code below in a less messy way

% These are MxN matrices.  M rows for M hours, N columns for N animals.  The values are the number of minutes per hour spent 
% in that state by that animal.   
% - WT M -
Wake24hWTBL_M      = TimeInState_struct.WT.(Sexes{1}).BL.Wake;
Wake24hWTSD_M      = TimeInState_struct.WT.(Sexes{1}).SD.Wake;
Wake24hWTSDexcSD_M = TimeInState_struct.WT.(Sexes{1}).SDexcSD.Wake;
NREM24hWTBL_M      = TimeInState_struct.WT.(Sexes{1}).BL.NREM;
NREM24hWTSD_M      = TimeInState_struct.WT.(Sexes{1}).SD.NREM;
NREM24hWTSDexcSD_M = TimeInState_struct.WT.(Sexes{1}).SDexcSD.NREM;
REM24hWTBL_M       = TimeInState_struct.WT.(Sexes{1}).BL.REM;
REM24hWTSD_M       = TimeInState_struct.WT.(Sexes{1}).SD.REM;
REM24hWTSDexcSD_M  = TimeInState_struct.WT.(Sexes{1}).SDexcSD.REM;

% - Mut M - 
Wake24hMutBL_M      = TimeInState_struct.Mut.(Sexes{1}).BL.Wake;
Wake24hMutSD_M      = TimeInState_struct.Mut.(Sexes{1}).SD.Wake;
Wake24hMutSDexcSD_M = TimeInState_struct.Mut.(Sexes{1}).SDexcSD.Wake;
NREM24hMutBL_M      = TimeInState_struct.Mut.(Sexes{1}).BL.NREM;
NREM24hMutSD_M      = TimeInState_struct.Mut.(Sexes{1}).SD.NREM;
NREM24hMutSDexcSD_M = TimeInState_struct.Mut.(Sexes{1}).SDexcSD.NREM; 
REM24hMutBL_M       = TimeInState_struct.Mut.(Sexes{1}).BL.REM;
REM24hMutSD_M       = TimeInState_struct.Mut.(Sexes{1}).SD.REM;
REM24hMutSDexcSD_M  = TimeInState_struct.Mut.(Sexes{1}).SDexcSD.REM;

% - WT F -
Wake24hWTBL_F      = TimeInState_struct.WT.(Sexes{2}).BL.Wake;
Wake24hWTSD_F      = TimeInState_struct.WT.(Sexes{2}).SD.Wake;
Wake24hWTSDexcSD_F = TimeInState_struct.WT.(Sexes{2}).SDexcSD.Wake;
NREM24hWTBL_F      = TimeInState_struct.WT.(Sexes{2}).BL.NREM;
NREM24hWTSD_F      = TimeInState_struct.WT.(Sexes{2}).SD.NREM;
NREM24hWTSDexcSD_F = TimeInState_struct.WT.(Sexes{2}).SDexcSD.NREM;
REM24hWTBL_F       = TimeInState_struct.WT.(Sexes{2}).BL.REM;
REM24hWTSD_F       = TimeInState_struct.WT.(Sexes{2}).SD.REM;
REM24hWTSDexcSD_F  = TimeInState_struct.WT.(Sexes{2}).SDexcSD.REM;

% - Mut F - 
Wake24hMutBL_F      = TimeInState_struct.Mut.(Sexes{2}).BL.Wake;
Wake24hMutSD_F      = TimeInState_struct.Mut.(Sexes{2}).SD.Wake;
Wake24hMutSDexcSD_F = TimeInState_struct.Mut.(Sexes{2}).SDexcSD.Wake;
NREM24hMutBL_F      = TimeInState_struct.Mut.(Sexes{2}).BL.NREM;
NREM24hMutSD_F      = TimeInState_struct.Mut.(Sexes{2}).SD.NREM;
NREM24hMutSDexcSD_F = TimeInState_struct.Mut.(Sexes{2}).SDexcSD.NREM; 
REM24hMutBL_F       = TimeInState_struct.Mut.(Sexes{2}).BL.REM;
REM24hMutSD_F       = TimeInState_struct.Mut.(Sexes{2}).SD.REM;
REM24hMutSDexcSD_F  = TimeInState_struct.Mut.(Sexes{2}).SDexcSD.REM;

%% time in stage per hour normalized to recording time (aka, relative to 60min) Male
% --- WT Baseline time in stage per hour Male
WakeWTBLRelmean_M = mean(Wake24hWTBL_M./60*100,2,'omitnan')';
WakeWTBLRelstd_M  = std((Wake24hWTBL_M./60*100)','omitnan')/sqrt(size(Wake24hWTBL_M,2));
NREMWTBLRelmean_M = mean(NREM24hWTBL_M./60*100,2,'omitnan')';
NREMWTBLRelstd_M  = std((NREM24hWTBL_M./60*100)','omitnan')/sqrt(size(NREM24hWTBL_M,2));
REMWTBLRelmean_M  = mean(REM24hWTBL_M./60*100,2,'omitnan')';
REMWTBLRelstd_M   = std((REM24hWTBL_M./60*100)','omitnan')/sqrt(size(REM24hWTBL_M,2));

% --- WT Baseline time in stage per hour Female
WakeWTBLRelmean_F = mean(Wake24hWTBL_F./60*100,2,'omitnan')';
WakeWTBLRelstd_F  = std((Wake24hWTBL_F./60*100)','omitnan')/sqrt(size(Wake24hWTBL_F,2));
NREMWTBLRelmean_F = mean(NREM24hWTBL_F./60*100,2,'omitnan')';
NREMWTBLRelstd_F  = std((NREM24hWTBL_F./60*100)','omitnan')/sqrt(size(NREM24hWTBL_F,2));
REMWTBLRelmean_F  = mean(REM24hWTBL_F./60*100,2,'omitnan')';
REMWTBLRelstd_F   = std((REM24hWTBL_F./60*100)','omitnan')/sqrt(size(REM24hWTBL_F,2));

% --- Mutant Baseline time in stage per hour normalized Male
WakeMutBLRelmean_M = mean(Wake24hMutBL_M./60*100,2,'omitnan')';
WakeMutBLRelstd_M  = std((Wake24hMutBL_M./60*100)','omitnan')/sqrt(size(Wake24hMutBL_M,2));
NREMMutBLRelmean_M = mean(NREM24hMutBL_M./60*100,2,'omitnan')';
NREMMutBLRelstd_M  = std((NREM24hMutBL_M./60*100)','omitnan')/sqrt(size(NREM24hMutBL_M,2));
REMMutBLRelmean_M  = mean(REM24hMutBL_M./60*100,2,'omitnan')';
REMMutBLRelstd_M   = std((REM24hMutBL_M./60*100)','omitnan')/sqrt(size(REM24hMutBL_M,2));

% --- Mutant Baseline time in stage per hour normalized Female
WakeMutBLRelmean_F = mean(Wake24hMutBL_F./60*100,2,'omitnan')';
WakeMutBLRelstd_F  = std((Wake24hMutBL_F./60*100)','omitnan')/sqrt(size(Wake24hMutBL_F,2));
NREMMutBLRelmean_F = mean(NREM24hMutBL_F./60*100,2,'omitnan')';
NREMMutBLRelstd_F  = std((NREM24hMutBL_F./60*100)','omitnan')/sqrt(size(NREM24hMutBL_F,2));
REMMutBLRelmean_F  = mean(REM24hMutBL_F./60*100,2,'omitnan')';
REMMutBLRelstd_F   = std((REM24hMutBL_F./60*100)','omitnan')/sqrt(size(REM24hMutBL_F,2));

% --- WT Sleep Dep time in stage per hour Male
WakeWTSDRelmean_M = mean(Wake24hWTSD_M./60*100,2,'omitnan')';
WakeWTSDRelstd_M  = std((Wake24hWTSD_M./60*100)','omitnan')/sqrt(size(Wake24hWTSD_M,2));
NREMWTSDRelmean_M = mean(NREM24hWTSD_M./60*100,2,'omitnan')';
NREMWTSDRelstd_M  = std((NREM24hWTSD_M./60*100)','omitnan')/sqrt(size(NREM24hWTSD_M,2));
REMWTSDRelmean_M  = mean(REM24hWTSD_M./60*100,2,'omitnan')';
REMWTSDRelstd_M   = std((REM24hWTSD_M./60*100)','omitnan')/sqrt(size(REM24hWTSD_M,2));

% --- WT Sleep Dep time in stage per hour Female
WakeWTSDRelmean_F = mean(Wake24hWTSD_F./60*100,2,'omitnan')';
WakeWTSDRelstd_F  = std((Wake24hWTSD_F./60*100)','omitnan')/sqrt(size(Wake24hWTSD_F,2));
NREMWTSDRelmean_F = mean(NREM24hWTSD_F./60*100,2,'omitnan')';
NREMWTSDRelstd_F  = std((NREM24hWTSD_F./60*100)','omitnan')/sqrt(size(NREM24hWTSD_F,2));
REMWTSDRelmean_F  = mean(REM24hWTSD_F./60*100,2,'omitnan')';
REMWTSDRelstd_F   = std((REM24hWTSD_F./60*100)','omitnan')/sqrt(size(REM24hWTSD_F,2));

% --- Mutant Sleep Dep time in stage per hour Male
WakeMutSDRelmean_M = mean(Wake24hMutSD_M./60*100,2,'omitnan')';
WakeMutSDRelstd_M  = std((Wake24hMutSD_M./60*100)','omitnan')/sqrt(size(Wake24hMutSD_M,2));
NREMMutSDRelmean_M = mean(NREM24hMutSD_M./60*100,2,'omitnan')';
NREMMutSDRelstd_M  = std((NREM24hMutSD_M./60*100)','omitnan')/sqrt(size(NREM24hMutSD_M,2));
REMMutSDRelmean_M  = mean(REM24hMutSD_M./60*100,2,'omitnan')';
REMMutSDRelstd_M   = std((REM24hMutSD_M./60*100)','omitnan')/sqrt(size(REM24hMutSD_M,2));

% --- Mutant Sleep Dep time in stage per hour Female
WakeMutSDRelmean_F = mean(Wake24hMutSD_F./60*100,2,'omitnan')';
WakeMutSDRelstd_F  = std((Wake24hMutSD_F./60*100)','omitnan')/sqrt(size(Wake24hMutSD_F,2));
NREMMutSDRelmean_F = mean(NREM24hMutSD_F./60*100,2,'omitnan')';
NREMMutSDRelstd_F  = std((NREM24hMutSD_F./60*100)','omitnan')/sqrt(size(NREM24hMutSD_F,2));
REMMutSDRelmean_F  = mean(REM24hMutSD_F./60*100,2,'omitnan')';
REMMutSDRelstd_F   = std((REM24hMutSD_F./60*100)','omitnan')/sqrt(size(REM24hMutSD_F,2));
%%
% ---------------------------------------------------------------------------------




% Figure out which ANOVAs need a post-hoc done, so we know where 
% to add asterisks to the plot. 
% ------ BL first 12 hours ------
% -- Wake --
% - Genotype 
p_vals.BL.LP.Wake.Posthocs.([Sexes{2} '_WTvsMut']) = table;
p_vals.BL.LP.Wake.Posthocs.([Sexes{1} '_WTvsMut']) = table;
interaction_term_text = ['Genotype:' SexVarName];
if RepMeas_ANOVA_tables.BL.First12hr.Wake.pValue('Genotype') < 0.05 | RepMeas_ANOVA_tables.BL.First12hr.Wake.pValue(interaction_term_text) < 0.05 | RepMeas_ANOVA_tables.BL.First12hr.Wake.pValue('Genotype:Time') < 0.05
	PostHoc_bool.Genotype.BLFirst12hr.W = 1;
	
	mc_table.Genotype.BL.First12hr.Wake 	   = multcompare(RM_Models.BL.LP.Wake,'Genotype','By','Sex','ComparisonType','lsd');  % For * symbols on graphs
	[~, ia] = unique(mc_table.Genotype.BL.First12hr.Wake.Sex(:), 'stable');  % remove duplicate rows
	mc_table.Genotype.BL.First12hr.Wake = mc_table.Genotype.BL.First12hr.Wake(ia,:);

	
	%unadjusted
	p_vals.BL.LP.Wake.Posthocs.([Sexes{2} '_WTvsMut']).Unadjusted = mc_table.Genotype.BL.First12hr.Wake.pValue(find(mc_table.Genotype.BL.First12hr.Wake.Sex==Sexes{2}));
	p_vals.BL.LP.Wake.Posthocs.([Sexes{1} '_WTvsMut']).Unadjusted = mc_table.Genotype.BL.First12hr.Wake.pValue(find(mc_table.Genotype.BL.First12hr.Wake.Sex==Sexes{1}));

	% adjusting
	mc_table.Genotype.BL.First12hr.Wake.pValue = mafdr(mc_table.Genotype.BL.First12hr.Wake.pValue,'BHFDR',true);
	% mc_table.Genotype.BL.First12hr.Wake = mc_table.Genotype.BL.First12hr.Wake(ia,:);

	p_vals.BL.LP.Wake.Posthocs.([Sexes{2} '_WTvsMut']).Adjusted = mc_table.Genotype.BL.First12hr.Wake.pValue(find(mc_table.Genotype.BL.First12hr.Wake.Sex==Sexes{2}));
	p_vals.BL.LP.Wake.Posthocs.([Sexes{1} '_WTvsMut']).Adjusted = mc_table.Genotype.BL.First12hr.Wake.pValue(find(mc_table.Genotype.BL.First12hr.Wake.Sex==Sexes{1}));


else 
	PostHoc_bool.Genotype.BLFirst12hr.W  = 0;
	mc_table.Genotype.BL.First12hr.Wake  = [];
	p_vals.BL.LP.Wake.Posthocs.([Sexes{2} '_WTvsMut']).Adjusted = zeros(0);
	p_vals.BL.LP.Wake.Posthocs.([Sexes{1} '_WTvsMut']).Adjusted = zeros(0);
end 

% - Sex - 
p_vals.BL.LP.Wake.Posthocs.(['WT_'  Sexes{1} 'vs' Sexes{2}])  = table;
p_vals.BL.LP.Wake.Posthocs.(['Mut_' Sexes{1} 'vs' Sexes{2}]) = table;
if RepMeas_ANOVA_tables.BL.First12hr.Wake.pValue(SexVarName) < 0.05 | RepMeas_ANOVA_tables.BL.First12hr.Wake.pValue(interaction_term_text) < 0.05 | RepMeas_ANOVA_tables.BL.First12hr.Wake.pValue([SexVarName ':Time']) < 0.05
	PostHoc_bool.(SexVarName).BLFirst12hr.W = 1;
	
	mc_table.(SexVarName).BL.First12hr.Wake = multcompare(RM_Models.BL.LP.Wake,'Sex','By','Genotype','ComparisonType','lsd');  % For # symbols on graphs
	[~, ia] = unique(mc_table.(SexVarName).BL.First12hr.Wake.Genotype(:), 'stable');  % remove duplicate rows
	mc_table.(SexVarName).BL.First12hr.Wake = mc_table.(SexVarName).BL.First12hr.Wake(ia,:);

	%unadjusted
	p_vals.BL.LP.Wake.Posthocs.(['WT_'  Sexes{1} 'vs' Sexes{2}]).Unadjusted = mc_table.(SexVarName).BL.First12hr.Wake.pValue(find(mc_table.(SexVarName).BL.First12hr.Wake.Genotype=='WT'));
	p_vals.BL.LP.Wake.Posthocs.(['Mut_' Sexes{1} 'vs' Sexes{2}]).Unadjusted = mc_table.(SexVarName).BL.First12hr.Wake.pValue(find(mc_table.(SexVarName).BL.First12hr.Wake.Genotype=='Mut'));
	
	% adjusting
	mc_table.(SexVarName).BL.First12hr.Wake.pValue = mafdr(mc_table.(SexVarName).BL.First12hr.Wake.pValue,'BHFDR',true);
	
	p_vals.BL.LP.Wake.Posthocs.(['WT_'  Sexes{1} 'vs' Sexes{2}]).Adjusted = mc_table.(SexVarName).BL.First12hr.Wake.pValue(find(mc_table.(SexVarName).BL.First12hr.Wake.Genotype=='WT'));
	p_vals.BL.LP.Wake.Posthocs.(['Mut_' Sexes{1} 'vs' Sexes{2}]).Adjusted = mc_table.(SexVarName).BL.First12hr.Wake.pValue(find(mc_table.(SexVarName).BL.First12hr.Wake.Genotype=='Mut'));
else 
	PostHoc_bool.(SexVarName).BLFirst12hr.W = 0;
	mc_table.(SexVarName).BL.First12hr.Wake = [];
	p_vals.BL.LP.Wake.Posthocs.(['WT_'  Sexes{1} 'vs' Sexes{2}]).Adjusted  = zeros(0);
	p_vals.BL.LP.Wake.Posthocs.(['Mut_' Sexes{1} 'vs' Sexes{2}]).Adjusted = zeros(0);
end 

% -- NREM --
% - Genotype -
p_vals.BL.LP.NREM.Posthocs.([Sexes{2} '_WTvsMut']) = table;
p_vals.BL.LP.NREM.Posthocs.([Sexes{1} '_WTvsMut']) = table;
if RepMeas_ANOVA_tables.BL.First12hr.NREM.pValue('Genotype') < 0.05 | RepMeas_ANOVA_tables.BL.First12hr.NREM.pValue(interaction_term_text) < 0.05 | RepMeas_ANOVA_tables.BL.First12hr.NREM.pValue('Genotype:Time') < 0.05
	PostHoc_bool.Genotype.BLFirst12hr.N = 1; 
	
	mc_table.Genotype.BL.First12hr.NREM 	   = multcompare(RM_Models.BL.LP.NREM,'Genotype','By','Sex','ComparisonType','lsd'); % For * symbols on graphs
	[~, ia] = unique(mc_table.Genotype.BL.First12hr.NREM.Sex(:), 'stable');  % remove duplicate rows
	mc_table.Genotype.BL.First12hr.NREM = mc_table.Genotype.BL.First12hr.NREM(ia,:);

	%unadjusted
	p_vals.BL.LP.NREM.Posthocs.([Sexes{2} '_WTvsMut']).Unadjusted = mc_table.Genotype.BL.First12hr.NREM.pValue(find(mc_table.Genotype.BL.First12hr.NREM.Sex==Sexes{2}));
	p_vals.BL.LP.NREM.Posthocs.([Sexes{1} '_WTvsMut']).Unadjusted = mc_table.Genotype.BL.First12hr.NREM.pValue(find(mc_table.Genotype.BL.First12hr.NREM.Sex==Sexes{1}));

	% adjusting
	mc_table.Genotype.BL.First12hr.NREM.pValue = mafdr(mc_table.Genotype.BL.First12hr.NREM.pValue,'BHFDR',true);
	
	p_vals.BL.LP.NREM.Posthocs.([Sexes{2} '_WTvsMut']).Adjusted = mc_table.Genotype.BL.First12hr.NREM.pValue(find(mc_table.Genotype.BL.First12hr.NREM.Sex==Sexes{2}));
	p_vals.BL.LP.NREM.Posthocs.([Sexes{1} '_WTvsMut']).Adjusted = mc_table.Genotype.BL.First12hr.NREM.pValue(find(mc_table.Genotype.BL.First12hr.NREM.Sex==Sexes{1}));


else 
	PostHoc_bool.Genotype.BLFirst12hr.N  = 0; 
	mc_table.Genotype.BL.First12hr.NREM  = [];
	p_vals.BL.LP.NREM.Posthocs.([Sexes{2} '_WTvsMut']).Adjusted = zeros(0);
	p_vals.BL.LP.NREM.Posthocs.([Sexes{1} 'M_WTvsMut']).Adjusted = zeros(0);
end 

% - Sex - 
p_vals.BL.LP.NREM.Posthocs.(['WT_'  Sexes{1} 'vs' Sexes{2}]) = table;
p_vals.BL.LP.NREM.Posthocs.(['Mut_' Sexes{1} 'vs' Sexes{2}]) = table;
if RepMeas_ANOVA_tables.BL.First12hr.NREM.pValue(SexVarName) < 0.05 | RepMeas_ANOVA_tables.BL.First12hr.NREM.pValue(interaction_term_text) < 0.05 | RepMeas_ANOVA_tables.BL.First12hr.NREM.pValue([ SexVarName ':Time']) < 0.05
	PostHoc_bool.(SexVarName).BLFirst12hr.N = 1;
	
	mc_table.(SexVarName).BL.First12hr.NREM = multcompare(RM_Models.BL.LP.NREM,'Sex','By','Genotype','ComparisonType','lsd'); % For * symbols on graphs
	[~, ia] = unique(mc_table.(SexVarName).BL.First12hr.NREM.Genotype(:), 'stable');  % remove duplicate rows
	mc_table.(SexVarName).BL.First12hr.NREM = mc_table.(SexVarName).BL.First12hr.NREM(ia,:);

	% unadjusted
	p_vals.BL.LP.NREM.Posthocs.(['WT_'  Sexes{1} 'vs' Sexes{2}]).Unadjusted = mc_table.(SexVarName).BL.First12hr.NREM.pValue(find(mc_table.(SexVarName).BL.First12hr.NREM.Genotype=='WT'));
	p_vals.BL.LP.NREM.Posthocs.(['Mut_' Sexes{1} 'vs' Sexes{2}]).Unadjusted = mc_table.(SexVarName).BL.First12hr.NREM.pValue(find(mc_table.(SexVarName).BL.First12hr.NREM.Genotype=='Mut'));

	% adjusting
	mc_table.(SexVarName).BL.First12hr.NREM.pValue = mafdr(mc_table.(SexVarName).BL.First12hr.NREM.pValue,'BHFDR',true);
	
	p_vals.BL.LP.NREM.Posthocs.(['WT_'  Sexes{1} 'vs' Sexes{2}]).Adjusted = mc_table.(SexVarName).BL.First12hr.NREM.pValue(find(mc_table.(SexVarName).BL.First12hr.NREM.Genotype=='WT'));
	p_vals.BL.LP.NREM.Posthocs.(['Mut_' Sexes{1} 'vs' Sexes{2}]).Adjusted = mc_table.(SexVarName).BL.First12hr.NREM.pValue(find(mc_table.(SexVarName).BL.First12hr.NREM.Genotype=='Mut'));

else 
	PostHoc_bool.(SexVarName).BLFirst12hr.N     = 0; 
	mc_table.(SexVarName).BL.First12hr.NREM     = [];
	p_vals.BL.LP.NREM.Posthocs.(['WT_'  Sexes{1} 'vs' Sexes{2}]).Adjusted = zeros(0);
	p_vals.BL.LP.NREM.Posthocs.(['Mut_' Sexes{1} 'vs' Sexes{2}]).Adjusted = zeros(0);
end 

% - REM - 
% - Genotype - 
p_vals.BL.LP.REM.Posthocs.([Sexes{2} '_WTvsMut']) = table;
p_vals.BL.LP.REM.Posthocs.([Sexes{1} '_WTvsMut']) = table;
if RepMeas_ANOVA_tables.BL.First12hr.REM.pValue('Genotype') < 0.05 | RepMeas_ANOVA_tables.BL.First12hr.REM.pValue(interaction_term_text) < 0.05 | RepMeas_ANOVA_tables.BL.First12hr.REM.pValue('Genotype:Time') < 0.05
	PostHoc_bool.Genotype.BLFirst12hr.R = 1;
	
	mc_table.Genotype.BL.First12hr.REM 		  = multcompare(RM_Models.BL.LP.REM,'Genotype','By','Sex','ComparisonType','lsd');
	[~, ia] = unique(mc_table.Genotype.BL.First12hr.REM.Sex(:), 'stable');  % remove duplicate rows
	mc_table.Genotype.BL.First12hr.REM = mc_table.Genotype.BL.First12hr.REM(ia,:);

	% unadjusted
	p_vals.BL.LP.REM.Posthocs.([Sexes{2} '_WTvsMut']).Unadjusted = mc_table.Genotype.BL.First12hr.REM.pValue(find(mc_table.Genotype.BL.First12hr.REM.Sex==Sexes{2}));
	p_vals.BL.LP.REM.Posthocs.([Sexes{1} '_WTvsMut']).Unadjusted = mc_table.Genotype.BL.First12hr.REM.pValue(find(mc_table.Genotype.BL.First12hr.REM.Sex==Sexes{1}));
	
	% adjusting
	mc_table.Genotype.BL.First12hr.REM.pValue = mafdr(mc_table.Genotype.BL.First12hr.REM.pValue,'BHFDR',true);
	
	p_vals.BL.LP.REM.Posthocs.([Sexes{2} '_WTvsMut']).Adjusted = mc_table.Genotype.BL.First12hr.REM.pValue(find(mc_table.Genotype.BL.First12hr.REM.Sex==Sexes{2}));
	p_vals.BL.LP.REM.Posthocs.([Sexes{1} '_WTvsMut']).Adjusted = mc_table.Genotype.BL.First12hr.REM.pValue(find(mc_table.Genotype.BL.First12hr.REM.Sex==Sexes{1}));

else 
	PostHoc_bool.Genotype.BLFirst12hr.R = 0;
	mc_table.Genotype.BL.First12hr.REM  = [];
	p_vals.BL.LP.REM.Posthocs.([Sexes{2} '_WTvsMut']).Adjusted = zeros(0);
	p_vals.BL.LP.REM.Posthocs.([Sexes{1} '_WTvsMut']).Adjusted = zeros(0);
end 

% - Sex - 
p_vals.BL.LP.REM.Posthocs.(['WT_'  Sexes{1} 'vs' Sexes{2}]) = table;
p_vals.BL.LP.REM.Posthocs.(['Mut_' Sexes{1} 'vs' Sexes{2}]) = table;
if RepMeas_ANOVA_tables.BL.First12hr.REM.pValue(SexVarName) < 0.05 | RepMeas_ANOVA_tables.BL.First12hr.REM.pValue(interaction_term_text) < 0.05 | RepMeas_ANOVA_tables.BL.First12hr.REM.pValue([SexVarName ':Time']) < 0.05
	PostHoc_bool.(SexVarName).BLFirst12hr.R = 1;
	
	mc_table.(SexVarName).BL.First12hr.REM = multcompare(RM_Models.BL.LP.REM,'Sex','By','Genotype','ComparisonType','lsd');
	[~, ia] = unique(mc_table.(SexVarName).BL.First12hr.REM.Genotype(:), 'stable');  % remove duplicate rows
	mc_table.(SexVarName).BL.First12hr.REM = mc_table.(SexVarName).BL.First12hr.REM(ia,:);

	%unadjusted
	p_vals.BL.LP.REM.Posthocs.(['WT_'  Sexes{1} 'vs' Sexes{2}]).Unadjusted = mc_table.(SexVarName).BL.First12hr.REM.pValue(find(mc_table.(SexVarName).BL.First12hr.REM.Genotype=='WT'));
	p_vals.BL.LP.REM.Posthocs.(['Mut_' Sexes{1} 'vs' Sexes{2}]).Unadjusted = mc_table.(SexVarName).BL.First12hr.REM.pValue(find(mc_table.(SexVarName).BL.First12hr.REM.Genotype=='Mut'));

	%adjusted
	mc_table.(SexVarName).BL.First12hr.REM.pValue = mafdr(mc_table.(SexVarName).BL.First12hr.REM.pValue,'BHFDR',true);
	
	p_vals.BL.LP.REM.Posthocs.(['WT_'  Sexes{1} 'vs' Sexes{2}]).Adjusted = mc_table.(SexVarName).BL.First12hr.REM.pValue(find(mc_table.(SexVarName).BL.First12hr.REM.Genotype=='WT'));
	p_vals.BL.LP.REM.Posthocs.(['Mut_' Sexes{1} 'vs' Sexes{2}]).Adjusted = mc_table.(SexVarName).BL.First12hr.REM.pValue(find(mc_table.(SexVarName).BL.First12hr.REM.Genotype=='Mut'));

else 
	PostHoc_bool.(SexVarName).BLFirst12hr.R     = 0;
	mc_table.(SexVarName).BL.First12hr.REM      = [];
	p_vals.BL.LP.REM.Posthocs.(['WT_'  Sexes{1} 'vs' Sexes{2}]).Adjusted = zeros(0);
	p_vals.BL.LP.REM.Posthocs.(['Mut_' Sexes{1} 'vs' Sexes{2}]).Adjusted = zeros(0);
end 

% ------ BL last 12 hours ------
% -- Wake --
% - Genotype 
p_vals.BL.DP.Wake.Posthocs.([Sexes{2} '_WTvsMut']) = table;
p_vals.BL.DP.Wake.Posthocs.([Sexes{1} '_WTvsMut']) = table;
if RepMeas_ANOVA_tables.BL.Last12hr.Wake.pValue('Genotype') < 0.05 | RepMeas_ANOVA_tables.BL.Last12hr.Wake.pValue(interaction_term_text) < 0.05 | RepMeas_ANOVA_tables.BL.Last12hr.Wake.pValue('Genotype:Time') < 0.05
	PostHoc_bool.Genotype.BLLast12hr.W = 1;
	
	mc_table.Genotype.BL.Last12hr.Wake = multcompare(RM_Models.BL.DP.Wake,'Genotype','By','Sex','ComparisonType','lsd');  % For * symbols on graphs
	[~, ia] = unique(mc_table.Genotype.BL.Last12hr.Wake.Sex(:), 'stable');  % remove duplicate rows
	mc_table.Genotype.BL.Last12hr.Wake = mc_table.Genotype.BL.Last12hr.Wake(ia,:);

	% unadjusted
	p_vals.BL.DP.Wake.Posthocs.([Sexes{2} '_WTvsMut']).Unadjusted = mc_table.Genotype.BL.Last12hr.Wake.pValue(find(mc_table.Genotype.BL.Last12hr.Wake.Sex==Sexes{2}));
	p_vals.BL.DP.Wake.Posthocs.([Sexes{1} '_WTvsMut']).Unadjusted = mc_table.Genotype.BL.Last12hr.Wake.pValue(find(mc_table.Genotype.BL.Last12hr.Wake.Sex==Sexes{1}));

	% adjusted
	mc_table.Genotype.BL.Last12hr.Wake.pValue = mafdr(mc_table.Genotype.BL.Last12hr.Wake.pValue,'BHFDR',true);
	
	p_vals.BL.DP.Wake.Posthocs.([Sexes{2} '_WTvsMut']).Adjusted = mc_table.Genotype.BL.Last12hr.Wake.pValue(find(mc_table.Genotype.BL.Last12hr.Wake.Sex==Sexes{2}));
	p_vals.BL.DP.Wake.Posthocs.([Sexes{1} '_WTvsMut']).Adjusted = mc_table.Genotype.BL.Last12hr.Wake.pValue(find(mc_table.Genotype.BL.Last12hr.Wake.Sex==Sexes{1}));


else 
	PostHoc_bool.Genotype.BLLast12hr.W   = 0;
	mc_table.Genotype.BL.Last12hr.Wake   = [];
	p_vals.BL.DP.Wake.Posthocs.([Sexes{2} '_WTvsMut']).Adjusted = zeros(0);
	p_vals.BL.DP.Wake.Posthocs.([Sexes{1} '_WTvsMut']).Adjusted = zeros(0);
end 

% - Sex - 
p_vals.BL.DP.Wake.Posthocs.(['WT_'  Sexes{1} 'vs' Sexes{2}]) = table;
p_vals.BL.DP.Wake.Posthocs.(['Mut_' Sexes{1} 'vs' Sexes{2}]) = table;
if RepMeas_ANOVA_tables.BL.Last12hr.Wake.pValue(SexVarName) < 0.05 | RepMeas_ANOVA_tables.BL.Last12hr.Wake.pValue(interaction_term_text) < 0.05 | RepMeas_ANOVA_tables.BL.Last12hr.Wake.pValue([SexVarName ':Time']) < 0.05
	PostHoc_bool.(SexVarName).BLLast12hr.W = 1;
	
	mc_table.(SexVarName).BL.Last12hr.Wake = multcompare(RM_Models.BL.DP.Wake,'Sex','By','Genotype','ComparisonType','lsd');  % For # symbols on graphs
	[~, ia] = unique(mc_table.(SexVarName).BL.Last12hr.Wake.Genotype(:), 'stable');  % remove duplicate rows
	mc_table.(SexVarName).BL.Last12hr.Wake = mc_table.(SexVarName).BL.Last12hr.Wake(ia,:);
	
	%unadjusted
	p_vals.BL.DP.Wake.Posthocs.(['WT_'  Sexes{1} 'vs' Sexes{2}]).Unadjusted = mc_table.(SexVarName).BL.Last12hr.Wake.pValue(find(mc_table.(SexVarName).BL.Last12hr.Wake.Genotype=='WT'));
	p_vals.BL.DP.Wake.Posthocs.(['Mut_' Sexes{1} 'vs' Sexes{2}]).Unadjusted = mc_table.(SexVarName).BL.Last12hr.Wake.pValue(find(mc_table.(SexVarName).BL.Last12hr.Wake.Genotype=='Mut'));

	%adjusted
	mc_table.(SexVarName).BL.Last12hr.Wake.pValue = mafdr(mc_table.(SexVarName).BL.Last12hr.Wake.pValue,'BHFDR',true);
	
	p_vals.BL.DP.Wake.Posthocs.(['WT_'  Sexes{1} 'vs' Sexes{2}]).Adjusted = mc_table.(SexVarName).BL.Last12hr.Wake.pValue(find(mc_table.(SexVarName).BL.Last12hr.Wake.Genotype=='WT'));
	p_vals.BL.DP.Wake.Posthocs.(['Mut_' Sexes{1} 'vs' Sexes{2}]).Adjusted = mc_table.(SexVarName).BL.Last12hr.Wake.pValue(find(mc_table.(SexVarName).BL.Last12hr.Wake.Genotype=='Mut'));

else 
	PostHoc_bool.(SexVarName).BLLast12hr.W = 0;
	mc_table.(SexVarName).BL.Last12hr.Wake = [];
	p_vals.BL.DP.Wake.Posthocs.(['WT_'  Sexes{1} 'vs' Sexes{2}]).Adjusted = zeros(0);
	p_vals.BL.DP.Wake.Posthocs.(['Mut_' Sexes{1} 'vs' Sexes{2}]).Adjusted = zeros(0);
end 

% -- NREM --
% - Genotype -
p_vals.BL.DP.NREM.Posthocs.([Sexes{2} '_WTvsMut']) = table;
p_vals.BL.DP.NREM.Posthocs.([Sexes{1} '_WTvsMut']) = table;
if RepMeas_ANOVA_tables.BL.Last12hr.NREM.pValue('Genotype') < 0.05 | RepMeas_ANOVA_tables.BL.Last12hr.NREM.pValue(interaction_term_text) < 0.05 | RepMeas_ANOVA_tables.BL.Last12hr.NREM.pValue('Genotype:Time') < 0.05
	PostHoc_bool.Genotype.BLLast12hr.N = 1; 
	
	mc_table.Genotype.BL.Last12hr.NREM = multcompare(RM_Models.BL.DP.NREM,'Genotype','By','Sex','ComparisonType','lsd'); % For * symbols on graphs
	[~, ia] = unique(mc_table.Genotype.BL.Last12hr.NREM.Sex(:), 'stable');  % remove duplicate rows
	mc_table.Genotype.BL.Last12hr.NREM = mc_table.Genotype.BL.Last12hr.NREM(ia,:);

	%unadjusted
	p_vals.BL.DP.NREM.Posthocs.([Sexes{2} '_WTvsMut']).Unadjusted = mc_table.Genotype.BL.Last12hr.NREM.pValue(find(mc_table.Genotype.BL.Last12hr.NREM.Sex==Sexes{2}));
	p_vals.BL.DP.NREM.Posthocs.([Sexes{1} '_WTvsMut']).Unadjusted = mc_table.Genotype.BL.Last12hr.NREM.pValue(find(mc_table.Genotype.BL.Last12hr.NREM.Sex==Sexes{1}));

	%adjusted
	mc_table.Genotype.BL.Last12hr.NREM.pValue = mafdr(mc_table.Genotype.BL.Last12hr.NREM.pValue,'BHFDR',true);
	
	p_vals.BL.DP.NREM.Posthocs.([Sexes{2} '_WTvsMut']).Adjusted = mc_table.Genotype.BL.Last12hr.NREM.pValue(find(mc_table.Genotype.BL.Last12hr.NREM.Sex==Sexes{2}));
	p_vals.BL.DP.NREM.Posthocs.([Sexes{1} '_WTvsMut']).Adjusted = mc_table.Genotype.BL.Last12hr.NREM.pValue(find(mc_table.Genotype.BL.Last12hr.NREM.Sex==Sexes{1}));

else 
	PostHoc_bool.Genotype.BLLast12hr.N = 0; 
	mc_table.Genotype.BL.Last12hr.NREM = [];
	p_vals.BL.DP.NREM.Posthocs.([Sexes{2} '_WTvsMut']).Adjusted = zeros(0);
	p_vals.BL.DP.NREM.Posthocs.([Sexes{1} '_WTvsMut']).Adjusted = zeros(0);
end 

% - Sex - 
p_vals.BL.DP.NREM.Posthocs.(['WT_'  Sexes{1} 'vs' Sexes{2}]) = table;
p_vals.BL.DP.NREM.Posthocs.(['Mut_' Sexes{1} 'vs' Sexes{2}]) = table;
if RepMeas_ANOVA_tables.BL.Last12hr.NREM.pValue(SexVarName) < 0.05 | RepMeas_ANOVA_tables.BL.Last12hr.NREM.pValue(interaction_term_text) < 0.05 | RepMeas_ANOVA_tables.BL.Last12hr.NREM.pValue([SexVarName ':Time']) < 0.05
	PostHoc_bool.(SexVarName).BLLast12hr.N = 1; 
	
	mc_table.(SexVarName).BL.Last12hr.NREM = multcompare(RM_Models.BL.DP.NREM,'Sex','By','Genotype','ComparisonType','lsd'); % For * symbols on graphs
	[~, ia] = unique(mc_table.(SexVarName).BL.Last12hr.NREM.Genotype(:), 'stable');  % remove duplicate rows
	mc_table.(SexVarName).BL.Last12hr.NREM = mc_table.(SexVarName).BL.Last12hr.NREM(ia,:);

	% unadjusted
	p_vals.BL.DP.NREM.Posthocs.(['WT_'  Sexes{1} 'vs' Sexes{2}]).Unadjusted = mc_table.(SexVarName).BL.Last12hr.NREM.pValue(find(mc_table.(SexVarName).BL.Last12hr.NREM.Genotype=='WT'));
	p_vals.BL.DP.NREM.Posthocs.(['Mut_' Sexes{1} 'vs' Sexes{2}]).Unadjusted = mc_table.(SexVarName).BL.Last12hr.NREM.pValue(find(mc_table.(SexVarName).BL.Last12hr.NREM.Genotype=='Mut'));

	% adjusted
	mc_table.(SexVarName).BL.Last12hr.NREM.pValue = mafdr(mc_table.(SexVarName).BL.Last12hr.NREM.pValue,'BHFDR',true);
	
	p_vals.BL.DP.NREM.Posthocs.(['WT_'  Sexes{1} 'vs' Sexes{2}]).Adjusted = mc_table.(SexVarName).BL.Last12hr.NREM.pValue(find(mc_table.(SexVarName).BL.Last12hr.NREM.Genotype=='WT'));
	p_vals.BL.DP.NREM.Posthocs.(['Mut_' Sexes{1} 'vs' Sexes{2}]).Adjusted = mc_table.(SexVarName).BL.Last12hr.NREM.pValue(find(mc_table.(SexVarName).BL.Last12hr.NREM.Genotype=='Mut'));

else 
	PostHoc_bool.(SexVarName).BLLast12hr.N = 0; 
	mc_table.(SexVarName).BL.Last12hr.NREM = [];
	p_vals.BL.DP.NREM.Posthocs.(['WT_'  Sexes{1} 'vs' Sexes{2}]).Adjusted = zeros(0);
	p_vals.BL.DP.NREM.Posthocs.(['Mut_' Sexes{1} 'vs' Sexes{2}]).Adjusted = zeros(0);
end 

% - REM - 
% - Genotype - 
p_vals.BL.DP.REM.Posthocs.([Sexes{2} '_WTvsMut']) = table;
p_vals.BL.DP.REM.Posthocs.([Sexes{1} '_WTvsMut']) = table;
if RepMeas_ANOVA_tables.BL.Last12hr.REM.pValue('Genotype') < 0.05 | RepMeas_ANOVA_tables.BL.Last12hr.REM.pValue(interaction_term_text) < 0.05 | RepMeas_ANOVA_tables.BL.Last12hr.REM.pValue('Genotype:Time') < 0.05
	PostHoc_bool.Genotype.BLLast12hr.R = 1;
	
	mc_table.Genotype.BL.Last12hr.REM = multcompare(RM_Models.BL.DP.REM,'Genotype','By','Sex','ComparisonType','lsd');
	[~, ia] = unique(mc_table.Genotype.BL.Last12hr.REM.Sex(:), 'stable');  % remove duplicate rows
	mc_table.Genotype.BL.Last12hr.REM = mc_table.Genotype.BL.Last12hr.REM(ia,:);

	%unadjusted
	p_vals.BL.DP.REM.Posthocs.([Sexes{2} '_WTvsMut']).Unadjusted = mc_table.Genotype.BL.Last12hr.REM.pValue(find(mc_table.Genotype.BL.Last12hr.REM.Sex==Sexes{2}));
	p_vals.BL.DP.REM.Posthocs.([Sexes{1} '_WTvsMut']).Unadjusted = mc_table.Genotype.BL.Last12hr.REM.pValue(find(mc_table.Genotype.BL.Last12hr.REM.Sex==Sexes{1}));

	% adjusted
	mc_table.Genotype.BL.Last12hr.REM.pValue = mafdr(mc_table.Genotype.BL.Last12hr.REM.pValue,'BHFDR',true);
	
	p_vals.BL.DP.REM.Posthocs.([Sexes{2} '_WTvsMut']).Adjusted = mc_table.Genotype.BL.Last12hr.REM.pValue(find(mc_table.Genotype.BL.Last12hr.REM.Sex==Sexes{2}));
	p_vals.BL.DP.REM.Posthocs.([Sexes{1} '_WTvsMut']).Adjusted = mc_table.Genotype.BL.Last12hr.REM.pValue(find(mc_table.Genotype.BL.Last12hr.REM.Sex==Sexes{1}));

else 
	PostHoc_bool.Genotype.BLLast12hr.R = 0;
	mc_table.Genotype.BL.Last12hr.REM  = [];
	p_vals.BL.DP.REM.Posthocs.([Sexes{2} '_WTvsMut']).Adjusted = zeros(0);
	p_vals.BL.DP.REM.Posthocs.([Sexes{1} '_WTvsMut']).Adjusted = zeros(0);
end 

% - Sex - 
p_vals.BL.DP.REM.Posthocs.(['WT_'  Sexes{1} 'vs' Sexes{2}]) = table;
p_vals.BL.DP.REM.Posthocs.(['Mut_' Sexes{1} 'vs' Sexes{2}]) = table;
if RepMeas_ANOVA_tables.BL.Last12hr.REM.pValue(SexVarName) < 0.05 | RepMeas_ANOVA_tables.BL.Last12hr.REM.pValue(interaction_term_text) < 0.05 | RepMeas_ANOVA_tables.BL.Last12hr.REM.pValue([SexVarName ':Time']) < 0.05
	PostHoc_bool.(SexVarName).BLLast12hr.R = 1;
	
	mc_table.(SexVarName).BL.Last12hr.REM = multcompare(RM_Models.BL.DP.REM,'Sex','By','Genotype','ComparisonType','lsd');
	[~, ia] = unique(mc_table.(SexVarName).BL.Last12hr.REM.Genotype(:), 'stable');  % remove duplicate rows
	mc_table.(SexVarName).BL.Last12hr.REM = mc_table.(SexVarName).BL.Last12hr.REM(ia,:);

	% unadjusted
	p_vals.BL.DP.REM.Posthocs.(['WT_'  Sexes{1} 'vs' Sexes{2}]).Unadjusted  = mc_table.(SexVarName).BL.Last12hr.REM.pValue(find(mc_table.(SexVarName).BL.Last12hr.REM.Genotype=='WT'));
	p_vals.BL.DP.REM.Posthocs.(['Mut_' Sexes{1} 'vs' Sexes{2}]).Unadjusted = mc_table.(SexVarName).BL.Last12hr.REM.pValue(find(mc_table.(SexVarName).BL.Last12hr.REM.Genotype=='Mut'));

	% adjust
	mc_table.(SexVarName).BL.Last12hr.REM.pValue = mafdr(mc_table.(SexVarName).BL.Last12hr.REM.pValue,'BHFDR',true);
	
	p_vals.BL.DP.REM.Posthocs.(['WT_'  Sexes{1} 'vs' Sexes{2}]).Adjusted = mc_table.(SexVarName).BL.Last12hr.REM.pValue(find(mc_table.(SexVarName).BL.Last12hr.REM.Genotype=='WT'));
	p_vals.BL.DP.REM.Posthocs.(['Mut_' Sexes{1} 'vs' Sexes{2}]).Adjusted = mc_table.(SexVarName).BL.Last12hr.REM.pValue(find(mc_table.(SexVarName).BL.Last12hr.REM.Genotype=='Mut'));

else 
	PostHoc_bool.(SexVarName).BLLast12hr.R = 0;
	mc_table.(SexVarName).BL.Last12hr.REM  = [];
	p_vals.BL.DP.REM.Posthocs.(['WT_'  Sexes{1} 'vs' Sexes{2}]).Adjusted = zeros(0);
	p_vals.BL.DP.REM.Posthocs.(['Mut_' Sexes{1} 'vs' Sexes{2}]).Adjusted = zeros(0);
end 


% ------ SD first 12 hours ------
% -- Wake --
% - Genotype 
p_vals.SDexcSD.LP.Wake.Posthocs.([Sexes{2} '_WTvsMut']) = table;
p_vals.SDexcSD.LP.Wake.Posthocs.([Sexes{1} '_WTvsMut']) = table;
if RepMeas_ANOVA_tables.SD.First12hr.Wake.pValue('Genotype') < 0.05 | RepMeas_ANOVA_tables.SD.First12hr.Wake.pValue(interaction_term_text) < 0.05 | RepMeas_ANOVA_tables.SD.First12hr.Wake.pValue('Genotype:Time') < 0.05
	PostHoc_bool.Genotype.SDFirst12hr.W = 1;
	
	mc_table.Genotype.SD.First12hr.Wake = multcompare(RM_Models.SD.LP.Wake,'Genotype','By','Sex','ComparisonType','lsd');  % For * symbols on graphs
	[~, ia] = unique(mc_table.Genotype.SD.First12hr.Wake.Sex(:), 'stable');  % remove duplicate rows
	mc_table.Genotype.SD.First12hr.Wake = mc_table.Genotype.SD.First12hr.Wake(ia,:);

	%unadjusted
	p_vals.SDexcSD.LP.Wake.Posthocs.([Sexes{2} '_WTvsMut']).Unadjusted = mc_table.Genotype.SD.First12hr.Wake.pValue(find(mc_table.Genotype.SD.First12hr.Wake.Sex==Sexes{2}));
	p_vals.SDexcSD.LP.Wake.Posthocs.([Sexes{1} '_WTvsMut']).Unadjusted = mc_table.Genotype.SD.First12hr.Wake.pValue(find(mc_table.Genotype.SD.First12hr.Wake.Sex==Sexes{1}));

	% adjusting
	mc_table.Genotype.SD.First12hr.Wake.pValue = mafdr(mc_table.Genotype.SD.First12hr.Wake.pValue,'BHFDR',true);
	
	p_vals.SDexcSD.LP.Wake.Posthocs.([Sexes{2} '_WTvsMut']).Adjusted = mc_table.Genotype.SD.First12hr.Wake.pValue(find(mc_table.Genotype.SD.First12hr.Wake.Sex==Sexes{2}));
	p_vals.SDexcSD.LP.Wake.Posthocs.([Sexes{1} '_WTvsMut']).Adjusted = mc_table.Genotype.SD.First12hr.Wake.pValue(find(mc_table.Genotype.SD.First12hr.Wake.Sex==Sexes{1}));


else 
	PostHoc_bool.Genotype.SDFirst12hr.W = 0;
	mc_table.Genotype.SD.First12hr.Wake = [];
	p_vals.SDexcSD.LP.Wake.Posthocs.([Sexes{2} '_WTvsMut']).Adjusted = zeros(0);
	p_vals.SDexcSD.LP.Wake.Posthocs.([Sexes{1} '_WTvsMut']).Adjusted = zeros(0);
end 

% - Sex - 
p_vals.SDexcSD.LP.Wake.Posthocs.(['WT_'  Sexes{1} 'vs' Sexes{2}]) = table;
p_vals.SDexcSD.LP.Wake.Posthocs.(['Mut_' Sexes{1} 'vs' Sexes{2}]) = table;
if RepMeas_ANOVA_tables.SD.First12hr.Wake.pValue(SexVarName) < 0.05 | RepMeas_ANOVA_tables.SD.First12hr.Wake.pValue(interaction_term_text) < 0.05 | RepMeas_ANOVA_tables.SD.First12hr.Wake.pValue([SexVarName ':Time']) < 0.05
	PostHoc_bool.(SexVarName).SDFirst12hr.W = 1;
	
	mc_table.(SexVarName).SD.First12hr.Wake = multcompare(RM_Models.SD.LP.Wake,'Sex','By','Genotype','ComparisonType','lsd');  % For # symbols on graphs
	[~, ia] = unique(mc_table.(SexVarName).SD.First12hr.Wake.Genotype(:), 'stable');  % remove duplicate rows
	mc_table.(SexVarName).SD.First12hr.Wake = mc_table.(SexVarName).SD.First12hr.Wake(ia,:);

	% unadjusted 
	p_vals.SDexcSD.LP.Wake.Posthocs.(['WT_'  Sexes{1} 'vs' Sexes{2}]).Unadjusted = mc_table.(SexVarName).SD.First12hr.Wake.pValue(find(mc_table.(SexVarName).SD.First12hr.Wake.Genotype=='WT'));
	p_vals.SDexcSD.LP.Wake.Posthocs.(['Mut_' Sexes{1} 'vs' Sexes{2}]).Unadjusted = mc_table.(SexVarName).SD.First12hr.Wake.pValue(find(mc_table.(SexVarName).SD.First12hr.Wake.Genotype=='Mut'));

	% adjust
	mc_table.(SexVarName).SD.First12hr.Wake.pValue = mafdr(mc_table.(SexVarName).SD.First12hr.Wake.pValue,'BHFDR',true);
	
	p_vals.SDexcSD.LP.Wake.Posthocs.(['WT_'  Sexes{1} 'vs' Sexes{2}]).Adjusted = mc_table.(SexVarName).SD.First12hr.Wake.pValue(find(mc_table.(SexVarName).SD.First12hr.Wake.Genotype=='WT'));
	p_vals.SDexcSD.LP.Wake.Posthocs.(['Mut_' Sexes{1} 'vs' Sexes{2}]).Adjusted = mc_table.(SexVarName).SD.First12hr.Wake.pValue(find(mc_table.(SexVarName).SD.First12hr.Wake.Genotype=='Mut'));


else 
	PostHoc_bool.(SexVarName).SDFirst12hr.W = 0;
	mc_table.(SexVarName).SD.First12hr.Wake = [];
	p_vals.SDexcSD.LP.Wake.Posthocs.(['WT_'  Sexes{1} 'vs' Sexes{2}]).Adjusted = zeros(0);
	p_vals.SDexcSD.LP.Wake.Posthocs.(['Mut_' Sexes{1} 'vs' Sexes{2}]).Adjusted = zeros(0);

end 

% -- NREM --
% - Genotype -
p_vals.SDexcSD.LP.NREM.Posthocs.([Sexes{2} '_WTvsMut']) = table;
p_vals.SDexcSD.LP.NREM.Posthocs.([Sexes{1} '_WTvsMut']) = table;
if RepMeas_ANOVA_tables.SD.First12hr.NREM.pValue('Genotype') < 0.05 | RepMeas_ANOVA_tables.SD.First12hr.NREM.pValue(interaction_term_text) < 0.05 | RepMeas_ANOVA_tables.SD.First12hr.NREM.pValue('Genotype:Time') < 0.05
	PostHoc_bool.Genotype.SDFirst12hr.N = 1;
	
	mc_table.Genotype.SD.First12hr.NREM = multcompare(RM_Models.SD.LP.NREM,'Genotype','By','Sex','ComparisonType','lsd'); % For * symbols on graphs
	[~, ia] = unique(mc_table.Genotype.SD.First12hr.NREM.Sex(:), 'stable');  % remove duplicate rows
	mc_table.Genotype.SD.First12hr.NREM = mc_table.Genotype.SD.First12hr.NREM(ia,:);

	% unadjusted
	p_vals.SDexcSD.LP.NREM.Posthocs.([Sexes{2} '_WTvsMut']).Unadjusted = mc_table.Genotype.SD.First12hr.NREM.pValue(find(mc_table.Genotype.SD.First12hr.NREM.Sex==Sexes{2}));
	p_vals.SDexcSD.LP.NREM.Posthocs.([Sexes{1} '_WTvsMut']).Unadjusted = mc_table.Genotype.SD.First12hr.NREM.pValue(find(mc_table.Genotype.SD.First12hr.NREM.Sex==Sexes{1}));

	% adjust
	mc_table.Genotype.SD.First12hr.NREM.pValue = mafdr(mc_table.Genotype.SD.First12hr.NREM.pValue,'BHFDR',true);
	
	p_vals.SDexcSD.LP.NREM.Posthocs.([Sexes{2} '_WTvsMut']).Adjusted = mc_table.Genotype.SD.First12hr.NREM.pValue(find(mc_table.Genotype.SD.First12hr.NREM.Sex==Sexes{2}));
	p_vals.SDexcSD.LP.NREM.Posthocs.([Sexes{1} '_WTvsMut']).Adjusted = mc_table.Genotype.SD.First12hr.NREM.pValue(find(mc_table.Genotype.SD.First12hr.NREM.Sex==Sexes{1}));

else 
	PostHoc_bool.Genotype.SDFirst12hr.N = 0; 
	mc_table.Genotype.SD.First12hr.NREM = [];
	p_vals.SDexcSD.LP.NREM.Posthocs.([Sexes{2} '_WTvsMut']).Adjusted = zeros(0);
	p_vals.SDexcSD.LP.NREM.Posthocs.([Sexes{1} '_WTvsMut']).Adjusted = zeros(0);
end 

% - Sex - 
p_vals.SDexcSD.LP.NREM.Posthocs.(['WT_'  Sexes{1} 'vs' Sexes{2}]) = table;
p_vals.SDexcSD.LP.NREM.Posthocs.(['Mut_' Sexes{1} 'vs' Sexes{2}]) = table;
if RepMeas_ANOVA_tables.SD.First12hr.NREM.pValue(SexVarName) < 0.05 | RepMeas_ANOVA_tables.SD.First12hr.NREM.pValue(interaction_term_text) < 0.05 | RepMeas_ANOVA_tables.SD.First12hr.NREM.pValue([SexVarName ':Time']) < 0.05
	PostHoc_bool.(SexVarName).SDFirst12hr.N = 1;
	
	mc_table.(SexVarName).SD.First12hr.NREM = multcompare(RM_Models.SD.LP.NREM,'Sex','By','Genotype','ComparisonType','lsd'); % For * symbols on graphs
	[~, ia] = unique(mc_table.(SexVarName).SD.First12hr.NREM.Genotype(:), 'stable');  % remove duplicate rows
	mc_table.(SexVarName).SD.First12hr.NREM = mc_table.(SexVarName).SD.First12hr.NREM(ia,:);

	%unadjusted
	p_vals.SDexcSD.LP.NREM.Posthocs.(['WT_'  Sexes{1} 'vs' Sexes{2}]).Unadjusted = mc_table.(SexVarName).SD.First12hr.NREM.pValue(find(mc_table.(SexVarName).SD.First12hr.NREM.Genotype=='WT'));
	p_vals.SDexcSD.LP.NREM.Posthocs.(['Mut_' Sexes{1} 'vs' Sexes{2}]).Unadjusted = mc_table.(SexVarName).SD.First12hr.NREM.pValue(find(mc_table.(SexVarName).SD.First12hr.NREM.Genotype=='Mut'));

	% adjust
	mc_table.(SexVarName).SD.First12hr.NREM.pValue = mafdr(mc_table.(SexVarName).SD.First12hr.NREM.pValue,'BHFDR',true);
	
	p_vals.SDexcSD.LP.NREM.Posthocs.(['WT_'  Sexes{1} 'vs' Sexes{2}]).Adjusted = mc_table.(SexVarName).SD.First12hr.NREM.pValue(find(mc_table.(SexVarName).SD.First12hr.NREM.Genotype=='WT'));
	p_vals.SDexcSD.LP.NREM.Posthocs.(['Mut_' Sexes{1} 'vs' Sexes{2}]).Adjusted = mc_table.(SexVarName).SD.First12hr.NREM.pValue(find(mc_table.(SexVarName).SD.First12hr.NREM.Genotype=='Mut'));

else 
	PostHoc_bool.(SexVarName).SDFirst12hr.N = 0; 
	mc_table.(SexVarName).SD.First12hr.NREM = [];
	p_vals.SDexcSD.LP.NREM.Posthocs.(['WT_'  Sexes{1} 'vs' Sexes{2}]).Adjusted = zeros(0);
	p_vals.SDexcSD.LP.NREM.Posthocs.(['Mut_' Sexes{1} 'vs' Sexes{2}]).Adjusted = zeros(0);

end 

% - REM - 
% - Genotype - 
p_vals.SDexcSD.LP.REM.Posthocs.([Sexes{2} '_WTvsMut']) = table;
p_vals.SDexcSD.LP.REM.Posthocs.([Sexes{1} '_WTvsMut']) = table;
if RepMeas_ANOVA_tables.SD.First12hr.REM.pValue('Genotype') < 0.05 | RepMeas_ANOVA_tables.SD.First12hr.REM.pValue(interaction_term_text) < 0.05 | RepMeas_ANOVA_tables.SD.First12hr.REM.pValue('Genotype:Time') < 0.05
	PostHoc_bool.Genotype.SDFirst12hr.R = 1;
	
	mc_table.Genotype.SD.First12hr.REM = multcompare(RM_Models.SD.LP.REM,'Genotype','By','Sex','ComparisonType','lsd');
	[~, ia] = unique(mc_table.Genotype.SD.First12hr.REM.Sex(:), 'stable');  % remove duplicate rows
	mc_table.Genotype.SD.First12hr.REM = mc_table.Genotype.SD.First12hr.REM(ia,:);

	% unadjusted
	p_vals.SDexcSD.LP.REM.Posthocs.([Sexes{2} '_WTvsMut']).Unadjusted = mc_table.Genotype.SD.First12hr.REM.pValue(find(mc_table.Genotype.SD.First12hr.REM.Sex==Sexes{2}));
	p_vals.SDexcSD.LP.REM.Posthocs.([Sexes{1} '_WTvsMut']).Unadjusted = mc_table.Genotype.SD.First12hr.REM.pValue(find(mc_table.Genotype.SD.First12hr.REM.Sex==Sexes{1}));

	% adjust
	mc_table.Genotype.SD.First12hr.REM.pValue = mafdr(mc_table.Genotype.SD.First12hr.REM.pValue,'BHFDR',true);
	
	p_vals.SDexcSD.LP.REM.Posthocs.([Sexes{2} '_WTvsMut']).Adjusted = mc_table.Genotype.SD.First12hr.REM.pValue(find(mc_table.Genotype.SD.First12hr.REM.Sex==Sexes{2}));
	p_vals.SDexcSD.LP.REM.Posthocs.([Sexes{1} '_WTvsMut']).Adjusted = mc_table.Genotype.SD.First12hr.REM.pValue(find(mc_table.Genotype.SD.First12hr.REM.Sex==Sexes{1}));

else 
	PostHoc_bool.Genotype.SDFirst12hr.R = 0;
	mc_table.Genotype.SD.First12hr.REM  = [];
	p_vals.SDexcSD.LP.REM.Posthocs.([Sexes{2} '_WTvsMut']).Adjusted = zeros(0);
	p_vals.SDexcSD.LP.REM.Posthocs.([Sexes{1} '_WTvsMut']).Adjusted = zeros(0);
end 

% - Sex - 
p_vals.SDexcSD.LP.REM.Posthocs.(['WT_'  Sexes{1} 'vs' Sexes{2}]) = table;
p_vals.SDexcSD.LP.REM.Posthocs.(['Mut_' Sexes{1} 'vs' Sexes{2}]) = table;
if RepMeas_ANOVA_tables.SD.First12hr.REM.pValue(SexVarName) < 0.05 | RepMeas_ANOVA_tables.SD.First12hr.REM.pValue(interaction_term_text) < 0.05 | RepMeas_ANOVA_tables.SD.First12hr.REM.pValue([SexVarName ':Time']) < 0.05
	PostHoc_bool.(SexVarName).SDFirst12hr.R = 1;
	
	mc_table.(SexVarName).SD.First12hr.REM = multcompare(RM_Models.SD.LP.REM,'Sex','By','Genotype','ComparisonType','lsd');
	[~, ia] = unique(mc_table.(SexVarName).SD.First12hr.REM.Genotype(:), 'stable');  % remove duplicate rows
	mc_table.(SexVarName).SD.First12hr.REM = mc_table.(SexVarName).SD.First12hr.REM(ia,:);

	%unadjusted
	p_vals.SDexcSD.LP.REM.Posthocs.(['WT_'  Sexes{1} 'vs' Sexes{2}]).Unadjusted = mc_table.(SexVarName).SD.First12hr.REM.pValue(find(mc_table.(SexVarName).SD.First12hr.REM.Genotype=='WT'));
	p_vals.SDexcSD.LP.REM.Posthocs.(['Mut_' Sexes{1} 'vs' Sexes{2}]).Unadjusted = mc_table.(SexVarName).SD.First12hr.REM.pValue(find(mc_table.(SexVarName).SD.First12hr.REM.Genotype=='Mut'));

	% adjust
	mc_table.(SexVarName).SD.First12hr.REM.pValue = mafdr(mc_table.(SexVarName).SD.First12hr.REM.pValue,'BHFDR',true);
	
	p_vals.SDexcSD.LP.REM.Posthocs.(['WT_'  Sexes{1} 'vs' Sexes{2}]).Adjusted = mc_table.(SexVarName).SD.First12hr.REM.pValue(find(mc_table.(SexVarName).SD.First12hr.REM.Genotype=='WT'));
	p_vals.SDexcSD.LP.REM.Posthocs.(['Mut_' Sexes{1} 'vs' Sexes{2}]).Adjusted = mc_table.(SexVarName).SD.First12hr.REM.pValue(find(mc_table.(SexVarName).SD.First12hr.REM.Genotype=='Mut'));

else 
	PostHoc_bool.(SexVarName).SDFirst12hr.R = 0;
	mc_table.(SexVarName).SD.First12hr.REM  = [];
	p_vals.SDexcSD.LP.REM.Posthocs.(['WT_'  Sexes{1} 'vs' Sexes{2}]).Adjusted = zeros(0);
	p_vals.SDexcSD.LP.REM.Posthocs.(['Mut_' Sexes{1} 'vs' Sexes{2}]).Adjusted = zeros(0);


end 

% ------ SD last 12 hours ------
% -- Wake --
% - Genotype 
p_vals.SD.DP.Wake.Posthocs.([Sexes{2} '_WTvsMut']) = table;
p_vals.SD.DP.Wake.Posthocs.([Sexes{1} '_WTvsMut']) = table;
if RepMeas_ANOVA_tables.SD.Last12hr.Wake.pValue('Genotype') < 0.05 | RepMeas_ANOVA_tables.SD.Last12hr.Wake.pValue(interaction_term_text) < 0.05 | RepMeas_ANOVA_tables.SD.Last12hr.Wake.pValue('Genotype:Time') < 0.05
	PostHoc_bool.Genotype.SDLast12hr.W = 1;
	
	mc_table.Genotype.SD.Last12hr.Wake = multcompare(RM_Models.SD.DP.Wake,'Genotype','By','Sex','ComparisonType','lsd');  % For * symbols on graphs
	[~, ia] = unique(mc_table.Genotype.SD.Last12hr.Wake.Sex(:), 'stable');  % remove duplicate rows
	mc_table.Genotype.SD.Last12hr.Wake = mc_table.Genotype.SD.Last12hr.Wake(ia,:);

	% unadjusted
	p_vals.SD.DP.Wake.Posthocs.([Sexes{2} '_WTvsMut']).Unadjusted = mc_table.Genotype.SD.Last12hr.Wake.pValue(find(mc_table.Genotype.SD.Last12hr.Wake.Sex==Sexes{2}));
	p_vals.SD.DP.Wake.Posthocs.([Sexes{1} '_WTvsMut']).Unadjusted = mc_table.Genotype.SD.Last12hr.Wake.pValue(find(mc_table.Genotype.SD.Last12hr.Wake.Sex==Sexes{1}));

	% adjust
	mc_table.Genotype.SD.Last12hr.Wake.pValue = mafdr(mc_table.Genotype.SD.Last12hr.Wake.pValue,'BHFDR',true);
	
	p_vals.SD.DP.Wake.Posthocs.([Sexes{2} '_WTvsMut']).Adjusted = mc_table.Genotype.SD.Last12hr.Wake.pValue(find(mc_table.Genotype.SD.Last12hr.Wake.Sex==Sexes{2}));
	p_vals.SD.DP.Wake.Posthocs.([Sexes{1} '_WTvsMut']).Adjusted = mc_table.Genotype.SD.Last12hr.Wake.pValue(find(mc_table.Genotype.SD.Last12hr.Wake.Sex==Sexes{1}));


else 
	PostHoc_bool.Genotype.SDLast12hr.W = 0;
	mc_table.Genotype.SD.Last12hr.Wake = [];
	p_vals.SD.DP.Wake.Posthocs.([Sexes{2} '_WTvsMut']).Adjusted = zeros(0);
	p_vals.SD.DP.Wake.Posthocs.([Sexes{1} '_WTvsMut']).Adjusted = zeros(0);

end 

% - Sex - 
p_vals.SD.DP.Wake.Posthocs.(['WT_'  Sexes{1} 'vs' Sexes{2}]) = table;
p_vals.SD.DP.Wake.Posthocs.(['Mut_' Sexes{1} 'vs' Sexes{2}]) = table;
if RepMeas_ANOVA_tables.SD.Last12hr.Wake.pValue(SexVarName) < 0.05 | RepMeas_ANOVA_tables.SD.Last12hr.Wake.pValue(interaction_term_text) < 0.05 | RepMeas_ANOVA_tables.SD.Last12hr.Wake.pValue([SexVarName ':Time']) < 0.05
	PostHoc_bool.(SexVarName).SDLast12hr.W = 1;
	
	mc_table.(SexVarName).SD.Last12hr.Wake = multcompare(RM_Models.SD.DP.Wake,'Sex','By','Genotype','ComparisonType','lsd');  % For # symbols on graphs
	[~, ia] = unique(mc_table.(SexVarName).SD.Last12hr.Wake.Genotype(:), 'stable');  % remove duplicate rows
	mc_table.(SexVarName).SD.Last12hr.Wake = mc_table.(SexVarName).SD.Last12hr.Wake(ia,:);

	% unadjusted
	p_vals.SD.DP.Wake.Posthocs.(['WT_'  Sexes{1} 'vs' Sexes{2}]).Unadjusted = mc_table.(SexVarName).SD.Last12hr.Wake.pValue(find(mc_table.(SexVarName).SD.Last12hr.Wake.Genotype=='WT'));
	p_vals.SD.DP.Wake.Posthocs.(['Mut_' Sexes{1} 'vs' Sexes{2}]).Unadjusted = mc_table.(SexVarName).SD.Last12hr.Wake.pValue(find(mc_table.(SexVarName).SD.Last12hr.Wake.Genotype=='Mut'));

	% adjust
	mc_table.(SexVarName).SD.Last12hr.Wake.pValue = mafdr(mc_table.(SexVarName).SD.Last12hr.Wake.pValue,'BHFDR',true);
	
	p_vals.SD.DP.Wake.Posthocs.(['WT_'  Sexes{1} 'vs' Sexes{2}]).Adjusted = mc_table.(SexVarName).SD.Last12hr.Wake.pValue(find(mc_table.(SexVarName).SD.Last12hr.Wake.Genotype=='WT'));
	p_vals.SD.DP.Wake.Posthocs.(['Mut_' Sexes{1} 'vs' Sexes{2}]).Adjusted = mc_table.(SexVarName).SD.Last12hr.Wake.pValue(find(mc_table.(SexVarName).SD.Last12hr.Wake.Genotype=='Mut'));


else 
	PostHoc_bool.(SexVarName).SDLast12hr.W = 0;
	mc_table.(SexVarName).SD.Last12hr.Wake = [];
	p_vals.SD.DP.Wake.Posthocs.(['WT_'  Sexes{1} 'vs' Sexes{2}]).Adjusted = zeros(0);
	p_vals.SD.DP.Wake.Posthocs.(['Mut_' Sexes{1} 'vs' Sexes{2}]).Adjusted = zeros(0);

end 

% -- NREM --
% - Genotype -
p_vals.SD.DP.NREM.Posthocs.([Sexes{2} '_WTvsMut']) = table;
p_vals.SD.DP.NREM.Posthocs.([Sexes{1} '_WTvsMut']) = table;
if RepMeas_ANOVA_tables.SD.Last12hr.NREM.pValue('Genotype') < 0.05 | RepMeas_ANOVA_tables.SD.Last12hr.NREM.pValue(interaction_term_text) < 0.05 | RepMeas_ANOVA_tables.SD.Last12hr.NREM.pValue('Genotype:Time') < 0.05
	PostHoc_bool.Genotype.SDLast12hr.N = 1; 
	
	mc_table.Genotype.SD.Last12hr.NREM = multcompare(RM_Models.SD.DP.NREM,'Genotype','By','Sex','ComparisonType','lsd'); % For * symbols on graphs
	[~, ia] = unique(mc_table.Genotype.SD.Last12hr.NREM.Sex(:), 'stable');  % remove duplicate rows
	mc_table.Genotype.SD.Last12hr.NREM = mc_table.Genotype.SD.Last12hr.NREM(ia,:);

	% unadjusted
	p_vals.SD.DP.NREM.Posthocs.([Sexes{2} '_WTvsMut']).Unadjusted = mc_table.Genotype.SD.Last12hr.NREM.pValue(find(mc_table.Genotype.SD.Last12hr.NREM.Sex==Sexes{2}));
	p_vals.SD.DP.NREM.Posthocs.([Sexes{1} '_WTvsMut']).Unadjusted = mc_table.Genotype.SD.Last12hr.NREM.pValue(find(mc_table.Genotype.SD.Last12hr.NREM.Sex==Sexes{1}));

	% adjust
	mc_table.Genotype.SD.Last12hr.NREM.pValue = mafdr(mc_table.Genotype.SD.Last12hr.NREM.pValue,'BHFDR',true);
	
	p_vals.SD.DP.NREM.Posthocs.([Sexes{2} '_WTvsMut']).Adjusted = mc_table.Genotype.SD.Last12hr.NREM.pValue(find(mc_table.Genotype.SD.Last12hr.NREM.Sex==Sexes{2}));
	p_vals.SD.DP.NREM.Posthocs.([Sexes{1} '_WTvsMut']).Adjusted = mc_table.Genotype.SD.Last12hr.NREM.pValue(find(mc_table.Genotype.SD.Last12hr.NREM.Sex==Sexes{1}));

else 
	PostHoc_bool.Genotype.SDLast12hr.N = 0; 
	mc_table.Genotype.SD.Last12hr.NREM = [];
	p_vals.SD.DP.NREM.Posthocs.([Sexes{2} '_WTvsMut']).Adjusted = zeros(0);
	p_vals.SD.DP.NREM.Posthocs.([Sexes{1} '_WTvsMut']).Adjusted = zeros(0);
end 

% - Sex - 
p_vals.SD.DP.NREM.Posthocs.(['WT_'  Sexes{1} 'vs' Sexes{2}]) = table;
p_vals.SD.DP.NREM.Posthocs.(['Mut_' Sexes{1} 'vs' Sexes{2}]) = table;
if RepMeas_ANOVA_tables.SD.Last12hr.NREM.pValue(SexVarName) < 0.05 | RepMeas_ANOVA_tables.SD.Last12hr.NREM.pValue(interaction_term_text) < 0.05 | RepMeas_ANOVA_tables.SD.Last12hr.NREM.pValue([SexVarName ':Time']) < 0.05
	PostHoc_bool.(SexVarName).SDLast12hr.N = 1; 
	

	mc_table.(SexVarName).SD.Last12hr.NREM = multcompare(RM_Models.SD.DP.NREM,'Sex','By','Genotype','ComparisonType','lsd'); % For * symbols on graphs
	[~, ia] = unique(mc_table.(SexVarName).SD.Last12hr.NREM.Genotype(:), 'stable');  % remove duplicate rows
	mc_table.(SexVarName).SD.Last12hr.NREM = mc_table.(SexVarName).SD.Last12hr.NREM(ia,:);

	% unadjustd
	p_vals.SD.DP.NREM.Posthocs.(['WT_'  Sexes{1} 'vs' Sexes{2}]).Unadjusted = mc_table.(SexVarName).SD.Last12hr.NREM.pValue(find(mc_table.(SexVarName).SD.Last12hr.NREM.Genotype=='WT'));
	p_vals.SD.DP.NREM.Posthocs.(['Mut_' Sexes{1} 'vs' Sexes{2}]).Unadjusted = mc_table.(SexVarName).SD.Last12hr.NREM.pValue(find(mc_table.(SexVarName).SD.Last12hr.NREM.Genotype=='Mut'));

	% adjust
	mc_table.(SexVarName).SD.Last12hr.NREM.pValue = mafdr(mc_table.(SexVarName).SD.Last12hr.NREM.pValue,'BHFDR',true);
	
	p_vals.SD.DP.NREM.Posthocs.(['WT_'  Sexes{1} 'vs' Sexes{2}]).Adjusted = mc_table.(SexVarName).SD.Last12hr.NREM.pValue(find(mc_table.(SexVarName).SD.Last12hr.NREM.Genotype=='WT'));
	p_vals.SD.DP.NREM.Posthocs.(['Mut_' Sexes{1} 'vs' Sexes{2}]).Adjusted = mc_table.(SexVarName).SD.Last12hr.NREM.pValue(find(mc_table.(SexVarName).SD.Last12hr.NREM.Genotype=='Mut'));

else 
	PostHoc_bool.(SexVarName).SDLast12hr.N = 0; 
	mc_table.(SexVarName).SD.Last12hr.NREM = [];
	p_vals.SD.DP.NREM.Posthocs.(['WT_'  Sexes{1} 'vs' Sexes{2}]).Adjusted = zeros(0);
	p_vals.SD.DP.NREM.Posthocs.(['Mut_' Sexes{1} 'vs' Sexes{2}]).Adjusted = zeros(0);

end 

% - REM - 
% - Genotype - 
p_vals.SD.DP.REM.Posthocs.([Sexes{2} '_WTvsMut']) = table;
p_vals.SD.DP.REM.Posthocs.([Sexes{1} '_WTvsMut']) = table;
if RepMeas_ANOVA_tables.SD.Last12hr.REM.pValue('Genotype') < 0.05 | RepMeas_ANOVA_tables.SD.Last12hr.REM.pValue(interaction_term_text) < 0.05 | RepMeas_ANOVA_tables.SD.Last12hr.REM.pValue('Genotype:Time') < 0.05
	PostHoc_bool.Genotype.SDLast12hr.R = 1;
	
	mc_table.Genotype.SD.Last12hr.REM = multcompare(RM_Models.SD.DP.REM,'Genotype','By','Sex','ComparisonType','lsd');
	[~, ia] = unique(mc_table.Genotype.SD.Last12hr.REM.Sex(:), 'stable');  % remove duplicate rows
	mc_table.Genotype.SD.Last12hr.REM = mc_table.Genotype.SD.Last12hr.REM(ia,:);

	% unadjusted
	p_vals.SD.DP.REM.Posthocs.([Sexes{2} '_WTvsMut']).Unadjusted = mc_table.Genotype.SD.Last12hr.REM.pValue(find(mc_table.Genotype.SD.Last12hr.REM.Sex==Sexes{2}));
	p_vals.SD.DP.REM.Posthocs.([Sexes{1} '_WTvsMut']).Unadjusted = mc_table.Genotype.SD.Last12hr.REM.pValue(find(mc_table.Genotype.SD.Last12hr.REM.Sex==Sexes{1}));

	% adjust
	mc_table.Genotype.SD.Last12hr.REM.pValue = mafdr(mc_table.Genotype.SD.Last12hr.REM.pValue,'BHFDR',true);
	
	p_vals.SD.DP.REM.Posthocs.([Sexes{2} '_WTvsMut']).Adjusted = mc_table.Genotype.SD.Last12hr.REM.pValue(find(mc_table.Genotype.SD.Last12hr.REM.Sex==Sexes{2}));
	p_vals.SD.DP.REM.Posthocs.([Sexes{1} '_WTvsMut']).Adjusted = mc_table.Genotype.SD.Last12hr.REM.pValue(find(mc_table.Genotype.SD.Last12hr.REM.Sex==Sexes{1}));



else 
	PostHoc_bool.Genotype.SDLast12hr.R = 0;
	mc_table.Genotype.SD.Last12hr.REM  = [];
	p_vals.SD.DP.REM.Posthocs.([Sexes{2} '_WTvsMut']).Adjusted = zeros(0);
	p_vals.SD.DP.REM.Posthocs.([Sexes{1} '_WTvsMut']).Adjusted = zeros(0);

end 

% - Sex - 
p_vals.SD.DP.REM.Posthocs.(['WT_'  Sexes{1} 'vs' Sexes{2}]) = table;
p_vals.SD.DP.REM.Posthocs.(['Mut_' Sexes{1} 'vs' Sexes{2}]) = table;
if RepMeas_ANOVA_tables.SD.Last12hr.REM.pValue(SexVarName) < 0.05 | RepMeas_ANOVA_tables.SD.Last12hr.REM.pValue(interaction_term_text) < 0.05 | RepMeas_ANOVA_tables.SD.Last12hr.REM.pValue([SexVarName ':Time']) < 0.05
	PostHoc_bool.(SexVarName).SDLast12hr.R = 1;
	
	mc_table.(SexVarName).SD.Last12hr.REM = multcompare(RM_Models.SD.DP.REM,'Sex','By','Genotype','ComparisonType','lsd');
	[~, ia] = unique(mc_table.(SexVarName).SD.Last12hr.REM.Genotype(:), 'stable');  % remove duplicate rows
	mc_table.(SexVarName).SD.Last12hr.REM = mc_table.(SexVarName).SD.Last12hr.REM(ia,:);

	% unadjusted
	p_vals.SD.DP.REM.Posthocs.(['WT_'  Sexes{1} 'vs' Sexes{2}]).Unadjusted = mc_table.(SexVarName).SD.Last12hr.REM.pValue(find(mc_table.(SexVarName).SD.Last12hr.REM.Genotype=='WT'));
	p_vals.SD.DP.REM.Posthocs.(['Mut_' Sexes{1} 'vs' Sexes{2}]).Unadjusted = mc_table.(SexVarName).SD.Last12hr.REM.pValue(find(mc_table.(SexVarName).SD.Last12hr.REM.Genotype=='Mut'));

	% adjust
	mc_table.(SexVarName).SD.Last12hr.REM.pValue = mafdr(mc_table.(SexVarName).SD.Last12hr.REM.pValue,'BHFDR',true);
	
	p_vals.SD.DP.REM.Posthocs.(['WT_'  Sexes{1} 'vs' Sexes{2}]).Adjusted = mc_table.(SexVarName).SD.Last12hr.REM.pValue(find(mc_table.(SexVarName).SD.Last12hr.REM.Genotype=='WT'));
	p_vals.SD.DP.REM.Posthocs.(['Mut_' Sexes{1} 'vs' Sexes{2}]).Adjusted = mc_table.(SexVarName).SD.Last12hr.REM.pValue(find(mc_table.(SexVarName).SD.Last12hr.REM.Genotype=='Mut'));

else 
	PostHoc_bool.(SexVarName).SDLast12hr.R = 0;
	mc_table.(SexVarName).SD.Last12hr.REM  = [];
	p_vals.SD.DP.REM.Posthocs.(['WT_'  Sexes{1} 'vs' Sexes{2}]).Adjusted = zeros(0);
	p_vals.SD.DP.REM.Posthocs.(['Mut_' Sexes{1} 'vs' Sexes{2}]).Adjusted = zeros(0);
end 






% ------- SD First 6hr Dark Period -----------
if isfield(RepMeas_ANOVA_tables.SD,'First6hrDP')
	% - Wake - 
	if RepMeas_ANOVA_tables.SD.First6hrDP.Wake.pValue('Genotype') < 0.05 | RepMeas_ANOVA_tables.SD.First6hrDP.Wake.pValue(interaction_term_text) < 0.05 RepMeas_ANOVA_tables.SD.First6hrDP.Wake.pValue('Genotype:Time') < 0.05
		PostHoc_bool.SDFirst6hrDP.W = 1;
		mc_table.SD.First6hrDP.Wake 	   = multcompare(RM_Models.SD.First6hrDP.Wake,'Genotype','By','Time' ,'ComparisonType','lsd');
		mc_table.SD.First6hrDP.Wake.pValue = mafdr(mc_table.SD.First6hrDP.Wake.pValue,'BHFDR',true);
		[~, ia] = unique(mc_table.SD.First6hrDP.Wake.Time(:), 'stable');  % remove duplicate rows
		mc_table.SD.First6hrDP.Wake = mc_table.SD.First6hrDP.Wake(ia,:);
		xlocs_sig_effects.SD.WT_vs_Mut_First6hrDP_W = double(string(mc_table.SD.First6hrDP.Wake.Time(mc_table.SD.First6hrDP.Wake.pValue<0.05))); % these are times (h), not indices
	else 
		PostHoc_bool.Genotype.SDFirst6hrDP.W = 0;
		mc_table.Genotype.SD.First6hrDP.Wake = [];
	end 

	% add post-hocs for sex here 

	% - NREM - 
	if RepMeas_ANOVA_tables.SD.First6hrDP.NREM.pValue('Genotype') < 0.05 | RepMeas_ANOVA_tables.SD.First6hrDP.NREM.pValue(interaction_term_text) < 0.05 | RepMeas_ANOVA_tables.SD.First6hrDP.NREM.pValue('Genotype:Time') < 0.05
		PostHoc_bool.SDFirst6hrDP.N = 1;
		mc_table.SD.First6hrDP.NREM 	   = multcompare(RM_Models.SD.First6hrDP.NREM,'Genotype','By','Time','ComparisonType','lsd');
		mc_table.SD.First6hrDP.NREM.pValue = mafdr(mc_table.SD.First6hrDP.NREM.pValue,'BHFDR',true);
		[~, ia] = unique(mc_table.SD.First6hrDP.NREM.Time(:), 'stable');  % remove duplicate rows
		mc_table.SD.First6hrDP.NREM = mc_table.SD.First6hrDP.NREM(ia,:);
		xlocs_sig_effects.SD.WT_vs_Mut_First6hrDP_N = double(string(mc_table.SD.First6hrDP.NREM.Time(mc_table.SD.First6hrDP.NREM.pValue<0.05))); % these are times (h), not indices
	else 
		PostHoc_bool.SDFirst6hrDP.N = 0;
		mc_table.SD.First6hrDP.NREM = [];
	end 

	% add post-hocs for sex here

	% - REM - 
	if RepMeas_ANOVA_tables.SD.First6hrDP.REM.pValue('Genotype') < 0.05 | RepMeas_ANOVA_tables.SD.First6hrDP.REM.pValue(interaction_term_text) < 0.05 | RepMeas_ANOVA_tables.SD.First6hrDP.REM.pValue('Genotype:Time') < 0.05
		PostHoc_bool.SDFirst6hrDP.R = 1;
		mc_table.SD.First6hrDP.REM 		  = multcompare(RM_Models.SD.First6hrDP.REM,'Genotype','By','Time' ,'ComparisonType','lsd');
		mc_table.SD.First6hrDP.REM.pValue = mafdr(mc_table.SD.First6hrDP.REM.pValue,'BHFDR',true);
		[~, ia] = unique(mc_table.SD.First6hrDP.REM.Time(:), 'stable');  % remove duplicate rows
		mc_table.SD.First6hrDP.REM = mc_table.SD.First6hrDP.REM(ia,:);
		xlocs_sig_effects.SD.WT_vs_Mut_First6hrDP_R = double(string(mc_table.SD.First6hrDP.REM.Time(mc_table.SD.First6hrDP.REM.pValue<0.05))); % these are times (h), not indices
	else 
		PostHoc_bool.SDFirst6hrDP.R = 0;
		mc_table.SD.First6hrDP.REM  = [];
	end 
	% add post-hocs for sex here

end % end of SD first 6hr Dark Period



% ------- SD Last 6hr Dark Period -----------
if isfield(RepMeas_ANOVA_tables.SD,'Last6hrDP')
	% - Wake - 
	if RepMeas_ANOVA_tables.SD.Last6hrDP.Wake.pValue('Genotype') < 0.05 | RepMeas_ANOVA_tables.SD.Last6hrDP.Wake.pValue(interaction_term_text) < 0.05 | RepMeas_ANOVA_tables.SD.Last6hrDP.Wake.pValue('Genotype:Time') < 0.05
		PostHoc_bool.SDLast6hrDP.W = 1;
		mc_table.SD.Last6hrDP.Wake 		  = multcompare(RM_Models.SD.Last6hrDP.Wake,'Genotype','By','Time' ,'ComparisonType','lsd');
		mc_table.SD.Last6hrDP.Wake.pValue = mafdr(mc_table.SD.Last6hrDP.Wake.pValue,'BHFDR',true);
		[~, ia] = unique(mc_table.SD.Last6hrDP.Wake.Time(:), 'stable');  % remove duplicate rows
		mc_table.SD.Last6hrDP.Wake = mc_table.SD.Last6hrDP.Wake(ia,:);
		xlocs_sig_effects.SD.WT_vs_Mut_Last6hrDP_W = double(string(mc_table.SD.Last6hrDP.Wake.Time(mc_table.SD.Last6hrDP.Wake.pValue<0.05))); % these are times (h), not indices
	else 
		PostHoc_bool.SDLast6hrDP.W = 0;
		mc_table.SD.Last6hrDP.Wake = [];
	end 
	% add post-hocs for sex here

	% - NREM - 
	if RepMeas_ANOVA_tables.SD.Last6hrDP.NREM.pValue('Genotype') < 0.05 | RepMeas_ANOVA_tables.SD.Last6hrDP.NREM.pValue(interaction_term_text) < 0.05 | RepMeas_ANOVA_tables.SD.Last6hrDP.NREM.pValue('Genotype:Time') < 0.05
		PostHoc_bool.SDLast6hrDP.N = 1;
		mc_table.SD.Last6hrDP.NREM 		  = multcompare(RM_Models.SD.Last6hrDP.NREM,'Genotype','By','Time' ,'ComparisonType','lsd');
		mc_table.SD.Last6hrDP.NREM.pValue = mafdr(mc_table.SD.Last6hrDP.NREM.pValue,'BHFDR',true);
		[~, ia] = unique(mc_table.SD.Last6hrDP.NREM.Time(:), 'stable');  % remove duplicate rows
		mc_table.SD.Last6hrDP.NREM = mc_table.SD.Last6hrDP.NREM(ia,:);
		xlocs_sig_effects.SD.WT_vs_Mut_Last6hrDP_N = double(string(mc_table.SD.Last6hrDP.NREM.Time(mc_table.SD.Last6hrDP.NREM.pValue<0.05))); % these are times (h), not indices
	else 
		PostHoc_bool.SDLast6hrDP.N = 0;
		mc_table.SD.Last6hrDP.NREM = [];
	end 
	% add post-hocs for sex here
	
	% - REM - 
	if RepMeas_ANOVA_tables.SD.Last6hrDP.REM.pValue('Genotype') < 0.05 | RepMeas_ANOVA_tables.SD.Last6hrDP.REM.pValue(interaction_term_text) < 0.05 | RepMeas_ANOVA_tables.SD.Last6hrDP.REM.pValue('Genotype:Time') < 0.05
		PostHoc_bool.SDLast6hrDP.R = 1;
		mc_table.SD.Last6hrDP.REM 		 = multcompare(RM_Models.SD.Last6hrDP.REM,'Genotype','By','Time' ,'ComparisonType','lsd');
		mc_table.SD.Last6hrDP.REM.pValue = mafdr(mc_table.SD.Last6hrDP.REM.pValue,'BHFDR',true);
		[~, ia] = unique(mc_table.SD.Last6hrDP.REM.Time(:), 'stable');  % remove duplicate rows
		mc_table.SD.Last6hrDP.REM = mc_table.SD.Last6hrDP.REM(ia,:);
		xlocs_sig_effects.SD.WT_vs_Mut_Last6hrDP_R = double(string(mc_table.SD.Last6hrDP.REM.Time(mc_table.SD.Last6hrDP.REM.pValue<0.05))); % these are times (h), not indices
	else 
		PostHoc_bool.SDLast6hrDP.R = 0;
		mc_table.SD.Last6hrDP.REM  = [];
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


% ------ 3-panel figure for all of BL Males---------
H=figure;
H.Name = ['Three-Panel BL Time in State Figure ' Sexes{1}];
H.Renderer = 'painters';
H.Position = [429 118 678 698];

t= tiledlayout(3,1,'TileSpacing','compact','Padding','compact','TileIndexing', 'columnmajor');
nexttile

% -- BL W -- 
rectangle('Position',LPboxpos,'FaceColor','w','EdgeColor','k','LineWidth',2);     % LP rectangle
hold on
rectangle('Position',DPboxpos,'FaceColor','k','EdgeColor','k','LineWidth',2);    % DP rectangle
e1=errorbar(WakeWTBLRelmean_M,WakeWTBLRelstd_M,'Color','k', 'LineWidth',mylinewidth,'Marker','o','MarkerFaceColor','k','MarkerSize',mymarkersize,'CapSize',0);
e2=errorbar(WakeMutBLRelmean_M,WakeMutBLRelstd_M,'Color','r','LineWidth',mylinewidth,'Marker','o','MarkerFaceColor','r','MarkerSize',mymarkersize,'CapSize',0);
hold off 
ax=H.Children.Children;
set(ax, 'Layer','top')
%title('Wake, Baseline');
ylabel({'WAKE','Time in Wake (% TRT)'});
%xlabel('Time (hours)');
l = legend([e1 e2],LegendLabels);
%ax=gca;
%l=legend(ax,LegendLabels);
l.Box = 'off';
l.Location = 'None';
%l.Position = [0.7326 0.1649 0.1634 0.0798];
l.Location = 'SouthEast';
l.FontSize = 12;
%set(gca, 'XTick',[1 4 7 10 13 16 19 22])
set(ax,'XTick',[])
%set(gca, 'XTickLabels',{'1' '4' '7' '10' '13' '16' '19' '22'})
xlim([0 24])
ylim([-5 100])
ax.FontSize = myfontsize;
ax.LineWidth = mylinewidth;
% % -- Add stars if sig diffs (LP)
% if PostHoc_bool.BLFirst12hr.W
% 	hold on 
% 	bool_xlocsBL_LP_sig_values = ismember(e1.XData,xlocs_sig_effects.BL.WT_vs_Mut_LP_W);
% 	%bool_xlocsMutBL_LP_sig_values = ismember(e2.XData,xlocs_sig_effects.WTBL_vs_WTSD_LP_W);
% 	largest_y1_vals = e1.YData(bool_xlocsBL_LP_sig_values) + e1.YPositiveDelta(bool_xlocsBL_LP_sig_values);
% 	largest_y2_vals = e2.YData(bool_xlocsBL_LP_sig_values) + e2.YPositiveDelta(bool_xlocsBL_LP_sig_values);
% 	ylocs_sig_effects = max([largest_y1_vals; largest_y2_vals]) + 10*ones(1,length(largest_y1_vals));
% 	ylocs_sig_effects(ylocs_sig_effects>ax.YLim(2))=ax.YLim(2);
% 	text(xlocs_sig_effects.BL.WT_vs_Mut_LP_W, ylocs_sig_effects, '*', 'FontSize',astr_size,'HorizontalAlignment','center','VerticalAlignment','middle');
% 	hold off 
% end 
% % -- Add stars if sig diffs (DP)
% if PostHoc_bool.BLLast12hr.W
% 	hold on 
% 	bool_xlocsBL_DP_sig_values = ismember(e1.XData,xlocs_sig_effects.BL.WT_vs_Mut_DP_W);
% 	%bool_xlocsMutBL_LP_sig_values = ismember(e2.XData,xlocs_sig_effects.WTBL_vs_WTSD_LP_W);
% 	largest_y1_vals = e1.YData(bool_xlocsBL_DP_sig_values) + e1.YPositiveDelta(bool_xlocsBL_DP_sig_values);
% 	largest_y2_vals = e2.YData(bool_xlocsBL_DP_sig_values) + e2.YPositiveDelta(bool_xlocsBL_DP_sig_values);
% 	ylocs_sig_effects = max([largest_y1_vals; largest_y2_vals]) + 10*ones(1,length(largest_y1_vals));
% 	ylocs_sig_effects(ylocs_sig_effects>ax.YLim(2))=ax.YLim(2);
% 	text(xlocs_sig_effects.BL.WT_vs_Mut_DP_W, ylocs_sig_effects, '*', 'FontSize',astr_size,'HorizontalAlignment','center','VerticalAlignment','middle');
% 	hold off 
% end 



% -- BL NREM Males--
nexttile
rectangle('Position',LPboxpos,'FaceColor','w','EdgeColor','k','LineWidth',2);     % LP rectangle
hold on
rectangle('Position',DPboxpos,'FaceColor','k','EdgeColor','k','LineWidth',2);    % DP rectangle

e1=errorbar(NREMMutBLRelmean_M,NREMMutBLRelstd_M,'Color','r','LineWidth',mylinewidth,'Marker','o','MarkerFaceColor','r','MarkerSize',mymarkersize,'CapSize',0);
e2=errorbar(NREMWTBLRelmean_M,NREMWTBLRelstd_M,'Color','k','LineWidth',mylinewidth,'Marker','o','MarkerFaceColor','k','MarkerSize',mymarkersize,'CapSize',0);
hold off
ax=H.Children.Children(1);
set(ax, 'Layer','top')
%title('NREM Baseline');
% legend('Mutant','WT');
% legend 'boxoff';
%set(gca, 'XTick',[1 4 7 10 13 16 19 22])
set(ax,'XTick',[])
%set(gca, 'XTickLabels',{'1' '4' '7' '10' '13' '16' '19' '22'})
xlim([0 24])
ylim([-5 100])
ylabel({'NREM','Time in NREM (% TRT)'});
%xlabel('Time (hours)')

ax.FontSize = myfontsize;
ax.LineWidth = mylinewidth;
% -- Add stars if sig diffs (LP)
% if PostHoc_bool.BLFirst12hr.N
% 	hold on 
% 	bool_xlocsWTBL_LP_sig_values = ismember(e1.XData,xlocs_sig_effects.BL.WT_vs_Mut_LP_N);
% 	%bool_xlocsMutBL_LP_sig_values = ismember(e2.XData,xlocs_sig_effects.WTBL_vs_WTSD_LP_W);
% 	largest_y1_vals = e1.YData(bool_xlocsWTBL_LP_sig_values) + e1.YPositiveDelta(bool_xlocsWTBL_LP_sig_values);
% 	largest_y2_vals = e2.YData(bool_xlocsWTBL_LP_sig_values) + e2.YPositiveDelta(bool_xlocsWTBL_LP_sig_values);
% 	ylocs_sig_effects = max([largest_y1_vals; largest_y2_vals]) + 10*ones(1,length(largest_y1_vals));
% 	ylocs_sig_effects(ylocs_sig_effects>ax.YLim(2))=ax.YLim(2);
% 	text(xlocs_sig_effects.BL.WT_vs_Mut_LP_N, ylocs_sig_effects, '*', 'FontSize',astr_size,'HorizontalAlignment','center','VerticalAlignment','middle');
% 	hold off 
% end 
% % -- Add stars if sig diffs (DP)
% if PostHoc_bool.BLLast12hr.N
% 	hold on 
% 	bool_xlocsWTBL_DP_sig_values = ismember(e1.XData,xlocs_sig_effects.BL.WT_vs_Mut_DP_N);
% 	%bool_xlocsMutBL_LP_sig_values = ismember(e2.XData,xlocs_sig_effects.WTBL_vs_WTSD_LP_W);
% 	largest_y1_vals = e1.YData(bool_xlocsWTBL_DP_sig_values) + e1.YPositiveDelta(bool_xlocsWTBL_DP_sig_values);
% 	largest_y2_vals = e2.YData(bool_xlocsWTBL_DP_sig_values) + e2.YPositiveDelta(bool_xlocsWTBL_DP_sig_values);
% 	ylocs_sig_effects = max([largest_y1_vals; largest_y2_vals]) + 10*ones(1,length(largest_y1_vals));
% 	ylocs_sig_effects(ylocs_sig_effects>ax.YLim(2))=ax.YLim(2);
% 	text(xlocs_sig_effects.BL.WT_vs_Mut_DP_N, ylocs_sig_effects, '*', 'FontSize',astr_size,'HorizontalAlignment','center','VerticalAlignment','middle');
% 	hold off 
% end 
text(0,95, ...
		{'NOTE: These figures do not have *, # or bars added automatically.';'See RepMeas_ANOVA_tables_BothSexes_TIS and'; 'RepMeas_ANOVA_TIS_BothSexes_p_values and PostHocs_TIS_MandF.'}, ...
		'Interpreter','none','fontsize',12)

% -- BL REM Males--
nexttile
rectangle('Position',LPboxposREM,'FaceColor','w','EdgeColor','k','LineWidth',2);     % LP rectangle
hold on
rectangle('Position',DPboxposREM,'FaceColor','k','EdgeColor','k','LineWidth',2);    % DP rectangle

e1=errorbar(REMMutBLRelmean_M,REMMutBLRelstd_M,'Color','r','LineWidth',mylinewidth,'Marker','o','MarkerFaceColor','r','MarkerSize',mymarkersize,'CapSize',0);
e2=errorbar(REMWTBLRelmean_M,REMWTBLRelstd_M,'Color','k','LineWidth',mylinewidth,'Marker','o','MarkerFaceColor','k','MarkerSize',mymarkersize,'CapSize',0);
hold off
ax=H.Children.Children(1);
set(ax, 'Layer','top')
%title('REM Baseline');
% legend('Mutant','WT');
% legend 'boxoff';
set(ax, 'XTick',[1 4 7 10 13 16 19 22])
set(ax, 'XTickLabels',{'1' '4' '7' '10' '13' '16' '19' '22'})
xlim([0 24])
ylim ([-1 30])
xlabel('Time (hours)')
ylabel({'REM','Time in REM (% TRT)'});

ax.FontSize = myfontsize;
ax.LineWidth = mylinewidth;
% -- Add stars if sig diffs (LP)
% if PostHoc_bool.BLFirst12hr.R
% 	hold on 
% 	bool_xlocsWTBL_LP_sig_values = ismember(e1.XData,xlocs_sig_effects.BL.WT_vs_Mut_LP_R);
% 	%bool_xlocsMutBL_LP_sig_values = ismember(e2.XData,xlocs_sig_effects.WTBL_vs_WTSD_LP_W);
% 	largest_y1_vals = e1.YData(bool_xlocsWTBL_LP_sig_values) + e1.YPositiveDelta(bool_xlocsWTBL_LP_sig_values);
% 	largest_y2_vals = e2.YData(bool_xlocsWTBL_LP_sig_values) + e2.YPositiveDelta(bool_xlocsWTBL_LP_sig_values);
% 	ylocs_sig_effects = max([largest_y1_vals; largest_y2_vals]) + 10*ones(1,length(largest_y1_vals));
% 	ylocs_sig_effects(ylocs_sig_effects>ax.YLim(2))=ax.YLim(2);
% 	text(xlocs_sig_effects.BL.WT_vs_Mut_LP_R, ylocs_sig_effects, '*', 'FontSize',astr_size,'HorizontalAlignment','center','VerticalAlignment','middle');
% 	hold off 
% end 
% % -- Add stars if sig diffs (DP)
% if PostHoc_bool.BLLast12hr.R
% 	hold on 
% 	bool_xlocsWTBL_DP_sig_values = ismember(e1.XData,xlocs_sig_effects.BL.WT_vs_Mut_DP_R);
% 	%bool_xlocsMutBL_LP_sig_values = ismember(e2.XData,xlocs_sig_effects.WTBL_vs_WTSD_LP_W);
% 	largest_y1_vals = e1.YData(bool_xlocsWTBL_DP_sig_values) + e1.YPositiveDelta(bool_xlocsWTBL_DP_sig_values);
% 	largest_y2_vals = e2.YData(bool_xlocsWTBL_DP_sig_values) + e2.YPositiveDelta(bool_xlocsWTBL_DP_sig_values);
% 	ylocs_sig_effects = max([largest_y1_vals; largest_y2_vals]) + 10*ones(1,length(largest_y1_vals));
% 	ylocs_sig_effects(ylocs_sig_effects>ax.YLim(2))=ax.YLim(2);
% 	text(xlocs_sig_effects.BL.WT_vs_Mut_DP_R, ylocs_sig_effects, '*', 'FontSize',astr_size,'HorizontalAlignment','center','VerticalAlignment','middle');
% 	hold off 
% end 


% ------ 3-panel figure for all of BL Females---------
H=figure;
H.Name = ['Three-Panel BL Time in State Figure ' Sexes{2}];
H.Renderer = 'painters';
H.Position = [429 118 678 698];

t= tiledlayout(3,1,'TileSpacing','compact','Padding','compact','TileIndexing', 'columnmajor');
nexttile

% -- BL W Females-- 
rectangle('Position',LPboxpos,'FaceColor','w','EdgeColor','k','LineWidth',2);     % LP rectangle
hold on
rectangle('Position',DPboxpos,'FaceColor','k','EdgeColor','k','LineWidth',2);    % DP rectangle
e1=errorbar(WakeWTBLRelmean_F,WakeWTBLRelstd_F,'Color','k', 'LineWidth',mylinewidth,'Marker','o','MarkerFaceColor','k','MarkerSize',mymarkersize,'CapSize',0);
e2=errorbar(WakeMutBLRelmean_F,WakeMutBLRelstd_F,'Color','r','LineWidth',mylinewidth,'Marker','o','MarkerFaceColor','r','MarkerSize',mymarkersize,'CapSize',0);
hold off 
ax=H.Children.Children;
set(ax, 'Layer','top')
%title('Wake, Baseline');
ylabel({'WAKE','Time in Wake (% TRT)'});
%xlabel('Time (hours)');
l = legend([e1 e2],LegendLabels);

%l=legend(ax,LegendLabels);
l.Box = 'off';
l.Location = 'None';
%l.Position = [0.7326 0.1649 0.1634 0.0798];
l.Location = 'SouthEast';
l.FontSize = 12;
%set(gca, 'XTick',[1 4 7 10 13 16 19 22])
set(ax,'XTick',[])
%set(gca, 'XTickLabels',{'1' '4' '7' '10' '13' '16' '19' '22'})
xlim([0 24])
ylim([-5 100])
ax.FontSize = myfontsize;
ax.LineWidth = mylinewidth;
% % -- Add stars if sig diffs (LP)
% if PostHoc_bool.BLFirst12hr.W
% 	hold on 
% 	bool_xlocsBL_LP_sig_values = ismember(e1.XData,xlocs_sig_effects.BL.WT_vs_Mut_LP_W);
% 	%bool_xlocsMutBL_LP_sig_values = ismember(e2.XData,xlocs_sig_effects.WTBL_vs_WTSD_LP_W);
% 	largest_y1_vals = e1.YData(bool_xlocsBL_LP_sig_values) + e1.YPositiveDelta(bool_xlocsBL_LP_sig_values);
% 	largest_y2_vals = e2.YData(bool_xlocsBL_LP_sig_values) + e2.YPositiveDelta(bool_xlocsBL_LP_sig_values);
% 	ylocs_sig_effects = max([largest_y1_vals; largest_y2_vals]) + 10*ones(1,length(largest_y1_vals));
% 	ylocs_sig_effects(ylocs_sig_effects>ax.YLim(2))=ax.YLim(2);
% 	text(xlocs_sig_effects.BL.WT_vs_Mut_LP_W, ylocs_sig_effects, '*', 'FontSize',astr_size,'HorizontalAlignment','center','VerticalAlignment','middle');
% 	hold off 
% end 
% % -- Add stars if sig diffs (DP)
% if PostHoc_bool.BLLast12hr.W
% 	hold on 
% 	bool_xlocsBL_DP_sig_values = ismember(e1.XData,xlocs_sig_effects.BL.WT_vs_Mut_DP_W);
% 	%bool_xlocsMutBL_LP_sig_values = ismember(e2.XData,xlocs_sig_effects.WTBL_vs_WTSD_LP_W);
% 	largest_y1_vals = e1.YData(bool_xlocsBL_DP_sig_values) + e1.YPositiveDelta(bool_xlocsBL_DP_sig_values);
% 	largest_y2_vals = e2.YData(bool_xlocsBL_DP_sig_values) + e2.YPositiveDelta(bool_xlocsBL_DP_sig_values);
% 	ylocs_sig_effects = max([largest_y1_vals; largest_y2_vals]) + 10*ones(1,length(largest_y1_vals));
% 	ylocs_sig_effects(ylocs_sig_effects>ax.YLim(2))=ax.YLim(2);
% 	text(xlocs_sig_effects.BL.WT_vs_Mut_DP_W, ylocs_sig_effects, '*', 'FontSize',astr_size,'HorizontalAlignment','center','VerticalAlignment','middle');
% 	hold off 
% end 



% -- BL NREM Females--
nexttile
rectangle('Position',LPboxpos,'FaceColor','w','EdgeColor','k','LineWidth',2);     % LP rectangle
hold on
rectangle('Position',DPboxpos,'FaceColor','k','EdgeColor','k','LineWidth',2);    % DP rectangle

e1=errorbar(NREMMutBLRelmean_F,NREMMutBLRelstd_F,'Color','r','LineWidth',mylinewidth,'Marker','o','MarkerFaceColor','r','MarkerSize',mymarkersize,'CapSize',0);
e2=errorbar(NREMWTBLRelmean_F,NREMWTBLRelstd_F,'Color','k','LineWidth',mylinewidth,'Marker','o','MarkerFaceColor','k','MarkerSize',mymarkersize,'CapSize',0);
hold off
ax=H.Children.Children(1);
set(ax, 'Layer','top')
%title('NREM Baseline');
% legend('Mutant','WT');
% legend 'boxoff';
%set(gca, 'XTick',[1 4 7 10 13 16 19 22])
set(ax,'XTick',[])
%set(gca, 'XTickLabels',{'1' '4' '7' '10' '13' '16' '19' '22'})
xlim([0 24])
ylim([-5 100])
ylabel({'NREM','Time in NREM (% TRT)'});
%xlabel('Time (hours)')

ax.FontSize = myfontsize;
ax.LineWidth = mylinewidth;
% -- Add stars if sig diffs (LP)
% if PostHoc_bool.BLFirst12hr.N
% 	hold on 
% 	bool_xlocsWTBL_LP_sig_values = ismember(e1.XData,xlocs_sig_effects.BL.WT_vs_Mut_LP_N);
% 	%bool_xlocsMutBL_LP_sig_values = ismember(e2.XData,xlocs_sig_effects.WTBL_vs_WTSD_LP_W);
% 	largest_y1_vals = e1.YData(bool_xlocsWTBL_LP_sig_values) + e1.YPositiveDelta(bool_xlocsWTBL_LP_sig_values);
% 	largest_y2_vals = e2.YData(bool_xlocsWTBL_LP_sig_values) + e2.YPositiveDelta(bool_xlocsWTBL_LP_sig_values);
% 	ylocs_sig_effects = max([largest_y1_vals; largest_y2_vals]) + 10*ones(1,length(largest_y1_vals));
% 	ylocs_sig_effects(ylocs_sig_effects>ax.YLim(2))=ax.YLim(2);
% 	text(xlocs_sig_effects.BL.WT_vs_Mut_LP_N, ylocs_sig_effects, '*', 'FontSize',astr_size,'HorizontalAlignment','center','VerticalAlignment','middle');
% 	hold off 
% end 
% % -- Add stars if sig diffs (DP)
% if PostHoc_bool.BLLast12hr.N
% 	hold on 
% 	bool_xlocsWTBL_DP_sig_values = ismember(e1.XData,xlocs_sig_effects.BL.WT_vs_Mut_DP_N);
% 	%bool_xlocsMutBL_LP_sig_values = ismember(e2.XData,xlocs_sig_effects.WTBL_vs_WTSD_LP_W);
% 	largest_y1_vals = e1.YData(bool_xlocsWTBL_DP_sig_values) + e1.YPositiveDelta(bool_xlocsWTBL_DP_sig_values);
% 	largest_y2_vals = e2.YData(bool_xlocsWTBL_DP_sig_values) + e2.YPositiveDelta(bool_xlocsWTBL_DP_sig_values);
% 	ylocs_sig_effects = max([largest_y1_vals; largest_y2_vals]) + 10*ones(1,length(largest_y1_vals));
% 	ylocs_sig_effects(ylocs_sig_effects>ax.YLim(2))=ax.YLim(2);
% 	text(xlocs_sig_effects.BL.WT_vs_Mut_DP_N, ylocs_sig_effects, '*', 'FontSize',astr_size,'HorizontalAlignment','center','VerticalAlignment','middle');
% 	hold off 
% end 
text(0,95, ...
		{'NOTE: These figures do not have *, # or bars added automatically.';'See RepMeas_ANOVA_tables_BothSexes_TIS and'; 'RepMeas_ANOVA_TIS_BothSexes_p_values and PostHocs_TIS_MandF.'}, ...
		'Interpreter','none','fontsize',12)

% -- BL REM Females--
nexttile
rectangle('Position',LPboxposREM,'FaceColor','w','EdgeColor','k','LineWidth',2);     % LP rectangle
hold on
rectangle('Position',DPboxposREM,'FaceColor','k','EdgeColor','k','LineWidth',2);    % DP rectangle

e1=errorbar(REMMutBLRelmean_F,REMMutBLRelstd_F,'Color','r','LineWidth',mylinewidth,'Marker','o','MarkerFaceColor','r','MarkerSize',mymarkersize,'CapSize',0);
e2=errorbar(REMWTBLRelmean_F,REMWTBLRelstd_F,'Color','k','LineWidth',mylinewidth,'Marker','o','MarkerFaceColor','k','MarkerSize',mymarkersize,'CapSize',0);
hold off
ax=H.Children.Children(1);
set(ax, 'Layer','top')
%title('REM Baseline');
% legend('Mutant','WT');
% legend 'boxoff';
set(ax, 'XTick',[1 4 7 10 13 16 19 22])
set(ax, 'XTickLabels',{'1' '4' '7' '10' '13' '16' '19' '22'})
xlim([0 24])
ylim ([-1 30])
xlabel('Time (hours)')
ylabel({'REM','Time in REM (% TRT)'});

ax.FontSize = myfontsize;
ax.LineWidth = mylinewidth;
% -- Add stars if sig diffs (LP)
% if PostHoc_bool.BLFirst12hr.R
% 	hold on 
% 	bool_xlocsWTBL_LP_sig_values = ismember(e1.XData,xlocs_sig_effects.BL.WT_vs_Mut_LP_R);
% 	%bool_xlocsMutBL_LP_sig_values = ismember(e2.XData,xlocs_sig_effects.WTBL_vs_WTSD_LP_W);
% 	largest_y1_vals = e1.YData(bool_xlocsWTBL_LP_sig_values) + e1.YPositiveDelta(bool_xlocsWTBL_LP_sig_values);
% 	largest_y2_vals = e2.YData(bool_xlocsWTBL_LP_sig_values) + e2.YPositiveDelta(bool_xlocsWTBL_LP_sig_values);
% 	ylocs_sig_effects = max([largest_y1_vals; largest_y2_vals]) + 10*ones(1,length(largest_y1_vals));
% 	ylocs_sig_effects(ylocs_sig_effects>ax.YLim(2))=ax.YLim(2);
% 	text(xlocs_sig_effects.BL.WT_vs_Mut_LP_R, ylocs_sig_effects, '*', 'FontSize',astr_size,'HorizontalAlignment','center','VerticalAlignment','middle');
% 	hold off 
% end 
% % -- Add stars if sig diffs (DP)
% if PostHoc_bool.BLLast12hr.R
% 	hold on 
% 	bool_xlocsWTBL_DP_sig_values = ismember(e1.XData,xlocs_sig_effects.BL.WT_vs_Mut_DP_R);
% 	%bool_xlocsMutBL_LP_sig_values = ismember(e2.XData,xlocs_sig_effects.WTBL_vs_WTSD_LP_W);
% 	largest_y1_vals = e1.YData(bool_xlocsWTBL_DP_sig_values) + e1.YPositiveDelta(bool_xlocsWTBL_DP_sig_values);
% 	largest_y2_vals = e2.YData(bool_xlocsWTBL_DP_sig_values) + e2.YPositiveDelta(bool_xlocsWTBL_DP_sig_values);
% 	ylocs_sig_effects = max([largest_y1_vals; largest_y2_vals]) + 10*ones(1,length(largest_y1_vals));
% 	ylocs_sig_effects(ylocs_sig_effects>ax.YLim(2))=ax.YLim(2);
% 	text(xlocs_sig_effects.BL.WT_vs_Mut_DP_R, ylocs_sig_effects, '*', 'FontSize',astr_size,'HorizontalAlignment','center','VerticalAlignment','middle');
% 	hold off 
% end 


% ------ SD Figures ----------------------------------------------------------------

% ------ 3-panel figure for all of SD Male---------
H=figure;
H.Name = ['Three-Panel SD Time in State Figure ' Sexes{1}];
H.Renderer = 'painters';
H.Position = [429 118 678 698];

t= tiledlayout(3,1,'TileSpacing','compact','Padding','compact','TileIndexing', 'columnmajor');
nexttile
%-- SD W as part of 3-panel --
rectangle('Position',LPboxpos,'FaceColor','w','EdgeColor','k','LineWidth',2);     % LP rectangle
hold on
rectangle('Position',DPboxpos,'FaceColor','k','EdgeColor','k','LineWidth',2);    % DP rectangle
rectangle('Position',[0 -5 SD_length_hrs 5],'FaceColor',[0.5 0.5 0.5],'EdgeColor','k','LineWidth',2);  % SD rectangle

% p2 = patch([0 0 SD_length_hrs SD_length_hrs],[-5 0 0 -5],[1 1 1]);  % SD hatching
% hatchfill2(p2,'single','HatchAngle',45,'HatchDensity',35,'HatchLineWidth',1.5);

% -- Wake Males --
%rectangle('Position',posB,'FaceColor',[0.5 0.5 0.5],'EdgeColor','none');
e1=errorbar(WakeWTSDRelmean_M,WakeWTSDRelstd_M,'Color','k','LineWidth',mylinewidth,'Marker','o','MarkerFaceColor','k','MarkerSize',mymarkersize,'CapSize',0);
e2=errorbar(WakeMutSDRelmean_M,WakeMutSDRelstd_M,'Color','r','LineWidth',mylinewidth,'Marker','o','MarkerFaceColor','r','MarkerSize',mymarkersize,'CapSize',0);
hold off 
ax=H.Children.Children;
set(ax, 'Layer','top')
%title('RS Wake')
ylabel({'WAKE','Time in WAKE (% TRT)'});
xlabel('Time (hours)')
xlabel('')
l = legend([e1 e2],LegendLabels);
l.Box = 'off';
l.Location = 'SouthEast';
%l.Position = [0.7326 0.1649 0.1634 0.0798];
l.FontSize = 12;
%set(gca, 'XTick',[1 4 7 10 13 16 19 22])
set(ax,'XTick',[])
%set(gca, 'XTickLabels',{'1' '4' '7' '10' '13' '16' '19' '22'})
set(ax, 'XTickLabels',{'' '' '' '' '' '' '' ''})

xlim([0 24])
ylim([-5 100])

ax.FontSize = myfontsize;
ax.LineWidth = mylinewidth;
% -- Add stars if sig diffs (LP)
% if PostHoc_bool.SDFirst12hr.W
% 	hold on 
% 	bool_xlocsSD_LP_sig_values = ismember(e1.XData,xlocs_sig_effects.SD.WT_vs_Mut_LP_W);
% 	%bool_xlocsMutBL_LP_sig_values = ismember(e2.XData,xlocs_sig_effects.WTBL_vs_WTSD_LP_W);
% 	largest_y1_vals = e1.YData(bool_xlocsSD_LP_sig_values) + e1.YPositiveDelta(bool_xlocsSD_LP_sig_values);
% 	largest_y2_vals = e2.YData(bool_xlocsSD_LP_sig_values) + e2.YPositiveDelta(bool_xlocsSD_LP_sig_values);
% 	ylocs_sig_effects = max([largest_y1_vals; largest_y2_vals]) + 10*ones(1,length(largest_y1_vals));
% 	ylocs_sig_effects(ylocs_sig_effects>ax.YLim(2))=ax.YLim(2);
% 	text(xlocs_sig_effects.SD.WT_vs_Mut_LP_W, ylocs_sig_effects, '*', 'FontSize',astr_size,'HorizontalAlignment','center','VerticalAlignment','middle');
% 	hold off 
% end 
% % -- Add stars if sig diffs (DP)
% if PostHoc_bool.SDLast12hr.W
% 	hold on 
% 	bool_xlocsSD_DP_sig_values = ismember(e1.XData,xlocs_sig_effects.SD.WT_vs_Mut_DP_W);
% 	%bool_xlocsMutBL_LP_sig_values = ismember(e2.XData,xlocs_sig_effects.WTBL_vs_WTSD_LP_W);
% 	largest_y1_vals = e1.YData(bool_xlocsSD_DP_sig_values) + e1.YPositiveDelta(bool_xlocsSD_DP_sig_values);
% 	largest_y2_vals = e2.YData(bool_xlocsSD_DP_sig_values) + e2.YPositiveDelta(bool_xlocsSD_DP_sig_values);
% 	ylocs_sig_effects = max([largest_y1_vals; largest_y2_vals]) + 10*ones(1,length(largest_y1_vals));
% 	ylocs_sig_effects(ylocs_sig_effects>ax.YLim(2))=ax.YLim(2);
% 	text(xlocs_sig_effects.SD.WT_vs_Mut_DP_W, ylocs_sig_effects, '*', 'FontSize',astr_size,'HorizontalAlignment','center','VerticalAlignment','middle');
% 	hold off 
% end 

% % -- Add stars if sig diffs (First 6 hours of DP)
% if isfield(PostHoc_bool,'SDFirst6hrDP') & PostHoc_bool.SDFirst6hrDP.W
% 	hold on 
% 	bool_xlocsSD_First6hrDP_sig_values = ismember(e1.XData,xlocs_sig_effects.SD.WT_vs_Mut_First6hrDP_W);
% 	largest_y1_vals = e1.YData(bool_xlocsSD_First6hrDP_sig_values) + e1.YPositiveDelta(bool_xlocsSD_First6hrDP_sig_values);
% 	largest_y2_vals = e2.YData(bool_xlocsSD_First6hrDP_sig_values) + e2.YPositiveDelta(bool_xlocsSD_First6hrDP_sig_values);
% 	ylocs_sig_effects = max([largest_y1_vals; largest_y2_vals]) + 10*ones(1,length(largest_y1_vals));
% 	ylocs_sig_effects(ylocs_sig_effects>ax.YLim(2))=ax.YLim(2);
% 	text(xlocs_sig_effects.SD.WT_vs_Mut_First6hrDP_W, ylocs_sig_effects, '*', 'FontSize',astr_size,'HorizontalAlignment','center','VerticalAlignment','middle','Color',[0.5 0.5 0.5]);
% 	hold off 
% end 
% % -- Add stars if sig diffs (Last 6 hours of DP)
% if isfield(PostHoc_bool,'SDLast6hrDP') & PostHoc_bool.SDLast6hrDP.W
% 	hold on 
% 	bool_xlocsSD_Last6hrDP_sig_values = ismember(e1.XData,xlocs_sig_effects.SD.WT_vs_Mut_Last6hrDP_W);
% 	largest_y1_vals = e1.YData(bool_xlocsSD_Last6hrDP_sig_values) + e1.YPositiveDelta(bool_xlocsSD_Last6hrDP_sig_values);
% 	largest_y2_vals = e2.YData(bool_xlocsSD_Last6hrDP_sig_values) + e2.YPositiveDelta(bool_xlocsSD_Last6hrDP_sig_values);
% 	ylocs_sig_effects = max([largest_y1_vals; largest_y2_vals]) + 10*ones(1,length(largest_y1_vals));
% 	ylocs_sig_effects(ylocs_sig_effects>ax.YLim(2))=ax.YLim(2);
% 	text(xlocs_sig_effects.SD.WT_vs_Mut_Last6hrDP_W, ylocs_sig_effects, '*', 'FontSize',astr_size,'HorizontalAlignment','center','VerticalAlignment','middle','Color',[0.5 0.5 0.5]);
% 	hold off 
% end 

%-- SD NREM as part of 3-panel Males--
nexttile 
rectangle('Position',LPboxpos,'FaceColor','w','EdgeColor','k','LineWidth',2);     % LP rectangle
hold on
rectangle('Position',DPboxpos,'FaceColor','k','EdgeColor','k','LineWidth',2);    % DP rectangle
rectangle('Position',[0 -5 SD_length_hrs 5],'FaceColor',[0.5 0.5 0.5],'EdgeColor','k','LineWidth',2);  % SD rectangle

e1=errorbar(NREMMutSDRelmean_M,NREMMutSDRelstd_M,'Color','r','LineWidth',mylinewidth,'Marker','o','MarkerFaceColor','r','MarkerSize',mymarkersize,'CapSize',0);
e2=errorbar(NREMWTSDRelmean_M,NREMWTSDRelstd_M,'Color','k','LineWidth',mylinewidth,'Marker','o','MarkerFaceColor','k','MarkerSize',mymarkersize,'CapSize',0);
hold off
ax=H.Children.Children(1);
set(ax, 'Layer','top')
%title('RS NREM');
% legend('Mutant','WT');
% legend 'boxoff';
%set(gca, 'XTick',[1 4 7 10 13 16 19 22])
set(ax,'XTick',[])
%set(gca, 'XTickLabels',{'1' '4' '7' '10' '13' '16' '19' '22'})

xlim([0 24])
ylim([-5 100])
%xlabel('Time (hours)')
xlabel('')
ylabel({'NREM','Time in NREM (% TRT)'});

ax.FontSize = myfontsize;
ax.LineWidth = mylinewidth;

% -- Add stars if sig diffs (LP)
% if PostHoc_bool.SDFirst12hr.N
% 	hold on 
% 	bool_xlocsSD_LP_sig_values = ismember(e1.XData,xlocs_sig_effects.SD.WT_vs_Mut_LP_N);
% 	%bool_xlocsMutBL_LP_sig_values = ismember(e2.XData,xlocs_sig_effects.WTBL_vs_WTSD_LP_W);
% 	largest_y1_vals = e1.YData(bool_xlocsSD_LP_sig_values) + e1.YPositiveDelta(bool_xlocsSD_LP_sig_values);
% 	largest_y2_vals = e2.YData(bool_xlocsSD_LP_sig_values) + e2.YPositiveDelta(bool_xlocsSD_LP_sig_values);
% 	ylocs_sig_effects = max([largest_y1_vals; largest_y2_vals]) + 10*ones(1,length(largest_y1_vals));
% 	ylocs_sig_effects(ylocs_sig_effects>ax.YLim(2))=ax.YLim(2);
% 	text(xlocs_sig_effects.SD.WT_vs_Mut_LP_N, ylocs_sig_effects, '*', 'FontSize',astr_size,'HorizontalAlignment','center','VerticalAlignment','middle');
% 	hold off 
% end 
% % -- Add stars if sig diffs (DP)
% if PostHoc_bool.SDLast12hr.N
% 	hold on 
% 	bool_xlocsSD_DP_sig_values = ismember(e1.XData,xlocs_sig_effects.SD.WT_vs_Mut_DP_N);
% 	%bool_xlocsMutBL_LP_sig_values = ismember(e2.XData,xlocs_sig_effects.WTBL_vs_WTSD_LP_W);
% 	largest_y1_vals = e1.YData(bool_xlocsSD_DP_sig_values) + e1.YPositiveDelta(bool_xlocsSD_DP_sig_values);
% 	largest_y2_vals = e2.YData(bool_xlocsSD_DP_sig_values) + e2.YPositiveDelta(bool_xlocsSD_DP_sig_values);
% 	ylocs_sig_effects = max([largest_y1_vals; largest_y2_vals]) + 10*ones(1,length(largest_y1_vals));
% 	ylocs_sig_effects(ylocs_sig_effects>ax.YLim(2))=ax.YLim(2);
% 	text(xlocs_sig_effects.SD.WT_vs_Mut_DP_N, ylocs_sig_effects, '*', 'FontSize',astr_size,'HorizontalAlignment','center','VerticalAlignment','middle');
% 	hold off 
% end 
% % -- Add stars if sig diffs (First 6 hours of DP)
% if isfield(PostHoc_bool,'SDFirst6hrDP') & PostHoc_bool.SDFirst6hrDP.N
% 	hold on 
% 	bool_xlocsSD_First6hrDP_sig_values = ismember(e1.XData,xlocs_sig_effects.SD.WT_vs_Mut_First6hrDP_N);
% 	largest_y1_vals = e1.YData(bool_xlocsSD_First6hrDP_sig_values) + e1.YPositiveDelta(bool_xlocsSD_First6hrDP_sig_values);
% 	largest_y2_vals = e2.YData(bool_xlocsSD_First6hrDP_sig_values) + e2.YPositiveDelta(bool_xlocsSD_First6hrDP_sig_values);
% 	ylocs_sig_effects = max([largest_y1_vals; largest_y2_vals]) + 10*ones(1,length(largest_y1_vals));
% 	ylocs_sig_effects(ylocs_sig_effects>ax.YLim(2))=ax.YLim(2);
% 	text(xlocs_sig_effects.SD.WT_vs_Mut_First6hrDP_N, ylocs_sig_effects, '*', 'FontSize',astr_size,'HorizontalAlignment','center','VerticalAlignment','middle','Color',[0.5 0.5 0.5]);
% 	hold off 
% end 
% % -- Add stars if sig diffs (Last 6 hours of DP)
% if isfield(PostHoc_bool,'SDLast6hrDP') & PostHoc_bool.SDLast6hrDP.N
% 	hold on 
% 	bool_xlocsSD_Last6hrDP_sig_values = ismember(e1.XData,xlocs_sig_effects.SD.WT_vs_Mut_Last6hrDP_N);
% 	largest_y1_vals = e1.YData(bool_xlocsSD_Last6hrDP_sig_values) + e1.YPositiveDelta(bool_xlocsSD_Last6hrDP_sig_values);
% 	largest_y2_vals = e2.YData(bool_xlocsSD_Last6hrDP_sig_values) + e2.YPositiveDelta(bool_xlocsSD_Last6hrDP_sig_values);
% 	ylocs_sig_effects = max([largest_y1_vals; largest_y2_vals]) + 10*ones(1,length(largest_y1_vals));
% 	ylocs_sig_effects(ylocs_sig_effects>ax.YLim(2))=ax.YLim(2);
% 	text(xlocs_sig_effects.SD.WT_vs_Mut_Last6hrDP_N, ylocs_sig_effects, '*', 'FontSize',astr_size,'HorizontalAlignment','center','VerticalAlignment','middle','Color',[0.5 0.5 0.5]);
% 	hold off 
% end 
text(0,95, ...
		{'NOTE: These figures do not have *, # or bars added automatically.';'See RepMeas_ANOVA_tables_BothSexes_TIS and'; 'RepMeas_ANOVA_TIS_BothSexes_p_values and PostHocs_TIS_MandF.'}, ...
		'Interpreter','none','fontsize',12)

% -- SD REM as part of 3-panel Males--
nexttile
rectangle('Position',LPboxposREM,'FaceColor','w','EdgeColor','k','LineWidth',2);     % LP rectangle
hold on
rectangle('Position',DPboxposREM,'FaceColor','k','EdgeColor','k','LineWidth',2);    % DP rectangle
rectangle('Position',[0 -1 SD_length_hrs 1],'FaceColor',[0.5 0.5 0.5],'EdgeColor','k','LineWidth',2);  % SD rectangle

e1=errorbar(REMMutSDRelmean_M,REMMutSDRelstd_M,'Color','r','LineWidth',mylinewidth,'Marker','o','MarkerFaceColor','r','MarkerSize',mymarkersize,'CapSize',0);
e2=errorbar(REMWTSDRelmean_M,REMWTSDRelstd_M,'Color','k','LineWidth',mylinewidth,'Marker','o','MarkerFaceColor','k','MarkerSize',mymarkersize,'CapSize',0);
hold off 
ax=H.Children.Children(1);
set(ax, 'Layer','top')
% legend('Mutant','WT')
% legend 'boxoff';
%title('RS REM');
set(ax, 'XTick',[1 4 7 10 13 16 19 22])
set(ax, 'XTickLabels',{'1' '4' '7' '10' '13' '16' '19' '22'})
xlim([0 24])
ylim([-1 30])
xlabel('Time (hours)')
ylabel({'REM','Time in REM (% TRT)'});

ax.FontSize = myfontsize;
ax.LineWidth = mylinewidth;

% -- Add stars if sig diffs (LP)
% if PostHoc_bool.SDFirst12hr.R
% 	hold on 
% 	bool_xlocsSD_LP_sig_values = ismember(e1.XData,xlocs_sig_effects.SD.WT_vs_Mut_LP_R);
% 	%bool_xlocsMutBL_LP_sig_values = ismember(e2.XData,xlocs_sig_effects.WTBL_vs_WTSD_LP_W);
% 	largest_y1_vals = e1.YData(bool_xlocsSD_LP_sig_values) + e1.YPositiveDelta(bool_xlocsSD_LP_sig_values);
% 	largest_y2_vals = e2.YData(bool_xlocsSD_LP_sig_values) + e2.YPositiveDelta(bool_xlocsSD_LP_sig_values);
% 	ylocs_sig_effects = max([largest_y1_vals; largest_y2_vals]) + 10*ones(1,length(largest_y1_vals));
% 	ylocs_sig_effects(ylocs_sig_effects>ax.YLim(2))=ax.YLim(2);
% 	text(xlocs_sig_effects.SD.WT_vs_Mut_LP_R, ylocs_sig_effects, '*', 'FontSize',astr_size,'HorizontalAlignment','center','VerticalAlignment','middle');
% 	hold off 
% end 
% % -- Add stars if sig diffs (DP)
% if PostHoc_bool.SDLast12hr.R
% 	hold on 
% 	bool_xlocsSD_DP_sig_values = ismember(e1.XData,xlocs_sig_effects.SD.WT_vs_Mut_DP_R);
% 	%bool_xlocsMutBL_LP_sig_values = ismember(e2.XData,xlocs_sig_effects.WTBL_vs_WTSD_LP_W);
% 	largest_y1_vals = e1.YData(bool_xlocsSD_DP_sig_values) + e1.YPositiveDelta(bool_xlocsSD_DP_sig_values);
% 	largest_y2_vals = e2.YData(bool_xlocsSD_DP_sig_values) + e2.YPositiveDelta(bool_xlocsSD_DP_sig_values);
% 	ylocs_sig_effects = max([largest_y1_vals; largest_y2_vals]) + 10*ones(1,length(largest_y1_vals));
% 	ylocs_sig_effects(ylocs_sig_effects>ax.YLim(2))=ax.YLim(2);
% 	text(xlocs_sig_effects.SD.WT_vs_Mut_DP_R, ylocs_sig_effects, '*', 'FontSize',astr_size,'HorizontalAlignment','center','VerticalAlignment','middle');
% 	hold off 
% end 
% % -- Add stars if sig diffs (First 6 hours of DP)
% if isfield(PostHoc_bool,'SDFirst6hrDP') & PostHoc_bool.SDFirst6hrDP.R
% 	hold on 
% 	bool_xlocsSD_First6hrDP_sig_values = ismember(e1.XData,xlocs_sig_effects.SD.WT_vs_Mut_First6hrDP_R);
% 	largest_y1_vals = e1.YData(bool_xlocsSD_First6hrDP_sig_values) + e1.YPositiveDelta(bool_xlocsSD_First6hrDP_sig_values);
% 	largest_y2_vals = e2.YData(bool_xlocsSD_First6hrDP_sig_values) + e2.YPositiveDelta(bool_xlocsSD_First6hrDP_sig_values);
% 	ylocs_sig_effects = max([largest_y1_vals; largest_y2_vals]) + 10*ones(1,length(largest_y1_vals));
% 	ylocs_sig_effects(ylocs_sig_effects>ax.YLim(2))=ax.YLim(2);
% 	text(xlocs_sig_effects.SD.WT_vs_Mut_First6hrDP_R, ylocs_sig_effects, '*', 'FontSize',astr_size,'HorizontalAlignment','center','VerticalAlignment','middle','Color',[0.5 0.5 0.5]);
% 	hold off 
% end 
% % -- Add stars if sig diffs (Last 6 hours of DP)
% if isfield(PostHoc_bool,'SDLast6hrDP') & PostHoc_bool.SDLast6hrDP.R
% 	hold on 
% 	bool_xlocsSD_Last6hrDP_sig_values = ismember(e1.XData,xlocs_sig_effects.SD.WT_vs_Mut_Last6hrDP_R);
% 	largest_y1_vals = e1.YData(bool_xlocsSD_Last6hrDP_sig_values) + e1.YPositiveDelta(bool_xlocsSD_Last6hrDP_sig_values);
% 	largest_y2_vals = e2.YData(bool_xlocsSD_Last6hrDP_sig_values) + e2.YPositiveDelta(bool_xlocsSD_Last6hrDP_sig_values);
% 	ylocs_sig_effects = max([largest_y1_vals; largest_y2_vals]) + 10*ones(1,length(largest_y1_vals));
% 	ylocs_sig_effects(ylocs_sig_effects>ax.YLim(2))=ax.YLim(2);
% 	text(xlocs_sig_effects.SD.WT_vs_Mut_Last6hrDP_R, ylocs_sig_effects, '*', 'FontSize',astr_size,'HorizontalAlignment','center','VerticalAlignment','middle','Color',[0.5 0.5 0.5]);
% 	hold off 
% end 



% ------ 3-panel figure for all of SD Female---------
H=figure;
H.Name = ['Three-Panel SD Time in State Figure ' Sexes{2}];
H.Renderer = 'painters';
H.Position = [429 118 678 698];

t= tiledlayout(3,1,'TileSpacing','compact','Padding','compact','TileIndexing', 'columnmajor');
nexttile
%-- SD W as part of 3-panel --
rectangle('Position',LPboxpos,'FaceColor','w','EdgeColor','k','LineWidth',2);     % LP rectangle
hold on
rectangle('Position',DPboxpos,'FaceColor','k','EdgeColor','k','LineWidth',2);    % DP rectangle
rectangle('Position',[0 -5 SD_length_hrs 5],'FaceColor',[0.5 0.5 0.5],'EdgeColor','k','LineWidth',2);  % SD rectangle

% p2 = patch([0 0 SD_length_hrs SD_length_hrs],[-5 0 0 -5],[1 1 1]);  % SD hatching
% hatchfill2(p2,'single','HatchAngle',45,'HatchDensity',35,'HatchLineWidth',1.5);

% -- Wake Females --
%rectangle('Position',posB,'FaceColor',[0.5 0.5 0.5],'EdgeColor','none');
e1=errorbar(WakeWTSDRelmean_F,WakeWTSDRelstd_F,'Color','k','LineWidth',mylinewidth,'Marker','o','MarkerFaceColor','k','MarkerSize',mymarkersize,'CapSize',0);
e2=errorbar(WakeMutSDRelmean_F,WakeMutSDRelstd_F,'Color','r','LineWidth',mylinewidth,'Marker','o','MarkerFaceColor','r','MarkerSize',mymarkersize,'CapSize',0);
hold off 
ax=H.Children.Children;
set(ax, 'Layer','top')
%title('RS Wake')
ylabel({'WAKE','Time in WAKE (% TRT)'});
xlabel('Time (hours)')
xlabel('')
l = legend([e1 e2],LegendLabels);
l.Box = 'off';
l.Location = 'SouthEast';
%l.Position = [0.7326 0.1649 0.1634 0.0798];
l.FontSize = 12;
%set(gca, 'XTick',[1 4 7 10 13 16 19 22])
set(ax,'XTick',[])
%set(gca, 'XTickLabels',{'1' '4' '7' '10' '13' '16' '19' '22'})
set(ax, 'XTickLabels',{'' '' '' '' '' '' '' ''})

xlim([0 24])
ylim([-5 100])
ax.FontSize = myfontsize;
ax.LineWidth = mylinewidth;
% -- Add stars if sig diffs (LP)
% if PostHoc_bool.SDFirst12hr.W
% 	hold on 
% 	bool_xlocsSD_LP_sig_values = ismember(e1.XData,xlocs_sig_effects.SD.WT_vs_Mut_LP_W);
% 	%bool_xlocsMutBL_LP_sig_values = ismember(e2.XData,xlocs_sig_effects.WTBL_vs_WTSD_LP_W);
% 	largest_y1_vals = e1.YData(bool_xlocsSD_LP_sig_values) + e1.YPositiveDelta(bool_xlocsSD_LP_sig_values);
% 	largest_y2_vals = e2.YData(bool_xlocsSD_LP_sig_values) + e2.YPositiveDelta(bool_xlocsSD_LP_sig_values);
% 	ylocs_sig_effects = max([largest_y1_vals; largest_y2_vals]) + 10*ones(1,length(largest_y1_vals));
% 	ylocs_sig_effects(ylocs_sig_effects>ax.YLim(2))=ax.YLim(2);
% 	text(xlocs_sig_effects.SD.WT_vs_Mut_LP_W, ylocs_sig_effects, '*', 'FontSize',astr_size,'HorizontalAlignment','center','VerticalAlignment','middle');
% 	hold off 
% end 
% % -- Add stars if sig diffs (DP)
% if PostHoc_bool.SDLast12hr.W
% 	hold on 
% 	bool_xlocsSD_DP_sig_values = ismember(e1.XData,xlocs_sig_effects.SD.WT_vs_Mut_DP_W);
% 	%bool_xlocsMutBL_LP_sig_values = ismember(e2.XData,xlocs_sig_effects.WTBL_vs_WTSD_LP_W);
% 	largest_y1_vals = e1.YData(bool_xlocsSD_DP_sig_values) + e1.YPositiveDelta(bool_xlocsSD_DP_sig_values);
% 	largest_y2_vals = e2.YData(bool_xlocsSD_DP_sig_values) + e2.YPositiveDelta(bool_xlocsSD_DP_sig_values);
% 	ylocs_sig_effects = max([largest_y1_vals; largest_y2_vals]) + 10*ones(1,length(largest_y1_vals));
% 	ylocs_sig_effects(ylocs_sig_effects>ax.YLim(2))=ax.YLim(2);
% 	text(xlocs_sig_effects.SD.WT_vs_Mut_DP_W, ylocs_sig_effects, '*', 'FontSize',astr_size,'HorizontalAlignment','center','VerticalAlignment','middle');
% 	hold off 
% end 

% % -- Add stars if sig diffs (First 6 hours of DP)
% if isfield(PostHoc_bool,'SDFirst6hrDP') & PostHoc_bool.SDFirst6hrDP.W
% 	hold on 
% 	bool_xlocsSD_First6hrDP_sig_values = ismember(e1.XData,xlocs_sig_effects.SD.WT_vs_Mut_First6hrDP_W);
% 	largest_y1_vals = e1.YData(bool_xlocsSD_First6hrDP_sig_values) + e1.YPositiveDelta(bool_xlocsSD_First6hrDP_sig_values);
% 	largest_y2_vals = e2.YData(bool_xlocsSD_First6hrDP_sig_values) + e2.YPositiveDelta(bool_xlocsSD_First6hrDP_sig_values);
% 	ylocs_sig_effects = max([largest_y1_vals; largest_y2_vals]) + 10*ones(1,length(largest_y1_vals));
% 	ylocs_sig_effects(ylocs_sig_effects>ax.YLim(2))=ax.YLim(2);
% 	text(xlocs_sig_effects.SD.WT_vs_Mut_First6hrDP_W, ylocs_sig_effects, '*', 'FontSize',astr_size,'HorizontalAlignment','center','VerticalAlignment','middle','Color',[0.5 0.5 0.5]);
% 	hold off 
% end 
% % -- Add stars if sig diffs (Last 6 hours of DP)
% if isfield(PostHoc_bool,'SDLast6hrDP') & PostHoc_bool.SDLast6hrDP.W
% 	hold on 
% 	bool_xlocsSD_Last6hrDP_sig_values = ismember(e1.XData,xlocs_sig_effects.SD.WT_vs_Mut_Last6hrDP_W);
% 	largest_y1_vals = e1.YData(bool_xlocsSD_Last6hrDP_sig_values) + e1.YPositiveDelta(bool_xlocsSD_Last6hrDP_sig_values);
% 	largest_y2_vals = e2.YData(bool_xlocsSD_Last6hrDP_sig_values) + e2.YPositiveDelta(bool_xlocsSD_Last6hrDP_sig_values);
% 	ylocs_sig_effects = max([largest_y1_vals; largest_y2_vals]) + 10*ones(1,length(largest_y1_vals));
% 	ylocs_sig_effects(ylocs_sig_effects>ax.YLim(2))=ax.YLim(2);
% 	text(xlocs_sig_effects.SD.WT_vs_Mut_Last6hrDP_W, ylocs_sig_effects, '*', 'FontSize',astr_size,'HorizontalAlignment','center','VerticalAlignment','middle','Color',[0.5 0.5 0.5]);
% 	hold off 
% end 

%-- SD NREM as part of 3-panel Females--
nexttile 
rectangle('Position',LPboxpos,'FaceColor','w','EdgeColor','k','LineWidth',2);     % LP rectangle
hold on
rectangle('Position',DPboxpos,'FaceColor','k','EdgeColor','k','LineWidth',2);    % DP rectangle
rectangle('Position',[0 -5 SD_length_hrs 5],'FaceColor',[0.5 0.5 0.5],'EdgeColor','k','LineWidth',2);  % SD rectangle

e1=errorbar(NREMMutSDRelmean_F,NREMMutSDRelstd_F,'Color','r','LineWidth',mylinewidth,'Marker','o','MarkerFaceColor','r','MarkerSize',mymarkersize,'CapSize',0);
e2=errorbar(NREMWTSDRelmean_F,NREMWTSDRelstd_F,'Color','k','LineWidth',mylinewidth,'Marker','o','MarkerFaceColor','k','MarkerSize',mymarkersize,'CapSize',0);
hold off
ax=H.Children.Children(1);
set(ax, 'Layer','top')
%title('RS NREM');
% legend('Mutant','WT');
% legend 'boxoff';
%set(gca, 'XTick',[1 4 7 10 13 16 19 22])
set(ax,'XTick',[])
%set(gca, 'XTickLabels',{'1' '4' '7' '10' '13' '16' '19' '22'})

xlim([0 24])
ylim([-5 100])
%xlabel('Time (hours)')
xlabel('')
ylabel({'NREM','Time in NREM (% TRT)'});
ax.FontSize = myfontsize;
ax.LineWidth = mylinewidth;

% -- Add stars if sig diffs (LP)
% if PostHoc_bool.SDFirst12hr.N
% 	hold on 
% 	bool_xlocsSD_LP_sig_values = ismember(e1.XData,xlocs_sig_effects.SD.WT_vs_Mut_LP_N);
% 	%bool_xlocsMutBL_LP_sig_values = ismember(e2.XData,xlocs_sig_effects.WTBL_vs_WTSD_LP_W);
% 	largest_y1_vals = e1.YData(bool_xlocsSD_LP_sig_values) + e1.YPositiveDelta(bool_xlocsSD_LP_sig_values);
% 	largest_y2_vals = e2.YData(bool_xlocsSD_LP_sig_values) + e2.YPositiveDelta(bool_xlocsSD_LP_sig_values);
% 	ylocs_sig_effects = max([largest_y1_vals; largest_y2_vals]) + 10*ones(1,length(largest_y1_vals));
% 	ylocs_sig_effects(ylocs_sig_effects>ax.YLim(2))=ax.YLim(2);
% 	text(xlocs_sig_effects.SD.WT_vs_Mut_LP_N, ylocs_sig_effects, '*', 'FontSize',astr_size,'HorizontalAlignment','center','VerticalAlignment','middle');
% 	hold off 
% end 
% % -- Add stars if sig diffs (DP)
% if PostHoc_bool.SDLast12hr.N
% 	hold on 
% 	bool_xlocsSD_DP_sig_values = ismember(e1.XData,xlocs_sig_effects.SD.WT_vs_Mut_DP_N);
% 	%bool_xlocsMutBL_LP_sig_values = ismember(e2.XData,xlocs_sig_effects.WTBL_vs_WTSD_LP_W);
% 	largest_y1_vals = e1.YData(bool_xlocsSD_DP_sig_values) + e1.YPositiveDelta(bool_xlocsSD_DP_sig_values);
% 	largest_y2_vals = e2.YData(bool_xlocsSD_DP_sig_values) + e2.YPositiveDelta(bool_xlocsSD_DP_sig_values);
% 	ylocs_sig_effects = max([largest_y1_vals; largest_y2_vals]) + 10*ones(1,length(largest_y1_vals));
% 	ylocs_sig_effects(ylocs_sig_effects>ax.YLim(2))=ax.YLim(2);
% 	text(xlocs_sig_effects.SD.WT_vs_Mut_DP_N, ylocs_sig_effects, '*', 'FontSize',astr_size,'HorizontalAlignment','center','VerticalAlignment','middle');
% 	hold off 
% end 
% % -- Add stars if sig diffs (First 6 hours of DP)
% if isfield(PostHoc_bool,'SDFirst6hrDP') & PostHoc_bool.SDFirst6hrDP.N
% 	hold on 
% 	bool_xlocsSD_First6hrDP_sig_values = ismember(e1.XData,xlocs_sig_effects.SD.WT_vs_Mut_First6hrDP_N);
% 	largest_y1_vals = e1.YData(bool_xlocsSD_First6hrDP_sig_values) + e1.YPositiveDelta(bool_xlocsSD_First6hrDP_sig_values);
% 	largest_y2_vals = e2.YData(bool_xlocsSD_First6hrDP_sig_values) + e2.YPositiveDelta(bool_xlocsSD_First6hrDP_sig_values);
% 	ylocs_sig_effects = max([largest_y1_vals; largest_y2_vals]) + 10*ones(1,length(largest_y1_vals));
% 	ylocs_sig_effects(ylocs_sig_effects>ax.YLim(2))=ax.YLim(2);
% 	text(xlocs_sig_effects.SD.WT_vs_Mut_First6hrDP_N, ylocs_sig_effects, '*', 'FontSize',astr_size,'HorizontalAlignment','center','VerticalAlignment','middle','Color',[0.5 0.5 0.5]);
% 	hold off 
% end 
% % -- Add stars if sig diffs (Last 6 hours of DP)
% if isfield(PostHoc_bool,'SDLast6hrDP') & PostHoc_bool.SDLast6hrDP.N
% 	hold on 
% 	bool_xlocsSD_Last6hrDP_sig_values = ismember(e1.XData,xlocs_sig_effects.SD.WT_vs_Mut_Last6hrDP_N);
% 	largest_y1_vals = e1.YData(bool_xlocsSD_Last6hrDP_sig_values) + e1.YPositiveDelta(bool_xlocsSD_Last6hrDP_sig_values);
% 	largest_y2_vals = e2.YData(bool_xlocsSD_Last6hrDP_sig_values) + e2.YPositiveDelta(bool_xlocsSD_Last6hrDP_sig_values);
% 	ylocs_sig_effects = max([largest_y1_vals; largest_y2_vals]) + 10*ones(1,length(largest_y1_vals));
% 	ylocs_sig_effects(ylocs_sig_effects>ax.YLim(2))=ax.YLim(2);
% 	text(xlocs_sig_effects.SD.WT_vs_Mut_Last6hrDP_N, ylocs_sig_effects, '*', 'FontSize',astr_size,'HorizontalAlignment','center','VerticalAlignment','middle','Color',[0.5 0.5 0.5]);
% 	hold off 
% end 
text(0,95, ...
		{'NOTE: These figures do not have *, # or bars added automatically.';'See RepMeas_ANOVA_tables_BothSexes_TIS and'; 'RepMeas_ANOVA_TIS_BothSexes_p_values and PostHocs_TIS_MandF.'}, ...
		'Interpreter','none','fontsize',12)

% -- SD REM as part of 3-panel Females--
nexttile
rectangle('Position',LPboxposREM,'FaceColor','w','EdgeColor','k','LineWidth',2);     % LP rectangle
hold on
rectangle('Position',DPboxposREM,'FaceColor','k','EdgeColor','k','LineWidth',2);    % DP rectangle
rectangle('Position',[0 -1 SD_length_hrs 1],'FaceColor',[0.5 0.5 0.5],'EdgeColor','k','LineWidth',2);  % SD rectangle

e1=errorbar(REMMutSDRelmean_F,REMMutSDRelstd_F,'Color','r','LineWidth',mylinewidth,'Marker','o','MarkerFaceColor','r','MarkerSize',mymarkersize,'CapSize',0);
e2=errorbar(REMWTSDRelmean_F,REMWTSDRelstd_F,'Color','k','LineWidth',mylinewidth,'Marker','o','MarkerFaceColor','k','MarkerSize',mymarkersize,'CapSize',0);
hold off 
ax=H.Children.Children(1);
set(ax, 'Layer','top')
% legend('Mutant','WT')
% legend 'boxoff';
%title('RS REM');
set(ax, 'XTick',[1 4 7 10 13 16 19 22])
set(ax, 'XTickLabels',{'1' '4' '7' '10' '13' '16' '19' '22'})
xlim([0 24])
ylim([-1 30])
xlabel('Time (hours)')
ylabel({'REM','Time in REM (% TRT)'});

ax.FontSize = myfontsize;
ax.LineWidth = mylinewidth;

% -- Add stars if sig diffs (LP)
% if PostHoc_bool.SDFirst12hr.R
% 	hold on 
% 	bool_xlocsSD_LP_sig_values = ismember(e1.XData,xlocs_sig_effects.SD.WT_vs_Mut_LP_R);
% 	%bool_xlocsMutBL_LP_sig_values = ismember(e2.XData,xlocs_sig_effects.WTBL_vs_WTSD_LP_W);
% 	largest_y1_vals = e1.YData(bool_xlocsSD_LP_sig_values) + e1.YPositiveDelta(bool_xlocsSD_LP_sig_values);
% 	largest_y2_vals = e2.YData(bool_xlocsSD_LP_sig_values) + e2.YPositiveDelta(bool_xlocsSD_LP_sig_values);
% 	ylocs_sig_effects = max([largest_y1_vals; largest_y2_vals]) + 10*ones(1,length(largest_y1_vals));
% 	ylocs_sig_effects(ylocs_sig_effects>ax.YLim(2))=ax.YLim(2);
% 	text(xlocs_sig_effects.SD.WT_vs_Mut_LP_R, ylocs_sig_effects, '*', 'FontSize',astr_size,'HorizontalAlignment','center','VerticalAlignment','middle');
% 	hold off 
% end 
% % -- Add stars if sig diffs (DP)
% if PostHoc_bool.SDLast12hr.R
% 	hold on 
% 	bool_xlocsSD_DP_sig_values = ismember(e1.XData,xlocs_sig_effects.SD.WT_vs_Mut_DP_R);
% 	%bool_xlocsMutBL_LP_sig_values = ismember(e2.XData,xlocs_sig_effects.WTBL_vs_WTSD_LP_W);
% 	largest_y1_vals = e1.YData(bool_xlocsSD_DP_sig_values) + e1.YPositiveDelta(bool_xlocsSD_DP_sig_values);
% 	largest_y2_vals = e2.YData(bool_xlocsSD_DP_sig_values) + e2.YPositiveDelta(bool_xlocsSD_DP_sig_values);
% 	ylocs_sig_effects = max([largest_y1_vals; largest_y2_vals]) + 10*ones(1,length(largest_y1_vals));
% 	ylocs_sig_effects(ylocs_sig_effects>ax.YLim(2))=ax.YLim(2);
% 	text(xlocs_sig_effects.SD.WT_vs_Mut_DP_R, ylocs_sig_effects, '*', 'FontSize',astr_size,'HorizontalAlignment','center','VerticalAlignment','middle');
% 	hold off 
% end 
% % -- Add stars if sig diffs (First 6 hours of DP)
% if isfield(PostHoc_bool,'SDFirst6hrDP') & PostHoc_bool.SDFirst6hrDP.R
% 	hold on 
% 	bool_xlocsSD_First6hrDP_sig_values = ismember(e1.XData,xlocs_sig_effects.SD.WT_vs_Mut_First6hrDP_R);
% 	largest_y1_vals = e1.YData(bool_xlocsSD_First6hrDP_sig_values) + e1.YPositiveDelta(bool_xlocsSD_First6hrDP_sig_values);
% 	largest_y2_vals = e2.YData(bool_xlocsSD_First6hrDP_sig_values) + e2.YPositiveDelta(bool_xlocsSD_First6hrDP_sig_values);
% 	ylocs_sig_effects = max([largest_y1_vals; largest_y2_vals]) + 10*ones(1,length(largest_y1_vals));
% 	ylocs_sig_effects(ylocs_sig_effects>ax.YLim(2))=ax.YLim(2);
% 	text(xlocs_sig_effects.SD.WT_vs_Mut_First6hrDP_R, ylocs_sig_effects, '*', 'FontSize',astr_size,'HorizontalAlignment','center','VerticalAlignment','middle','Color',[0.5 0.5 0.5]);
% 	hold off 
% end 
% % -- Add stars if sig diffs (Last 6 hours of DP)
% if isfield(PostHoc_bool,'SDLast6hrDP') & PostHoc_bool.SDLast6hrDP.R
% 	hold on 
% 	bool_xlocsSD_Last6hrDP_sig_values = ismember(e1.XData,xlocs_sig_effects.SD.WT_vs_Mut_Last6hrDP_R);
% 	largest_y1_vals = e1.YData(bool_xlocsSD_Last6hrDP_sig_values) + e1.YPositiveDelta(bool_xlocsSD_Last6hrDP_sig_values);
% 	largest_y2_vals = e2.YData(bool_xlocsSD_Last6hrDP_sig_values) + e2.YPositiveDelta(bool_xlocsSD_Last6hrDP_sig_values);
% 	ylocs_sig_effects = max([largest_y1_vals; largest_y2_vals]) + 10*ones(1,length(largest_y1_vals));
% 	ylocs_sig_effects(ylocs_sig_effects>ax.YLim(2))=ax.YLim(2);
% 	text(xlocs_sig_effects.SD.WT_vs_Mut_Last6hrDP_R, ylocs_sig_effects, '*', 'FontSize',astr_size,'HorizontalAlignment','center','VerticalAlignment','middle','Color',[0.5 0.5 0.5]);
% 	hold off 
% end 
% ---------------------------------------------------------------------------------------
% ---------------------------------------------------------------------------------------


% ---- Set up the Post hoc tables to return ----
% -- Genotype --
if ~isempty(mc_table.Genotype.BL.First12hr.Wake)  PostHoc_tables.Genotype.BL.First12hr.Wake  = mc_table.Genotype.BL.First12hr.Wake;  else PostHoc_tables.Genotype.BL.First12hr.Wake  = []; end  
if ~isempty(mc_table.Genotype.BL.First12hr.NREM)  PostHoc_tables.Genotype.BL.First12hr.NREM  = mc_table.Genotype.BL.First12hr.NREM;  else PostHoc_tables.Genotype.BL.First12hr.NREM  = []; end  
if ~isempty(mc_table.Genotype.BL.First12hr.REM)   PostHoc_tables.Genotype.BL.First12hr.REM   = mc_table.Genotype.BL.First12hr.REM;   else PostHoc_tables.Genotype.BL.First12hr.REM   = []; end
if ~isempty(mc_table.Genotype.BL.Last12hr.Wake)   PostHoc_tables.Genotype.BL.Last12hr.Wake   = mc_table.Genotype.BL.Last12hr.Wake;   else PostHoc_tables.Genotype.BL.Last12hr.Wake   = []; end
if ~isempty(mc_table.Genotype.BL.Last12hr.NREM)   PostHoc_tables.Genotype.BL.Last12hr.NREM   = mc_table.Genotype.BL.Last12hr.NREM;   else PostHoc_tables.Genotype.BL.Last12hr.NREM   = []; end
if ~isempty(mc_table.Genotype.BL.Last12hr.REM)    PostHoc_tables.Genotype.BL.Last12hr.REM    = mc_table.Genotype.BL.Last12hr.REM;    else PostHoc_tables.Genotype.BL.Last12hr.REM    = []; end 
if ~isempty(mc_table.Genotype.SD.First12hr.Wake)  PostHoc_tables.Genotype.SD.First12hr.Wake  = mc_table.Genotype.SD.First12hr.Wake;  else PostHoc_tables.Genotype.SD.First12hr.Wake  = []; end
if ~isempty(mc_table.Genotype.SD.First12hr.NREM)  PostHoc_tables.Genotype.SD.First12hr.NREM  = mc_table.Genotype.SD.First12hr.NREM;  else PostHoc_tables.Genotype.SD.First12hr.NREM  = []; end
if ~isempty(mc_table.Genotype.SD.First12hr.REM)   PostHoc_tables.Genotype.SD.First12hr.REM   = mc_table.Genotype.SD.First12hr.REM;   else PostHoc_tables.Genotype.SD.First12hr.REM   = []; end
if ~isempty(mc_table.Genotype.SD.Last12hr.Wake)   PostHoc_tables.Genotype.SD.Last12hr.Wake   = mc_table.Genotype.SD.Last12hr.Wake;   else PostHoc_tables.Genotype.SD.Last12hr.Wake   = []; end
if ~isempty(mc_table.Genotype.SD.Last12hr.NREM)   PostHoc_tables.Genotype.SD.Last12hr.NREM   = mc_table.Genotype.SD.Last12hr.NREM;   else PostHoc_tables.Genotype.SD.Last12hr.NREM   = []; end
if ~isempty(mc_table.Genotype.SD.Last12hr.REM)    PostHoc_tables.Genotype.SD.Last12hr.REM    = mc_table.Genotype.SD.Last12hr.REM;    else PostHoc_tables.Genotype.SD.Last12hr.REM    = []; end 

% -- Genotype --
if ~isempty(mc_table.(SexVarName).BL.First12hr.Wake)  PostHoc_tables.(SexVarName).BL.First12hr.Wake  = mc_table.(SexVarName).BL.First12hr.Wake;  else PostHoc_tables.(SexVarName).BL.First12hr.Wake  = []; end  
if ~isempty(mc_table.(SexVarName).BL.First12hr.NREM)  PostHoc_tables.(SexVarName).BL.First12hr.NREM  = mc_table.(SexVarName).BL.First12hr.NREM;  else PostHoc_tables.(SexVarName).BL.First12hr.NREM  = []; end  
if ~isempty(mc_table.(SexVarName).BL.First12hr.REM)   PostHoc_tables.(SexVarName).BL.First12hr.REM   = mc_table.(SexVarName).BL.First12hr.REM;   else PostHoc_tables.(SexVarName).BL.First12hr.REM   = []; end
if ~isempty(mc_table.(SexVarName).BL.Last12hr.Wake)   PostHoc_tables.(SexVarName).BL.Last12hr.Wake   = mc_table.(SexVarName).BL.Last12hr.Wake;   else PostHoc_tables.(SexVarName).BL.Last12hr.Wake   = []; end
if ~isempty(mc_table.(SexVarName).BL.Last12hr.NREM)   PostHoc_tables.(SexVarName).BL.Last12hr.NREM   = mc_table.(SexVarName).BL.Last12hr.NREM;   else PostHoc_tables.(SexVarName).BL.Last12hr.NREM   = []; end
if ~isempty(mc_table.(SexVarName).BL.Last12hr.REM)    PostHoc_tables.(SexVarName).BL.Last12hr.REM    = mc_table.(SexVarName).BL.Last12hr.REM;    else PostHoc_tables.(SexVarName).BL.Last12hr.REM    = []; end 
if ~isempty(mc_table.(SexVarName).SD.First12hr.Wake)  PostHoc_tables.(SexVarName).SD.First12hr.Wake  = mc_table.(SexVarName).SD.First12hr.Wake;  else PostHoc_tables.(SexVarName).SD.First12hr.Wake  = []; end
if ~isempty(mc_table.(SexVarName).SD.First12hr.NREM)  PostHoc_tables.(SexVarName).SD.First12hr.NREM  = mc_table.(SexVarName).SD.First12hr.NREM;  else PostHoc_tables.(SexVarName).SD.First12hr.NREM  = []; end
if ~isempty(mc_table.(SexVarName).SD.First12hr.REM)   PostHoc_tables.(SexVarName).SD.First12hr.REM   = mc_table.(SexVarName).SD.First12hr.REM;   else PostHoc_tables.(SexVarName).SD.First12hr.REM   = []; end
if ~isempty(mc_table.(SexVarName).SD.Last12hr.Wake)   PostHoc_tables.(SexVarName).SD.Last12hr.Wake   = mc_table.(SexVarName).SD.Last12hr.Wake;   else PostHoc_tables.(SexVarName).SD.Last12hr.Wake   = []; end
if ~isempty(mc_table.(SexVarName).SD.Last12hr.NREM)   PostHoc_tables.(SexVarName).SD.Last12hr.NREM   = mc_table.(SexVarName).SD.Last12hr.NREM;   else PostHoc_tables.(SexVarName).SD.Last12hr.NREM   = []; end
if ~isempty(mc_table.(SexVarName).SD.Last12hr.REM)    PostHoc_tables.(SexVarName).SD.Last12hr.REM    = mc_table.(SexVarName).SD.Last12hr.REM;    else PostHoc_tables.(SexVarName).SD.Last12hr.REM    = []; end 


% % If you are doing 6-hour segments of the DP
% if isfield(mc_table.SD,'First6hrDP')
% 	if ~isempty(mc_table.SD.First6hrDP.Wake) PostHoc_tables.SD.First6hrDP.Wake = mc_table.SD.First6hrDP.Wake; else PostHoc_tables.SD.First6hrDP.Wake = []; end
% 	if ~isempty(mc_table.SD.First6hrDP.NREM) PostHoc_tables.SD.First6hrDP.NREM = mc_table.SD.First6hrDP.NREM; else PostHoc_tables.SD.First6hrDP.NREM = []; end
% 	if ~isempty(mc_table.SD.First6hrDP.REM)  PostHoc_tables.SD.First6hrDP.REM  = mc_table.SD.First6hrDP.REM;  else PostHoc_tables.SD.First6hrDP.REM  = []; end 
% end
% if isfield(mc_table.SD,'Last6hrDP')
% 	if ~isempty(mc_table.SD.Last6hrDP.Wake)  PostHoc_tables.SD.Last6hrDP.Wake  = mc_table.SD.Last6hrDP.Wake;  else PostHoc_tables.SD.Last6hrDP.Wake  = []; end
% 	if ~isempty(mc_table.SD.Last6hrDP.NREM)  PostHoc_tables.SD.Last6hrDP.NREM  = mc_table.SD.Last6hrDP.NREM;  else PostHoc_tables.SD.Last6hrDP.NREM  = []; end
% 	if ~isempty(mc_table.SD.Last6hrDP.REM)   PostHoc_tables.SD.Last6hrDP.REM   = mc_table.SD.Last6hrDP.REM;   else PostHoc_tables.SD.Last6hrDP.REM   = []; end 
% end 




























