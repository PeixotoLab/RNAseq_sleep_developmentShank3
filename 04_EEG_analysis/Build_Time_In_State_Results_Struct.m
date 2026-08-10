function TimeInStateResults = Build_Time_In_State_Results_Struct(RepMeas_ANOVA_tables_Males_hourly_percentageFirst12hrBL,RepMeas_ANOVA_tables_Males_hourly_percentageLast12hrBL, ...
                                                        RepMeas_ANOVA_tables_Males_hourly_percentageFirst12hrSD,RepMeas_ANOVA_tables_Males_hourly_percentageLast12hrSD, ...
                                                        RepMeas_ANOVA_tables_Females_hourly_percentageFirst12hrBL,RepMeas_ANOVA_tables_Females_hourly_percentageLast12hrBL, ...
                                                        RepMeas_ANOVA_tables_Females_hourly_percentageFirst12hrSD,RepMeas_ANOVA_tables_Females_hourly_percentageLast12hrSD);

% This helper function simply puts the Repeated Measures ANOVA results for the time in state analysis into a more convienient output

% ------------------------------------------  Genotype -------------------------------------------------

% ------------------------ P-values ----------------------------

% ------------ Males ---------------------
% --- BL LP ---
TimeInStateResults.Male.RepMeasANOVAs.p_vals.Genotype.BL.LP.Wake =  RepMeas_ANOVA_tables_Males_hourly_percentageFirst12hrBL.BL_custom.Wake{'Genotype','pValue'};
TimeInStateResults.Male.RepMeasANOVAs.p_vals.Genotype.BL.LP.NREM =  RepMeas_ANOVA_tables_Males_hourly_percentageFirst12hrBL.BL_custom.NREM{'Genotype','pValue'};
TimeInStateResults.Male.RepMeasANOVAs.p_vals.Genotype.BL.LP.REM  =  RepMeas_ANOVA_tables_Males_hourly_percentageFirst12hrBL.BL_custom.REM{'Genotype','pValue'};

% --- BL DP ---
TimeInStateResults.Male.RepMeasANOVAs.p_vals.Genotype.BL.DP.Wake =  RepMeas_ANOVA_tables_Males_hourly_percentageLast12hrBL.BL_custom.Wake{'Genotype','pValue'};
TimeInStateResults.Male.RepMeasANOVAs.p_vals.Genotype.BL.DP.NREM =  RepMeas_ANOVA_tables_Males_hourly_percentageLast12hrBL.BL_custom.NREM{'Genotype','pValue'};
TimeInStateResults.Male.RepMeasANOVAs.p_vals.Genotype.BL.DP.REM  =  RepMeas_ANOVA_tables_Males_hourly_percentageLast12hrBL.BL_custom.REM{'Genotype','pValue'};

% --- SD LP ---
TimeInStateResults.Male.RepMeasANOVAs.p_vals.Genotype.SD.LP.Wake =  RepMeas_ANOVA_tables_Males_hourly_percentageFirst12hrSD.SD_custom.Wake{'Genotype','pValue'};
TimeInStateResults.Male.RepMeasANOVAs.p_vals.Genotype.SD.LP.NREM =  RepMeas_ANOVA_tables_Males_hourly_percentageFirst12hrSD.SD_custom.NREM{'Genotype','pValue'};
TimeInStateResults.Male.RepMeasANOVAs.p_vals.Genotype.SD.LP.REM  =  RepMeas_ANOVA_tables_Males_hourly_percentageFirst12hrSD.SD_custom.REM{'Genotype','pValue'};

% --- SD DP ---
TimeInStateResults.Male.RepMeasANOVAs.p_vals.Genotype.SD.DP.Wake =  RepMeas_ANOVA_tables_Males_hourly_percentageLast12hrSD.SD_custom.Wake{'Genotype','pValue'};
TimeInStateResults.Male.RepMeasANOVAs.p_vals.Genotype.SD.DP.NREM =  RepMeas_ANOVA_tables_Males_hourly_percentageLast12hrSD.SD_custom.NREM{'Genotype','pValue'};
TimeInStateResults.Male.RepMeasANOVAs.p_vals.Genotype.SD.DP.REM  =  RepMeas_ANOVA_tables_Males_hourly_percentageLast12hrSD.SD_custom.REM{'Genotype','pValue'};


% ------------ Females -------------------
if ~isempty(RepMeas_ANOVA_tables_Females_hourly_percentageFirst12hrBL)
	% --- BL LP ---
	TimeInStateResults.Female.RepMeasANOVAs.p_vals.Genotype.BL.LP.Wake =  RepMeas_ANOVA_tables_Females_hourly_percentageFirst12hrBL.BL_custom.Wake{'Genotype','pValue'};
	TimeInStateResults.Female.RepMeasANOVAs.p_vals.Genotype.BL.LP.NREM =  RepMeas_ANOVA_tables_Females_hourly_percentageFirst12hrBL.BL_custom.NREM{'Genotype','pValue'};
	TimeInStateResults.Female.RepMeasANOVAs.p_vals.Genotype.BL.LP.REM  =  RepMeas_ANOVA_tables_Females_hourly_percentageFirst12hrBL.BL_custom.REM{'Genotype','pValue'};

	% --- BL DP ---
	TimeInStateResults.Female.RepMeasANOVAs.p_vals.Genotype.BL.DP.Wake =  RepMeas_ANOVA_tables_Females_hourly_percentageLast12hrBL.BL_custom.Wake{'Genotype','pValue'};
	TimeInStateResults.Female.RepMeasANOVAs.p_vals.Genotype.BL.DP.NREM =  RepMeas_ANOVA_tables_Females_hourly_percentageLast12hrBL.BL_custom.NREM{'Genotype','pValue'};
	TimeInStateResults.Female.RepMeasANOVAs.p_vals.Genotype.BL.DP.REM  =  RepMeas_ANOVA_tables_Females_hourly_percentageLast12hrBL.BL_custom.REM{'Genotype','pValue'};

	% --- SD LP ---
	TimeInStateResults.Female.RepMeasANOVAs.p_vals.Genotype.SD.LP.Wake =  RepMeas_ANOVA_tables_Females_hourly_percentageFirst12hrSD.SD_custom.Wake{'Genotype','pValue'};
	TimeInStateResults.Female.RepMeasANOVAs.p_vals.Genotype.SD.LP.NREM =  RepMeas_ANOVA_tables_Females_hourly_percentageFirst12hrSD.SD_custom.NREM{'Genotype','pValue'};
	TimeInStateResults.Female.RepMeasANOVAs.p_vals.Genotype.SD.LP.REM  =  RepMeas_ANOVA_tables_Females_hourly_percentageFirst12hrSD.SD_custom.REM{'Genotype','pValue'};

	% --- SD DP ---
	TimeInStateResults.Female.RepMeasANOVAs.p_vals.Genotype.SD.DP.Wake =  RepMeas_ANOVA_tables_Females_hourly_percentageLast12hrSD.SD_custom.Wake{'Genotype','pValue'};
	TimeInStateResults.Female.RepMeasANOVAs.p_vals.Genotype.SD.DP.NREM =  RepMeas_ANOVA_tables_Females_hourly_percentageLast12hrSD.SD_custom.NREM{'Genotype','pValue'};
	TimeInStateResults.Female.RepMeasANOVAs.p_vals.Genotype.SD.DP.REM  =  RepMeas_ANOVA_tables_Females_hourly_percentageLast12hrSD.SD_custom.REM{'Genotype','pValue'};
else 
	% --- BL LP ---
	TimeInStateResults.Female.RepMeasANOVAs.p_vals.Genotype.BL.LP.Wake =  {};
	TimeInStateResults.Female.RepMeasANOVAs.p_vals.Genotype.BL.LP.NREM =  {};
	TimeInStateResults.Female.RepMeasANOVAs.p_vals.Genotype.BL.LP.REM  =  {};

	% --- BL DP ---
	TimeInStateResults.Female.RepMeasANOVAs.p_vals.Genotype.BL.DP.Wake = {};
	TimeInStateResults.Female.RepMeasANOVAs.p_vals.Genotype.BL.DP.NREM = {};
	TimeInStateResults.Female.RepMeasANOVAs.p_vals.Genotype.BL.DP.REM  = {};

	% --- SD LP ---
	TimeInStateResults.Female.RepMeasANOVAs.p_vals.Genotype.SD.LP.Wake =  {};
	TimeInStateResults.Female.RepMeasANOVAs.p_vals.Genotype.SD.LP.NREM =  {};
	TimeInStateResults.Female.RepMeasANOVAs.p_vals.Genotype.SD.LP.REM  =  {};

	% --- SD DP ---
	TimeInStateResults.Female.RepMeasANOVAs.p_vals.Genotype.SD.DP.Wake = {};
	TimeInStateResults.Female.RepMeasANOVAs.p_vals.Genotype.SD.DP.NREM = {};
	TimeInStateResults.Female.RepMeasANOVAs.p_vals.Genotype.SD.DP.REM  = {};
end	

% ------------------------ ANOVA Tables ----------------------------

% ------------ Males ---------------------
% --- BL LP ---
TimeInStateResults.Male.RepMeasANOVAs.Tables.Genotype.BL.LP.Wake =  RepMeas_ANOVA_tables_Males_hourly_percentageFirst12hrBL.BL_custom.Wake;
TimeInStateResults.Male.RepMeasANOVAs.Tables.Genotype.BL.LP.NREM =  RepMeas_ANOVA_tables_Males_hourly_percentageFirst12hrBL.BL_custom.NREM;
TimeInStateResults.Male.RepMeasANOVAs.Tables.Genotype.BL.LP.REM  =  RepMeas_ANOVA_tables_Males_hourly_percentageFirst12hrBL.BL_custom.REM;

% --- BL DP ---
TimeInStateResults.Male.RepMeasANOVAs.Tables.Genotype.BL.DP.Wake =  RepMeas_ANOVA_tables_Males_hourly_percentageLast12hrBL.BL_custom.Wake;
TimeInStateResults.Male.RepMeasANOVAs.Tables.Genotype.BL.DP.NREM =  RepMeas_ANOVA_tables_Males_hourly_percentageLast12hrBL.BL_custom.NREM;
TimeInStateResults.Male.RepMeasANOVAs.Tables.Genotype.BL.DP.REM  =  RepMeas_ANOVA_tables_Males_hourly_percentageLast12hrBL.BL_custom.REM;

% --- SD LP ---
TimeInStateResults.Male.RepMeasANOVAs.Tables.Genotype.SD.LP.Wake =  RepMeas_ANOVA_tables_Males_hourly_percentageFirst12hrSD.SD_custom.Wake;
TimeInStateResults.Male.RepMeasANOVAs.Tables.Genotype.SD.LP.NREM =  RepMeas_ANOVA_tables_Males_hourly_percentageFirst12hrSD.SD_custom.NREM;
TimeInStateResults.Male.RepMeasANOVAs.Tables.Genotype.SD.LP.REM  =  RepMeas_ANOVA_tables_Males_hourly_percentageFirst12hrSD.SD_custom.REM;

% --- SD DP ---
TimeInStateResults.Male.RepMeasANOVAs.Tables.Genotype.SD.DP.Wake =  RepMeas_ANOVA_tables_Males_hourly_percentageLast12hrSD.SD_custom.Wake;
TimeInStateResults.Male.RepMeasANOVAs.Tables.Genotype.SD.DP.NREM =  RepMeas_ANOVA_tables_Males_hourly_percentageLast12hrSD.SD_custom.NREM;
TimeInStateResults.Male.RepMeasANOVAs.Tables.Genotype.SD.DP.REM  =  RepMeas_ANOVA_tables_Males_hourly_percentageLast12hrSD.SD_custom.REM;


% ------------ Females ---------------------
if ~isempty(RepMeas_ANOVA_tables_Females_hourly_percentageFirst12hrBL)
	% --- BL LP ---
	TimeInStateResults.Female.RepMeasANOVAs.Tables.Genotype.BL.LP.Wake =  RepMeas_ANOVA_tables_Females_hourly_percentageFirst12hrBL.BL_custom.Wake;
	TimeInStateResults.Female.RepMeasANOVAs.Tables.Genotype.BL.LP.NREM =  RepMeas_ANOVA_tables_Females_hourly_percentageFirst12hrBL.BL_custom.NREM;
	TimeInStateResults.Female.RepMeasANOVAs.Tables.Genotype.BL.LP.REM  =  RepMeas_ANOVA_tables_Females_hourly_percentageFirst12hrBL.BL_custom.REM;

	% --- BL DP ---
	TimeInStateResults.Female.RepMeasANOVAs.Tables.Genotype.BL.DP.Wake =  RepMeas_ANOVA_tables_Females_hourly_percentageLast12hrBL.BL_custom.Wake;
	TimeInStateResults.Female.RepMeasANOVAs.Tables.Genotype.BL.DP.NREM =  RepMeas_ANOVA_tables_Females_hourly_percentageLast12hrBL.BL_custom.NREM;
	TimeInStateResults.Female.RepMeasANOVAs.Tables.Genotype.BL.DP.REM  =  RepMeas_ANOVA_tables_Females_hourly_percentageLast12hrBL.BL_custom.REM;

	% --- SD LP ---
	TimeInStateResults.Female.RepMeasANOVAs.Tables.Genotype.SD.LP.Wake =  RepMeas_ANOVA_tables_Females_hourly_percentageFirst12hrSD.SD_custom.Wake;
	TimeInStateResults.Female.RepMeasANOVAs.Tables.Genotype.SD.LP.NREM =  RepMeas_ANOVA_tables_Females_hourly_percentageFirst12hrSD.SD_custom.NREM;
	TimeInStateResults.Female.RepMeasANOVAs.Tables.Genotype.SD.LP.REM  =  RepMeas_ANOVA_tables_Females_hourly_percentageFirst12hrSD.SD_custom.REM;

	% --- SD DP ---
	TimeInStateResults.Female.RepMeasANOVAs.Tables.Genotype.SD.DP.Wake =  RepMeas_ANOVA_tables_Females_hourly_percentageLast12hrSD.SD_custom.Wake;
	TimeInStateResults.Female.RepMeasANOVAs.Tables.Genotype.SD.DP.NREM =  RepMeas_ANOVA_tables_Females_hourly_percentageLast12hrSD.SD_custom.NREM;
	TimeInStateResults.Female.RepMeasANOVAs.Tables.Genotype.SD.DP.REM  =  RepMeas_ANOVA_tables_Females_hourly_percentageLast12hrSD.SD_custom.REM;
else 
	% --- BL LP ---
	TimeInStateResults.Female.RepMeasANOVAs.Tables.Genotype.BL.LP.Wake =  {};
	TimeInStateResults.Female.RepMeasANOVAs.Tables.Genotype.BL.LP.NREM =  {};
	TimeInStateResults.Female.RepMeasANOVAs.Tables.Genotype.BL.LP.REM  =  {};

	% --- BL DP ---
	TimeInStateResults.Female.RepMeasANOVAs.Tables.Genotype.BL.DP.Wake = {};
	TimeInStateResults.Female.RepMeasANOVAs.Tables.Genotype.BL.DP.NREM = {};
	TimeInStateResults.Female.RepMeasANOVAs.Tables.Genotype.BL.DP.REM  = {};

	% --- SD LP ---
	TimeInStateResults.Female.RepMeasANOVAs.Tables.Genotype.SD.LP.Wake =  {};
	TimeInStateResults.Female.RepMeasANOVAs.Tables.Genotype.SD.LP.NREM =  {};
	TimeInStateResults.Female.RepMeasANOVAs.Tables.Genotype.SD.LP.REM  =  {};

	% --- SD DP ---
	TimeInStateResults.Female.RepMeasANOVAs.Tables.Genotype.SD.DP.Wake = {};
	TimeInStateResults.Female.RepMeasANOVAs.Tables.Genotype.SD.DP.NREM = {};
	TimeInStateResults.Female.RepMeasANOVAs.Tables.Genotype.SD.DP.REM  = {};
end 
% ------------------------------------------------------------------------------------------------------
% ------------------------------------------------------------------------------------------------------


  
% ------------------------------------------  BL vs SD -------------------------------------------------

% ------------------------ P-values ----------------------------

% ------------ Males ---------------------
% --- BL LP ---
TimeInStateResults.Male.RepMeasANOVAs.p_vals.BLvsSD.LP.Wake =  RepMeas_ANOVA_tables_Males_hourly_percentageFirst12hrBL.BL_custom.Wake{'Genotype','pValue'};
TimeInStateResults.Male.RepMeasANOVAs.p_vals.BLvsSD.LP.NREM =  RepMeas_ANOVA_tables_Males_hourly_percentageFirst12hrBL.BL_custom.NREM{'Genotype','pValue'};
TimeInStateResults.Male.RepMeasANOVAs.p_vals.BLvsSD.LP.REM  =  RepMeas_ANOVA_tables_Males_hourly_percentageFirst12hrBL.BL_custom.REM{'Genotype','pValue'};

% --- BL DP ---
TimeInStateResults.Male.RepMeasANOVAs.p_vals.BLvsSD.DP.Wake =  RepMeas_ANOVA_tables_Males_hourly_percentageLast12hrBL.BL_custom.Wake{'Genotype','pValue'};
TimeInStateResults.Male.RepMeasANOVAs.p_vals.BLvsSD.DP.NREM =  RepMeas_ANOVA_tables_Males_hourly_percentageLast12hrBL.BL_custom.NREM{'Genotype','pValue'};
TimeInStateResults.Male.RepMeasANOVAs.p_vals.BLvsSD.DP.REM  =  RepMeas_ANOVA_tables_Males_hourly_percentageLast12hrBL.BL_custom.REM{'Genotype','pValue'};

% --- SD LP ---
TimeInStateResults.Male.RepMeasANOVAs.p_vals.BLvsSD.LP.Wake =  RepMeas_ANOVA_tables_Males_hourly_percentageFirst12hrSD.SD_custom.Wake{'Genotype','pValue'};
TimeInStateResults.Male.RepMeasANOVAs.p_vals.BLvsSD.LP.NREM =  RepMeas_ANOVA_tables_Males_hourly_percentageFirst12hrSD.SD_custom.NREM{'Genotype','pValue'};
TimeInStateResults.Male.RepMeasANOVAs.p_vals.BLvsSD.LP.REM  =  RepMeas_ANOVA_tables_Males_hourly_percentageFirst12hrSD.SD_custom.REM{'Genotype','pValue'};

% --- SD DP ---
TimeInStateResults.Male.RepMeasANOVAs.p_vals.BLvsSD.DP.Wake =  RepMeas_ANOVA_tables_Males_hourly_percentageLast12hrSD.SD_custom.Wake{'Genotype','pValue'};
TimeInStateResults.Male.RepMeasANOVAs.p_vals.BLvsSD.DP.NREM =  RepMeas_ANOVA_tables_Males_hourly_percentageLast12hrSD.SD_custom.NREM{'Genotype','pValue'};
TimeInStateResults.Male.RepMeasANOVAs.p_vals.BLvsSD.DP.REM  =  RepMeas_ANOVA_tables_Males_hourly_percentageLast12hrSD.SD_custom.REM{'Genotype','pValue'};


% ------------ Females -------------------
if ~isempty(RepMeas_ANOVA_tables_Females_hourly_percentageFirst12hrBL)
	% --- BL LP ---
	TimeInStateResults.Female.RepMeasANOVAs.p_vals.BLvsSD.LP.Wake =  RepMeas_ANOVA_tables_Females_hourly_percentageFirst12hrBL.BL_custom.Wake{'Genotype','pValue'};
	TimeInStateResults.Female.RepMeasANOVAs.p_vals.BLvsSD.LP.NREM =  RepMeas_ANOVA_tables_Females_hourly_percentageFirst12hrBL.BL_custom.NREM{'Genotype','pValue'};
	TimeInStateResults.Female.RepMeasANOVAs.p_vals.BLvsSD.LP.REM  =  RepMeas_ANOVA_tables_Females_hourly_percentageFirst12hrBL.BL_custom.REM{'Genotype','pValue'};

	% --- BL DP ---
	TimeInStateResults.Female.RepMeasANOVAs.p_vals.BLvsSD.DP.Wake =  RepMeas_ANOVA_tables_Females_hourly_percentageLast12hrBL.BL_custom.Wake{'Genotype','pValue'};
	TimeInStateResults.Female.RepMeasANOVAs.p_vals.BLvsSD.DP.NREM =  RepMeas_ANOVA_tables_Females_hourly_percentageLast12hrBL.BL_custom.NREM{'Genotype','pValue'};
	TimeInStateResults.Female.RepMeasANOVAs.p_vals.BLvsSD.DP.REM  =  RepMeas_ANOVA_tables_Females_hourly_percentageLast12hrBL.BL_custom.REM{'Genotype','pValue'};

	% --- SD LP ---
	TimeInStateResults.Female.RepMeasANOVAs.p_vals.BLvsSD.LP.Wake =  RepMeas_ANOVA_tables_Females_hourly_percentageFirst12hrSD.SD_custom.Wake{'Genotype','pValue'};
	TimeInStateResults.Female.RepMeasANOVAs.p_vals.BLvsSD.LP.NREM =  RepMeas_ANOVA_tables_Females_hourly_percentageFirst12hrSD.SD_custom.NREM{'Genotype','pValue'};
	TimeInStateResults.Female.RepMeasANOVAs.p_vals.BLvsSD.LP.REM  =  RepMeas_ANOVA_tables_Females_hourly_percentageFirst12hrSD.SD_custom.REM{'Genotype','pValue'};

	% --- SD DP ---
	TimeInStateResults.Female.RepMeasANOVAs.p_vals.BLvsSD.DP.Wake =  RepMeas_ANOVA_tables_Females_hourly_percentageLast12hrSD.SD_custom.Wake{'Genotype','pValue'};
	TimeInStateResults.Female.RepMeasANOVAs.p_vals.BLvsSD.DP.NREM =  RepMeas_ANOVA_tables_Females_hourly_percentageLast12hrSD.SD_custom.NREM{'Genotype','pValue'};
	TimeInStateResults.Female.RepMeasANOVAs.p_vals.BLvsSD.DP.REM  =  RepMeas_ANOVA_tables_Females_hourly_percentageLast12hrSD.SD_custom.REM{'Genotype','pValue'};
else 
	% --- BL LP ---
	TimeInStateResults.Female.RepMeasANOVAs.p_vals.BLvsSD.LP.Wake =  {};
	TimeInStateResults.Female.RepMeasANOVAs.p_vals.BLvsSD.LP.NREM =  {};
	TimeInStateResults.Female.RepMeasANOVAs.p_vals.BLvsSD.LP.REM  =  {};

	% --- BL DP ---
	TimeInStateResults.Female.RepMeasANOVAs.p_vals.BLvsSD.DP.Wake = {};
	TimeInStateResults.Female.RepMeasANOVAs.p_vals.BLvsSD.DP.NREM = {};
	TimeInStateResults.Female.RepMeasANOVAs.p_vals.BLvsSD.DP.REM  = {};

	% --- SD LP ---
	TimeInStateResults.Female.RepMeasANOVAs.p_vals.BLvsSD.LP.Wake =  {};
	TimeInStateResults.Female.RepMeasANOVAs.p_vals.BLvsSD.LP.NREM =  {};
	TimeInStateResults.Female.RepMeasANOVAs.p_vals.BLvsSD.LP.REM  =  {};

	% --- SD DP ---
	TimeInStateResults.Female.RepMeasANOVAs.p_vals.BLvsSD.DP.Wake = {};
	TimeInStateResults.Female.RepMeasANOVAs.p_vals.BLvsSD.DP.NREM = {};
	TimeInStateResults.Female.RepMeasANOVAs.p_vals.BLvsSD.DP.REM  = {};
end	

% ------------------------ ANOVA Tables ----------------------------

% ------------ Males ---------------------
% --- BL LP ---
TimeInStateResults.Male.RepMeasANOVAs.Tables.BLvsSD.LP.Wake =  RepMeas_ANOVA_tables_Males_hourly_percentageFirst12hrBL.BL_custom.Wake;
TimeInStateResults.Male.RepMeasANOVAs.Tables.BLvsSD.LP.NREM =  RepMeas_ANOVA_tables_Males_hourly_percentageFirst12hrBL.BL_custom.NREM;
TimeInStateResults.Male.RepMeasANOVAs.Tables.BLvsSD.LP.REM  =  RepMeas_ANOVA_tables_Males_hourly_percentageFirst12hrBL.BL_custom.REM;

% --- BL DP ---
TimeInStateResults.Male.RepMeasANOVAs.Tables.BLvsSD.DP.Wake =  RepMeas_ANOVA_tables_Males_hourly_percentageLast12hrBL.BL_custom.Wake;
TimeInStateResults.Male.RepMeasANOVAs.Tables.BLvsSD.DP.NREM =  RepMeas_ANOVA_tables_Males_hourly_percentageLast12hrBL.BL_custom.NREM;
TimeInStateResults.Male.RepMeasANOVAs.Tables.BLvsSD.DP.REM  =  RepMeas_ANOVA_tables_Males_hourly_percentageLast12hrBL.BL_custom.REM;

% --- SD LP ---
TimeInStateResults.Male.RepMeasANOVAs.Tables.BLvsSD.LP.Wake =  RepMeas_ANOVA_tables_Males_hourly_percentageFirst12hrSD.SD_custom.Wake;
TimeInStateResults.Male.RepMeasANOVAs.Tables.BLvsSD.LP.NREM =  RepMeas_ANOVA_tables_Males_hourly_percentageFirst12hrSD.SD_custom.NREM;
TimeInStateResults.Male.RepMeasANOVAs.Tables.BLvsSD.LP.REM  =  RepMeas_ANOVA_tables_Males_hourly_percentageFirst12hrSD.SD_custom.REM;

% --- SD DP ---
TimeInStateResults.Male.RepMeasANOVAs.Tables.BLvsSD.DP.Wake =  RepMeas_ANOVA_tables_Males_hourly_percentageLast12hrSD.SD_custom.Wake;
TimeInStateResults.Male.RepMeasANOVAs.Tables.BLvsSD.DP.NREM =  RepMeas_ANOVA_tables_Males_hourly_percentageLast12hrSD.SD_custom.NREM;
TimeInStateResults.Male.RepMeasANOVAs.Tables.BLvsSD.DP.REM  =  RepMeas_ANOVA_tables_Males_hourly_percentageLast12hrSD.SD_custom.REM;


% ------------ Females ---------------------
if ~isempty(RepMeas_ANOVA_tables_Females_hourly_percentageFirst12hrBL)
	% --- BL LP ---
	TimeInStateResults.Female.RepMeasANOVAs.Tables.BLvsSD.LP.Wake =  RepMeas_ANOVA_tables_Females_hourly_percentageFirst12hrBL.BL_custom.Wake;
	TimeInStateResults.Female.RepMeasANOVAs.Tables.BLvsSD.LP.NREM =  RepMeas_ANOVA_tables_Females_hourly_percentageFirst12hrBL.BL_custom.NREM;
	TimeInStateResults.Female.RepMeasANOVAs.Tables.BLvsSD.LP.REM  =  RepMeas_ANOVA_tables_Females_hourly_percentageFirst12hrBL.BL_custom.REM;

	% --- BL DP ---
	TimeInStateResults.Female.RepMeasANOVAs.Tables.BLvsSD.DP.Wake =  RepMeas_ANOVA_tables_Females_hourly_percentageLast12hrBL.BL_custom.Wake;
	TimeInStateResults.Female.RepMeasANOVAs.Tables.BLvsSD.DP.NREM =  RepMeas_ANOVA_tables_Females_hourly_percentageLast12hrBL.BL_custom.NREM;
	TimeInStateResults.Female.RepMeasANOVAs.Tables.BLvsSD.DP.REM  =  RepMeas_ANOVA_tables_Females_hourly_percentageLast12hrBL.BL_custom.REM;

	% --- SD LP ---
	TimeInStateResults.Female.RepMeasANOVAs.Tables.BLvsSD.LP.Wake =  RepMeas_ANOVA_tables_Females_hourly_percentageFirst12hrSD.SD_custom.Wake;
	TimeInStateResults.Female.RepMeasANOVAs.Tables.BLvsSD.LP.NREM =  RepMeas_ANOVA_tables_Females_hourly_percentageFirst12hrSD.SD_custom.NREM;
	TimeInStateResults.Female.RepMeasANOVAs.Tables.BLvsSD.LP.REM  =  RepMeas_ANOVA_tables_Females_hourly_percentageFirst12hrSD.SD_custom.REM;

	% --- SD DP ---
	TimeInStateResults.Female.RepMeasANOVAs.Tables.BLvsSD.DP.Wake =  RepMeas_ANOVA_tables_Females_hourly_percentageLast12hrSD.SD_custom.Wake;
	TimeInStateResults.Female.RepMeasANOVAs.Tables.BLvsSD.DP.NREM =  RepMeas_ANOVA_tables_Females_hourly_percentageLast12hrSD.SD_custom.NREM;
	TimeInStateResults.Female.RepMeasANOVAs.Tables.BLvsSD.DP.REM  =  RepMeas_ANOVA_tables_Females_hourly_percentageLast12hrSD.SD_custom.REM;
else 
	% --- BL LP ---
	TimeInStateResults.Female.RepMeasANOVAs.Tables.BLvsSD.LP.Wake =  {};
	TimeInStateResults.Female.RepMeasANOVAs.Tables.BLvsSD.LP.NREM =  {};
	TimeInStateResults.Female.RepMeasANOVAs.Tables.BLvsSD.LP.REM  =  {};

	% --- BL DP ---
	TimeInStateResults.Female.RepMeasANOVAs.Tables.BLvsSD.DP.Wake = {};
	TimeInStateResults.Female.RepMeasANOVAs.Tables.BLvsSD.DP.NREM = {};
	TimeInStateResults.Female.RepMeasANOVAs.Tables.BLvsSD.DP.REM  = {};

	% --- SD LP ---
	TimeInStateResults.Female.RepMeasANOVAs.Tables.BLvsSD.LP.Wake =  {};
	TimeInStateResults.Female.RepMeasANOVAs.Tables.BLvsSD.LP.NREM =  {};
	TimeInStateResults.Female.RepMeasANOVAs.Tables.BLvsSD.LP.REM  =  {};

	% --- SD DP ---
	TimeInStateResults.Female.RepMeasANOVAs.Tables.BLvsSD.DP.Wake = {};
	TimeInStateResults.Female.RepMeasANOVAs.Tables.BLvsSD.DP.NREM = {};
	TimeInStateResults.Female.RepMeasANOVAs.Tables.BLvsSD.DP.REM  = {};
end 














