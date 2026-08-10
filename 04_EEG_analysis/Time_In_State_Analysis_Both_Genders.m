function [RepMeasuresANOVA_tables_hourly_percentage, RM_model] = Time_In_State_Analysis_Both_Genders(options)
%
% 	USAGE: out = Time_In_State_Analysis_Both_Genders(WT_24_hr_data_struct,Mut_24_hr_data_struct)
%
% 	This function carries out the analysis of time in each arousal state and makes corresponding figure.
%   
% 	INPUTS: 	WT_24_hr_data_struct:  A struct containing the  


arguments
    options.TimeInState
    options.TimeFrameRMANOVA = []
    options.inBLorSD         = 'Both'
    options.Sexes            = []
    options.SexVarName       = []
end 

TimeInState_struct = options.TimeInState;
TimeFrameRMANOVA   = options.TimeFrameRMANOVA;
inBLorSD           = options.inBLorSD;
Sexes 			   = options.Sexes;
SexVarName  	   = options.SexVarName;

%Male data WT
Wake24hWTBL_Males = TimeInState_struct.WT.(Sexes{1}).BL.Wake;
Wake24hWTSD_Males = TimeInState_struct.WT.(Sexes{1}).SD.Wake;
NREM24hWTBL_Males = TimeInState_struct.WT.(Sexes{1}).BL.NREM;
NREM24hWTSD_Males = TimeInState_struct.WT.(Sexes{1}).SD.NREM;
REM24hWTBL_Males  = TimeInState_struct.WT.(Sexes{1}).BL.REM;
REM24hWTSD_Males  = TimeInState_struct.WT.(Sexes{1}).SD.REM;

%Male data Mut
Wake24hMutBL_Males = TimeInState_struct.Mut.(Sexes{1}).BL.Wake;
Wake24hMutSD_Males = TimeInState_struct.Mut.(Sexes{1}).SD.Wake;
NREM24hMutBL_Males = TimeInState_struct.Mut.(Sexes{1}).BL.NREM;
NREM24hMutSD_Males = TimeInState_struct.Mut.(Sexes{1}).SD.NREM;
REM24hMutBL_Males  = TimeInState_struct.Mut.(Sexes{1}).BL.REM;
REM24hMutSD_Males  = TimeInState_struct.Mut.(Sexes{1}).SD.REM;

%Female data WT
Wake24hWTBL_Females = TimeInState_struct.WT.(Sexes{2}).BL.Wake;
Wake24hWTSD_Females = TimeInState_struct.WT.(Sexes{2}).SD.Wake;
NREM24hWTBL_Females = TimeInState_struct.WT.(Sexes{2}).BL.NREM;
NREM24hWTSD_Females = TimeInState_struct.WT.(Sexes{2}).SD.NREM;
REM24hWTBL_Females  = TimeInState_struct.WT.(Sexes{2}).BL.REM;
REM24hWTSD_Females  = TimeInState_struct.WT.(Sexes{2}).SD.REM;

%Female data Mut
Wake24hMutBL_Females = TimeInState_struct.Mut.(Sexes{2}).BL.Wake;
Wake24hMutSD_Females = TimeInState_struct.Mut.(Sexes{2}).SD.Wake;
NREM24hMutBL_Females = TimeInState_struct.Mut.(Sexes{2}).BL.NREM;
NREM24hMutSD_Females = TimeInState_struct.Mut.(Sexes{2}).SD.NREM;
REM24hMutBL_Females  = TimeInState_struct.Mut.(Sexes{2}).BL.REM;
REM24hMutSD_Females  = TimeInState_struct.Mut.(Sexes{2}).SD.REM;

% Convert to percentages by hour (with nicer names)
% Wake BL (WT+Mut, M+F)
WakePercByHour_WT_BL_Males   = Wake24hWTBL_Males/60;
WakePercByHour_Mut_BL_Males  = Wake24hMutBL_Males/60;
WakePercByHour_WT_BL_Females = Wake24hWTBL_Females/60;
WakePercByHour_Mut_BL_Females= Wake24hMutBL_Females/60;

% NREM BL (WT+Mut, M+F)
NREM_PercByHour_WT_BL_Males    = NREM24hWTBL_Males/60;
NREM_PercByHour_Mut_BL_Males   = NREM24hMutBL_Males/60;
NREM_PercByHour_WT_BL_Females  = NREM24hWTBL_Females/60;
NREM_PercByHour_Mut_BL_Females = NREM24hMutBL_Females/60;

% REM BL (WT+Mut, M+F)
REM_PercByHour_WT_BL_Males    = REM24hWTBL_Males/60;
REM_PercByHour_Mut_BL_Males   = REM24hMutBL_Males/60;
REM_PercByHour_WT_BL_Females  = REM24hWTBL_Females/60;
REM_PercByHour_Mut_BL_Females = REM24hMutBL_Females/60;

% Wake SD (WT+Mut, M+F)
WakePercByHour_WT_SD_Males   = Wake24hWTSD_Males/60;
WakePercByHour_Mut_SD_Males  = Wake24hMutSD_Males/60;
WakePercByHour_WT_SD_Females = Wake24hWTSD_Females/60;
WakePercByHour_Mut_SD_Females= Wake24hMutSD_Females/60;

% NREM SD (WT+Mut, M+F)
NREMPercByHour_WT_SD_Males   = NREM24hWTSD_Males/60;
NREMPercByHour_Mut_SD_Males  = NREM24hMutSD_Males/60;
NREMPercByHour_WT_SD_Females = NREM24hWTSD_Females/60;
NREMPercByHour_Mut_SD_Females= NREM24hMutSD_Females/60;

% REM SD (WT+Mut, M+F)
REMPercByHour_WT_SD_Males   = REM24hWTSD_Males/60;
REMPercByHour_Mut_SD_Males  = REM24hMutSD_Males/60;
REMPercByHour_WT_SD_Females = REM24hWTSD_Females/60;
REMPercByHour_Mut_SD_Females= REM24hMutSD_Females/60;

% ---------------------------------------------------------------------------------------
% --- RUN THE STATS ---------------------------------------------------------------------

% % Baseline 
% % Repeated Measures Wake Percentage by Hour (Baseline)
% WakePercByHour_WT_BL_Males   = Wake24hWTBL_Males/60;
% WakePercByHour_Mut_BL_Males  = Wake24hMutBL_Males/60;
% WakePercByHour_WT_BL_Females = Wake24hWTBL_Females/60;
% WakePercByHour_Mut_BL_Females= Wake24hMutBL_Females/60;
% %Genotype = categorical([0 0 0 0 0 1 1 1 1 1 1]');  % Categorical is important!  otherwise ranova produces incorrect F values
% %Genotype = categorical([zeros(1,size(WakePercByHour_Mut_BL,2)) ones(1,size(WakePercByHour_WT_BL,2))]');  % Categorical is important!  otherwise ranova produces incorrect F values
% Genotype = categorical([zeros(1,size(WakePercByHour_Mut_BL_Males,2)) ones(1,size(WakePercByHour_WT_BL_Males,2)) zeros(1,size(WakePercByHour_Mut_BL_Females,2)) ones(1,size(WakePercByHour_WT_BL_Females,2))]');  % Categorical is important!  otherwise ranova produces incorrect F values 0 for mutants, 1 for WT
% Sex      = categorical([zeros(1,size(WakePercByHour_Mut_BL_Males,2)) zeros(1,size(WakePercByHour_WT_BL_Males,2)) ones(1,size(WakePercByHour_Mut_BL_Females,2)) ones(1,size(WakePercByHour_WT_BL_Females,2))]');

% Y = [WakePercByHour_Mut_BL_Males'; WakePercByHour_WT_BL_Males';WakePercByHour_Mut_BL_Females'; WakePercByHour_WT_BL_Females'];
% t = table(Genotype,Sex,Y(:,1),Y(:,2),Y(:,3),Y(:,4),Y(:,5),Y(:,6),Y(:,7),Y(:,8),Y(:,9),Y(:,10),Y(:,11),Y(:,12), ...
%                    Y(:,13),Y(:,14),Y(:,15),Y(:,16),Y(:,17),Y(:,18),Y(:,19),Y(:,20),Y(:,21),Y(:,22),Y(:,23),Y(:,24), ...
%                    'VariableNames',{'Genotype','Sex','Hour1','Hour2','Hour3','Hour4','Hour5','Hour6','Hour7','Hour8',...
%                     'Hour9','Hour10','Hour11','Hour12','Hour13','Hour14','Hour15','Hour16','Hour17','Hour18','Hour19','Hour20','Hour21','Hour22','Hour23','Hour24'});

% Time = 1:24;   % Within-subjects variable
% WithinStructure = table(Time','VariableNames',{'Time'});
% WithinStructure.Time = categorical(WithinStructure.Time);

% rm = fitrm(t,'Hour1-Hour24 ~ Genotype*Sex','WithinDesign',WithinStructure,'WithinModel','Time');
% RManovatbl_Wake_Percentage_Hourly_BL_interaction = ranova(rm,'WithinModel','Time');

% % Repeated Measures SWS Percentage by Hour (Baseline)
% NREM_PercByHour_WT_BL_Males    = NREM24hWTBL_Males/60;
% NREM_PercByHour_Mut_BL_Males   = NREM24hMutBL_Males/60;
% NREM_PercByHour_WT_BL_Females  = NREM24hWTBL_Females/60;
% NREM_PercByHour_Mut_BL_Females = NREM24hMutBL_Females/60;
% Y = [NREM_PercByHour_Mut_BL_Males'; NREM_PercByHour_WT_BL_Males'; NREM_PercByHour_Mut_BL_Females'; NREM_PercByHour_WT_BL_Females'];

% t = table(Genotype,Sex,Y(:,1),Y(:,2),Y(:,3),Y(:,4),Y(:,5),Y(:,6),Y(:,7),Y(:,8),Y(:,9),Y(:,10),Y(:,11),Y(:,12), ...
%                    Y(:,13),Y(:,14),Y(:,15),Y(:,16),Y(:,17),Y(:,18),Y(:,19),Y(:,20),Y(:,21),Y(:,22),Y(:,23),Y(:,24), ...
%                    'VariableNames',{'Genotype','Sex','Hour1','Hour2','Hour3','Hour4','Hour5','Hour6','Hour7','Hour8',...
%                     'Hour9','Hour10','Hour11','Hour12','Hour13','Hour14','Hour15','Hour16','Hour17','Hour18','Hour19','Hour20','Hour21','Hour22','Hour23','Hour24'});


% rm = fitrm(t,'Hour1-Hour24 ~ Genotype*Sex','WithinDesign',WithinStructure,'WithinModel','Time');
% RManovatbl_NREM_Percentage_Hourly_BL_interaction = ranova(rm,'WithinModel','Time');


% % Repeated Measures REM Percentage by Hour (Baseline)
% REM_PercByHour_WT_BL_Males    = REM24hWTBL_Males/60;
% REM_PercByHour_Mut_BL_Males   = REM24hMutBL_Males/60;
% REM_PercByHour_WT_BL_Females  = REM24hWTBL_Females/60;
% REM_PercByHour_Mut_BL_Females = REM24hMutBL_Females/60;
% Y = [REM_PercByHour_Mut_BL_Males'; REM_PercByHour_WT_BL_Males';REM_PercByHour_Mut_BL_Females'; REM_PercByHour_WT_BL_Females'];

% t = table(Genotype,Sex,Y(:,1),Y(:,2),Y(:,3),Y(:,4),Y(:,5),Y(:,6),Y(:,7),Y(:,8),Y(:,9),Y(:,10),Y(:,11),Y(:,12), ...
%                    Y(:,13),Y(:,14),Y(:,15),Y(:,16),Y(:,17),Y(:,18),Y(:,19),Y(:,20),Y(:,21),Y(:,22),Y(:,23),Y(:,24), ...
%                    'VariableNames',{'Genotype','Sex','Hour1','Hour2','Hour3','Hour4','Hour5','Hour6','Hour7','Hour8',...
%                     'Hour9','Hour10','Hour11','Hour12','Hour13','Hour14','Hour15','Hour16','Hour17','Hour18','Hour19','Hour20','Hour21','Hour22','Hour23','Hour24'});


% rm = fitrm(t,'Hour1-Hour24 ~ Genotype*Sex','WithinDesign',WithinStructure,'WithinModel','Time');
% RManovatbl_REM_Percentage_Hourly_BL_interaction = ranova(rm,'WithinModel','Time');


% % ---- Sleep Dep ----
% % Repeated Measures Wake Percentage by Hour (Sleep Dep)
% WakePercByHour_WT_SD_Males   = Wake24hWTSD_Males/60;
% WakePercByHour_Mut_SD_Males  = Wake24hMutSD_Males/60;
% WakePercByHour_WT_SD_Females = Wake24hWTSD_Females/60;
% WakePercByHour_Mut_SD_Females= Wake24hMutSD_Females/60;
% %Genotype = categorical([zeros(1,size(WakePercByHour_Mut_SD,2)) ones(1,size(WakePercByHour_WT_SD,2))]');  % Categorical is important!  otherwise ranova produces incorrect F values

% Y = [WakePercByHour_Mut_SD_Males'; WakePercByHour_WT_SD_Males'; WakePercByHour_Mut_SD_Females'; WakePercByHour_WT_SD_Females'];
% Genotype = categorical([zeros(1,size(WakePercByHour_Mut_SD_Males,2)) ones(1,size(WakePercByHour_WT_SD_Males,2)) zeros(1,size(WakePercByHour_Mut_SD_Females,2)) ones(1,size(WakePercByHour_WT_SD_Females,2))]');
% Sex      = categorical([zeros(1,size(WakePercByHour_Mut_SD_Males,2)) zeros(1,size(WakePercByHour_WT_SD_Males,2)) ones(1,size(WakePercByHour_Mut_SD_Females,2)) ones(1,size(WakePercByHour_WT_SD_Females,2))]');


% t = table(Genotype,Sex,Y(:,1),Y(:,2),Y(:,3),Y(:,4),Y(:,5),Y(:,6),Y(:,7),Y(:,8),Y(:,9),Y(:,10),Y(:,11),Y(:,12), ...
%                    Y(:,13),Y(:,14),Y(:,15),Y(:,16),Y(:,17),Y(:,18),Y(:,19),Y(:,20),Y(:,21),Y(:,22),Y(:,23),Y(:,24), ...
%                    'VariableNames',{'Genotype','Sex','Hour1','Hour2','Hour3','Hour4','Hour5','Hour6','Hour7','Hour8',...
%                     'Hour9','Hour10','Hour11','Hour12','Hour13','Hour14','Hour15','Hour16','Hour17','Hour18','Hour19','Hour20','Hour21','Hour22','Hour23','Hour24'});

% rm = fitrm(t,'Hour1-Hour24 ~ Genotype*Sex','WithinDesign',WithinStructure,'WithinModel','Time');
% RManovatbl_Wake_Percentage_Hourly_SD_interaction = ranova(rm,'WithinModel','Time');

% % Repeated Measures NREM Percentage by Hour (Sleep Dep)
% NREMPercByHour_WT_SD_Males   = NREM24hWTSD_Males/60;
% NREMPercByHour_Mut_SD_Males  = NREM24hMutSD_Males/60;
% NREMPercByHour_WT_SD_Females = NREM24hWTSD_Females/60;
% NREMPercByHour_Mut_SD_Females= NREM24hMutSD_Females/60;

% Y = [NREMPercByHour_Mut_SD_Males'; NREMPercByHour_WT_SD_Males';NREMPercByHour_Mut_SD_Females'; NREMPercByHour_WT_SD_Females'];
% t = table(Genotype,Sex,Y(:,1),Y(:,2),Y(:,3),Y(:,4),Y(:,5),Y(:,6),Y(:,7),Y(:,8),Y(:,9),Y(:,10),Y(:,11),Y(:,12), ...
%                    Y(:,13),Y(:,14),Y(:,15),Y(:,16),Y(:,17),Y(:,18),Y(:,19),Y(:,20),Y(:,21),Y(:,22),Y(:,23),Y(:,24), ...
%                    'VariableNames',{'Genotype','Sex','Hour1','Hour2','Hour3','Hour4','Hour5','Hour6','Hour7','Hour8',...
%                     'Hour9','Hour10','Hour11','Hour12','Hour13','Hour14','Hour15','Hour16','Hour17','Hour18','Hour19','Hour20','Hour21','Hour22','Hour23','Hour24'});


% rm = fitrm(t,'Hour1-Hour24 ~ Genotype*Sex','WithinDesign',WithinStructure,'WithinModel','Time');
% RManovatbl_NREM_Percentage_Hourly_SD_interaction = ranova(rm,'WithinModel','Time');

% % Repeated Measures REM Percentage by Hour (Sleep Dep)
% REMPercByHour_WT_SD_Males   = REM24hWTSD_Males/60;
% REMPercByHour_Mut_SD_Males  = REM24hMutSD_Males/60;
% REMPercByHour_WT_SD_Females = REM24hWTSD_Females/60;
% REMPercByHour_Mut_SD_Females= REM24hMutSD_Females/60;

% Y = [REMPercByHour_Mut_SD_Males'; REMPercByHour_WT_SD_Males';REMPercByHour_Mut_SD_Females'; REMPercByHour_WT_SD_Females'];
% t = table(Genotype,Sex,Y(:,1),Y(:,2),Y(:,3),Y(:,4),Y(:,5),Y(:,6),Y(:,7),Y(:,8),Y(:,9),Y(:,10),Y(:,11),Y(:,12), ...
%                    Y(:,13),Y(:,14),Y(:,15),Y(:,16),Y(:,17),Y(:,18),Y(:,19),Y(:,20),Y(:,21),Y(:,22),Y(:,23),Y(:,24), ...
%                    'VariableNames',{'Genotype','Sex','Hour1','Hour2','Hour3','Hour4','Hour5','Hour6','Hour7','Hour8',...
%                     'Hour9','Hour10','Hour11','Hour12','Hour13','Hour14','Hour15','Hour16','Hour17','Hour18','Hour19','Hour20','Hour21','Hour22','Hour23','Hour24'});


% rm = fitrm(t,'Hour1-Hour24 ~ Genotype*Sex','WithinDesign',WithinStructure,'WithinModel','Time');
% RManovatbl_REM_Percentage_Hourly_SD_interaction = ranova(rm,'WithinModel','Time');
% % ------------------------------------------------------------------------------------------------------


% % Baseline
% RepMeasuresANOVA_tables_hourly_percentage.BL.Wake_interaction = RManovatbl_Wake_Percentage_Hourly_BL_interaction;
% RepMeasuresANOVA_tables_hourly_percentage.BL.NREM_interaction = RManovatbl_NREM_Percentage_Hourly_BL_interaction;
% RepMeasuresANOVA_tables_hourly_percentage.BL.REM_interaction  = RManovatbl_REM_Percentage_Hourly_BL_interaction;

% % Sleep Dep
% RepMeasuresANOVA_tables_hourly_percentage.SD.Wake_interaction = RManovatbl_Wake_Percentage_Hourly_SD_interaction;
% RepMeasuresANOVA_tables_hourly_percentage.SD.NREM_interaction = RManovatbl_NREM_Percentage_Hourly_SD_interaction;
% RepMeasuresANOVA_tables_hourly_percentage.SD.REM_interaction  = RManovatbl_REM_Percentage_Hourly_SD_interaction;





% % Compare BL to Sleep Dep?  






% --- If requested, run a repeated measures ANOVA on a custom timeframe, given as hours from start of recording. [1 12] means first 12 hours. don't start with 0
% if nargin==5
% 	error('You called Time_In_State_Analysis.m with 5 arguments. Either 4 or 6 are acceptable.')
% end 
if ~isempty(TimeFrameRMANOVA)
	if strcmp(inBLorSD,'BL')
		% Initialize variables that will be used for all 3 RMANOVA calls.  
		%Genotype = categorical([0 0 0 0 0 1 1 1 1 1 1]');  % Categorical is important!  otherwise ranova produces incorrect F values
		%Genotype = categorical([zeros(1,size(WakePercByHour_Mut_BL_Males,2)) ones(1,size(WakePercByHour_WT_BL_Males,2)) zeros(1,size(WakePercByHour_Mut_BL_Females,2)) ones(1,size(WakePercByHour_WT_BL_Females,2))]');  % Categorical is important!  otherwise ranova produces incorrect F values
		Genotype = categorical([repmat({'Mut'},1,size(WakePercByHour_Mut_BL_Males,2)), repmat({'WT'},1,size(WakePercByHour_WT_BL_Males,2)), repmat({'Mut'},1,size(WakePercByHour_Mut_BL_Females,2)), repmat({'WT'},1,size(WakePercByHour_WT_BL_Females,2))]');  
		%Sex      = categorical([zeros(1,size(WakePercByHour_Mut_BL_Males,2)) zeros(1,size(WakePercByHour_WT_BL_Males,2)) ones(1,size(WakePercByHour_Mut_BL_Females,2)) ones(1,size(WakePercByHour_WT_BL_Females,2))]');
		Sex      = categorical([repmat({Sexes{1}},1,size(WakePercByHour_Mut_BL_Males,2)),repmat({Sexes{1}},1, size(WakePercByHour_WT_BL_Males,2)),repmat({Sexes{2}},1,size(WakePercByHour_Mut_BL_Females,2)),repmat({Sexes{2}},1,size(WakePercByHour_WT_BL_Females,2))]');

		Time = TimeFrameRMANOVA(1):TimeFrameRMANOVA(2);   % Within-subjects variable
		WithinStructure = table(Time','VariableNames',{'Time'});
		WithinStructure.Time = categorical(WithinStructure.Time);

		rm_model_string = ['Hour',num2str(TimeFrameRMANOVA(1)),'-Hour',num2str(TimeFrameRMANOVA(2)),' ~ Genotype*Sex']; 


		% Repeated Measures Wake Percentage by Hour (Baseline)
		WakePercByHour_WT_BL_Males   = Wake24hWTBL_Males/60;
		WakePercByHour_Mut_BL_Males  = Wake24hMutBL_Males/60;
		WakePercByHour_WT_BL_Females = Wake24hWTBL_Females/60;
		WakePercByHour_Mut_BL_Females= Wake24hMutBL_Females/60;
		
		Y = [WakePercByHour_Mut_BL_Males'; WakePercByHour_WT_BL_Males';WakePercByHour_Mut_BL_Females'; WakePercByHour_WT_BL_Females'];
		t_W_custom_BL = table(Genotype,Sex,Y(:,1),Y(:,2),Y(:,3),Y(:,4),Y(:,5),Y(:,6),Y(:,7),Y(:,8),Y(:,9),Y(:,10),Y(:,11),Y(:,12), ...
                   	   				   Y(:,13),Y(:,14),Y(:,15),Y(:,16),Y(:,17),Y(:,18),Y(:,19),Y(:,20),Y(:,21),Y(:,22),Y(:,23),Y(:,24), ...
                   	   				   'VariableNames',{'Genotype','Sex','Hour1','Hour2','Hour3','Hour4','Hour5','Hour6','Hour7','Hour8',...
                       				   'Hour9','Hour10','Hour11','Hour12','Hour13','Hour14','Hour15','Hour16','Hour17','Hour18','Hour19','Hour20','Hour21','Hour22','Hour23','Hour24'});

		
		rm_W_custom_BL = fitrm(t_W_custom_BL,rm_model_string,'WithinDesign',WithinStructure,'WithinModel','Time');
		RManovatbl_Wake_Percentage_Hourly_BL_custom_interaction = ranova(rm_W_custom_BL,'WithinModel','Time');

		% Repeated Measures SWS Percentage by Hour (Baseline)
		NREM_PercByHour_WT_BL_Males    = NREM24hWTBL_Males/60;
		NREM_PercByHour_Mut_BL_Males   = NREM24hMutBL_Males/60;
		NREM_PercByHour_WT_BL_Females  = NREM24hWTBL_Females/60;
		NREM_PercByHour_Mut_BL_Females = NREM24hMutBL_Females/60;
		Y = [NREM_PercByHour_Mut_BL_Males'; NREM_PercByHour_WT_BL_Males';NREM_PercByHour_Mut_BL_Females'; NREM_PercByHour_WT_BL_Females'];

		t_N_custom_BL = table(Genotype,Sex,Y(:,1),Y(:,2),Y(:,3),Y(:,4),Y(:,5),Y(:,6),Y(:,7),Y(:,8),Y(:,9),Y(:,10),Y(:,11),Y(:,12), ...
                   	   					Y(:,13),Y(:,14),Y(:,15),Y(:,16),Y(:,17),Y(:,18),Y(:,19),Y(:,20),Y(:,21),Y(:,22),Y(:,23),Y(:,24), ...
                       					'VariableNames',{'Genotype','Sex','Hour1','Hour2','Hour3','Hour4','Hour5','Hour6','Hour7','Hour8',...
                       					'Hour9','Hour10','Hour11','Hour12','Hour13','Hour14','Hour15','Hour16','Hour17','Hour18','Hour19','Hour20','Hour21','Hour22','Hour23','Hour24'});


		rm_N_custom_BL = fitrm(t_N_custom_BL,rm_model_string,'WithinDesign',WithinStructure,'WithinModel','Time');
		RManovatbl_NREM_Percentage_Hourly_BL_custom_interaction = ranova(rm_N_custom_BL,'WithinModel','Time');


		% Repeated Measures REM Percentage by Hour (Baseline)
		REM_PercByHour_WT_BL_Males    = REM24hWTBL_Males/60;
		REM_PercByHour_Mut_BL_Males   = REM24hMutBL_Males/60;
		REM_PercByHour_WT_BL_Females  = REM24hWTBL_Females/60;
		REM_PercByHour_Mut_BL_Females = REM24hMutBL_Females/60;
		Y = [REM_PercByHour_Mut_BL_Males'; REM_PercByHour_WT_BL_Males';REM_PercByHour_Mut_BL_Females'; REM_PercByHour_WT_BL_Females'];

		t_R_custom_BL = table(Genotype,Sex,Y(:,1),Y(:,2),Y(:,3),Y(:,4),Y(:,5),Y(:,6),Y(:,7),Y(:,8),Y(:,9),Y(:,10),Y(:,11),Y(:,12), ...
                   	   	   			   Y(:,13),Y(:,14),Y(:,15),Y(:,16),Y(:,17),Y(:,18),Y(:,19),Y(:,20),Y(:,21),Y(:,22),Y(:,23),Y(:,24), ...
                           			   'VariableNames',{'Genotype','Sex','Hour1','Hour2','Hour3','Hour4','Hour5','Hour6','Hour7','Hour8',...
                           			   'Hour9','Hour10','Hour11','Hour12','Hour13','Hour14','Hour15','Hour16','Hour17','Hour18','Hour19','Hour20','Hour21','Hour22','Hour23','Hour24'});


		rm_R_custom_BL = fitrm(t_R_custom_BL,rm_model_string,'WithinDesign',WithinStructure,'WithinModel','Time');
		RManovatbl_REM_Percentage_Hourly_BL_custom_interaction = ranova(rm_R_custom_BL,'WithinModel','Time');

		% Baseline
		RepMeasuresANOVA_tables_hourly_percentage.BL_custom_interaction.Wake = RManovatbl_Wake_Percentage_Hourly_BL_custom_interaction;
		RepMeasuresANOVA_tables_hourly_percentage.BL_custom_interaction.NREM = RManovatbl_NREM_Percentage_Hourly_BL_custom_interaction;
		RepMeasuresANOVA_tables_hourly_percentage.BL_custom_interaction.REM  = RManovatbl_REM_Percentage_Hourly_BL_custom_interaction;

		% fix the names (ie if Sex is really Age)
		if(~strcmp(SexVarName,'Sex'))
			RepMeasuresANOVA_tables_hourly_percentage.BL_custom_interaction.Wake = change_ANOVA_table_names(RepMeasuresANOVA_tables_hourly_percentage.BL_custom_interaction.Wake,SexVarName);
			RepMeasuresANOVA_tables_hourly_percentage.BL_custom_interaction.NREM = change_ANOVA_table_names(RepMeasuresANOVA_tables_hourly_percentage.BL_custom_interaction.NREM,SexVarName);
			RepMeasuresANOVA_tables_hourly_percentage.BL_custom_interaction.REM  = change_ANOVA_table_names(RepMeasuresANOVA_tables_hourly_percentage.BL_custom_interaction.REM,SexVarName);
		end

		RM_model.Wake = rm_W_custom_BL;
        RM_model.NREM = rm_N_custom_BL;
        RM_model.REM  = rm_R_custom_BL;

	elseif strcmp(inBLorSD,'SD')
	% --- Sleep Dep day (custom time frame) ----

		% Sleep Dep day (custom time frame)
		% Genotype = categorical([zeros(1,size(WakePercByHour_Mut_SD_Males,2)) ones(1,size(WakePercByHour_WT_SD_Males,2)) zeros(1,size(WakePercByHour_Mut_SD_Females,2)) ones(1,size(WakePercByHour_WT_SD_Females,2))]');  % Categorical is important!  otherwise ranova produces incorrect F values
		% Sex      = categorical([zeros(1,size(WakePercByHour_Mut_SD_Males,2)) zeros(1,size(WakePercByHour_WT_SD_Males,2)) ones(1,size(WakePercByHour_Mut_SD_Females,2)) ones(1,size(WakePercByHour_WT_SD_Females,2))]');
		Genotype = categorical([repmat({'Mut'},1,size(WakePercByHour_Mut_SD_Males,2)), repmat({'WT'},1,size(WakePercByHour_WT_SD_Males,2)), repmat({'Mut'},1,size(WakePercByHour_Mut_SD_Females,2)), repmat({'WT'},1,size(WakePercByHour_WT_SD_Females,2))]');  
		Sex      = categorical([repmat({Sexes{1}},1,size(WakePercByHour_Mut_SD_Males,2)),   repmat({Sexes{1}},1, size(WakePercByHour_WT_SD_Males,2)), repmat({Sexes{2}},1,size(WakePercByHour_Mut_SD_Females,2)),   repmat({Sexes{2}},1,size(WakePercByHour_WT_SD_Females,2))]');



		Time = TimeFrameRMANOVA(1):TimeFrameRMANOVA(2);   % Within-subjects variable
		WithinStructure = table(Time','VariableNames',{'Time'});
		WithinStructure.Time = categorical(WithinStructure.Time);

		rm_model_string = ['Hour',num2str(TimeFrameRMANOVA(1)),'-Hour',num2str(TimeFrameRMANOVA(2)),' ~ Genotype*Sex']; 

		
		% Wake percentage
		WakePercByHour_WT_SD_Males    = Wake24hWTSD_Males/60;
		WakePercByHour_Mut_SD_Males   = Wake24hMutSD_Males/60;
		WakePercByHour_WT_SD_Females  = Wake24hWTSD_Females/60;
		WakePercByHour_Mut_SD_Females = Wake24hMutSD_Females/60;
		
		Y = [WakePercByHour_Mut_SD_Males'; WakePercByHour_WT_SD_Males';WakePercByHour_Mut_SD_Females'; WakePercByHour_WT_SD_Females'];
		t_W_custom_SD = table(Genotype,Sex,Y(:,1),Y(:,2),Y(:,3),Y(:,4),Y(:,5),Y(:,6),Y(:,7),Y(:,8),Y(:,9),Y(:,10),Y(:,11),Y(:,12), ...
                   		   			   Y(:,13),Y(:,14),Y(:,15),Y(:,16),Y(:,17),Y(:,18),Y(:,19),Y(:,20),Y(:,21),Y(:,22),Y(:,23),Y(:,24), ...
                   		   			   'VariableNames',{'Genotype','Sex','Hour1','Hour2','Hour3','Hour4','Hour5','Hour6','Hour7','Hour8',...
                    	   			   'Hour9','Hour10','Hour11','Hour12','Hour13','Hour14','Hour15','Hour16','Hour17','Hour18','Hour19','Hour20','Hour21','Hour22','Hour23','Hour24'});

		rm_W_custom_SD = fitrm(t_W_custom_SD,rm_model_string,'WithinDesign',WithinStructure,'WithinModel','Time');
		RManovatbl_Wake_Percentage_Hourly_SD_custom_interaction = ranova(rm_W_custom_SD,'WithinModel','Time');

		
		% Repeated Measures NREM Percentage by Hour (Sleep Dep entire 24 hours of sleep dep day)
		NREMPercByHour_WT_SD_Males   = NREM24hWTSD_Males/60;
		NREMPercByHour_Mut_SD_Males  = NREM24hMutSD_Males/60;
		NREMPercByHour_WT_SD_Females = NREM24hWTSD_Females/60;
		NREMPercByHour_Mut_SD_Females= NREM24hMutSD_Females/60;

		Y = [NREMPercByHour_Mut_SD_Males'; NREMPercByHour_WT_SD_Males';NREMPercByHour_Mut_SD_Females'; NREMPercByHour_WT_SD_Females'];
		t_N_custom_SD = table(Genotype,Sex,Y(:,1),Y(:,2),Y(:,3),Y(:,4),Y(:,5),Y(:,6),Y(:,7),Y(:,8),Y(:,9),Y(:,10),Y(:,11),Y(:,12), ...
                   		   			   Y(:,13),Y(:,14),Y(:,15),Y(:,16),Y(:,17),Y(:,18),Y(:,19),Y(:,20),Y(:,21),Y(:,22),Y(:,23),Y(:,24), ...
                   		   			   'VariableNames',{'Genotype','Sex','Hour1','Hour2','Hour3','Hour4','Hour5','Hour6','Hour7','Hour8',...
                    	   			   'Hour9','Hour10','Hour11','Hour12','Hour13','Hour14','Hour15','Hour16','Hour17','Hour18','Hour19','Hour20','Hour21','Hour22','Hour23','Hour24'});


		rm_N_custom_SD = fitrm(t_N_custom_SD,rm_model_string,'WithinDesign',WithinStructure,'WithinModel','Time');
		RManovatbl_NREM_Percentage_Hourly_SD_custom_interaction = ranova(rm_N_custom_SD,'WithinModel','Time');

		% Repeated Measures REM Percentage by Hour (Sleep Dep entire 24 hours of sleep dep day)
		REMPercByHour_WT_SD_Males   = REM24hWTSD_Males/60;
		REMPercByHour_Mut_SD_Males  = REM24hMutSD_Males/60;
		REMPercByHour_WT_SD_Females = REM24hWTSD_Females/60;
		REMPercByHour_Mut_SD_Females= REM24hMutSD_Females/60;

		Y = [REMPercByHour_Mut_SD_Males'; REMPercByHour_WT_SD_Males';REMPercByHour_Mut_SD_Females'; REMPercByHour_WT_SD_Females'];
		t_R_custom_SD = table(Genotype,Sex,Y(:,1),Y(:,2),Y(:,3),Y(:,4),Y(:,5),Y(:,6),Y(:,7),Y(:,8),Y(:,9),Y(:,10),Y(:,11),Y(:,12), ...
                   		   			   Y(:,13),Y(:,14),Y(:,15),Y(:,16),Y(:,17),Y(:,18),Y(:,19),Y(:,20),Y(:,21),Y(:,22),Y(:,23),Y(:,24), ...
                   		   			   'VariableNames',{'Genotype','Sex','Hour1','Hour2','Hour3','Hour4','Hour5','Hour6','Hour7','Hour8',...
                    	   			   'Hour9','Hour10','Hour11','Hour12','Hour13','Hour14','Hour15','Hour16','Hour17','Hour18','Hour19','Hour20','Hour21','Hour22','Hour23','Hour24'});


		rm_R_custom_SD = fitrm(t_R_custom_SD,rm_model_string,'WithinDesign',WithinStructure,'WithinModel','Time');
		RManovatbl_REM_Percentage_Hourly_SD_custom_interaction = ranova(rm_R_custom_SD,'WithinModel','Time');
		% ------------------------------------------------------------------------------------------------------

		
		RepMeasuresANOVA_tables_hourly_percentage.SD_custom_interaction.Wake = RManovatbl_Wake_Percentage_Hourly_SD_custom_interaction;
		RepMeasuresANOVA_tables_hourly_percentage.SD_custom_interaction.NREM = RManovatbl_NREM_Percentage_Hourly_SD_custom_interaction;
		RepMeasuresANOVA_tables_hourly_percentage.SD_custom_interaction.REM  = RManovatbl_REM_Percentage_Hourly_SD_custom_interaction;

		% fix the names (ie if Sex is really Age)
		if(~strcmp(SexVarName,'Sex'))
			RepMeasuresANOVA_tables_hourly_percentage.SD_custom_interaction.Wake = change_ANOVA_table_names(RepMeasuresANOVA_tables_hourly_percentage.SD_custom_interaction.Wake,SexVarName);
			RepMeasuresANOVA_tables_hourly_percentage.SD_custom_interaction.NREM = change_ANOVA_table_names(RepMeasuresANOVA_tables_hourly_percentage.SD_custom_interaction.NREM,SexVarName);
			RepMeasuresANOVA_tables_hourly_percentage.SD_custom_interaction.REM  = change_ANOVA_table_names(RepMeasuresANOVA_tables_hourly_percentage.SD_custom_interaction.REM,SexVarName);
		end

		RM_model.Wake = rm_W_custom_SD;
        RM_model.NREM = rm_N_custom_SD;
        RM_model.REM  = rm_R_custom_SD;

	elseif strcmp(inBLorSD,'Both')
		
		error('You selected ''Both'' for inBlorSD.  This is not implemented yet.')

		% % --- Baseline --- 
		% Genotype = categorical([zeros(1,size(WakePercByHour_Mut_BL_Males,2)) ones(1,size(WakePercByHour_WT_BL_Males,2)) zeros(1,size(WakePercByHour_Mut_BL_Females,2)) ones(1,size(WakePercByHour_WT_BL_Females,2))]');  % Categorical is important!  otherwise ranova produces incorrect F values
		% Sex      = categorical([zeros(1,size(WakePercByHour_Mut_BL_Males,2)) zeros(1,size(WakePercByHour_WT_BL_Males,2)) ones(1,size(WakePercByHour_Mut_BL_Females,2)) ones(1,size(WakePercByHour_WT_BL_Females,2))]');

		% Time = TimeFrameRMANOVA(1):TimeFrameRMANOVA(2);   % Within-subjects variable
		% WithinStructure = table(Time','VariableNames',{'Time'});
		% WithinStructure.Time = categorical(WithinStructure.Time);

		% rm_model_string = ['Hour',num2str(TimeFrameRMANOVA(1)),'-Hour',num2str(TimeFrameRMANOVA(2)),' ~ Genotype*Sex']; 




		% % Repeated Measures Wake Percentage by Hour (Baseline)
		% WakePercByHour_WT_BL_Males   = Wake24hWTBL_Males/60;
		% WakePercByHour_Mut_BL_Males  = Wake24hMutBL_Males/60;
		% WakePercByHour_WT_BL_Females = Wake24hWTBL_Females/60;
		% WakePercByHour_Mut_BL_Females= Wake24hMutBL_Females/60;
		
		% Y = [WakePercByHour_Mut_BL_Males'; WakePercByHour_WT_BL_Males';WakePercByHour_Mut_BL_Females'; WakePercByHour_WT_BL_Females'];
		% t_W_custom_BL = table(Genotype,Sex,Y(:,1),Y(:,2),Y(:,3),Y(:,4),Y(:,5),Y(:,6),Y(:,7),Y(:,8),Y(:,9),Y(:,10),Y(:,11),Y(:,12), ...
        %            	   				   Y(:,13),Y(:,14),Y(:,15),Y(:,16),Y(:,17),Y(:,18),Y(:,19),Y(:,20),Y(:,21),Y(:,22),Y(:,23),Y(:,24), ...
        %            	   				   'VariableNames',{'Genotype','Sex','Hour1','Hour2','Hour3','Hour4','Hour5','Hour6','Hour7','Hour8',...
        %                				   'Hour9','Hour10','Hour11','Hour12','Hour13','Hour14','Hour15','Hour16','Hour17','Hour18','Hour19','Hour20','Hour21','Hour22','Hour23','Hour24'});

		
		% rm_W_custom_BL = fitrm(t_W_custom_BL,rm_model_string,'WithinDesign',WithinStructure,'WithinModel','Time');
		% RManovatbl_Wake_Percentage_Hourly_BL_custom_interaction = ranova(rm_W_custom_BL,'WithinModel','Time');

		% % Repeated Measures SWS Percentage by Hour (Baseline)
		% NREM_PercByHour_WT_BL_Males    = NREM24hWTBL_Males/60;
		% NREM_PercByHour_Mut_BL_Males   = NREM24hMutBL_Males/60;
		% NREM_PercByHour_WT_BL_Females  = NREM24hWTBL_Females/60;
		% NREM_PercByHour_Mut_BL_Females = NREM24hMutBL_Females/60;
		% Y = [NREM_PercByHour_Mut_BL_Males'; NREM_PercByHour_WT_BL_Males';NREM_PercByHour_Mut_BL_Females'; NREM_PercByHour_WT_BL_Females'];

		% t_N_custom_BL = table(Genotype,Sex,Y(:,1),Y(:,2),Y(:,3),Y(:,4),Y(:,5),Y(:,6),Y(:,7),Y(:,8),Y(:,9),Y(:,10),Y(:,11),Y(:,12), ...
        %            	   					Y(:,13),Y(:,14),Y(:,15),Y(:,16),Y(:,17),Y(:,18),Y(:,19),Y(:,20),Y(:,21),Y(:,22),Y(:,23),Y(:,24), ...
        %                					'VariableNames',{'Genotype','Sex','Hour1','Hour2','Hour3','Hour4','Hour5','Hour6','Hour7','Hour8',...
        %                					'Hour9','Hour10','Hour11','Hour12','Hour13','Hour14','Hour15','Hour16','Hour17','Hour18','Hour19','Hour20','Hour21','Hour22','Hour23','Hour24'});


		% rm_N_custom_BL = fitrm(t_N_custom_BL,rm_model_string,'WithinDesign',WithinStructure,'WithinModel','Time');
		% RManovatbl_NREM_Percentage_Hourly_BL_custom_interaction = ranova(rm_N_custom_BL,'WithinModel','Time');


		% % Repeated Measures REM Percentage by Hour (Baseline)
		% REM_PercByHour_WT_BL_Males    = REM24hWTBL_Males/60;
		% REM_PercByHour_Mut_BL_Males   = REM24hMutBL_Males/60;
		% REM_PercByHour_WT_BL_Females  = REM24hWTBL_Females/60;
		% REM_PercByHour_Mut_BL_Females = REM24hMutBL_Females/60;
		% Y = [REM_PercByHour_Mut_BL_Males'; REM_PercByHour_WT_BL_Males';REM_PercByHour_Mut_BL_Females'; REM_PercByHour_WT_BL_Females'];

		% t_R_custom_BL = table(Genotype,Sex,Y(:,1),Y(:,2),Y(:,3),Y(:,4),Y(:,5),Y(:,6),Y(:,7),Y(:,8),Y(:,9),Y(:,10),Y(:,11),Y(:,12), ...
        %            	   	   			   Y(:,13),Y(:,14),Y(:,15),Y(:,16),Y(:,17),Y(:,18),Y(:,19),Y(:,20),Y(:,21),Y(:,22),Y(:,23),Y(:,24), ...
        %                    			   'VariableNames',{'Genotype','Sex','Hour1','Hour2','Hour3','Hour4','Hour5','Hour6','Hour7','Hour8',...
        %                    			   'Hour9','Hour10','Hour11','Hour12','Hour13','Hour14','Hour15','Hour16','Hour17','Hour18','Hour19','Hour20','Hour21','Hour22','Hour23','Hour24'});


		% rm_R_custom_BL = fitrm(t_R_custom_BL,rm_model_string,'WithinDesign',WithinStructure,'WithinModel','Time');
		% RManovatbl_REM_Percentage_Hourly_BL_custom_interaction = ranova(rm_R_custom_BL,'WithinModel','Time');

		
		% % --- Sleep Dep (entire 24 hours of sleep dep day ----
		% Genotype = categorical([zeros(1,size(WakePercByHour_Mut_SD_Males,2)) ones(1,size(WakePercByHour_WT_SD_Males,2)) zeros(1,size(WakePercByHour_Mut_SD_Females,2)) ones(1,size(WakePercByHour_WT_SD_Females,2))]');  % Categorical is important!  otherwise ranova produces incorrect F values
		% Sex      = categorical([zeros(1,size(WakePercByHour_Mut_SD_Males,2)) zeros(1,size(WakePercByHour_WT_SD_Males,2)) ones(1,size(WakePercByHour_Mut_SD_Females,2)) ones(1,size(WakePercByHour_WT_SD_Females,2))]');

		% Time = TimeFrameRMANOVA(1):TimeFrameRMANOVA(2);   % Within-subjects variable
		% WithinStructure = table(Time','VariableNames',{'Time'});
		% WithinStructure.Time = categorical(WithinStructure.Time);

		% rm_model_string = ['Hour',num2str(TimeFrameRMANOVA(1)),'-Hour',num2str(TimeFrameRMANOVA(2)),' ~ Genotype*Sex']; 

		
		% % Wake 
		% WakePercByHour_WT_SD_Males    = Wake24hWTSD_Males/60;
		% WakePercByHour_Mut_SD_Males   = Wake24hMutSD_Males/60;
		% WakePercByHour_WT_SD_Females  = Wake24hWTSD_Females/60;
		% WakePercByHour_Mut_SD_Females = Wake24hMutSD_Females/60;
		
		% Y = [WakePercByHour_Mut_SD_Males'; WakePercByHour_WT_SD_Males';WakePercByHour_Mut_SD_Females'; WakePercByHour_WT_SD_Females'];
		% t_W_custom_SD = table(Genotype,Sex,Y(:,1),Y(:,2),Y(:,3),Y(:,4),Y(:,5),Y(:,6),Y(:,7),Y(:,8),Y(:,9),Y(:,10),Y(:,11),Y(:,12), ...
        %            		   			   Y(:,13),Y(:,14),Y(:,15),Y(:,16),Y(:,17),Y(:,18),Y(:,19),Y(:,20),Y(:,21),Y(:,22),Y(:,23),Y(:,24), ...
        %            		   			   'VariableNames',{'Genotype','Sex','Hour1','Hour2','Hour3','Hour4','Hour5','Hour6','Hour7','Hour8',...
        %             	   			   'Hour9','Hour10','Hour11','Hour12','Hour13','Hour14','Hour15','Hour16','Hour17','Hour18','Hour19','Hour20','Hour21','Hour22','Hour23','Hour24'});

		% rm_W_custom_SD = fitrm(t_W_custom_SD,rm_model_string,'WithinDesign',WithinStructure,'WithinModel','Time');
		% RManovatbl_Wake_Percentage_Hourly_SD_custom_interaction = ranova(rm_W_custom_SD,'WithinModel','Time');

		
		% % Repeated Measures NREM Percentage by Hour (Sleep Dep entire 24 hours)
		% NREMPercByHour_WT_SD_Males   = NREM24hWTSD_Males/60;
		% NREMPercByHour_Mut_SD_Males  = NREM24hMutSD_Males/60;
		% NREMPercByHour_WT_SD_Females = NREM24hWTSD_Females/60;
		% NREMPercByHour_Mut_SD_Females= NREM24hMutSD_Females/60;

		% Y = [NREMPercByHour_Mut_SD_Males'; NREMPercByHour_WT_SD_Males';NREMPercByHour_Mut_SD_Females'; NREMPercByHour_WT_SD_Females'];
		% t_N_custom_SD = table(Genotype,Sex,Y(:,1),Y(:,2),Y(:,3),Y(:,4),Y(:,5),Y(:,6),Y(:,7),Y(:,8),Y(:,9),Y(:,10),Y(:,11),Y(:,12), ...
        %            		   			   Y(:,13),Y(:,14),Y(:,15),Y(:,16),Y(:,17),Y(:,18),Y(:,19),Y(:,20),Y(:,21),Y(:,22),Y(:,23),Y(:,24), ...
        %            		   			   'VariableNames',{'Genotype','Sex','Hour1','Hour2','Hour3','Hour4','Hour5','Hour6','Hour7','Hour8',...
        %             	   			   'Hour9','Hour10','Hour11','Hour12','Hour13','Hour14','Hour15','Hour16','Hour17','Hour18','Hour19','Hour20','Hour21','Hour22','Hour23','Hour24'});


		% rm_N_custom_SD = fitrm(t_N_custom_SD,rm_model_string,'WithinDesign',WithinStructure,'WithinModel','Time');
		% RManovatbl_NREM_Percentage_Hourly_SD_custom_interaction = ranova(rm_N_custom_SD,'WithinModel','Time');

		% % Repeated Measures REM Percentage by Hour (Sleep Dep entire 24 hours)
		% REMPercByHour_WT_SD_Males   = REM24hWTSD_Males/60;
		% REMPercByHour_Mut_SD_Males  = REM24hMutSD_Males/60;
		% REMPercByHour_WT_SD_Females = REM24hWTSD_Females/60;
		% REMPercByHour_Mut_SD_Females= REM24hMutSD_Females/60;

		% Y = [REMPercByHour_Mut_SD_Males'; REMPercByHour_WT_SD_Males';REMPercByHour_Mut_SD_Females'; REMPercByHour_WT_SD_Females'];
		% t_R_custom_SD = table(Genotype,Sex,Y(:,1),Y(:,2),Y(:,3),Y(:,4),Y(:,5),Y(:,6),Y(:,7),Y(:,8),Y(:,9),Y(:,10),Y(:,11),Y(:,12), ...
        %            		   			   Y(:,13),Y(:,14),Y(:,15),Y(:,16),Y(:,17),Y(:,18),Y(:,19),Y(:,20),Y(:,21),Y(:,22),Y(:,23),Y(:,24), ...
        %            		   			   'VariableNames',{'Genotype','Sex','Hour1','Hour2','Hour3','Hour4','Hour5','Hour6','Hour7','Hour8',...
        %             	   			   'Hour9','Hour10','Hour11','Hour12','Hour13','Hour14','Hour15','Hour16','Hour17','Hour18','Hour19','Hour20','Hour21','Hour22','Hour23','Hour24'});


		% rm_R_custom_SD = fitrm(t_R_custom_SD,rm_model_string,'WithinDesign',WithinStructure,'WithinModel','Time');
		% RManovatbl_REM_Percentage_Hourly_SD_custom_interaction = ranova(rm_R_custom_SD,'WithinModel','Time');
		% % ------------------------------------------------------------------------------------------------------

		% % Baseline
		% RepMeasuresANOVA_tables_hourly_percentage.BL_custom_interaction.Wake = RManovatbl_Wake_Percentage_Hourly_BL_custom_interaction;
		% RepMeasuresANOVA_tables_hourly_percentage.BL_custom_interaction.NREM = RManovatbl_NREM_Percentage_Hourly_BL_custom_interaction;
		% RepMeasuresANOVA_tables_hourly_percentage.BL_custom_interaction.REM  = RManovatbl_REM_Percentage_Hourly_BL_custom_interaction;

		% % SD
		% RepMeasuresANOVA_tables_hourly_percentage.SD_custom_interaction.Wake = RManovatbl_Wake_Percentage_Hourly_SD_custom_interaction;
		% RepMeasuresANOVA_tables_hourly_percentage.SD_custom_interaction.NREM = RManovatbl_NREM_Percentage_Hourly_SD_custom_interaction;
		% RepMeasuresANOVA_tables_hourly_percentage.SD_custom_interaction.REM  = RManovatbl_REM_Percentage_Hourly_SD_custom_interaction;

	end 
end % end of if ~isempty(TimeFrameRMANOVA)


