function [RM_Anova_tables,PostHoc_tables] = Make_9_Panel_TimeInState_Figure(HourlyStatePercentages,TimeInState,SD_length,LegendLabels,MorF)
%
% USAGE:  Make_9_Panel_TimeInState_Figure(HourlyStatePercentages,SD_length,LegendLabels)
%
% This function makes a figure like Figure 3 figure supplement 1 in Ingiosi et al eLife 2019
% 
% Each row of the figure corresponds to states (W, NREM, REM) and each panel has two curves. 
%
% You may need to adjust the parameters in the next few lines.  




% Set some parameters to use for all panels (line thickness, font sizes, etc. )

LPboxpos = ([0 -1 12 5]);  %[x y w h] Light period box along bottom of plot
DPboxpos = ([12 -1 12 5]); %[x y w h] Dark period box along bottom of plot
SDboxpos = ([0 -1 SD_length 100]);	% [x y w h] Sleep Dep hatching

LPboxposREM = ([0 -1 12 1]);  %[x y w h] Light period box along bottom of plot REM VERSION: SHORTER
DPboxposREM = ([12 -1 12 1]); %[x y w h] Dark period box along bottom of plot  REM VERSION: SHORTER

Sig12hrColor = [0.8 0.8 0.8];  % RGB for color of bar that indicates significant RM ANova over 12 hours

myfontsize  = 14;
mylinewidth = 2; 
myaxislinewidth = 2;
EBMarkersize = 10;
SDbox_color = [0.5 0.5 ];

MutantText = LegendLabels{2};

% Handle the case where no data are present (no Female data for instance)
if isempty(HourlyStatePercentages)
	RM_Anova_tables = [];
	PostHoc_tables  = [];
	return
end 

% ------------------------------------------------------------------------------------------------------------------------
% --- Stats --------------------------------------------------------------------------------------------------------------
TimeFrameLP = [SD_length+1 12];
TimeFrameDP = [13 24];


sig_effect.WTBL_vs_WTSD_LP_W = 0;
sig_effect.WTBL_vs_WTSD_DP_W = 0;
sig_effect.WTBL_vs_WTSD_LP_N = 0;
sig_effect.WTBL_vs_WTSD_DP_N = 0;
sig_effect.WTBL_vs_WTSD_LP_R = 0;
sig_effect.WTBL_vs_WTSD_DP_R = 0;
sig_effect.MutBL_vs_MutSD_LP_W = 0;
sig_effect.MutBL_vs_MutSD_DP_W = 0;
sig_effect.MutBL_vs_MutSD_LP_N = 0;
sig_effect.MutBL_vs_MutSD_DP_N = 0;
sig_effect.MutBL_vs_MutSD_LP_R = 0;
sig_effect.MutBL_vs_MutSD_DP_R = 0;
sig_effect.WTSD_vs_MutSD_LP_W = 0;
sig_effect.WTSD_vs_MutSD_DP_W = 0;
sig_effect.WTSD_vs_MutSD_LP_N = 0;
sig_effect.WTSD_vs_MutSD_DP_N = 0;
sig_effect.WTSD_vs_MutSD_LP_R = 0;
sig_effect.WTSD_vs_MutSD_DP_R = 0;
% --- Panel 1: RM Anova WT BL vs WT SD (LP exc SD Wake) --------------------------------------------------------
WakePercByHour_WT_BL = TimeInState.WT.(MorF).BL.Wake/60;  % Divide by 60 to get percentages, rather than minutes  
WakePercByHour_WT_SD = TimeInState.WT.(MorF).SD.Wake/60;
BLorSD = categorical([zeros(1,size(WakePercByHour_WT_BL,2)) ones(1,size(WakePercByHour_WT_SD,2))]');  % Categorical is important!  otherwise ranova produces incorrect F values

Time1LP = TimeFrameLP(1):TimeFrameLP(2);   % Within-subjects variable
WithinStructure1LP = table(Time1LP','VariableNames',{'Time'});
WithinStructure1LP.Time = categorical(WithinStructure1LP.Time);

rm_model_string_BLorSD_LP = ['Hour',num2str(Time1LP(1)),'-Hour',num2str(Time1LP(end)),' ~ BLorSD']; 


YW1 = [WakePercByHour_WT_BL'; WakePercByHour_WT_SD'];
t_panel1 = table(BLorSD,YW1(:,1),YW1(:,2),YW1(:,3),YW1(:,4),YW1(:,5),YW1(:,6),YW1(:,7),YW1(:,8),YW1(:,9),YW1(:,10),YW1(:,11),YW1(:,12), ...
           	   				   YW1(:,13),YW1(:,14),YW1(:,15),YW1(:,16),YW1(:,17),YW1(:,18),YW1(:,19),YW1(:,20),YW1(:,21),YW1(:,22),YW1(:,23),YW1(:,24), ...
           	   				     'VariableNames',{'BLorSD','Hour1','Hour2','Hour3','Hour4','Hour5','Hour6','Hour7','Hour8',...
               				     'Hour9','Hour10','Hour11','Hour12','Hour13','Hour14','Hour15','Hour16','Hour17','Hour18','Hour19','Hour20','Hour21','Hour22','Hour23','Hour24'});


rm_panel1_LP = fitrm(t_panel1,rm_model_string_BLorSD_LP,'WithinDesign',WithinStructure1LP,'WithinModel','Time');
RManovatbl_WTBL_vs_WTSD_LP_W = ranova(rm_panel1_LP,'WithinModel','Time');
% -- Post-hocs --
if RManovatbl_WTBL_vs_WTSD_LP_W{'BLorSD','pValue'} < 0.05 | RManovatbl_WTBL_vs_WTSD_LP_W{'BLorSD:Time','pValue'} < 0.05
	sig_effect.WTBL_vs_WTSD_LP_W = 1;  % set boolean for significant effect 
	mc_table.panel1_LP 		  = multcompare(rm_panel1_LP,'BLorSD','By','Time','ComparisonType','lsd');
	mc_table.panel1_LP.pValue = mafdr(mc_table.panel1_LP.pValue,'BHFDR',true);
	[~, ia] = unique(mc_table.panel1_LP.Time(:), 'stable');  % remove duplicate rows
	mc_table.panel1_LP = mc_table.panel1_LP(ia,:);
	xlocs_sig_effects.WTBL_vs_WTSD_LP_W = double(string(mc_table.panel1_LP.Time(mc_table.panel1_LP.pValue<0.05))); % these are times (h), not indices
end
% --------------------------------------------------------------------------------------------------------------


% --- Panel 1: RM Anova WT BL vs WT SD (DP Wake) ---------------------------------------------------------------
Time1DP = TimeFrameDP(1):TimeFrameDP(2);   % Within-subjects variable
WithinStructure1DP = table(Time1DP','VariableNames',{'Time'});
WithinStructure1DP.Time = categorical(WithinStructure1DP.Time);
rm_model_string_BLorSD_DP = ['Hour',num2str(Time1DP(1)),'-Hour',num2str(Time1DP(end)),' ~ BLorSD']; 

rm_panel1_DP = fitrm(t_panel1,rm_model_string_BLorSD_DP,'WithinDesign',WithinStructure1DP,'WithinModel','Time');
RManovatbl_WTBL_vs_WTSD_DP_W = ranova(rm_panel1_DP,'WithinModel','Time');
% -- Post-hocs --
if RManovatbl_WTBL_vs_WTSD_DP_W{'BLorSD','pValue'} < 0.05 | RManovatbl_WTBL_vs_WTSD_DP_W{'BLorSD:Time','pValue'} < 0.05
	sig_effect.WTBL_vs_WTSD_DP_W = 1;  % set boolean for significant effect 
	mc_table.panel1_DP 		  = multcompare(rm_panel1_DP,'BLorSD','By','Time','ComparisonType','lsd');
	mc_table.panel1_DP.pValue = mafdr(mc_table.panel1_DP.pValue,'BHFDR',true);
	[~, ia] = unique(mc_table.panel1_DP.Time(:), 'stable');  % remove duplicate rows
	mc_table.panel1_DP = mc_table.panel1_DP(ia,:);
	xlocs_sig_effects.WTBL_vs_WTSD_DP_W = double(string(mc_table.panel1_DP.Time(mc_table.panel1_DP.pValue<0.05)));
end
% --------------------------------------------------------------------------------------------------------------


% --- Panel 2: RM Anova WT BL vs WT SD (LP exc SD NREM) --------------------------------------------------------
NREMPercByHour_WT_BL = TimeInState.WT.(MorF).BL.NREM/60;  % Divide by 60 to get percentages, rather than minutes  
NREMPercByHour_WT_SD = TimeInState.WT.(MorF).SD.NREM/60;

YN1 = [NREMPercByHour_WT_BL'; NREMPercByHour_WT_SD'];
t_panel2 = table(BLorSD,YN1(:,1),YN1(:,2),YN1(:,3),YN1(:,4),YN1(:,5),YN1(:,6),YN1(:,7),YN1(:,8),YN1(:,9),YN1(:,10),YN1(:,11),YN1(:,12), ...
           	   				   YN1(:,13),YN1(:,14),YN1(:,15),YN1(:,16),YN1(:,17),YN1(:,18),YN1(:,19),YN1(:,20),YN1(:,21),YN1(:,22),YN1(:,23),YN1(:,24), ...
           	   				     'VariableNames',{'BLorSD','Hour1','Hour2','Hour3','Hour4','Hour5','Hour6','Hour7','Hour8',...
               				     'Hour9','Hour10','Hour11','Hour12','Hour13','Hour14','Hour15','Hour16','Hour17','Hour18','Hour19','Hour20','Hour21','Hour22','Hour23','Hour24'});


rm_panel2_LP = fitrm(t_panel2,rm_model_string_BLorSD_LP,'WithinDesign',WithinStructure1LP,'WithinModel','Time');
RManovatbl_WTBL_vs_WTSD_LP_N = ranova(rm_panel2_LP,'WithinModel','Time');
% -- Post-hocs --
if RManovatbl_WTBL_vs_WTSD_LP_N{'BLorSD','pValue'} < 0.05 | RManovatbl_WTBL_vs_WTSD_LP_N{'BLorSD:Time','pValue'} < 0.05
	sig_effect.WTBL_vs_WTSD_LP_N = 1;  % set boolean for significant effect 
	mc_table.panel2_LP 		  = multcompare(rm_panel2_LP,'BLorSD','By','Time','ComparisonType','lsd');
	mc_table.panel2_LP.pValue = mafdr(mc_table.panel2_LP.pValue,'BHFDR',true);
	[~, ia] = unique(mc_table.panel2_LP.Time(:), 'stable');  % remove duplicate rows
	mc_table.panel2_LP = mc_table.panel2_LP(ia,:);
	xlocs_sig_effects.WTBL_vs_WTSD_LP_N = double(string(mc_table.panel2_LP.Time(mc_table.panel2_LP.pValue<0.05)));
end
% --------------------------------------------------------------------------------------------------------------


% --- Panel 2: RM Anova WT BL vs WT SD (DP NREM) ---------------------------------------------------------------
rm_panel2_DP = fitrm(t_panel2,rm_model_string_BLorSD_DP,'WithinDesign',WithinStructure1DP,'WithinModel','Time');
RManovatbl_WTBL_vs_WTSD_DP_N = ranova(rm_panel2_DP,'WithinModel','Time');
% -- Post-hocs --
if RManovatbl_WTBL_vs_WTSD_DP_N{'BLorSD','pValue'} < 0.05 | RManovatbl_WTBL_vs_WTSD_DP_N{'BLorSD:Time','pValue'} < 0.05
	sig_effect.WTBL_vs_WTSD_DP_N = 1;  % set boolean for significant effect 
	mc_table.panel2_DP 		  = multcompare(rm_panel2_DP,'BLorSD','By','Time','ComparisonType','lsd');
	mc_table.panel2_DP.pValue = mafdr(mc_table.panel2_DP.pValue,'BHFDR',true);
	[~, ia] = unique(mc_table.panel2_DP.Time(:), 'stable');  % remove duplicate rows
	mc_table.panel2_DP = mc_table.panel2_DP(ia,:);
	xlocs_sig_effects.WTBL_vs_WTSD_DP_N = double(string(mc_table.panel2_DP.Time(mc_table.panel2_DP.pValue<0.05)));
end
% --------------------------------------------------------------------------------------------------------------


% --- Panel 3: RM Anova WT BL vs WT SD (LP exc SD REM) ---------------------------------------------------------
REMPercByHour_WT_BL = TimeInState.WT.(MorF).BL.REM/60;  % Divide by 60 to get percentages, rather than minutes  
REMPercByHour_WT_SD = TimeInState.WT.(MorF).SD.REM/60;

YR1 = [REMPercByHour_WT_BL'; REMPercByHour_WT_SD'];
t_panel3 = table(BLorSD,YR1(:,1),YR1(:,2),YR1(:,3),YR1(:,4),YR1(:,5),YR1(:,6),YR1(:,7),YR1(:,8),YR1(:,9),YR1(:,10),YR1(:,11),YR1(:,12), ...
           	   				   YR1(:,13),YR1(:,14),YR1(:,15),YR1(:,16),YR1(:,17),YR1(:,18),YR1(:,19),YR1(:,20),YR1(:,21),YR1(:,22),YR1(:,23),YR1(:,24), ...
           	   				     'VariableNames',{'BLorSD','Hour1','Hour2','Hour3','Hour4','Hour5','Hour6','Hour7','Hour8',...
               				     'Hour9','Hour10','Hour11','Hour12','Hour13','Hour14','Hour15','Hour16','Hour17','Hour18','Hour19','Hour20','Hour21','Hour22','Hour23','Hour24'});


rm_panel3_LP = fitrm(t_panel3,rm_model_string_BLorSD_LP,'WithinDesign',WithinStructure1LP,'WithinModel','Time');
RManovatbl_WTBL_vs_WTSD_LP_R = ranova(rm_panel3_LP,'WithinModel','Time');
% -- Post-hocs --
if RManovatbl_WTBL_vs_WTSD_LP_R{'BLorSD','pValue'} < 0.05 | RManovatbl_WTBL_vs_WTSD_LP_R{'BLorSD:Time','pValue'} < 0.05
	sig_effect.WTBL_vs_WTSD_LP_R = 1;  % set boolean for significant effect 
	mc_table.panel3_LP 		  = multcompare(rm_panel3_LP,'BLorSD','By','Time','ComparisonType','lsd');
	mc_table.panel3_LP.pValue = mafdr(mc_table.panel3_LP.pValue,'BHFDR',true);
	[~, ia] = unique(mc_table.panel3_LP.Time(:), 'stable');  % remove duplicate rows
	mc_table.panel3_LP = mc_table.panel3_LP(ia,:);
	xlocs_sig_effects.WTBL_vs_WTSD_LP_R = double(string(mc_table.panel3_LP.Time(mc_table.panel3_LP.pValue<0.05)));
end
% --------------------------------------------------------------------------------------------------------------


% --- Panel 3: RM Anova WT BL vs WT SD (DP REM) ----------------------------------------------------------------
rm_panel3_DP = fitrm(t_panel3,rm_model_string_BLorSD_DP,'WithinDesign',WithinStructure1DP,'WithinModel','Time');
RManovatbl_WTBL_vs_WTSD_DP_R = ranova(rm_panel3_DP,'WithinModel','Time');
% -- Post-hocs --
if RManovatbl_WTBL_vs_WTSD_DP_R{'BLorSD','pValue'} < 0.05 | RManovatbl_WTBL_vs_WTSD_DP_R{'BLorSD:Time','pValue'} < 0.05
	sig_effect.WTBL_vs_WTSD_DP_R = 1;  % set boolean for significant effect 
	mc_table.panel3_DP 		  = multcompare(rm_panel3_DP,'BLorSD','By','Time','ComparisonType','lsd');
	mc_table.panel3_DP.pValue = mafdr(mc_table.panel3_DP.pValue,'BHFDR',true);
	[~, ia] = unique(mc_table.panel3_DP.Time(:), 'stable');  % remove duplicate rows
	mc_table.panel3_DP = mc_table.panel3_DP(ia,:);
	xlocs_sig_effects.WTBL_vs_WTSD_DP_R = double(string(mc_table.panel3_DP.Time(mc_table.panel3_DP.pValue<0.05)));
end
% --------------------------------------------------------------------------------------------------------------


% --- Panel 4: RM Anova Mut BL vs Mut SD (LP exc SD) Wake --------------------------------------------------------
WakePercByHour_Mut_BL = TimeInState.Mut.(MorF).BL.Wake/60;  % Divide by 60 to get percentages, rather than minutes  
WakePercByHour_Mut_SD = TimeInState.Mut.(MorF).SD.Wake/60;
BLorSD = categorical([zeros(1,size(WakePercByHour_Mut_BL,2)) ones(1,size(WakePercByHour_Mut_SD,2))]');  % Categorical is important!  otherwise ranova produces incorrect F values


rm_model_string_BLorSD_LP = ['Hour',num2str(Time1LP(1)),'-Hour',num2str(Time1LP(end)),' ~ BLorSD']; 


YW4 = [WakePercByHour_Mut_BL'; WakePercByHour_Mut_SD'];
t_panel4 = table(BLorSD,YW4(:,1),YW4(:,2),YW4(:,3),YW4(:,4),YW4(:,5),YW4(:,6),YW4(:,7),YW4(:,8),YW4(:,9),YW4(:,10),YW4(:,11),YW4(:,12), ...
           	   				   YW4(:,13),YW4(:,14),YW4(:,15),YW4(:,16),YW4(:,17),YW4(:,18),YW4(:,19),YW4(:,20),YW4(:,21),YW4(:,22),YW4(:,23),YW4(:,24), ...
           	   				     'VariableNames',{'BLorSD','Hour1','Hour2','Hour3','Hour4','Hour5','Hour6','Hour7','Hour8',...
               				     'Hour9','Hour10','Hour11','Hour12','Hour13','Hour14','Hour15','Hour16','Hour17','Hour18','Hour19','Hour20','Hour21','Hour22','Hour23','Hour24'});


rm_panel4_LP = fitrm(t_panel4,rm_model_string_BLorSD_LP,'WithinDesign',WithinStructure1LP,'WithinModel','Time');
RManovatbl_MutBL_vs_MutSD_LP_W = ranova(rm_panel4_LP,'WithinModel','Time');
% -- Post-hocs --
if RManovatbl_MutBL_vs_MutSD_LP_W{'BLorSD','pValue'} < 0.05 | RManovatbl_MutBL_vs_MutSD_LP_W{'BLorSD:Time','pValue'} < 0.05
	sig_effect.MutBL_vs_MutSD_LP_W = 1;  % set boolean for significant effect 
	mc_table.panel4_LP 		  = multcompare(rm_panel4_LP,'BLorSD','By','Time','ComparisonType','lsd');
	mc_table.panel4_LP.pValue = mafdr(mc_table.panel4_LP.pValue,'BHFDR',true);
	[~, ia] = unique(mc_table.panel4_LP.Time(:), 'stable');  % remove duplicate rows
	mc_table.panel4_LP = mc_table.panel4_LP(ia,:);
	xlocs_sig_effects.MutBL_vs_MutSD_LP_W = double(string(mc_table.panel4_LP.Time(mc_table.panel4_LP.pValue<0.05)));
end
% --------------------------------------------------------------------------------------------------------------


% --- Panel 4: RM Anova Mut BL vs Mut SD (DP Wake) -------------------------------------------------------------
rm_model_string_BLorSD_DP = ['Hour',num2str(Time1DP(1)),'-Hour',num2str(Time1DP(end)),' ~ BLorSD']; 

rm_panel4_DP = fitrm(t_panel4,rm_model_string_BLorSD_DP,'WithinDesign',WithinStructure1DP,'WithinModel','Time');
RManovatbl_MutBL_vs_MutSD_DP_W = ranova(rm_panel4_DP,'WithinModel','Time');
% -- Post-hocs --
if RManovatbl_MutBL_vs_MutSD_DP_W{'BLorSD','pValue'} < 0.05 | RManovatbl_MutBL_vs_MutSD_DP_W{'BLorSD:Time','pValue'} < 0.05
	sig_effect.MutBL_vs_MutSD_DP_W = 1;  % set boolean for significant effect 
	mc_table.panel4_DP 		  = multcompare(rm_panel4_DP,'BLorSD','By','Time','ComparisonType','lsd');
	mc_table.panel4_DP.pValue = mafdr(mc_table.panel4_DP.pValue,'BHFDR',true);
	[~, ia] = unique(mc_table.panel4_DP.Time(:), 'stable');  % remove duplicate rows
	mc_table.panel4_DP = mc_table.panel4_DP(ia,:);
	xlocs_sig_effects.MutBL_vs_MutSD_DP_W = double(string(mc_table.panel4_DP.Time(mc_table.panel4_DP.pValue<0.05)));
end
% --------------------------------------------------------------------------------------------------------------


% --- Panel 5: RM Anova Mut BL vs Mut SD (LP exc SD NREM) --------------------------------------------------------
NREMPercByHour_Mut_BL = TimeInState.Mut.(MorF).BL.NREM/60;  % Divide by 60 to get percentages, rather than minutes  
NREMPercByHour_Mut_SD = TimeInState.Mut.(MorF).SD.NREM/60;

YN2 = [NREMPercByHour_Mut_BL'; NREMPercByHour_Mut_SD'];
t_panel5 = table(BLorSD,YN2(:,1),YN2(:,2),YN2(:,3),YN2(:,4),YN2(:,5),YN2(:,6),YN2(:,7),YN2(:,8),YN2(:,9),YN2(:,10),YN2(:,11),YN2(:,12), ...
           	   				   YN2(:,13),YN2(:,14),YN2(:,15),YN2(:,16),YN2(:,17),YN2(:,18),YN2(:,19),YN2(:,20),YN2(:,21),YN2(:,22),YN2(:,23),YN2(:,24), ...
           	   				     'VariableNames',{'BLorSD','Hour1','Hour2','Hour3','Hour4','Hour5','Hour6','Hour7','Hour8',...
               				     'Hour9','Hour10','Hour11','Hour12','Hour13','Hour14','Hour15','Hour16','Hour17','Hour18','Hour19','Hour20','Hour21','Hour22','Hour23','Hour24'});


rm_panel5_LP = fitrm(t_panel5,rm_model_string_BLorSD_LP,'WithinDesign',WithinStructure1LP,'WithinModel','Time');
RManovatbl_MutBL_vs_MutSD_LP_N = ranova(rm_panel5_LP,'WithinModel','Time');
% -- Post-hocs --
if RManovatbl_MutBL_vs_MutSD_LP_N{'BLorSD','pValue'} < 0.05 | RManovatbl_MutBL_vs_MutSD_LP_N{'BLorSD:Time','pValue'} < 0.05
	sig_effect.MutBL_vs_MutSD_LP_N = 1;  % set boolean for significant effect 
	mc_table.panel5_LP 		  = multcompare(rm_panel5_LP,'BLorSD','By','Time','ComparisonType','lsd');
	mc_table.panel5_LP.pValue = mafdr(mc_table.panel5_LP.pValue,'BHFDR',true);
	[~, ia] = unique(mc_table.panel5_LP.Time(:), 'stable');  % remove duplicate rows
	mc_table.panel5_LP = mc_table.panel5_LP(ia,:);
	xlocs_sig_effects.MutBL_vs_MutSD_LP_N = double(string(mc_table.panel5_LP.Time(mc_table.panel5_LP.pValue<0.05)));
end
% --------------------------------------------------------------------------------------------------------------


% --- Panel 5: RM Anova Mut BL vs Mut SD (DP NREM) -------------------------------------------------------------
rm_panel5_DP = fitrm(t_panel5,rm_model_string_BLorSD_DP,'WithinDesign',WithinStructure1DP,'WithinModel','Time');
RManovatbl_MutBL_vs_MutSD_DP_N = ranova(rm_panel5_DP,'WithinModel','Time');
% -- Post-hocs --
if RManovatbl_MutBL_vs_MutSD_DP_N{'BLorSD','pValue'} < 0.05 | RManovatbl_MutBL_vs_MutSD_DP_N{'BLorSD:Time','pValue'} < 0.05
	sig_effect.MutBL_vs_MutSD_DP_N = 1;  % set boolean for significant effect 
	mc_table.panel5_DP 		  = multcompare(rm_panel5_DP,'BLorSD','By','Time','ComparisonType','lsd');
	mc_table.panel5_DP.pValue = mafdr(mc_table.panel5_DP.pValue,'BHFDR',true);
	[~, ia] = unique(mc_table.panel5_DP.Time(:), 'stable');  % remove duplicate rows
	mc_table.panel5_DP = mc_table.panel5_DP(ia,:);
	xlocs_sig_effects.MutBL_vs_MutSD_DP_N = double(string(mc_table.panel5_DP.Time(mc_table.panel5_DP.pValue<0.05)));
end
% --------------------------------------------------------------------------------------------------------------


% --- Panel 6: RM Anova Mut BL vs Mut SD (LP exc SD REM) -------------------------------------------------------
REMPercByHour_Mut_BL = TimeInState.Mut.(MorF).BL.REM/60;  % Divide by 60 to get percentages, rather than minutes  
REMPercByHour_Mut_SD = TimeInState.Mut.(MorF).SD.REM/60;

YR2 = [REMPercByHour_Mut_BL'; REMPercByHour_Mut_SD'];
t_panel6 = table(BLorSD,YR2(:,1),YR2(:,2),YR2(:,3),YR2(:,4),YR2(:,5),YR2(:,6),YR2(:,7),YR2(:,8),YR2(:,9),YR2(:,10),YR2(:,11),YR2(:,12), ...
           	   				   YR2(:,13),YR2(:,14),YR2(:,15),YR2(:,16),YR2(:,17),YR2(:,18),YR2(:,19),YR2(:,20),YR2(:,21),YR2(:,22),YR2(:,23),YR2(:,24), ...
           	   				     'VariableNames',{'BLorSD','Hour1','Hour2','Hour3','Hour4','Hour5','Hour6','Hour7','Hour8',...
               				     'Hour9','Hour10','Hour11','Hour12','Hour13','Hour14','Hour15','Hour16','Hour17','Hour18','Hour19','Hour20','Hour21','Hour22','Hour23','Hour24'});


rm_panel6_LP = fitrm(t_panel6,rm_model_string_BLorSD_LP,'WithinDesign',WithinStructure1LP,'WithinModel','Time');
RManovatbl_MutBL_vs_MutSD_LP_R = ranova(rm_panel6_LP,'WithinModel','Time');
% -- Post-hocs --
if RManovatbl_MutBL_vs_MutSD_LP_R{'BLorSD','pValue'} < 0.05 | RManovatbl_MutBL_vs_MutSD_LP_R{'BLorSD:Time','pValue'} < 0.05
	sig_effect.MutBL_vs_MutSD_LP_R = 1;  % set boolean for significant effect 
	mc_table.panel6_LP 		  = multcompare(rm_panel6_LP,'BLorSD','By','Time','ComparisonType','lsd');
	mc_table.panel6_LP.pValue = mafdr(mc_table.panel6_LP.pValue,'BHFDR',true);
	[~, ia] = unique(mc_table.panel6_LP.Time(:), 'stable');  % remove duplicate rows
	mc_table.panel6_LP = mc_table.panel6_LP(ia,:);
	xlocs_sig_effects.MutBL_vs_MutSD_LP_R = double(string(mc_table.panel6_LP.Time(mc_table.panel6_LP.pValue<0.05)));
end

% --------------------------------------------------------------------------------------------------------------


% --- Panel 6: RM Anova Mut BL vs Mut SD (DP REM) --------------------------------------------------------------
rm_panel6_DP = fitrm(t_panel6,rm_model_string_BLorSD_DP,'WithinDesign',WithinStructure1DP,'WithinModel','Time');
RManovatbl_MutBL_vs_MutSD_DP_R = ranova(rm_panel6_DP,'WithinModel','Time');
% -- Post-hocs --
if RManovatbl_MutBL_vs_MutSD_DP_R{'BLorSD','pValue'} < 0.05 | RManovatbl_MutBL_vs_MutSD_DP_R{'BLorSD:Time','pValue'} < 0.05
	sig_effect.MutBL_vs_MutSD_DP_R = 1;  % set boolean for significant effect 
	mc_table.panel6_DP 		  = multcompare(rm_panel6_DP,'BLorSD','By','Time','ComparisonType','lsd');
	mc_table.panel6_DP.pValue = mafdr(mc_table.panel6_DP.pValue,'BHFDR',true);
	[~, ia] = unique(mc_table.panel6_DP.Time(:), 'stable');  % remove duplicate rows
	mc_table.panel6_DP = mc_table.panel6_DP(ia,:);
	xlocs_sig_effects.MutBL_vs_MutSD_DP_R = double(string(mc_table.panel6_DP.Time(mc_table.panel6_DP.pValue<0.05)));
end
% --------------------------------------------------------------------------------------------------------------


% --- Panel 7: RM Anova WT SD vs Mut SD (LP exc SD Wake) --------------------------------------------------------
WakePercByHour_WT_SD  = TimeInState.WT.(MorF).SD.Wake/60;  % Divide by 60 to get percentages, rather than minutes  
WakePercByHour_Mut_SD = TimeInState.Mut.(MorF).SD.Wake/60;
Genotype = categorical([zeros(1,size(WakePercByHour_WT_SD,2)) ones(1,size(WakePercByHour_Mut_SD,2))]');  % Categorical is important!  otherwise ranova produces incorrect F values


rm_model_string_WTorMut_LP = ['Hour',num2str(Time1LP(1)),'-Hour',num2str(Time1LP(end)),' ~ Genotype']; 


YW7 = [WakePercByHour_WT_SD'; WakePercByHour_Mut_SD'];
t_panel7 = table(Genotype,YW7(:,1),YW7(:,2),YW7(:,3),YW7(:,4),YW7(:,5),YW7(:,6),YW7(:,7),YW7(:,8),YW7(:,9),YW7(:,10),YW7(:,11),YW7(:,12), ...
           	   				   YW7(:,13),YW7(:,14),YW7(:,15),YW7(:,16),YW7(:,17),YW7(:,18),YW7(:,19),YW7(:,20),YW7(:,21),YW7(:,22),YW7(:,23),YW7(:,24), ...
           	   				     'VariableNames',{'Genotype','Hour1','Hour2','Hour3','Hour4','Hour5','Hour6','Hour7','Hour8',...
               				     'Hour9','Hour10','Hour11','Hour12','Hour13','Hour14','Hour15','Hour16','Hour17','Hour18','Hour19','Hour20','Hour21','Hour22','Hour23','Hour24'});


rm_panel7_LP = fitrm(t_panel7,rm_model_string_WTorMut_LP,'WithinDesign',WithinStructure1LP,'WithinModel','Time');
RManovatbl_WTSD_vs_MutSD_LP_W = ranova(rm_panel7_LP,'WithinModel','Time');
% -- Post-hocs --
if RManovatbl_WTSD_vs_MutSD_LP_W{'Genotype','pValue'} < 0.05 | RManovatbl_WTSD_vs_MutSD_LP_W{'Genotype:Time','pValue'} < 0.05
	sig_effect.WTSD_vs_MutSD_LP_W = 1;  % set boolean for significant effect 
	mc_table.panel7_LP 		  = multcompare(rm_panel7_LP,'Genotype','By','Time','ComparisonType','lsd');
	mc_table.panel7_LP.pValue = mafdr(mc_table.panel7_LP.pValue,'BHFDR',true);
	[~, ia] = unique(mc_table.panel7_LP.Time(:), 'stable');  % remove duplicate rows
	mc_table.panel7_LP = mc_table.panel7_LP(ia,:);
	xlocs_sig_effects.WTSD_vs_MutSD_LP_W = double(string(mc_table.panel7_LP.Time(mc_table.panel7_LP.pValue<0.05)));
end
% ---------------------------------------------------------------------------------------------------------------


% --- Panel 7: RM Anova WT SD vs Mut SD (DP Wake) ---------------------------------------------------------------
rm_model_string_WTorMut_DP = ['Hour',num2str(Time1DP(1)),'-Hour',num2str(Time1DP(end)),' ~ Genotype']; 

rm_panel7_DP = fitrm(t_panel7,rm_model_string_WTorMut_DP,'WithinDesign',WithinStructure1DP,'WithinModel','Time');
RManovatbl_WTSD_vs_MutSD_DP_W = ranova(rm_panel7_DP,'WithinModel','Time');
% -- Post-hocs --
if RManovatbl_WTSD_vs_MutSD_DP_W{'Genotype','pValue'} < 0.05 | RManovatbl_WTSD_vs_MutSD_DP_W{'Genotype:Time','pValue'} < 0.05
	sig_effect.WTSD_vs_MutSD_DP_W = 1;  % set boolean for significant effect 
	mc_table.panel7_DP 		  = multcompare(rm_panel7_DP,'Genotype','By','Time','ComparisonType','lsd');
	mc_table.panel7_DP.pValue = mafdr(mc_table.panel7_DP.pValue,'BHFDR',true);
	[~, ia] = unique(mc_table.panel7_DP.Time(:), 'stable');  % remove duplicate rows
	mc_table.panel7_DP = mc_table.panel7_DP(ia,:);
	xlocs_sig_effects.WTSD_vs_MutSD_DP_W = double(string(mc_table.panel7_DP.Time(mc_table.panel7_DP.pValue<0.05)));
end
% ---------------------------------------------------------------------------------------------------------------


% --- Panel 8: RM Anova WT SD vs Mut SD (LP exc SD NREM) --------------------------------------------------------
NREMPercByHour_WT_SD  = TimeInState.WT.(MorF).SD.NREM/60;  % Divide by 60 to get percentages, rather than minutes  
NREMPercByHour_Mut_SD = TimeInState.Mut.(MorF).SD.NREM/60;


YN8 = [NREMPercByHour_WT_SD'; NREMPercByHour_Mut_SD'];
t_panel8 = table(Genotype,YN8(:,1),YN8(:,2),YN8(:,3),YN8(:,4),YN8(:,5),YN8(:,6),YN8(:,7),YN8(:,8),YN8(:,9),YN8(:,10),YN8(:,11),YN8(:,12), ...
           	   				   YN8(:,13),YN8(:,14),YN8(:,15),YN8(:,16),YN8(:,17),YN8(:,18),YN8(:,19),YN8(:,20),YN8(:,21),YN8(:,22),YN8(:,23),YN8(:,24), ...
           	   				     'VariableNames',{'Genotype','Hour1','Hour2','Hour3','Hour4','Hour5','Hour6','Hour7','Hour8',...
               				     'Hour9','Hour10','Hour11','Hour12','Hour13','Hour14','Hour15','Hour16','Hour17','Hour18','Hour19','Hour20','Hour21','Hour22','Hour23','Hour24'});


rm_panel8_LP = fitrm(t_panel8,rm_model_string_WTorMut_LP,'WithinDesign',WithinStructure1LP,'WithinModel','Time');
RManovatbl_WTSD_vs_MutSD_LP_N = ranova(rm_panel8_LP,'WithinModel','Time');
% -- Post-hocs --
if RManovatbl_WTSD_vs_MutSD_LP_N{'Genotype','pValue'} < 0.05 | RManovatbl_WTSD_vs_MutSD_LP_N{'Genotype:Time','pValue'} < 0.05
	sig_effect.WTSD_vs_MutSD_LP_N = 1;  % set boolean for significant effect 
	mc_table.panel8_LP 		  = multcompare(rm_panel8_LP,'Genotype','By','Time','ComparisonType','lsd');
	mc_table.panel8_LP.pValue = mafdr(mc_table.panel8_LP.pValue,'BHFDR',true);
	[~, ia] = unique(mc_table.panel8_LP.Time(:), 'stable');  % remove duplicate rows
	mc_table.panel8_LP = mc_table.panel8_LP(ia,:);
	xlocs_sig_effects.WTSD_vs_MutSD_LP_N = double(string(mc_table.panel8_LP.Time(mc_table.panel8_LP.pValue<0.05)));
end
% ---------------------------------------------------------------------------------------------------------------


% --- Panel 8: RM Anova WT SD vs Mut SD (DP NREM) ---------------------------------------------------------------
rm_panel8_DP = fitrm(t_panel8,rm_model_string_WTorMut_DP,'WithinDesign',WithinStructure1DP,'WithinModel','Time');
RManovatbl_WTSD_vs_MutSD_DP_N = ranova(rm_panel8_DP,'WithinModel','Time');
% -- Post-hocs --
if RManovatbl_WTSD_vs_MutSD_DP_N{'Genotype','pValue'} < 0.05 | RManovatbl_WTSD_vs_MutSD_DP_N{'Genotype:Time','pValue'} < 0.05
	sig_effect.WTSD_vs_MutSD_DP_N = 1;  % set boolean for significant effect 
	mc_table.panel8_DP 		  = multcompare(rm_panel8_DP,'Genotype','By','Time','ComparisonType','lsd');
	mc_table.panel8_DP.pValue = mafdr(mc_table.panel8_DP.pValue,'BHFDR',true);
	[~, ia] = unique(mc_table.panel8_DP.Time(:), 'stable');  % remove duplicate rows
	mc_table.panel8_DP = mc_table.panel8_DP(ia,:);
	xlocs_sig_effects.WTSD_vs_MutSD_DP_N = double(string(mc_table.panel8_DP.Time(mc_table.panel8_DP.pValue<0.05)));
end
% ---------------------------------------------------------------------------------------------------------------


% --- Panel 9: RM Anova WT SD vs Mut SD (LP exc SD REM) ---------------------------------------------------------
REMPercByHour_WT_SD  = TimeInState.WT.(MorF).SD.REM/60;  % Divide by 60 to get percentages, rather than minutes  
REMPercByHour_Mut_SD = TimeInState.Mut.(MorF).SD.REM/60;


YR9 = [REMPercByHour_WT_SD'; REMPercByHour_Mut_SD'];
t_panel9 = table(Genotype,YR9(:,1),YR9(:,2),YR9(:,3),YR9(:,4),YR9(:,5),YR9(:,6),YR9(:,7),YR9(:,8),YR9(:,9),YR9(:,10),YR9(:,11),YR9(:,12), ...
           	   				   YR9(:,13),YR9(:,14),YR9(:,15),YR9(:,16),YR9(:,17),YR9(:,18),YR9(:,19),YR9(:,20),YR9(:,21),YR9(:,22),YR9(:,23),YR9(:,24), ...
           	   				     'VariableNames',{'Genotype','Hour1','Hour2','Hour3','Hour4','Hour5','Hour6','Hour7','Hour8',...
               				     'Hour9','Hour10','Hour11','Hour12','Hour13','Hour14','Hour15','Hour16','Hour17','Hour18','Hour19','Hour20','Hour21','Hour22','Hour23','Hour24'});


rm_panel9_LP = fitrm(t_panel9,rm_model_string_WTorMut_LP,'WithinDesign',WithinStructure1LP,'WithinModel','Time');
RManovatbl_WTSD_vs_MutSD_LP_R = ranova(rm_panel9_LP,'WithinModel','Time');
% -- Post-hocs --
if RManovatbl_WTSD_vs_MutSD_LP_R{'Genotype','pValue'} < 0.05 | RManovatbl_WTSD_vs_MutSD_LP_R{'Genotype:Time','pValue'} < 0.05
	sig_effect.WTSD_vs_MutSD_LP_R = 1;  % set boolean for significant effect 
	mc_table.panel9_LP 		  = multcompare(rm_panel9_LP,'Genotype','By','Time','ComparisonType','lsd');
	mc_table.panel9_LP.pValue = mafdr(mc_table.panel9_LP.pValue,'BHFDR',true);
	[~, ia] = unique(mc_table.panel9_LP.Time(:), 'stable');  % remove duplicate rows
	mc_table.panel9_LP = mc_table.panel9_LP(ia,:);
	xlocs_sig_effects.WTSD_vs_MutSD_LP_R = double(string(mc_table.panel9_LP.Time(mc_table.panel9_LP.pValue<0.05)));
end
% ---------------------------------------------------------------------------------------------------------------


% --- Panel 9: RM Anova WT SD vs Mut SD (DP REM) ----------------------------------------------------------------
rm_panel9_DP = fitrm(t_panel9,rm_model_string_WTorMut_DP,'WithinDesign',WithinStructure1DP,'WithinModel','Time');
RManovatbl_WTSD_vs_MutSD_DP_R = ranova(rm_panel9_DP,'WithinModel','Time');
% -- Post-hocs --
if RManovatbl_WTSD_vs_MutSD_DP_R{'Genotype','pValue'} < 0.05 | RManovatbl_WTSD_vs_MutSD_DP_R{'Genotype:Time','pValue'} < 0.05
	sig_effect.WTSD_vs_MutSD_DP_R = 1;  % set boolean for significant effect 
	mc_table.panel9_DP 		  = multcompare(rm_panel9_DP,'Genotype','By','Time','ComparisonType','lsd');
	mc_table.panel9_DP.pValue = mafdr(mc_table.panel9_DP.pValue,'BHFDR',true);
	[~, ia] = unique(mc_table.panel9_DP.Time(:), 'stable');  % remove duplicate rows
	mc_table.panel9_DP = mc_table.panel9_DP(ia,:);
	xlocs_sig_effects.WTSD_vs_MutSD_DP_R = double(string(mc_table.panel9_DP.Time(mc_table.panel9_DP.pValue<0.05)));
end
% ---------------------------------------------------------------------------------------------------------------


% --- End Stats ------------------------------------------------------------------------------------
% --------------------------------------------------------------------------------------------------









% --------------------------------------------------------------------------------------------------
% --- FIGURE CODE ----------------------------------------------------------------------------------
% Set the overall figure size and shape  
H9=figure;
H9.Name = 'Nine-Panel Time in State Figure';
H9.Renderer = 'painters';
H9.Position = [470 270 712 513];

t= tiledlayout(3,3,'TileSpacing','compact','Padding','compact','TileIndexing', 'columnmajor');


% --- Position (1,1): W    WT BL (black solid) vs WT SD (black dashed) -----------------------
%subplot(3,3,1)
nexttile 
ax = H9.Children.Children;
rectangle('Position',LPboxpos,'FaceColor','w','EdgeColor','k');  	% LP rectangle
hold on
rectangle('Position',DPboxpos,'FaceColor','k','EdgeColor','none');    % DP rectangle
p1 = patch([0 0 SD_length SD_length],[-1 100 100 -1],[1 1 1],'LineStyle','none'); 	% SD hatching
%p1.FaceColor = [1 1 1];  % white face color
hatchfill2(p1,'single','HatchAngle',45,'HatchDensity',35);
p2 = patch([0 0 SD_length SD_length],[-1 4 4 -1],[1 1 1]); 	% SD hatching
p2.Parent = ax;
hatchfill2(p2,'single','HatchAngle',45,'HatchDensity',35,'HatchLineWidth',1.5);
e1=errorbar([1:24],HourlyStatePercentages.Means.WT.BL.All24hrs.Wake,HourlyStatePercentages.SEMs.WT.BL.All24hrs.Wake,'Color','k','LineWidth',mylinewidth,'CapSize',0);
e2=errorbar([SD_length+1:24],HourlyStatePercentages.Means.WT.SD.All24hrexcSD.Wake,HourlyStatePercentages.SEMs.WT.SD.All24hrexcSD.Wake,'Color','k','LineWidth',mylinewidth,'CapSize',0,'LineStyle','--');
e1.Marker = '.';
e1.MarkerSize = EBMarkersize;
e2.Marker = '.';
e2.MarkerSize = EBMarkersize;
hold off 
set(ax, 'Layer','top')
title('')
ylabel({'{\bf Wakefulness}';'(% Recording Time)'});
xlabel('')
l = legend([e1 e2],'Wild Type BL','Wild Type SD');
l.Box = 'off';
l.Location = 'southeast';
l.FontSize = 8;
xlim([1 24])
ylim([-1 100])
ax.FontSize  = myfontsize;
ax.LineWidth = myaxislinewidth;
ax.XTick = 0:6:24;
ax.XTickLabel = {'' '' '' '' ''};
ax.YTick = 0:20:100;
ax.YTickLabel = {'0' '20' '40' '60' '80' '100'};
ax.TickDir = 'out';
% -- Add stars if sig diffs (LP)
if sig_effect.WTBL_vs_WTSD_LP_W
	hold on 
	
	rectangle('Position',[0 ax.YLim(2)-5 12 5],'FaceColor',Sig12hrColor,'EdgeColor','none'); 
	bool_xlocsWTBL_LP_sig_values = ismember(e1.XData,xlocs_sig_effects.WTBL_vs_WTSD_LP_W);
	bool_xlocsWTSD_LP_sig_values = ismember(e2.XData,xlocs_sig_effects.WTBL_vs_WTSD_LP_W);
	largest_y1_vals = e1.YData(bool_xlocsWTBL_LP_sig_values) + e1.YPositiveDelta(bool_xlocsWTBL_LP_sig_values);
	largest_y2_vals = e2.YData(bool_xlocsWTSD_LP_sig_values) + e2.YPositiveDelta(bool_xlocsWTSD_LP_sig_values);
	ylocs_sig_effects = max([largest_y1_vals; largest_y2_vals]) + 10*ones(1,length(largest_y1_vals));
	ylocs_sig_effects(ylocs_sig_effects>ax.YLim(2))=ax.YLim(2);
	text(xlocs_sig_effects.WTBL_vs_WTSD_LP_W, ylocs_sig_effects, '*', 'FontSize',18,'HorizontalAlignment','center','VerticalAlignment','middle');

	hold off 
end 
% -- Add stars if sig diffs (DP)
if sig_effect.WTBL_vs_WTSD_DP_W
	hold on 
	
	rectangle('Position',[12 ax.YLim(2)-5 12 5],'FaceColor',Sig12hrColor,'EdgeColor','none'); 
	bool_xlocsWTBL_DP_sig_values = ismember(e1.XData,xlocs_sig_effects.WTBL_vs_WTSD_DP_W);
	bool_xlocsWTSD_DP_sig_values = ismember(e2.XData,xlocs_sig_effects.WTBL_vs_WTSD_DP_W);
	largest_y1_vals = e1.YData(bool_xlocsWTBL_DP_sig_values) + e1.YPositiveDelta(bool_xlocsWTBL_DP_sig_values);
	largest_y2_vals = e2.YData(bool_xlocsWTSD_DP_sig_values) + e2.YPositiveDelta(bool_xlocsWTSD_DP_sig_values);
	ylocs_sig_effects = max([largest_y1_vals; largest_y2_vals]) + 10*ones(1,length(largest_y1_vals));
	ylocs_sig_effects(ylocs_sig_effects>ax.YLim(2))=ax.YLim(2);
	text(xlocs_sig_effects.WTBL_vs_WTSD_DP_W, ylocs_sig_effects, '*', 'FontSize',18,'HorizontalAlignment','center','VerticalAlignment','middle');

	hold off 
end 
% --- End Panel (1,1) -----------------------------------------------------------------------

% Position (2,1): NREM WT BL (black solid) vs WT SD (black dashed)
%subplot(3,3,4)
nexttile
ax=H9.Children.Children(1);
rectangle('Position',LPboxpos,'FaceColor','w','EdgeColor','k');  	% LP rectangle
hold on
rectangle('Position',DPboxpos,'FaceColor','k','EdgeColor','none');    % DP rectangle
p1 = patch([0 0 SD_length SD_length],[-1 100 100 -1],[1 1 1],'LineStyle','none'); 	% SD hatching
hatchfill2(p1,'single','HatchAngle',45,'HatchDensity',35);
p2 = patch([0 0 SD_length SD_length],[-1 4 4 -1],[1 1 1]); 	% SD hatching
p2.Parent = ax;
hatchfill2(p2,'single','HatchAngle',45,'HatchDensity',35,'HatchLineWidth',1.5);
e1=errorbar([1:24],HourlyStatePercentages.Means.WT.BL.All24hrs.NREM,HourlyStatePercentages.SEMs.WT.BL.All24hrs.NREM,'Color','k','LineWidth',mylinewidth,'CapSize',0);
e2=errorbar([SD_length+1:24],HourlyStatePercentages.Means.WT.SD.All24hrexcSD.NREM,HourlyStatePercentages.SEMs.WT.SD.All24hrexcSD.NREM,'Color','k','LineWidth',mylinewidth,'CapSize',0,'LineStyle','--');
e1.Marker = '.';
e1.MarkerSize = EBMarkersize;
e2.Marker = '.';
e2.MarkerSize = EBMarkersize;
hold off 
set(ax, 'Layer','top')
title('')
ylabel({'{\bf NREM}';'(% Recording Time)'});
xlabel('')
l = legend([e1 e2],'Wild Type BL','Wild Type SD');
l.Box = 'off';
l.Location = 'southeast';
l.Visible = 'off';
xlim([1 24])
ylim([-1 100])

ax.FontSize  = myfontsize;
ax.LineWidth = myaxislinewidth;
ax.XTick = 0:6:24;
ax.XTickLabel = {'' '' '' '' ''};
ax.YTick = 0:20:100;
ax.YTickLabel = {'0' '20' '40' '60' '80' '100'};
ax.TickDir = 'out';
% -- Add stars if sig diffs (LP)
if sig_effect.WTBL_vs_WTSD_LP_N
	hold on 
	
	rectangle('Position',[0 ax.YLim(2)-5 12 5],'FaceColor',Sig12hrColor,'EdgeColor','none'); 
	bool_xlocsWTBL_LP_sig_values = ismember(e1.XData,xlocs_sig_effects.WTBL_vs_WTSD_LP_N);
	bool_xlocsWTSD_LP_sig_values = ismember(e2.XData,xlocs_sig_effects.WTBL_vs_WTSD_LP_N);
	largest_y1_vals = e1.YData(bool_xlocsWTBL_LP_sig_values) + e1.YPositiveDelta(bool_xlocsWTBL_LP_sig_values);
	largest_y2_vals = e2.YData(bool_xlocsWTSD_LP_sig_values) + e2.YPositiveDelta(bool_xlocsWTSD_LP_sig_values);
	ylocs_sig_effects = max([largest_y1_vals; largest_y2_vals]) + 10*ones(1,length(largest_y1_vals));
	ylocs_sig_effects(ylocs_sig_effects>ax.YLim(2))=ax.YLim(2);
	text(xlocs_sig_effects.WTBL_vs_WTSD_LP_N, ylocs_sig_effects, '*', 'FontSize',18,'HorizontalAlignment','center','VerticalAlignment','middle');

	hold off 
end 
% -- Add stars if sig diffs (DP)
if sig_effect.WTBL_vs_WTSD_DP_N
	hold on 
	
	rectangle('Position',[12 ax.YLim(2)-5 12 5],'FaceColor',Sig12hrColor,'EdgeColor','none'); 
	bool_xlocsWTBL_DP_sig_values = ismember(e1.XData,xlocs_sig_effects.WTBL_vs_WTSD_DP_N);
	bool_xlocsWTSD_DP_sig_values = ismember(e2.XData,xlocs_sig_effects.WTBL_vs_WTSD_DP_N);
	largest_y1_vals = e1.YData(bool_xlocsWTBL_DP_sig_values) + e1.YPositiveDelta(bool_xlocsWTBL_DP_sig_values);
	largest_y2_vals = e2.YData(bool_xlocsWTSD_DP_sig_values) + e2.YPositiveDelta(bool_xlocsWTSD_DP_sig_values);
	ylocs_sig_effects = max([largest_y1_vals; largest_y2_vals]) + 10*ones(1,length(largest_y1_vals));
	ylocs_sig_effects(ylocs_sig_effects>ax.YLim(2))=ax.YLim(2);
	text(xlocs_sig_effects.WTBL_vs_WTSD_DP_N, ylocs_sig_effects, '*', 'FontSize',18,'HorizontalAlignment','center','VerticalAlignment','middle');

	hold off 
end 
% --- End panel (2,1) ---------------------------------------------------------------------------------------------

% Position (3,1): REM  WT BL (black solid) vs WT SD (black dashed)
%subplot(3,3,7)
nexttile 
ax=H9.Children.Children(1);
rectangle('Position',LPboxposREM,'FaceColor','w','EdgeColor','k');  	% LP rectangle
hold on
rectangle('Position',DPboxposREM,'FaceColor','k','EdgeColor','none');    % DP rectangle
p1 = patch([0 0 SD_length SD_length],[-1 100 100 -1],[1 1 1],'LineStyle','none'); 	% SD hatching
hatchfill2(p1,'single','HatchAngle',15,'HatchDensity',125);
p2 = patch([0 0 SD_length SD_length],[-1 0 0 -1],[1 1 1]); 	% SD hatching
p2.Parent = ax;
hatchfill2(p2,'single','HatchAngle',15,'HatchDensity',125,'HatchLineWidth',1.5);
e1=errorbar([1:24],HourlyStatePercentages.Means.WT.BL.All24hrs.REM,HourlyStatePercentages.SEMs.WT.BL.All24hrs.REM,'Color','k','LineWidth',mylinewidth,'CapSize',0);
e2=errorbar([SD_length+1:24],HourlyStatePercentages.Means.WT.SD.All24hrexcSD.REM,HourlyStatePercentages.SEMs.WT.SD.All24hrexcSD.REM,'Color','k','LineWidth',mylinewidth,'CapSize',0,'LineStyle','--');
e1.Marker = '.';
e1.MarkerSize = EBMarkersize;
e2.Marker = '.';
e2.MarkerSize = EBMarkersize;
hold off 
set(ax, 'Layer','top')
title('')
ylabel({'{\bf REM}';'(% Recording Time)'});
xlabel('')
l = legend([e1 e2],'Wild Type BL','Wild Type SD');
l.Box = 'off';
l.Location = 'southeast';
l.Visible = 'off';
ax.XLim=[1 24];
ax.YLim=[-1 20];
ax.FontSize  = myfontsize;
ax.LineWidth = myaxislinewidth;
ax.XTick = 0:6:24;
ax.XTickLabel = {'0' '6' '12' '18' '24'};
ax.YTick = 0:4:30;
%ax.YTickLabel = {'0' '4' '8' '12' '16' '30'};
ax.TickDir = 'out';
% -- Add stars if sig diffs (LP)
if sig_effect.WTBL_vs_WTSD_LP_R
	hold on 
	
	rectangle('Position',[0 ax.YLim(2)-1 12 1],'FaceColor',Sig12hrColor,'EdgeColor','none'); 
	bool_xlocsWTBL_LP_sig_values = ismember(e1.XData,xlocs_sig_effects.WTBL_vs_WTSD_LP_R);
	bool_xlocsWTSD_LP_sig_values = ismember(e2.XData,xlocs_sig_effects.WTBL_vs_WTSD_LP_R);
	largest_y1_vals = e1.YData(bool_xlocsWTBL_LP_sig_values) + e1.YPositiveDelta(bool_xlocsWTBL_LP_sig_values);
	largest_y2_vals = e2.YData(bool_xlocsWTSD_LP_sig_values) + e2.YPositiveDelta(bool_xlocsWTSD_LP_sig_values);
	ylocs_sig_effects = max([largest_y1_vals; largest_y2_vals]) + 4*ones(1,length(largest_y1_vals));
	ylocs_sig_effects(ylocs_sig_effects>ax.YLim(2))=ax.YLim(2);
	text(xlocs_sig_effects.WTBL_vs_WTSD_LP_R, ylocs_sig_effects, '*', 'FontSize',18,'HorizontalAlignment','center','VerticalAlignment','middle');

	hold off 
end 
% -- Add stars if sig diffs (DP)
if sig_effect.WTBL_vs_WTSD_DP_R
	hold on 
	
	rectangle('Position',[12 ax.YLim(2)-1 12 1],'FaceColor',Sig12hrColor,'EdgeColor','none'); 
	bool_xlocsWTBL_DP_sig_values = ismember(e1.XData,xlocs_sig_effects.WTBL_vs_WTSD_DP_R);
	bool_xlocsWTSD_DP_sig_values = ismember(e2.XData,xlocs_sig_effects.WTBL_vs_WTSD_DP_R);
	largest_y1_vals = e1.YData(bool_xlocsWTBL_DP_sig_values) + e1.YPositiveDelta(bool_xlocsWTBL_DP_sig_values);
	largest_y2_vals = e2.YData(bool_xlocsWTSD_DP_sig_values) + e2.YPositiveDelta(bool_xlocsWTSD_DP_sig_values);
	ylocs_sig_effects = max([largest_y1_vals; largest_y2_vals]) + 4*ones(1,length(largest_y1_vals));
	ylocs_sig_effects(ylocs_sig_effects>ax.YLim(2))=ax.YLim(2);
	text(xlocs_sig_effects.WTBL_vs_WTSD_DP_R, ylocs_sig_effects, '*', 'FontSize',18,'HorizontalAlignment','center','VerticalAlignment','middle');

	hold off 
end 
% --- End panel (3,1) ---------------------------------------------------------------------------------------------



% Position (1,2): W    Mut BL (red solid) vs Mut SD (red dashed)
%subplot(3,3,2)
nexttile 
ax=H9.Children.Children(1);

rectangle('Position',LPboxpos,'FaceColor','w','EdgeColor','k');  	% LP rectangle
hold on
rectangle('Position',DPboxpos,'FaceColor','k','EdgeColor','none');    % DP rectangle
p1 = patch([0 0 SD_length SD_length],[-1 100 100 -1],[1 1 1],'LineStyle','none'); 	% SD hatching
%p1.FaceColor = [1 1 1];  % white face color
hatchfill2(p1,'single','HatchAngle',45,'HatchDensity',35);
p2 = patch([0 0 SD_length SD_length],[-1 4 4 -1],[1 1 1]); 	% SD hatching
p2.Parent = ax;
hatchfill2(p2,'single','HatchAngle',45,'HatchDensity',35,'HatchLineWidth',1.5);
e1=errorbar([1:24],HourlyStatePercentages.Means.Mut.BL.All24hrs.Wake,HourlyStatePercentages.SEMs.Mut.BL.All24hrs.Wake,'Color','r','LineWidth',mylinewidth,'CapSize',0);
e2=errorbar([SD_length+1:24],HourlyStatePercentages.Means.Mut.SD.All24hrexcSD.Wake,HourlyStatePercentages.SEMs.Mut.SD.All24hrexcSD.Wake,'Color','r','LineWidth',mylinewidth,'CapSize',0,'LineStyle','--');
e1.Marker = '.';
e1.MarkerSize = EBMarkersize;
e2.Marker = '.';
e2.MarkerSize = EBMarkersize;
hold off 
set(ax, 'Layer','top')
title('')
ylabel({'{\bf Wakefulness}';'(% Recording Time)'});
ylabel('')
xlabel('')
l = legend([e1 e2],[MutantText ' BL'],[MutantText ' SD']);
l.Box = 'off';
l.Location = 'southeast';
l.FontSize = 8;
xlim([1 24])
ylim([-1 100])
ax.FontSize  = myfontsize;
ax.LineWidth = myaxislinewidth;
ax.XTick = 0:6:24;
ax.XTickLabel = {'' '' '' '' ''};
ax.YTick = 0:20:100;
ax.YTickLabel = {'' '' '' '' '' ''};
ax.TickDir = 'out';
% -- Add stars if sig diffs (LP)
if sig_effect.MutBL_vs_MutSD_LP_W
	hold on 

	rectangle('Position',[0 ax.YLim(2)-5 12 5],'FaceColor',Sig12hrColor,'EdgeColor','none'); 
	bool_xlocsMutBL_LP_sig_values = ismember(e1.XData,xlocs_sig_effects.MutBL_vs_MutSD_LP_W);
	bool_xlocsMutSD_LP_sig_values = ismember(e2.XData,xlocs_sig_effects.MutBL_vs_MutSD_LP_W);
	largest_y1_vals = e1.YData(bool_xlocsMutBL_LP_sig_values) + e1.YPositiveDelta(bool_xlocsMutBL_LP_sig_values);
	largest_y2_vals = e2.YData(bool_xlocsMutSD_LP_sig_values) + e2.YPositiveDelta(bool_xlocsMutSD_LP_sig_values);
	ylocs_sig_effects = max([largest_y1_vals; largest_y2_vals]) + 10*ones(1,length(largest_y1_vals));
	ylocs_sig_effects(ylocs_sig_effects>ax.YLim(2))=ax.YLim(2);
	text(xlocs_sig_effects.MutBL_vs_MutSD_LP_W, ylocs_sig_effects, '*', 'FontSize',18,'HorizontalAlignment','center','VerticalAlignment','middle');

	hold off 
end 
% -- Add stars if sig diffs (DP)
if sig_effect.MutBL_vs_MutSD_DP_W
	hold on 
	
	rectangle('Position',[12 ax.YLim(2)-5 12 5],'FaceColor',Sig12hrColor,'EdgeColor','none');
	bool_xlocsMutBL_DP_sig_values = ismember(e1.XData,xlocs_sig_effects.MutBL_vs_MutSD_DP_W);
	bool_xlocsMutSD_DP_sig_values = ismember(e2.XData,xlocs_sig_effects.MutBL_vs_MutSD_DP_W);
	largest_y1_vals = e1.YData(bool_xlocsMutBL_DP_sig_values) + e1.YPositiveDelta(bool_xlocsMutBL_DP_sig_values);
	largest_y2_vals = e2.YData(bool_xlocsMutSD_DP_sig_values) + e2.YPositiveDelta(bool_xlocsMutSD_DP_sig_values);
	ylocs_sig_effects = max([largest_y1_vals; largest_y2_vals]) + 10*ones(1,length(largest_y1_vals));
	ylocs_sig_effects(ylocs_sig_effects>ax.YLim(2))=ax.YLim(2);
	text(xlocs_sig_effects.MutBL_vs_MutSD_DP_W, ylocs_sig_effects, '*', 'FontSize',18,'HorizontalAlignment','center','VerticalAlignment','middle');

	hold off 
end 
% --- End panel (1,2) ---------------------------------------------------------------------------------------------- 

% Position (2,2): NREM Mut BL (red solid) vs Mut SD (red dashed)
%subplot(3,3,5)
nexttile
ax=H9.Children.Children(1);
rectangle('Position',LPboxpos,'FaceColor','w','EdgeColor','k');  	% LP rectangle
hold on
rectangle('Position',DPboxpos,'FaceColor','k','EdgeColor','none');    % DP rectangle
p1 = patch([0 0 SD_length SD_length],[-1 100 100 -1],[1 1 1],'LineStyle','none'); 	% SD hatching
hatchfill2(p1,'single','HatchAngle',45,'HatchDensity',35);
p2 = patch([0 0 SD_length SD_length],[-1 4 4 -1],[1 1 1]); 	% SD hatching
p2.Parent = ax;
hatchfill2(p2,'single','HatchAngle',45,'HatchDensity',35,'HatchLineWidth',1.5);
e1=errorbar([1:24],HourlyStatePercentages.Means.Mut.BL.All24hrs.NREM,HourlyStatePercentages.SEMs.Mut.BL.All24hrs.NREM,'Color','r','LineWidth',mylinewidth,'CapSize',0);
e2=errorbar([SD_length+1:24],HourlyStatePercentages.Means.Mut.SD.All24hrexcSD.NREM,HourlyStatePercentages.SEMs.Mut.SD.All24hrexcSD.NREM,'Color','r','LineWidth',mylinewidth,'CapSize',0,'LineStyle','--');
e1.Marker = '.';
e1.MarkerSize = EBMarkersize;
e2.Marker = '.';
e2.MarkerSize = EBMarkersize;
hold off 
set(ax, 'Layer','top')
title('')
ylabel({'{\bf NREM}';'(% Recording Time)'});
ylabel('')
xlabel('')
l = legend([e1 e2],'Wild Type BL','Wild Type SD');
l.Box = 'off';
l.Location = 'southeast';
l.Visible = 'off';
xlim([1 24])
ylim([-1 100])
ax.FontSize  = myfontsize;
ax.LineWidth = myaxislinewidth;
ax.XTick = 0:6:24;
ax.XTickLabel = {'' '' '' '' ''};
ax.YTick = 0:20:100;
ax.YTickLabel = {'' '' '' '' '' ''};
ax.TickDir = 'out';
% -- Add stars if sig diffs (LP)
if sig_effect.MutBL_vs_MutSD_LP_N
	hold on 
	
	rectangle('Position',[0 ax.YLim(2)-5 12 5],'FaceColor',Sig12hrColor,'EdgeColor','none'); 
	bool_xlocsMutBL_LP_sig_values = ismember(e1.XData,xlocs_sig_effects.MutBL_vs_MutSD_LP_N);
	bool_xlocsMutSD_LP_sig_values = ismember(e2.XData,xlocs_sig_effects.MutBL_vs_MutSD_LP_N);
	largest_y1_vals = e1.YData(bool_xlocsMutBL_LP_sig_values) + e1.YPositiveDelta(bool_xlocsMutBL_LP_sig_values);
	largest_y2_vals = e2.YData(bool_xlocsMutSD_LP_sig_values) + e2.YPositiveDelta(bool_xlocsMutSD_LP_sig_values);
	ylocs_sig_effects = max([largest_y1_vals; largest_y2_vals]) + 10*ones(1,length(largest_y1_vals));
	ylocs_sig_effects(ylocs_sig_effects>ax.YLim(2))=ax.YLim(2);
	text(xlocs_sig_effects.MutBL_vs_MutSD_LP_N, ylocs_sig_effects, '*', 'FontSize',18,'HorizontalAlignment','center','VerticalAlignment','middle');

	hold off 
end 
% -- Add stars if sig diffs (DP)
if sig_effect.MutBL_vs_MutSD_DP_N
	hold on 
	
	rectangle('Position',[12 ax.YLim(2)-5 12 5],'FaceColor',Sig12hrColor,'EdgeColor','none');
	bool_xlocsMutBL_DP_sig_values = ismember(e1.XData,xlocs_sig_effects.MutBL_vs_MutSD_DP_N);
	bool_xlocsMutSD_DP_sig_values = ismember(e2.XData,xlocs_sig_effects.MutBL_vs_MutSD_DP_N);
	largest_y1_vals = e1.YData(bool_xlocsMutBL_DP_sig_values) + e1.YPositiveDelta(bool_xlocsMutBL_DP_sig_values);
	largest_y2_vals = e2.YData(bool_xlocsMutSD_DP_sig_values) + e2.YPositiveDelta(bool_xlocsMutSD_DP_sig_values);
	ylocs_sig_effects = max([largest_y1_vals; largest_y2_vals]) + 10*ones(1,length(largest_y1_vals));
	ylocs_sig_effects(ylocs_sig_effects>ax.YLim(2))=ax.YLim(2);
	text(xlocs_sig_effects.MutBL_vs_MutSD_DP_N, ylocs_sig_effects, '*', 'FontSize',18,'HorizontalAlignment','center','VerticalAlignment','middle');

	hold off 
end 
% --- End Panel (2,2) ---------------------------------------------------------------------------------------------------



% Position (3,2): REM  Mut BL (red solid) vs Mut SD (red dashed)
%subplot(3,3,8)
nexttile
ax=H9.Children.Children(1);

rectangle('Position',LPboxposREM,'FaceColor','w','EdgeColor','k');  	% LP rectangle
hold on
rectangle('Position',DPboxposREM,'FaceColor','k','EdgeColor','none');    % DP rectangle
p1 = patch([0 0 SD_length SD_length],[-1 100 100 -1],[1 1 1],'LineStyle','none'); 	% SD hatching
hatchfill2(p1,'single','HatchAngle',15,'HatchDensity',125);
p2 = patch([0 0 SD_length SD_length],[-1 0 0 -1],[1 1 1]); 	% SD hatching
p2.Parent = ax;
hatchfill2(p2,'single','HatchAngle',15,'HatchDensity',125,'HatchLineWidth',1.5);
e1=errorbar([1:24],HourlyStatePercentages.Means.Mut.BL.All24hrs.REM,HourlyStatePercentages.SEMs.Mut.BL.All24hrs.REM,'Color','r','LineWidth',mylinewidth,'CapSize',0);
e2=errorbar([SD_length+1:24],HourlyStatePercentages.Means.Mut.SD.All24hrexcSD.REM,HourlyStatePercentages.SEMs.Mut.SD.All24hrexcSD.REM,'Color','r','LineWidth',mylinewidth,'CapSize',0,'LineStyle','--');
e1.Marker = '.';
e1.MarkerSize = EBMarkersize;
e2.Marker = '.';
e2.MarkerSize = EBMarkersize;
hold off 
set(ax, 'Layer','top')
title('')
ylabel({'{\bf REM}';'(% Recording Time)'});
ylabel('')
xlabel('Time (h)')
l = legend([e1 e2],'Wild Type BL','Wild Type SD');
l.Box = 'off';
l.Location = 'southeast';
l.Visible = 'off';
ax.XLim=[1 24];
ax.YLim=[-1 20];
ax.FontSize  = myfontsize;
ax.LineWidth = myaxislinewidth;
ax.XTick = 0:6:24;
ax.XTickLabel = {'0' '6' '12' '18' '24'};
ax.YTick = 0:4:30;
ax.YTickLabel = {'' '' '' '' '' ''};
ax.TickDir = 'out';
% -- Add stars if sig diffs (LP)
if sig_effect.MutBL_vs_MutSD_LP_R
	hold on 
	
	rectangle('Position',[0 ax.YLim(2)-1 12 1],'FaceColor',Sig12hrColor,'EdgeColor','none'); 
	bool_xlocsMutBL_LP_sig_values = ismember(e1.XData,xlocs_sig_effects.MutBL_vs_MutSD_LP_R);
	bool_xlocsMutSD_LP_sig_values = ismember(e2.XData,xlocs_sig_effects.MutBL_vs_MutSD_LP_R);
	largest_y1_vals = e1.YData(bool_xlocsMutBL_LP_sig_values) + e1.YPositiveDelta(bool_xlocsMutBL_LP_sig_values);
	largest_y2_vals = e2.YData(bool_xlocsMutSD_LP_sig_values) + e2.YPositiveDelta(bool_xlocsMutSD_LP_sig_values);
	ylocs_sig_effects = max([largest_y1_vals; largest_y2_vals]) + 4*ones(1,length(largest_y1_vals));
	ylocs_sig_effects(ylocs_sig_effects>ax.YLim(2))=ax.YLim(2);
	text(xlocs_sig_effects.MutBL_vs_MutSD_LP_R, ylocs_sig_effects, '*', 'FontSize',18,'HorizontalAlignment','center','VerticalAlignment','middle');

	hold off 
end 
% -- Add stars if sig diffs (DP)
if sig_effect.MutBL_vs_MutSD_DP_R
	hold on 
	
	rectangle('Position',[12 ax.YLim(2)-1 12 1],'FaceColor',Sig12hrColor,'EdgeColor','none');
	bool_xlocsMutBL_DP_sig_values = ismember(e1.XData,xlocs_sig_effects.MutBL_vs_MutSD_DP_R);
	bool_xlocsMutSD_DP_sig_values = ismember(e2.XData,xlocs_sig_effects.MutBL_vs_MutSD_DP_R);
	largest_y1_vals = e1.YData(bool_xlocsMutBL_DP_sig_values) + e1.YPositiveDelta(bool_xlocsMutBL_DP_sig_values);
	largest_y2_vals = e2.YData(bool_xlocsMutSD_DP_sig_values) + e2.YPositiveDelta(bool_xlocsMutSD_DP_sig_values);
	ylocs_sig_effects = max([largest_y1_vals; largest_y2_vals]) + 4*ones(1,length(largest_y1_vals));
	ylocs_sig_effects(ylocs_sig_effects>ax.YLim(2))=ax.YLim(2);
	text(xlocs_sig_effects.MutBL_vs_MutSD_DP_R, ylocs_sig_effects, '*', 'FontSize',18,'HorizontalAlignment','center','VerticalAlignment','middle');

	hold off 
end 
% --- End Panel (3,2) ---------------------------------------------------------------------------------------------------


% Position (1,3): W  WT SD (black solid) vs Mut SD (red solid)
%subplot(3,3,3)
nexttile
ax=H9.Children.Children(1);

rectangle('Position',LPboxpos,'FaceColor','w','EdgeColor','k');  	% LP rectangle
hold on
rectangle('Position',DPboxpos,'FaceColor','k','EdgeColor','none');    % DP rectangle
p1 = patch([0 0 SD_length SD_length],[-1 100 100 -1],[1 1 1],'LineStyle','none'); 	% SD hatching
hatchfill2(p1,'single','HatchAngle',45,'HatchDensity',35);
p2 = patch([0 0 SD_length SD_length],[-1 4 4 -1],[1 1 1]); 	% SD hatching
p2.Parent = ax;
hatchfill2(p2,'single','HatchAngle',45,'HatchDensity',35,'HatchLineWidth',1.5);
e1=errorbar([SD_length+1:24],HourlyStatePercentages.Means.WT.SD.All24hrexcSD.Wake,HourlyStatePercentages.SEMs.WT.SD.All24hrexcSD.Wake,'Color','k','LineWidth',mylinewidth,'CapSize',0,'LineStyle','--');
e2=errorbar([SD_length+1:24],HourlyStatePercentages.Means.Mut.SD.All24hrexcSD.Wake,HourlyStatePercentages.SEMs.Mut.SD.All24hrexcSD.Wake,'Color','r','LineWidth',mylinewidth,'CapSize',0,'LineStyle','--');
e1.Marker = '.';
e1.MarkerSize = EBMarkersize;
e2.Marker = '.';
e2.MarkerSize = EBMarkersize;
hold off 
set(ax, 'Layer','top')
title('')
ylabel({'{\bf Wakefulness}';'(% Recording Time)'});
ylabel('')
xlabel('')
l = legend([e1 e2],'Wild Type SD',[MutantText ' SD']);
l.Box = 'off';
l.Location = 'southeast';
l.FontSize = 8;
xlim([1 24])
ylim([-1 100])
ax.FontSize  = myfontsize;
ax.LineWidth = myaxislinewidth;
ax.XTick = 0:6:24;
ax.XTickLabel = {'' '' '' '' ''};
ax.YTick = 0:20:100;
ax.YTickLabel = {'' '' '' '' '' ''};
ax.TickDir = 'out';
% -- Add stars if sig diffs (LP)
if sig_effect.WTSD_vs_MutSD_LP_W
	hold on 
	
	rectangle('Position',[0 ax.YLim(2)-5 12 5],'FaceColor',Sig12hrColor,'EdgeColor','none'); 
	bool_xlocsWTSD_LP_sig_values  = ismember(e1.XData,xlocs_sig_effects.WTSD_vs_MutSD_LP_W);
	bool_xlocsMutSD_LP_sig_values = ismember(e2.XData,xlocs_sig_effects.WTSD_vs_MutSD_LP_W);
	largest_y1_vals = e1.YData(bool_xlocsWTSD_LP_sig_values) + e1.YPositiveDelta(bool_xlocsWTSD_LP_sig_values);
	largest_y2_vals = e2.YData(bool_xlocsMutSD_LP_sig_values) + e2.YPositiveDelta(bool_xlocsMutSD_LP_sig_values);
	ylocs_sig_effects = max([largest_y1_vals; largest_y2_vals]) + 10*ones(1,length(largest_y1_vals));
	ylocs_sig_effects(ylocs_sig_effects>ax.YLim(2))=ax.YLim(2);
	text(xlocs_sig_effects.WTSD_vs_MutSD_LP_W, ylocs_sig_effects, '*', 'FontSize',18,'HorizontalAlignment','center','VerticalAlignment','middle');

	hold off 
end 
% -- Add stars if sig diffs (DP)
if sig_effect.WTSD_vs_MutSD_DP_W
	hold on 
	
	rectangle('Position',[12 ax.YLim(2)-5 12 5],'FaceColor',Sig12hrColor,'EdgeColor','none');
	bool_xlocsWTSD_DP_sig_values  = ismember(e1.XData,xlocs_sig_effects.WTSD_vs_MutSD_DP_W);
	bool_xlocsMutSD_DP_sig_values = ismember(e2.XData,xlocs_sig_effects.WTSD_vs_MutSD_DP_W);
	largest_y1_vals = e1.YData(bool_xlocsWTSD_DP_sig_values) + e1.YPositiveDelta(bool_xlocsWTSD_DP_sig_values);
	largest_y2_vals = e2.YData(bool_xlocsMutSD_DP_sig_values) + e2.YPositiveDelta(bool_xlocsMutSD_DP_sig_values);
	ylocs_sig_effects = max([largest_y1_vals; largest_y2_vals]) + 10*ones(1,length(largest_y1_vals));
	ylocs_sig_effects(ylocs_sig_effects>ax.YLim(2))=ax.YLim(2);
	text(xlocs_sig_effects.WTSD_vs_MutSD_DP_W, ylocs_sig_effects, '*', 'FontSize',18,'HorizontalAlignment','center','VerticalAlignment','middle');

	hold off 
end 
% --- End Panel (1,3) ---------------------------------------------------------------------------------------------------



% Position (2,3): NREM  WT SD (black solid) vs Mut SD (red solid)
%subplot(3,3,6)
nexttile
ax=H9.Children.Children(1);
rectangle('Position',LPboxpos,'FaceColor','w','EdgeColor','k');  	% LP rectangle
hold on
rectangle('Position',DPboxpos,'FaceColor','k','EdgeColor','none');    % DP rectangle
p1 = patch([0 0 SD_length SD_length],[-1 100 100 -1],[1 1 1],'LineStyle','none'); 	% SD hatching
hatchfill2(p1,'single','HatchAngle',45,'HatchDensity',35);
p2 = patch([0 0 SD_length SD_length],[-1 4 4 -1],[1 1 1]); 	% SD hatching
p2.Parent = ax;
hatchfill2(p2,'single','HatchAngle',45,'HatchDensity',35,'HatchLineWidth',1.5);
e1=errorbar([SD_length+1:24],HourlyStatePercentages.Means.WT.SD.All24hrexcSD.NREM,HourlyStatePercentages.SEMs.WT.SD.All24hrexcSD.NREM,'Color','k','LineWidth',mylinewidth,'CapSize',0,'LineStyle','--');
e2=errorbar([SD_length+1:24],HourlyStatePercentages.Means.Mut.SD.All24hrexcSD.NREM,HourlyStatePercentages.SEMs.Mut.SD.All24hrexcSD.NREM,'Color','r','LineWidth',mylinewidth,'CapSize',0,'LineStyle','--');
e1.Marker = '.';
e1.MarkerSize = EBMarkersize;
e2.Marker = '.';
e2.MarkerSize = EBMarkersize;
hold off 
set(ax, 'Layer','top')
title('')
ylabel({'{\bf NREM}';'(% Recording Time)'});
ylabel('')
xlabel('')
l = legend([e1 e2],'Wild Type BL','Wild Type SD');
l.Box = 'off';
l.Location = 'southeast';
l.Visible = 'off';
xlim([1 24])
ylim([-1 100])
ax.FontSize  = myfontsize;
ax.LineWidth = myaxislinewidth;
ax.XTick = 0:6:24;
ax.XTickLabel = {'' '' '' '' ''};
ax.YTick = 0:20:100;
ax.YTickLabel = {'' '' '' '' '' ''};
ax.TickDir = 'out';
% -- Add stars if sig diffs (LP)
if sig_effect.WTSD_vs_MutSD_LP_N
	hold on 
	
	rectangle('Position',[0 ax.YLim(2)-5 12 5],'FaceColor',Sig12hrColor,'EdgeColor','none');
	bool_xlocsWTSD_LP_sig_values  = ismember(e1.XData,xlocs_sig_effects.WTSD_vs_MutSD_LP_N);
	bool_xlocsMutSD_LP_sig_values = ismember(e2.XData,xlocs_sig_effects.WTSD_vs_MutSD_LP_N);
	largest_y1_vals = e1.YData(bool_xlocsWTSD_LP_sig_values) + e1.YPositiveDelta(bool_xlocsWTSD_LP_sig_values);
	largest_y2_vals = e2.YData(bool_xlocsMutSD_LP_sig_values) + e2.YPositiveDelta(bool_xlocsMutSD_LP_sig_values);
	ylocs_sig_effects = max([largest_y1_vals; largest_y2_vals]) + 10*ones(1,length(largest_y1_vals));
	ylocs_sig_effects(ylocs_sig_effects>ax.YLim(2))=ax.YLim(2);
	text(xlocs_sig_effects.WTSD_vs_MutSD_LP_N, ylocs_sig_effects, '*', 'FontSize',18,'HorizontalAlignment','center','VerticalAlignment','middle');

	hold off 
end 
% -- Add stars if sig diffs (DP)
if sig_effect.WTSD_vs_MutSD_DP_N
	hold on 
	
	rectangle('Position',[12 ax.YLim(2)-5 12 5],'FaceColor',Sig12hrColor,'EdgeColor','none');
	bool_xlocsWTSD_DP_sig_values  = ismember(e1.XData,xlocs_sig_effects.WTSD_vs_MutSD_DP_N);
	bool_xlocsMutSD_DP_sig_values = ismember(e2.XData,xlocs_sig_effects.WTSD_vs_MutSD_DP_N);
	largest_y1_vals = e1.YData(bool_xlocsWTSD_DP_sig_values) + e1.YPositiveDelta(bool_xlocsWTSD_DP_sig_values);
	largest_y2_vals = e2.YData(bool_xlocsMutSD_DP_sig_values) + e2.YPositiveDelta(bool_xlocsMutSD_DP_sig_values);
	ylocs_sig_effects = max([largest_y1_vals; largest_y2_vals]) + 10*ones(1,length(largest_y1_vals));
	ylocs_sig_effects(ylocs_sig_effects>ax.YLim(2))=ax.YLim(2);
	text(xlocs_sig_effects.WTSD_vs_MutSD_DP_N, ylocs_sig_effects, '*', 'FontSize',18,'HorizontalAlignment','center','VerticalAlignment','middle');

	hold off 
end 
% --- End Panel (2,3) ---------------------------------------------------------------------------------------------------



% Position (3,3): REM  WT SD (black solid) vs Mut SD (red solid)
%subplot(3,3,9)
nexttile
ax=H9.Children.Children(1);
rectangle('Position',LPboxposREM,'FaceColor','w','EdgeColor','k');  	% LP rectangle
hold on
rectangle('Position',DPboxposREM,'FaceColor','k','EdgeColor','none');    % DP rectangle
p1 = patch([0 0 SD_length SD_length],[-1 100 100 -1],[1 1 1],'LineStyle','none'); 	% SD hatching
hatchfill2(p1,'single','HatchAngle',15,'HatchDensity',125);
p2 = patch([0 0 SD_length SD_length],[-1 0 0 -1],[1 1 1]); 	% SD hatching
p2.Parent = ax;
hatchfill2(p2,'single','HatchAngle',15,'HatchDensity',125,'HatchLineWidth',1.5);

e1=errorbar([SD_length+1:24],HourlyStatePercentages.Means.WT.SD.All24hrexcSD.REM,HourlyStatePercentages.SEMs.WT.SD.All24hrexcSD.REM,'Color','k','LineWidth',mylinewidth,'CapSize',0,'LineStyle','--');
e2=errorbar([SD_length+1:24],HourlyStatePercentages.Means.Mut.SD.All24hrexcSD.REM,HourlyStatePercentages.SEMs.Mut.SD.All24hrexcSD.REM,'Color','r','LineWidth',mylinewidth,'CapSize',0,'LineStyle','--');
e1.Marker = '.';
e1.MarkerSize = EBMarkersize;
e2.Marker = '.';
e2.MarkerSize = EBMarkersize;
hold off 
set(ax, 'Layer','top')
title('')
ylabel({'{\bf REM}';'(% Recording Time)'});
ylabel('')
xlabel('')
l = legend([e1 e2],'Wild Type BL','Wild Type SD');
l.Box = 'off';
l.Location = 'southeast';
l.Visible = 'off';

ax.XLim=[1 24];
ax.YLim=[-1 30];
ax.FontSize  = myfontsize;
ax.LineWidth = myaxislinewidth;
ax.XTick = 0:6:24;
ax.XTickLabel = {'0' '6' '12' '18' '24'};
ax.YTick = 0:4:30;
ax.YTickLabel = {'' '' '' '' '' ''};
ax.TickDir = 'out';
% -- Add stars if sig diffs (LP)
if sig_effect.WTSD_vs_MutSD_LP_R
	hold on 
	
	rectangle('Position',[0 ax.YLim(2)-1 12 1],'FaceColor',Sig12hrColor,'EdgeColor','none');
	bool_xlocsWTSD_LP_sig_values  = ismember(e1.XData,xlocs_sig_effects.WTSD_vs_MutSD_LP_R);
	bool_xlocsMutSD_LP_sig_values = ismember(e2.XData,xlocs_sig_effects.WTSD_vs_MutSD_LP_R);
	largest_y1_vals = e1.YData(bool_xlocsWTSD_LP_sig_values) + e1.YPositiveDelta(bool_xlocsWTSD_LP_sig_values);
	largest_y2_vals = e2.YData(bool_xlocsMutSD_LP_sig_values) + e2.YPositiveDelta(bool_xlocsMutSD_LP_sig_values);
	ylocs_sig_effects = max([largest_y1_vals; largest_y2_vals]) + 4*ones(1,length(largest_y1_vals));
	ylocs_sig_effects(ylocs_sig_effects>ax.YLim(2))=ax.YLim(2);
	text(xlocs_sig_effects.WTSD_vs_MutSD_LP_R, ylocs_sig_effects, '*', 'FontSize',18,'HorizontalAlignment','center','VerticalAlignment','middle');

	hold off 
end 
% -- Add stars if sig diffs (DP)
if sig_effect.WTSD_vs_MutSD_DP_R
	hold on 
	
	rectangle('Position',[12 ax.YLim(2)-1 12 1],'FaceColor',Sig12hrColor,'EdgeColor','none');
	bool_xlocsWTSD_DP_sig_values  = ismember(e1.XData,xlocs_sig_effects.WTSD_vs_MutSD_DP_R);
	bool_xlocsMutSD_DP_sig_values = ismember(e2.XData,xlocs_sig_effects.WTSD_vs_MutSD_DP_R);
	largest_y1_vals = e1.YData(bool_xlocsWTSD_DP_sig_values) + e1.YPositiveDelta(bool_xlocsWTSD_DP_sig_values);
	largest_y2_vals = e2.YData(bool_xlocsMutSD_DP_sig_values) + e2.YPositiveDelta(bool_xlocsMutSD_DP_sig_values);
	ylocs_sig_effects = max([largest_y1_vals; largest_y2_vals]) + 4*ones(1,length(largest_y1_vals));
	ylocs_sig_effects(ylocs_sig_effects>ax.YLim(2))=ax.YLim(2);
	text(xlocs_sig_effects.WTSD_vs_MutSD_DP_R, ylocs_sig_effects, '*', 'FontSize',18,'HorizontalAlignment','center','VerticalAlignment','middle');

	hold off 
end 
% --- End Panel (3,3) ---------------------------------------------------------------------------------------------------

title(t,MorF)





%  Return all RM anova tables and multcompare post-hoc tables? 
% -- Column 1: WT BL vs WT SD --
RM_Anova_tables.WTBLvsWTSD.Wake.LP = RManovatbl_WTBL_vs_WTSD_LP_W;  % Panel 1 LP W
RM_Anova_tables.WTBLvsWTSD.Wake.DP = RManovatbl_WTBL_vs_WTSD_DP_W;	% Panel 1 DP W

RM_Anova_tables.WTBLvsWTSD.NREM.LP = RManovatbl_WTBL_vs_WTSD_LP_N;	% Panel 2 LP NREM
RM_Anova_tables.WTBLvsWTSD.NREM.DP = RManovatbl_WTBL_vs_WTSD_DP_N;	% Panel 2 DP NREM

RM_Anova_tables.WTBLvsWTSD.REM.LP  = RManovatbl_WTBL_vs_WTSD_LP_R;	% Panel 3 LP REM
RM_Anova_tables.WTBLvsWTSD.REM.DP  = RManovatbl_WTBL_vs_WTSD_DP_R; 	% Panel 3 DP REM 

% -- Column 2: Mut BL vs Mut SD -- 
RM_Anova_tables.MutBLvsMutSD.Wake.LP = RManovatbl_MutBL_vs_MutSD_LP_W; 	% Panel 4 LP W
RM_Anova_tables.MutBLvsMutSD.Wake.DP = RManovatbl_MutBL_vs_MutSD_DP_W; 	% Panel 4 DP W

RM_Anova_tables.MutBLvsMutSD.NREM.LP = RManovatbl_MutBL_vs_MutSD_LP_N; 	% Panel 5 LP NREM
RM_Anova_tables.MutBLvsMutSD.NREM.DP = RManovatbl_MutBL_vs_MutSD_DP_N;	% Panel 5 DP NREM

RM_Anova_tables.MutBLvsMutSD.REM.LP  = RManovatbl_MutBL_vs_MutSD_LP_R; 	% Panel 6 LP REM
RM_Anova_tables.MutBLvsMutSD.REM.DP  = RManovatbl_MutBL_vs_MutSD_DP_R; 	% Panel 6 DP REM

% -- Column 3: WT SD vs Mut SD --
RM_Anova_tables.WTSDvsMutSD.Wake.LP = RManovatbl_WTSD_vs_MutSD_LP_W;	% Panel 7 LP W
RM_Anova_tables.WTSDvsMutSD.Wake.DP = RManovatbl_WTSD_vs_MutSD_DP_W; 	% Panel 7 DP W

RM_Anova_tables.WTSDvsMutSD.NREM.LP = RManovatbl_WTSD_vs_MutSD_LP_N; 	% Panel 8 LP NREM
RM_Anova_tables.WTSDvsMutSD.NREM.DP = RManovatbl_WTSD_vs_MutSD_DP_N; 	% Panel 8 DP NREM 

RM_Anova_tables.WTSDvsMutSD.REM.LP  = RManovatbl_WTSD_vs_MutSD_LP_R; 	% Panel 9 LP REM
RM_Anova_tables.WTSDvsMutSD.REM.DP  = RManovatbl_WTSD_vs_MutSD_DP_R; 	% Panel 9 DP REM
% --------------------------------------------------------------------------------------


% ---- Post-hoc tables -----------------------------------------------------------------
if ~exist('mc_table')
	return				% if the struct mc_table does not exist, exit the function
end 
% panel 1 W 
if isfield(mc_table,'panel1_LP')  PostHoc_tables.WTBLvsWTSD.Wake.LP  = mc_table.panel1_LP;  else PostHoc_tables.WTBLvsWTSD.Wake.LP = []; end  
if isfield(mc_table,'panel1_DP')  PostHoc_tables.WTBLvsWTSD.Wake.DP  = mc_table.panel1_DP;  else PostHoc_tables.WTBLvsWTSD.Wake.DP = []; end  

% panel 2 N
if isfield(mc_table,'panel2_LP')  PostHoc_tables.WTBLvsWTSD.NREM.LP  = mc_table.panel2_LP;  else PostHoc_tables.WTBLvsWTSD.NREM.LP = []; end  
if isfield(mc_table,'panel2_DP')  PostHoc_tables.WTBLvsWTSD.NREM.DP  = mc_table.panel2_DP;  else PostHoc_tables.WTBLvsWTSD.NREM.DP = []; end  

% panel 3 R
if isfield(mc_table,'panel3_LP')  PostHoc_tables.WTBLvsWTSD.REM.LP  = mc_table.panel3_LP;  else PostHoc_tables.WTBLvsWTSD.REM.LP = []; end  
if isfield(mc_table,'panel3_DP')  PostHoc_tables.WTBLvsWTSD.REM.DP  = mc_table.panel3_DP;  else PostHoc_tables.WTBLvsWTSD.REM.DP = []; end  

% panel 4 W MutBL_vs_MutSD_DP_W
if isfield(mc_table,'panel4_LP')  PostHoc_tables.MutBLvsMutSD.Wake.LP  = mc_table.panel4_LP;  else PostHoc_tables.MutBLvsMutSD.Wake.LP = []; end  
if isfield(mc_table,'panel4_DP')  PostHoc_tables.MutBLvsMutSD.Wake.DP  = mc_table.panel4_DP;  else PostHoc_tables.MutBLvsMutSD.Wake.DP = []; end  

% panel 5 N MutBL_vs_MutSD_DP_W
if isfield(mc_table,'panel5_LP')  PostHoc_tables.MutBLvsMutSD.NREM.LP  = mc_table.panel5_LP;  else PostHoc_tables.MutBLvsMutSD.NREM.LP = []; end  
if isfield(mc_table,'panel5_DP')  PostHoc_tables.MutBLvsMutSD.NREM.DP  = mc_table.panel5_DP;  else PostHoc_tables.MutBLvsMutSD.NREM.DP = []; end  

% panel 6 R MutBL_vs_MutSD_DP_W
if isfield(mc_table,'panel6_LP')  PostHoc_tables.MutBLvsMutSD.REM.LP  = mc_table.panel6_LP;  else PostHoc_tables.MutBLvsMutSD.REM.LP = []; end  
if isfield(mc_table,'panel6_DP')  PostHoc_tables.MutBLvsMutSD.REM.DP  = mc_table.panel6_DP;  else PostHoc_tables.MutBLvsMutSD.REM.DP = []; end  

% panel 7 W WTSDvsMutSD
if isfield(mc_table,'panel7_LP')  PostHoc_tables.WTSDvsMutSD.Wake.LP  = mc_table.panel7_LP;  else PostHoc_tables.WTSDvsMutSD.Wake.LP = []; end  
if isfield(mc_table,'panel7_DP')  PostHoc_tables.WTSDvsMutSD.Wake.DP  = mc_table.panel7_DP;  else PostHoc_tables.WTSDvsMutSD.Wake.DP = []; end  

% panel 8 N WTSDvsMutSD
if isfield(mc_table,'panel8_LP')  PostHoc_tables.WTSDvsMutSD.NREM.LP  = mc_table.panel8_LP;  else PostHoc_tables.WTSDvsMutSD.NREM.LP = []; end  
if isfield(mc_table,'panel8_DP')  PostHoc_tables.WTSDvsMutSD.NREM.DP  = mc_table.panel8_DP;  else PostHoc_tables.WTSDvsMutSD.NREM.DP = []; end  

% panel 9 R WTSDvsMutSD
if isfield(mc_table,'panel9_LP')  PostHoc_tables.WTSDvsMutSD.REM.LP  = mc_table.panel9_LP;  else PostHoc_tables.WTSDvsMutSD.REM.LP = []; end  
if isfield(mc_table,'panel9_DP')  PostHoc_tables.WTSDvsMutSD.REM.DP  = mc_table.panel9_DP;  else PostHoc_tables.WTSDvsMutSD.REM.DP = []; end  
































































































