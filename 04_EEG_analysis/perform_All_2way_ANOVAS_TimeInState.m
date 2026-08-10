function [p_vals,sig_difference_bool,ANOVA_tbls] = perform_All_2way_ANOVAS_TimeInState(Twelve_hour_avg_percentages_Male,Twelve_hour_avg_percentages_Female);
%
% This function simply calls perform_2way_Anova_SexGenotype_posthoc.m for Wake, NREM, REM in light phase dark phase and BL and sleep dep


% If you don't have female data, set up those structs as empty
if isempty(Twelve_hour_avg_percentages_Female)
	Twelve_hour_avg_percentages_Female.WT.BL.First12hrs.Wake  = [];
	Twelve_hour_avg_percentages_Female.Mut.BL.First12hrs.Wake = [];
	Twelve_hour_avg_percentages_Female.WT.BL.First12hrs.NREM  = [];
	Twelve_hour_avg_percentages_Female.Mut.BL.First12hrs.NREM = [];
	Twelve_hour_avg_percentages_Female.WT.BL.First12hrs.REM   = [];
	Twelve_hour_avg_percentages_Female.Mut.BL.First12hrs.REM  = [];
	Twelve_hour_avg_percentages_Female.WT.BL.Last12hrs.Wake   = [];
	Twelve_hour_avg_percentages_Female.Mut.BL.Last12hrs.Wake  = [];
	Twelve_hour_avg_percentages_Female.WT.BL.Last12hrs.NREM   = [];
	Twelve_hour_avg_percentages_Female.Mut.BL.Last12hrs.NREM  = [];
	Twelve_hour_avg_percentages_Female.WT.BL.Last12hrs.REM    = [];
	Twelve_hour_avg_percentages_Female.Mut.BL.Last12hrs.REM   = [];

	Twelve_hour_avg_percentages_Female.WT.SD.First12hrs.Wake  = [];
	Twelve_hour_avg_percentages_Female.Mut.SD.First12hrs.Wake = [];
	Twelve_hour_avg_percentages_Female.WT.SD.First12hrs.NREM  = [];
	Twelve_hour_avg_percentages_Female.Mut.SD.First12hrs.NREM = [];
	Twelve_hour_avg_percentages_Female.WT.SD.First12hrs.REM   = [];
	Twelve_hour_avg_percentages_Female.Mut.SD.First12hrs.REM  = [];
	Twelve_hour_avg_percentages_Female.WT.SD.Last12hrs.Wake   = [];
	Twelve_hour_avg_percentages_Female.Mut.SD.Last12hrs.Wake  = [];
	Twelve_hour_avg_percentages_Female.WT.SD.Last12hrs.NREM   = [];
	Twelve_hour_avg_percentages_Female.Mut.SD.Last12hrs.NREM  = [];
	Twelve_hour_avg_percentages_Female.WT.SD.Last12hrs.REM    = [];
	Twelve_hour_avg_percentages_Female.Mut.SD.Last12hrs.REM   = [];

	Twelve_hour_avg_percentages_Female.WT.SDexcSD.First12hrs.Wake  = [];
	Twelve_hour_avg_percentages_Female.Mut.SDexcSD.First12hrs.Wake = [];
	Twelve_hour_avg_percentages_Female.WT.SDexcSD.First12hrs.NREM  = [];
	Twelve_hour_avg_percentages_Female.Mut.SDexcSD.First12hrs.NREM = [];
	Twelve_hour_avg_percentages_Female.WT.SDexcSD.First12hrs.REM   = [];
	Twelve_hour_avg_percentages_Female.Mut.SDexcSD.First12hrs.REM  = [];
	Twelve_hour_avg_percentages_Female.WT.SDexcSD.Last12hrs.Wake   = [];
	Twelve_hour_avg_percentages_Female.Mut.SDexcSD.Last12hrs.Wake  = [];
	Twelve_hour_avg_percentages_Female.WT.SDexcSD.Last12hrs.NREM   = [];
	Twelve_hour_avg_percentages_Female.Mut.SDexcSD.Last12hrs.NREM  = [];
	Twelve_hour_avg_percentages_Female.WT.SDexcSD.Last12hrs.REM    = [];
	Twelve_hour_avg_percentages_Female.Mut.SDexcSD.Last12hrs.REM   = [];

end 


% Set up the matrices so the calls below are less messy
% BL LP Wake
Female_WT_BL_LP_Wake  = Twelve_hour_avg_percentages_Female.WT.BL.First12hrs.Wake;
Female_Mut_BL_LP_Wake = Twelve_hour_avg_percentages_Female.Mut.BL.First12hrs.Wake;
Male_WT_BL_LP_Wake    = Twelve_hour_avg_percentages_Male.WT.BL.First12hrs.Wake;
Male_Mut_BL_LP_Wake   = Twelve_hour_avg_percentages_Male.Mut.BL.First12hrs.Wake;

% BL LP NREM
Female_WT_BL_LP_NREM  = Twelve_hour_avg_percentages_Female.WT.BL.First12hrs.NREM;
Female_Mut_BL_LP_NREM = Twelve_hour_avg_percentages_Female.Mut.BL.First12hrs.NREM;
Male_WT_BL_LP_NREM    = Twelve_hour_avg_percentages_Male.WT.BL.First12hrs.NREM;
Male_Mut_BL_LP_NREM   = Twelve_hour_avg_percentages_Male.Mut.BL.First12hrs.NREM;

% BL LP REM
Female_WT_BL_LP_REM  = Twelve_hour_avg_percentages_Female.WT.BL.First12hrs.REM;
Female_Mut_BL_LP_REM = Twelve_hour_avg_percentages_Female.Mut.BL.First12hrs.REM;
Male_WT_BL_LP_REM    = Twelve_hour_avg_percentages_Male.WT.BL.First12hrs.REM;
Male_Mut_BL_LP_REM   = Twelve_hour_avg_percentages_Male.Mut.BL.First12hrs.REM;


% BL DP Wake
Female_WT_BL_DP_Wake  = Twelve_hour_avg_percentages_Female.WT.BL.Last12hrs.Wake;
Female_Mut_BL_DP_Wake = Twelve_hour_avg_percentages_Female.Mut.BL.Last12hrs.Wake;
Male_WT_BL_DP_Wake    = Twelve_hour_avg_percentages_Male.WT.BL.Last12hrs.Wake;
Male_Mut_BL_DP_Wake   = Twelve_hour_avg_percentages_Male.Mut.BL.Last12hrs.Wake;

% BL DP NREM
Female_WT_BL_DP_NREM  = Twelve_hour_avg_percentages_Female.WT.BL.Last12hrs.NREM;
Female_Mut_BL_DP_NREM = Twelve_hour_avg_percentages_Female.Mut.BL.Last12hrs.NREM;
Male_WT_BL_DP_NREM    = Twelve_hour_avg_percentages_Male.WT.BL.Last12hrs.NREM;
Male_Mut_BL_DP_NREM   = Twelve_hour_avg_percentages_Male.Mut.BL.Last12hrs.NREM;

% BL DP REM
Female_WT_BL_DP_REM  = Twelve_hour_avg_percentages_Female.WT.BL.Last12hrs.REM;
Female_Mut_BL_DP_REM = Twelve_hour_avg_percentages_Female.Mut.BL.Last12hrs.REM;
Male_WT_BL_DP_REM    = Twelve_hour_avg_percentages_Male.WT.BL.Last12hrs.REM;
Male_Mut_BL_DP_REM   = Twelve_hour_avg_percentages_Male.Mut.BL.Last12hrs.REM;

% Sleep Dep
% SD LP Wake 
Female_WT_SD_LP_Wake  = Twelve_hour_avg_percentages_Female.WT.SD.First12hrs.Wake;
Female_Mut_SD_LP_Wake = Twelve_hour_avg_percentages_Female.Mut.SD.First12hrs.Wake;
Male_WT_SD_LP_Wake    = Twelve_hour_avg_percentages_Male.WT.SD.First12hrs.Wake;
Male_Mut_SD_LP_Wake   = Twelve_hour_avg_percentages_Male.Mut.SD.First12hrs.Wake;

% SD LP NREM
Female_WT_SD_LP_NREM  = Twelve_hour_avg_percentages_Female.WT.SD.First12hrs.NREM;
Female_Mut_SD_LP_NREM = Twelve_hour_avg_percentages_Female.Mut.SD.First12hrs.NREM;
Male_WT_SD_LP_NREM    = Twelve_hour_avg_percentages_Male.WT.SD.First12hrs.NREM;
Male_Mut_SD_LP_NREM   = Twelve_hour_avg_percentages_Male.Mut.SD.First12hrs.NREM;

% SD LP REM
Female_WT_SD_LP_REM  = Twelve_hour_avg_percentages_Female.WT.SD.First12hrs.REM;
Female_Mut_SD_LP_REM = Twelve_hour_avg_percentages_Female.Mut.SD.First12hrs.REM;
Male_WT_SD_LP_REM    = Twelve_hour_avg_percentages_Male.WT.SD.First12hrs.REM;
Male_Mut_SD_LP_REM   = Twelve_hour_avg_percentages_Male.Mut.SD.First12hrs.REM;

% SDexcSD LP Wake 
Female_WT_SDexcSD_LP_Wake  = Twelve_hour_avg_percentages_Female.WT.SDexcSD.First12hrs.Wake;
Female_Mut_SDexcSD_LP_Wake = Twelve_hour_avg_percentages_Female.Mut.SDexcSD.First12hrs.Wake;
Male_WT_SDexcSD_LP_Wake    = Twelve_hour_avg_percentages_Male.WT.SDexcSD.First12hrs.Wake;
Male_Mut_SDexcSD_LP_Wake   = Twelve_hour_avg_percentages_Male.Mut.SDexcSD.First12hrs.Wake;

% SDexcSD LP NREM
Female_WT_SDexcSD_LP_NREM  = Twelve_hour_avg_percentages_Female.WT.SDexcSD.First12hrs.NREM;
Female_Mut_SDexcSD_LP_NREM = Twelve_hour_avg_percentages_Female.Mut.SDexcSD.First12hrs.NREM;
Male_WT_SDexcSD_LP_NREM    = Twelve_hour_avg_percentages_Male.WT.SDexcSD.First12hrs.NREM;
Male_Mut_SDexcSD_LP_NREM   = Twelve_hour_avg_percentages_Male.Mut.SDexcSD.First12hrs.NREM;

% SDexcSD LP REM
Female_WT_SDexcSD_LP_REM  = Twelve_hour_avg_percentages_Female.WT.SDexcSD.First12hrs.REM;
Female_Mut_SDexcSD_LP_REM = Twelve_hour_avg_percentages_Female.Mut.SDexcSD.First12hrs.REM;
Male_WT_SDexcSD_LP_REM    = Twelve_hour_avg_percentages_Male.WT.SDexcSD.First12hrs.REM;
Male_Mut_SDexcSD_LP_REM   = Twelve_hour_avg_percentages_Male.Mut.SDexcSD.First12hrs.REM;

% SD DP Wake
Female_WT_SD_DP_Wake  = Twelve_hour_avg_percentages_Female.WT.SD.Last12hrs.Wake;
Female_Mut_SD_DP_Wake = Twelve_hour_avg_percentages_Female.Mut.SD.Last12hrs.Wake;
Male_WT_SD_DP_Wake    = Twelve_hour_avg_percentages_Male.WT.SD.Last12hrs.Wake;
Male_Mut_SD_DP_Wake   = Twelve_hour_avg_percentages_Male.Mut.SD.Last12hrs.Wake;

% SD DP NREM
Female_WT_SD_DP_NREM  = Twelve_hour_avg_percentages_Female.WT.SD.Last12hrs.NREM;
Female_Mut_SD_DP_NREM = Twelve_hour_avg_percentages_Female.Mut.SD.Last12hrs.NREM;
Male_WT_SD_DP_NREM    = Twelve_hour_avg_percentages_Male.WT.SD.Last12hrs.NREM;
Male_Mut_SD_DP_NREM   = Twelve_hour_avg_percentages_Male.Mut.SD.Last12hrs.NREM;

% SD DP REM
Female_WT_SD_DP_REM  = Twelve_hour_avg_percentages_Female.WT.SD.Last12hrs.REM;
Female_Mut_SD_DP_REM = Twelve_hour_avg_percentages_Female.Mut.SD.Last12hrs.REM;
Male_WT_SD_DP_REM    = Twelve_hour_avg_percentages_Male.WT.SD.Last12hrs.REM;
Male_Mut_SD_DP_REM   = Twelve_hour_avg_percentages_Male.Mut.SD.Last12hrs.REM;


% --------------------------------------------------------------------------------
% ---- Done setting up matrices --------------------------------------------------

% ------------------------ Baseline ---------------------------------------
% ------- Light Phase --------------
% Wake 
[p_vals_BL_LP_Wake,sig_difference_BL_LP_Wake,ANOVA_tbl_BL_LP_Wake] = perform_2way_Anova_SexGenotype_posthoc(Female_WT_BL_LP_Wake,Female_Mut_BL_LP_Wake,Male_WT_BL_LP_Wake,Male_Mut_BL_LP_Wake);
% NREM
[p_vals_BL_LP_NREM,sig_difference_BL_LP_NREM,ANOVA_tbl_BL_LP_NREM] = perform_2way_Anova_SexGenotype_posthoc(Female_WT_BL_LP_NREM,Female_Mut_BL_LP_NREM,Male_WT_BL_LP_NREM,Male_Mut_BL_LP_NREM);
% REM 
[p_vals_BL_LP_REM,sig_difference_BL_LP_REM,ANOVA_tbl_BL_LP_REM]   = perform_2way_Anova_SexGenotype_posthoc(Female_WT_BL_LP_REM,Female_Mut_BL_LP_REM,Male_WT_BL_LP_REM,Male_Mut_BL_LP_REM);

% ------- Dark Phase --------------
% Wake 
[p_vals_BL_DP_Wake,sig_difference_BL_DP_Wake,ANOVA_tbl_BL_DP_Wake] = perform_2way_Anova_SexGenotype_posthoc(Female_WT_BL_DP_Wake,Female_Mut_BL_DP_Wake,Male_WT_BL_DP_Wake,Male_Mut_BL_DP_Wake);
% NREM
[p_vals_BL_DP_NREM,sig_difference_BL_DP_NREM,ANOVA_tbl_BL_DP_NREM] = perform_2way_Anova_SexGenotype_posthoc(Female_WT_BL_DP_NREM,Female_Mut_BL_DP_NREM,Male_WT_BL_DP_NREM,Male_Mut_BL_DP_NREM);
% REM 
[p_vals_BL_DP_REM,sig_difference_BL_DP_REM,ANOVA_tbl_BL_DP_REM]    = perform_2way_Anova_SexGenotype_posthoc(Female_WT_BL_DP_REM,Female_Mut_BL_DP_REM,Male_WT_BL_DP_REM,Male_Mut_BL_DP_REM);


% ------------------------ Sleep Dep ---------------------------------------
% ------- Light Phase --------------
% % Wake 
% [p_vals_SD_LP_Wake,sig_difference_SD_LP_Wake] = perform_2way_Anova_SexGenotype_posthoc(Female_WT_SD_LP_Wake,Female_Mut_SD_LP_Wake,Male_WT_SD_LP_Wake,Male_Mut_SD_LP_Wake);
% % NREM
% [p_vals_SD_LP_NREM,sig_difference_SD_LP_NREM] = perform_2way_Anova_SexGenotype_posthoc(Female_WT_SD_LP_NREM,Female_Mut_SD_LP_NREM,Male_WT_SD_LP_NREM,Male_Mut_SD_LP_NREM);
% % REM 
% [p_vals_SD_LP_REM,sig_difference_SD_LP_REM]   = perform_2way_Anova_SexGenotype_posthoc(Female_WT_SD_LP_REM,Female_Mut_SD_LP_REM,Male_WT_SD_LP_REM,Male_Mut_SD_LP_REM);

% ------- Dark Phase --------------
% Wake 
[p_vals_SD_DP_Wake,sig_difference_SD_DP_Wake,ANOVA_tbl_SD_DP_Wake] = perform_2way_Anova_SexGenotype_posthoc(Female_WT_SD_DP_Wake,Female_Mut_SD_DP_Wake,Male_WT_SD_DP_Wake,Male_Mut_SD_DP_Wake);
% NREM
[p_vals_SD_DP_NREM,sig_difference_SD_DP_NREM,ANOVA_tbl_SD_DP_NREM] = perform_2way_Anova_SexGenotype_posthoc(Female_WT_SD_DP_NREM,Female_Mut_SD_DP_NREM,Male_WT_SD_DP_NREM,Male_Mut_SD_DP_NREM);
% REM 
[p_vals_SD_DP_REM,sig_difference_SD_DP_REM,ANOVA_tbl_SD_DP_REM]   = perform_2way_Anova_SexGenotype_posthoc(Female_WT_SD_DP_REM,Female_Mut_SD_DP_REM,Male_WT_SD_DP_REM,Male_Mut_SD_DP_REM);

% ------- Light Phase (Excluding the actual sleep dep) --------------
% Wake 
[p_vals_SDexcSD_LP_Wake,sig_difference_SDexcSD_LP_Wake,ANOVA_tbl_SDexcSD_LP_Wake] = perform_2way_Anova_SexGenotype_posthoc(Female_WT_SDexcSD_LP_Wake,Female_Mut_SDexcSD_LP_Wake,Male_WT_SDexcSD_LP_Wake,Male_Mut_SDexcSD_LP_Wake);
% NREM
[p_vals_SDexcSD_LP_NREM,sig_difference_SDexcSD_LP_NREM,ANOVA_tbl_SDexcSD_LP_NREM] = perform_2way_Anova_SexGenotype_posthoc(Female_WT_SDexcSD_LP_NREM,Female_Mut_SDexcSD_LP_NREM,Male_WT_SDexcSD_LP_NREM,Male_Mut_SDexcSD_LP_NREM);
% REM 
[p_vals_SDexcSD_LP_REM,sig_difference_SDexcSD_LP_REM,ANOVA_tbl_SDexcSD_LP_REM]   = perform_2way_Anova_SexGenotype_posthoc(Female_WT_SDexcSD_LP_REM,Female_Mut_SDexcSD_LP_REM,Male_WT_SDexcSD_LP_REM,Male_Mut_SDexcSD_LP_REM);



% If you exclude the actual sleep dep, you don't get light phase and dark phase

% These are organized by panels of the figure (i.e. BL Light Period Wake, BL Dark Period NREM, etc)
% -- BL Male --
p_vals.Male.BL.LP.Wake = p_vals_BL_LP_Wake.Male;
p_vals.Male.BL.LP.NREM = p_vals_BL_LP_NREM.Male;
p_vals.Male.BL.LP.REM  = p_vals_BL_LP_REM.Male;
p_vals.Male.BL.DP.Wake = p_vals_BL_DP_Wake.Male;
p_vals.Male.BL.DP.NREM = p_vals_BL_DP_NREM.Male;
p_vals.Male.BL.DP.REM  = p_vals_BL_DP_REM.Male;

% -- BL Female --
p_vals.Female.BL.LP.Wake = p_vals_BL_LP_Wake.Female;
p_vals.Female.BL.LP.NREM = p_vals_BL_LP_NREM.Female;
p_vals.Female.BL.LP.REM  = p_vals_BL_LP_REM.Female;
p_vals.Female.BL.DP.Wake = p_vals_BL_DP_Wake.Female;
p_vals.Female.BL.DP.NREM = p_vals_BL_DP_NREM.Female;
p_vals.Female.BL.DP.REM  = p_vals_BL_DP_REM.Female;

% -- BL Both --
p_vals.BothMandF.BL.LP.Wake = p_vals_BL_LP_Wake.BothMandF;
p_vals.BothMandF.BL.LP.NREM = p_vals_BL_LP_NREM.BothMandF;
p_vals.BothMandF.BL.LP.REM  = p_vals_BL_LP_REM.BothMandF;
p_vals.BothMandF.BL.DP.Wake = p_vals_BL_DP_Wake.BothMandF;
p_vals.BothMandF.BL.DP.NREM = p_vals_BL_DP_NREM.BothMandF;
p_vals.BothMandF.BL.DP.REM  = p_vals_BL_DP_REM.BothMandF;


% Excluded because we never want this.  Removing to avoid confusion.  
% p_vals.SD.LP.Wake = p_vals_SD_LP_Wake;
% p_vals.SD.LP.NREM = p_vals_SD_LP_NREM;
% p_vals.SD.LP.REM  = p_vals_SD_LP_REM;

% -- SD Male -- 
p_vals.Male.SDexcSD.LP.Wake = p_vals_SDexcSD_LP_Wake.Male;
p_vals.Male.SDexcSD.LP.NREM = p_vals_SDexcSD_LP_NREM.Male;
p_vals.Male.SDexcSD.LP.REM  = p_vals_SDexcSD_LP_REM.Male;
p_vals.Male.SD.DP.Wake 		= p_vals_SD_DP_Wake.Male;
p_vals.Male.SD.DP.NREM 		= p_vals_SD_DP_NREM.Male;
p_vals.Male.SD.DP.REM  		= p_vals_SD_DP_REM.Male;

% -- SD Female -- 
p_vals.Female.SDexcSD.LP.Wake = p_vals_SDexcSD_LP_Wake.Female;
p_vals.Female.SDexcSD.LP.NREM = p_vals_SDexcSD_LP_NREM.Female;
p_vals.Female.SDexcSD.LP.REM  = p_vals_SDexcSD_LP_REM.Female;
p_vals.Female.SD.DP.Wake 	  = p_vals_SD_DP_Wake.Female;
p_vals.Female.SD.DP.NREM 	  = p_vals_SD_DP_NREM.Female;
p_vals.Female.SD.DP.REM  	  = p_vals_SD_DP_REM.Female;

% -- SD Both -- 
p_vals.BothMandF.SDexcSD.LP.Wake = p_vals_SDexcSD_LP_Wake.BothMandF;
p_vals.BothMandF.SDexcSD.LP.NREM = p_vals_SDexcSD_LP_NREM.BothMandF;
p_vals.BothMandF.SDexcSD.LP.REM  = p_vals_SDexcSD_LP_REM.BothMandF;
p_vals.BothMandF.SD.DP.Wake 	  = p_vals_SD_DP_Wake.BothMandF;
p_vals.BothMandF.SD.DP.NREM 	  = p_vals_SD_DP_NREM.BothMandF;
p_vals.BothMandF.SD.DP.REM  	  = p_vals_SD_DP_REM.BothMandF;

% -- post-hocs for M vs F WT --
% - BL LP -
p_vals.BothMandF.BL.LP.Wake.Posthocs.WT_MaleVsFemale = p_vals_BL_LP_Wake.Posthoc.WT.MalevsFemale;
p_vals.BothMandF.BL.LP.NREM.Posthocs.WT_MaleVsFemale = p_vals_BL_LP_NREM.Posthoc.WT.MalevsFemale;
p_vals.BothMandF.BL.LP.REM.Posthocs.WT_MaleVsFemale  = p_vals_BL_LP_REM.Posthoc.WT.MalevsFemale;

% - BL DP -
p_vals.BothMandF.BL.DP.Wake.Posthocs.WT_MaleVsFemale = p_vals_BL_DP_Wake.Posthoc.WT.MalevsFemale;
p_vals.BothMandF.BL.DP.NREM.Posthocs.WT_MaleVsFemale = p_vals_BL_DP_NREM.Posthoc.WT.MalevsFemale;
p_vals.BothMandF.BL.DP.REM.Posthocs.WT_MaleVsFemale  = p_vals_BL_DP_REM.Posthoc.WT.MalevsFemale;

% - SD LP -
p_vals.BothMandF.SDexcSD.LP.Wake.Posthocs.WT_MaleVsFemale = p_vals_SDexcSD_LP_Wake.Posthoc.WT.MalevsFemale;
p_vals.BothMandF.SDexcSD.LP.NREM.Posthocs.WT_MaleVsFemale = p_vals_SDexcSD_LP_NREM.Posthoc.WT.MalevsFemale;
p_vals.BothMandF.SDexcSD.LP.REM.Posthocs.WT_MaleVsFemale  = p_vals_SDexcSD_LP_REM.Posthoc.WT.MalevsFemale;

% - SD DP -
p_vals.BothMandF.SD.DP.Wake.Posthocs.WT_MaleVsFemale = p_vals_SD_DP_Wake.Posthoc.WT.MalevsFemale;
p_vals.BothMandF.SD.DP.NREM.Posthocs.WT_MaleVsFemale = p_vals_SD_DP_NREM.Posthoc.WT.MalevsFemale;
p_vals.BothMandF.SD.DP.REM.Posthocs.WT_MaleVsFemale  = p_vals_SD_DP_REM.Posthoc.WT.MalevsFemale;

% -- post-hocs for M vs F Mut --
% - BL LP -
p_vals.BothMandF.BL.LP.Wake.Posthocs.Mut_MaleVsFemale = p_vals_BL_LP_Wake.Posthoc.Mut.MalevsFemale;
p_vals.BothMandF.BL.LP.NREM.Posthocs.Mut_MaleVsFemale = p_vals_BL_LP_NREM.Posthoc.Mut.MalevsFemale;
p_vals.BothMandF.BL.LP.REM.Posthocs.Mut_MaleVsFemale  = p_vals_BL_LP_REM.Posthoc.Mut.MalevsFemale;

% - BL DP -
p_vals.BothMandF.BL.DP.Wake.Posthocs.Mut_MaleVsFemale = p_vals_BL_DP_Wake.Posthoc.Mut.MalevsFemale;
p_vals.BothMandF.BL.DP.NREM.Posthocs.Mut_MaleVsFemale = p_vals_BL_DP_NREM.Posthoc.Mut.MalevsFemale;
p_vals.BothMandF.BL.DP.REM.Posthocs.Mut_MaleVsFemale  = p_vals_BL_DP_REM.Posthoc.Mut.MalevsFemale;

% - SD LP -
p_vals.BothMandF.SDexcSD.LP.Wake.Posthocs.Mut_MaleVsFemale = p_vals_SDexcSD_LP_Wake.Posthoc.Mut.MalevsFemale;
p_vals.BothMandF.SDexcSD.LP.NREM.Posthocs.Mut_MaleVsFemale = p_vals_SDexcSD_LP_NREM.Posthoc.Mut.MalevsFemale;
p_vals.BothMandF.SDexcSD.LP.REM.Posthocs.Mut_MaleVsFemale  = p_vals_SDexcSD_LP_REM.Posthoc.Mut.MalevsFemale;

% - SD DP -
p_vals.BothMandF.SD.DP.Wake.Posthocs.Mut_MaleVsFemale = p_vals_SD_DP_Wake.Posthoc.Mut.MalevsFemale;
p_vals.BothMandF.SD.DP.NREM.Posthocs.Mut_MaleVsFemale = p_vals_SD_DP_NREM.Posthoc.Mut.MalevsFemale;
p_vals.BothMandF.SD.DP.REM.Posthocs.Mut_MaleVsFemale  = p_vals_SD_DP_REM.Posthoc.Mut.MalevsFemale;


% -- BL Male --
sig_difference_bool.Male.BL.LP.Wake = sig_difference_BL_LP_Wake.Male;
sig_difference_bool.Male.BL.LP.NREM = sig_difference_BL_LP_NREM.Male;
sig_difference_bool.Male.BL.LP.REM  = sig_difference_BL_LP_REM.Male;
sig_difference_bool.Male.BL.DP.Wake = sig_difference_BL_DP_Wake.Male;
sig_difference_bool.Male.BL.DP.NREM = sig_difference_BL_DP_NREM.Male;
sig_difference_bool.Male.BL.DP.REM  = sig_difference_BL_DP_REM.Male;

% -- BL Female --
sig_difference_bool.Female.BL.LP.Wake = sig_difference_BL_LP_Wake.Female;
sig_difference_bool.Female.BL.LP.NREM = sig_difference_BL_LP_NREM.Female;
sig_difference_bool.Female.BL.LP.REM  = sig_difference_BL_LP_REM.Female;
sig_difference_bool.Female.BL.DP.Wake = sig_difference_BL_DP_Wake.Female;
sig_difference_bool.Female.BL.DP.NREM = sig_difference_BL_DP_NREM.Female;
sig_difference_bool.Female.BL.DP.REM  = sig_difference_BL_DP_REM.Female;

% -- SD Male -- 
sig_difference_bool.Male.SDexcSD.LP.Wake = sig_difference_SDexcSD_LP_Wake.Male;
sig_difference_bool.Male.SDexcSD.LP.NREM = sig_difference_SDexcSD_LP_NREM.Male;
sig_difference_bool.Male.SDexcSD.LP.REM  = sig_difference_SDexcSD_LP_REM.Male;
sig_difference_bool.Male.SD.DP.Wake = sig_difference_SD_DP_Wake.Male;
sig_difference_bool.Male.SD.DP.NREM = sig_difference_SD_DP_NREM.Male;
sig_difference_bool.Male.SD.DP.REM  = sig_difference_SD_DP_REM.Male;

% -- SD Female -- 
sig_difference_bool.Female.SDexcSD.LP.Wake = sig_difference_SDexcSD_LP_Wake.Female;
sig_difference_bool.Female.SDexcSD.LP.NREM = sig_difference_SDexcSD_LP_NREM.Female;
sig_difference_bool.Female.SDexcSD.LP.REM  = sig_difference_SDexcSD_LP_REM.Female;
sig_difference_bool.Female.SD.DP.Wake = sig_difference_SD_DP_Wake.Female;
sig_difference_bool.Female.SD.DP.NREM = sig_difference_SD_DP_NREM.Female;
sig_difference_bool.Female.SD.DP.REM  = sig_difference_SD_DP_REM.Female;

% -- BL Male --
ANOVA_tbls.Male.BL.LP.Wake = ANOVA_tbl_BL_LP_Wake.Male;
ANOVA_tbls.Male.BL.LP.NREM = ANOVA_tbl_BL_LP_NREM.Male;
ANOVA_tbls.Male.BL.LP.REM  = ANOVA_tbl_BL_LP_REM.Male;
ANOVA_tbls.Male.BL.DP.Wake = ANOVA_tbl_BL_DP_Wake.Male;
ANOVA_tbls.Male.BL.DP.NREM = ANOVA_tbl_BL_DP_NREM.Male;
ANOVA_tbls.Male.BL.DP.REM  = ANOVA_tbl_BL_DP_REM.Male;

% -- SD Male -- 
ANOVA_tbls.Male.SDexcSD.LP.Wake = ANOVA_tbl_SDexcSD_LP_Wake.Male;
ANOVA_tbls.Male.SDexcSD.LP.NREM = ANOVA_tbl_SDexcSD_LP_NREM.Male;
ANOVA_tbls.Male.SDexcSD.LP.REM  = ANOVA_tbl_SDexcSD_LP_REM.Male;
ANOVA_tbls.Male.SD.DP.Wake      = ANOVA_tbl_SD_DP_Wake.Male;
ANOVA_tbls.Male.SD.DP.NREM      = ANOVA_tbl_SD_DP_NREM.Male;
ANOVA_tbls.Male.SD.DP.REM       = ANOVA_tbl_SD_DP_REM.Male;

% -- BL Female --
ANOVA_tbls.Female.BL.LP.Wake = ANOVA_tbl_BL_LP_Wake.Female;
ANOVA_tbls.Female.BL.LP.NREM = ANOVA_tbl_BL_LP_NREM.Female;
ANOVA_tbls.Female.BL.LP.REM  = ANOVA_tbl_BL_LP_REM.Female;
ANOVA_tbls.Female.BL.DP.Wake = ANOVA_tbl_BL_DP_Wake.Female;
ANOVA_tbls.Female.BL.DP.NREM = ANOVA_tbl_BL_DP_NREM.Female;
ANOVA_tbls.Female.BL.DP.REM  = ANOVA_tbl_BL_DP_REM.Female;

% -- SD Female -- 
ANOVA_tbls.Female.SDexcSD.LP.Wake = ANOVA_tbl_SDexcSD_LP_Wake.Female;
ANOVA_tbls.Female.SDexcSD.LP.NREM = ANOVA_tbl_SDexcSD_LP_NREM.Female;
ANOVA_tbls.Female.SDexcSD.LP.REM  = ANOVA_tbl_SDexcSD_LP_REM.Female;
ANOVA_tbls.Female.SD.DP.Wake      = ANOVA_tbl_SD_DP_Wake.Female;
ANOVA_tbls.Female.SD.DP.NREM      = ANOVA_tbl_SD_DP_NREM.Female;
ANOVA_tbls.Female.SD.DP.REM       = ANOVA_tbl_SD_DP_REM.Female;

















































