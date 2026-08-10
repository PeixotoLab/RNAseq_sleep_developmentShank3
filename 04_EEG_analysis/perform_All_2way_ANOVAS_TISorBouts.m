function [p_vals,sig_difference_bool,ANOVA_tbls] = perform_All_2way_ANOVAS_TISorBouts(Twelve_hour_avg_percentages_Male,Twelve_hour_avg_percentages_Female,WindowLength,Sexes,SexVarName);
%
% This function simply calls perform_2way_Anova_SexGenotype_posthoc.m for Wake, NREM, REM in light phase dark phase and BL and sleep dep
%
% INPUTS: 	Twelve_hour_avg_percentages_Male 	a struct with fields WT.BL.First12hrs for instance
% 			same for Female
% 			WindowLength 						do you want to analyze in 12-hour segments or 6-hour (only for dark period)? 

if WindowLength==12
	% If you don't have female data, set up those structs as empty
	if isempty(Twelve_hour_avg_percentages_Female)
		Female_Data_Present = false; 
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
	else 
		Female_Data_Present = true;
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
	[p_vals_BL_LP_Wake,sig_difference_BL_LP_Wake,ANOVA_tbl_BL_LP_Wake] = perform_2way_Anova_SexGenotype_posthoc(Female_WT_BL_LP_Wake,Female_Mut_BL_LP_Wake,Male_WT_BL_LP_Wake,Male_Mut_BL_LP_Wake,Sexes,SexVarName);
	% NREM
	[p_vals_BL_LP_NREM,sig_difference_BL_LP_NREM,ANOVA_tbl_BL_LP_NREM] = perform_2way_Anova_SexGenotype_posthoc(Female_WT_BL_LP_NREM,Female_Mut_BL_LP_NREM,Male_WT_BL_LP_NREM,Male_Mut_BL_LP_NREM,Sexes,SexVarName);
	% REM 
	[p_vals_BL_LP_REM,sig_difference_BL_LP_REM,ANOVA_tbl_BL_LP_REM]   = perform_2way_Anova_SexGenotype_posthoc(Female_WT_BL_LP_REM,Female_Mut_BL_LP_REM,Male_WT_BL_LP_REM,Male_Mut_BL_LP_REM,Sexes,SexVarName);

	% ------- Dark Phase --------------
	% Wake 
	[p_vals_BL_DP_Wake,sig_difference_BL_DP_Wake,ANOVA_tbl_BL_DP_Wake] = perform_2way_Anova_SexGenotype_posthoc(Female_WT_BL_DP_Wake,Female_Mut_BL_DP_Wake,Male_WT_BL_DP_Wake,Male_Mut_BL_DP_Wake,Sexes,SexVarName);
	% NREM
	[p_vals_BL_DP_NREM,sig_difference_BL_DP_NREM,ANOVA_tbl_BL_DP_NREM] = perform_2way_Anova_SexGenotype_posthoc(Female_WT_BL_DP_NREM,Female_Mut_BL_DP_NREM,Male_WT_BL_DP_NREM,Male_Mut_BL_DP_NREM,Sexes,SexVarName);
	% REM 
	[p_vals_BL_DP_REM,sig_difference_BL_DP_REM,ANOVA_tbl_BL_DP_REM]    = perform_2way_Anova_SexGenotype_posthoc(Female_WT_BL_DP_REM,Female_Mut_BL_DP_REM,Male_WT_BL_DP_REM,Male_Mut_BL_DP_REM,Sexes,SexVarName);


	% ------------------------ Sleep Dep ---------------------------------------
	% ------- Light Phase --------------
	% % Wake 
	% [p_vals_SD_LP_Wake,sig_difference_SD_LP_Wake] = perform_2way_Anova_SexGenotype_posthoc(Female_WT_SD_LP_Wake,Female_Mut_SD_LP_Wake,Male_WT_SD_LP_Wake,Male_Mut_SD_LP_Wake,Sexes,SexVarName);
	% % NREM
	% [p_vals_SD_LP_NREM,sig_difference_SD_LP_NREM] = perform_2way_Anova_SexGenotype_posthoc(Female_WT_SD_LP_NREM,Female_Mut_SD_LP_NREM,Male_WT_SD_LP_NREM,Male_Mut_SD_LP_NREM,Sexes,SexVarName);
	% % REM 
	% [p_vals_SD_LP_REM,sig_difference_SD_LP_REM]   = perform_2way_Anova_SexGenotype_posthoc(Female_WT_SD_LP_REM,Female_Mut_SD_LP_REM,Male_WT_SD_LP_REM,Male_Mut_SD_LP_REM,Sexes,SexVarName);

	% ------- Dark Phase --------------
	% Wake 
	[p_vals_SD_DP_Wake,sig_difference_SD_DP_Wake,ANOVA_tbl_SD_DP_Wake] = perform_2way_Anova_SexGenotype_posthoc(Female_WT_SD_DP_Wake,Female_Mut_SD_DP_Wake,Male_WT_SD_DP_Wake,Male_Mut_SD_DP_Wake,Sexes,SexVarName);
	% NREM
	[p_vals_SD_DP_NREM,sig_difference_SD_DP_NREM,ANOVA_tbl_SD_DP_NREM] = perform_2way_Anova_SexGenotype_posthoc(Female_WT_SD_DP_NREM,Female_Mut_SD_DP_NREM,Male_WT_SD_DP_NREM,Male_Mut_SD_DP_NREM,Sexes,SexVarName);
	% REM 
	[p_vals_SD_DP_REM,sig_difference_SD_DP_REM,ANOVA_tbl_SD_DP_REM]   = perform_2way_Anova_SexGenotype_posthoc(Female_WT_SD_DP_REM,Female_Mut_SD_DP_REM,Male_WT_SD_DP_REM,Male_Mut_SD_DP_REM,Sexes,SexVarName);

	% ------- Light Phase (Excluding the actual sleep dep) --------------
	% Wake 
	[p_vals_SDexcSD_LP_Wake,sig_difference_SDexcSD_LP_Wake,ANOVA_tbl_SDexcSD_LP_Wake] = perform_2way_Anova_SexGenotype_posthoc(Female_WT_SDexcSD_LP_Wake,Female_Mut_SDexcSD_LP_Wake,Male_WT_SDexcSD_LP_Wake,Male_Mut_SDexcSD_LP_Wake,Sexes,SexVarName);
	% NREM
	[p_vals_SDexcSD_LP_NREM,sig_difference_SDexcSD_LP_NREM,ANOVA_tbl_SDexcSD_LP_NREM] = perform_2way_Anova_SexGenotype_posthoc(Female_WT_SDexcSD_LP_NREM,Female_Mut_SDexcSD_LP_NREM,Male_WT_SDexcSD_LP_NREM,Male_Mut_SDexcSD_LP_NREM,Sexes,SexVarName);
	% REM 
	[p_vals_SDexcSD_LP_REM,sig_difference_SDexcSD_LP_REM,ANOVA_tbl_SDexcSD_LP_REM]   = perform_2way_Anova_SexGenotype_posthoc(Female_WT_SDexcSD_LP_REM,Female_Mut_SDexcSD_LP_REM,Male_WT_SDexcSD_LP_REM,Male_Mut_SDexcSD_LP_REM,Sexes,SexVarName);



	% If you exclude the actual sleep dep, you don't get light phase and dark phase

	% These are organized by panels of the figure (i.e. BL Light Period Wake, BL Dark Period NREM, etc)
	% -- BL Male --
	if Female_Data_Present p_vals.(Sexes{1}).BL.LP.Wake = p_vals_BL_LP_Wake.Posthoc.(Sexes{1}).WTvsMut; else p_vals.(Sexes{1}).BL.LP.Wake = p_vals_BL_LP_Wake.(Sexes{1}); end
	if Female_Data_Present p_vals.(Sexes{1}).BL.LP.NREM = p_vals_BL_LP_NREM.Posthoc.(Sexes{1}).WTvsMut; else p_vals.(Sexes{1}).BL.LP.NREM = p_vals_BL_LP_NREM.(Sexes{1}); end 
	if Female_Data_Present p_vals.(Sexes{1}).BL.LP.REM  = p_vals_BL_LP_REM.Posthoc.(Sexes{1}).WTvsMut;  else p_vals.(Sexes{1}).BL.LP.REM  = p_vals_BL_LP_REM.(Sexes{1});  end
	if Female_Data_Present p_vals.(Sexes{1}).BL.DP.Wake = p_vals_BL_DP_Wake.Posthoc.(Sexes{1}).WTvsMut; else p_vals.(Sexes{1}).BL.DP.Wake = p_vals_BL_DP_Wake.(Sexes{1}); end
	if Female_Data_Present p_vals.(Sexes{1}).BL.DP.NREM = p_vals_BL_DP_NREM.Posthoc.(Sexes{1}).WTvsMut; else p_vals.(Sexes{1}).BL.DP.NREM = p_vals_BL_DP_NREM.(Sexes{1}); end
	if Female_Data_Present p_vals.(Sexes{1}).BL.DP.REM  = p_vals_BL_DP_REM.Posthoc.(Sexes{1}).WTvsMut;  else p_vals.(Sexes{1}).BL.DP.REM  = p_vals_BL_DP_REM.(Sexes{1});  end

	% -- BL Female --
	p_vals.(Sexes{2}).BL.LP.Wake = p_vals_BL_LP_Wake.Posthoc.(Sexes{2});
	p_vals.(Sexes{2}).BL.LP.NREM = p_vals_BL_LP_NREM.Posthoc.(Sexes{2});
	p_vals.(Sexes{2}).BL.LP.REM  = p_vals_BL_LP_REM.Posthoc.(Sexes{2});
	p_vals.(Sexes{2}).BL.DP.Wake = p_vals_BL_DP_Wake.Posthoc.(Sexes{2});
	p_vals.(Sexes{2}).BL.DP.NREM = p_vals_BL_DP_NREM.Posthoc.(Sexes{2});
	p_vals.(Sexes{2}).BL.DP.REM  = p_vals_BL_DP_REM.Posthoc.(Sexes{2});

	% -- BL Both --
	p_vals.(['Both' Sexes{1} 'and' Sexes{2}]).BL.LP.Wake = p_vals_BL_LP_Wake.(['Both' Sexes{1} 'and' Sexes{2}]);
	p_vals.(['Both' Sexes{1} 'and' Sexes{2}]).BL.LP.NREM = p_vals_BL_LP_NREM.(['Both' Sexes{1} 'and' Sexes{2}]);
	p_vals.(['Both' Sexes{1} 'and' Sexes{2}]).BL.LP.REM  = p_vals_BL_LP_REM.(['Both' Sexes{1}  'and' Sexes{2}]);
	p_vals.(['Both' Sexes{1} 'and' Sexes{2}]).BL.DP.Wake = p_vals_BL_DP_Wake.(['Both' Sexes{1} 'and' Sexes{2}]);
	p_vals.(['Both' Sexes{1} 'and' Sexes{2}]).BL.DP.NREM = p_vals_BL_DP_NREM.(['Both' Sexes{1} 'and' Sexes{2}]);
	p_vals.(['Both' Sexes{1} 'and' Sexes{2}]).BL.DP.REM  = p_vals_BL_DP_REM.(['Both' Sexes{1}  'and' Sexes{2}]);


	% Excluded because we never want this.  Removing to avoid confusion.  
	% p_vals.SD.LP.Wake = p_vals_SD_LP_Wake;
	% p_vals.SD.LP.NREM = p_vals_SD_LP_NREM;
	% p_vals.SD.LP.REM  = p_vals_SD_LP_REM;

	% -- SD Male -- 
	if Female_Data_Present p_vals.(Sexes{1}).SDexcSD.LP.Wake = p_vals_SDexcSD_LP_Wake.Posthoc.(Sexes{1}).WTvsMut; else p_vals.(Sexes{1}).SDexcSD.LP.Wake = p_vals_SDexcSD_LP_Wake.(Sexes{1}); end 
	if Female_Data_Present p_vals.(Sexes{1}).SDexcSD.LP.NREM = p_vals_SDexcSD_LP_NREM.Posthoc.(Sexes{1}).WTvsMut; else p_vals.(Sexes{1}).SDexcSD.LP.NREM = p_vals_SDexcSD_LP_NREM.(Sexes{1}); end
	if Female_Data_Present p_vals.(Sexes{1}).SDexcSD.LP.REM  = p_vals_SDexcSD_LP_REM.Posthoc.(Sexes{1}).WTvsMut;  else p_vals.(Sexes{1}).SDexcSD.LP.REM  = p_vals_SDexcSD_LP_REM.(Sexes{1});  end
	if Female_Data_Present p_vals.(Sexes{1}).SD.DP.Wake      = p_vals_SD_DP_Wake.Posthoc.(Sexes{1}).WTvsMut;      else p_vals.(Sexes{1}).SD.DP.Wake 	 = p_vals_SD_DP_Wake.(Sexes{1});      end
	if Female_Data_Present p_vals.(Sexes{1}).SD.DP.NREM      = p_vals_SD_DP_NREM.Posthoc.(Sexes{1}).WTvsMut;      else p_vals.(Sexes{1}).SD.DP.NREM 	 = p_vals_SD_DP_NREM.(Sexes{1});      end
	if Female_Data_Present p_vals.(Sexes{1}).SD.DP.REM       = p_vals_SD_DP_REM.Posthoc.(Sexes{1}).WTvsMut;       else p_vals.(Sexes{1}).SD.DP.REM  	 = p_vals_SD_DP_REM.(Sexes{1});       end

	% -- SD Female -- 
	p_vals.(Sexes{2}).SDexcSD.LP.Wake = p_vals_SDexcSD_LP_Wake.Posthoc.(Sexes{2});
	p_vals.(Sexes{2}).SDexcSD.LP.NREM = p_vals_SDexcSD_LP_NREM.Posthoc.(Sexes{2});
	p_vals.(Sexes{2}).SDexcSD.LP.REM  = p_vals_SDexcSD_LP_REM.Posthoc.(Sexes{2});
	p_vals.(Sexes{2}).SD.DP.Wake 	  = p_vals_SD_DP_Wake.Posthoc.(Sexes{2});
	p_vals.(Sexes{2}).SD.DP.NREM 	  = p_vals_SD_DP_NREM.Posthoc.(Sexes{2});
	p_vals.(Sexes{2}).SD.DP.REM  	  = p_vals_SD_DP_REM.Posthoc.(Sexes{2});

	% -- SD Both -- 
	p_vals.(['Both' Sexes{1} 'and' Sexes{2}]).SDexcSD.LP.Wake = p_vals_SDexcSD_LP_Wake.(['Both' Sexes{1} 'and' Sexes{2}]);
	p_vals.(['Both' Sexes{1} 'and' Sexes{2}]).SDexcSD.LP.NREM = p_vals_SDexcSD_LP_NREM.(['Both' Sexes{1} 'and' Sexes{2}]);
	p_vals.(['Both' Sexes{1} 'and' Sexes{2}]).SDexcSD.LP.REM  = p_vals_SDexcSD_LP_REM.(['Both' Sexes{1} 'and' Sexes{2}]);
	p_vals.(['Both' Sexes{1} 'and' Sexes{2}]).SD.DP.Wake 	 = p_vals_SD_DP_Wake.(['Both' Sexes{1} 'and' Sexes{2}]);
	p_vals.(['Both' Sexes{1} 'and' Sexes{2}]).SD.DP.NREM 	 = p_vals_SD_DP_NREM.(['Both' Sexes{1} 'and' Sexes{2}]);
	p_vals.(['Both' Sexes{1} 'and' Sexes{2}]).SD.DP.REM  	 = p_vals_SD_DP_REM.(['Both' Sexes{1} 'and' Sexes{2}]);

	% -- post-hocs for M vs F WT --
	% - BL LP -
	p_vals.(['Both' Sexes{1} 'and' Sexes{2}]).BL.LP.Wake.Posthocs.(['WT_' Sexes{1} 'Vs' Sexes{2}]) = p_vals_BL_LP_Wake.Posthoc.WT.([Sexes{1} 'vs' Sexes{2}]);
	p_vals.(['Both' Sexes{1} 'and' Sexes{2}]).BL.LP.NREM.Posthocs.(['WT_' Sexes{1} 'Vs' Sexes{2}]) = p_vals_BL_LP_NREM.Posthoc.WT.([Sexes{1} 'vs' Sexes{2}]);
	p_vals.(['Both' Sexes{1} 'and' Sexes{2}]).BL.LP.REM.Posthocs.(['WT_' Sexes{1} 'Vs' Sexes{2}])  = p_vals_BL_LP_REM.Posthoc.WT.([Sexes{1} 'vs' Sexes{2}]);

	% - BL DP -
	p_vals.(['Both' Sexes{1} 'and' Sexes{2}]).BL.DP.Wake.Posthocs.(['WT_' Sexes{1} 'Vs' Sexes{2}]) = p_vals_BL_DP_Wake.Posthoc.WT.([Sexes{1} 'vs' Sexes{2}]);
	p_vals.(['Both' Sexes{1} 'and' Sexes{2}]).BL.DP.NREM.Posthocs.(['WT_' Sexes{1} 'Vs' Sexes{2}]) = p_vals_BL_DP_NREM.Posthoc.WT.([Sexes{1} 'vs' Sexes{2}]);
	p_vals.(['Both' Sexes{1} 'and' Sexes{2}]).BL.DP.REM.Posthocs.(['WT_' Sexes{1} 'Vs' Sexes{2}])  = p_vals_BL_DP_REM.Posthoc.WT.([Sexes{1} 'vs' Sexes{2}]);

	% - SD LP -
	p_vals.(['Both' Sexes{1} 'and' Sexes{2}]).SDexcSD.LP.Wake.Posthocs.(['WT_' Sexes{1} 'Vs' Sexes{2}]) = p_vals_SDexcSD_LP_Wake.Posthoc.WT.([Sexes{1} 'vs' Sexes{2}]);
	p_vals.(['Both' Sexes{1} 'and' Sexes{2}]).SDexcSD.LP.NREM.Posthocs.(['WT_' Sexes{1} 'Vs' Sexes{2}]) = p_vals_SDexcSD_LP_NREM.Posthoc.WT.([Sexes{1} 'vs' Sexes{2}]);
	p_vals.(['Both' Sexes{1} 'and' Sexes{2}]).SDexcSD.LP.REM.Posthocs.(['WT_' Sexes{1} 'Vs' Sexes{2}])  = p_vals_SDexcSD_LP_REM.Posthoc.WT.([Sexes{1} 'vs' Sexes{2}]);

	% - SD DP -
	p_vals.(['Both' Sexes{1} 'and' Sexes{2}]).SD.DP.Wake.Posthocs.(['WT_' Sexes{1} 'Vs' Sexes{2}]) = p_vals_SD_DP_Wake.Posthoc.WT.([Sexes{1} 'vs' Sexes{2}]);
	p_vals.(['Both' Sexes{1} 'and' Sexes{2}]).SD.DP.NREM.Posthocs.(['WT_' Sexes{1} 'Vs' Sexes{2}]) = p_vals_SD_DP_NREM.Posthoc.WT.([Sexes{1} 'vs' Sexes{2}]);
	p_vals.(['Both' Sexes{1} 'and' Sexes{2}]).SD.DP.REM.Posthocs.(['WT_' Sexes{1} 'Vs' Sexes{2}])  = p_vals_SD_DP_REM.Posthoc.WT.([Sexes{1} 'vs' Sexes{2}]);

	% -- post-hocs for M vs F Mut --
	% - BL LP -
	p_vals.(['Both' Sexes{1} 'and' Sexes{2}]).BL.LP.Wake.Posthocs.(['Mut_' Sexes{1} 'Vs' Sexes{2}]) = p_vals_BL_LP_Wake.Posthoc.Mut.([Sexes{1} 'vs' Sexes{2}]);
	p_vals.(['Both' Sexes{1} 'and' Sexes{2}]).BL.LP.NREM.Posthocs.(['Mut_' Sexes{1} 'Vs' Sexes{2}]) = p_vals_BL_LP_NREM.Posthoc.Mut.([Sexes{1} 'vs' Sexes{2}]);
	p_vals.(['Both' Sexes{1} 'and' Sexes{2}]).BL.LP.REM.Posthocs.(['Mut_' Sexes{1} 'Vs' Sexes{2}])  = p_vals_BL_LP_REM.Posthoc.Mut.([Sexes{1} 'vs' Sexes{2}]);

	% - BL DP -
	p_vals.(['Both' Sexes{1} 'and' Sexes{2}]).BL.DP.Wake.Posthocs.(['Mut_' Sexes{1} 'Vs' Sexes{2}]) = p_vals_BL_DP_Wake.Posthoc.Mut.([Sexes{1} 'vs' Sexes{2}]);
	p_vals.(['Both' Sexes{1} 'and' Sexes{2}]).BL.DP.NREM.Posthocs.(['Mut_' Sexes{1} 'Vs' Sexes{2}]) = p_vals_BL_DP_NREM.Posthoc.Mut.([Sexes{1} 'vs' Sexes{2}]);
	p_vals.(['Both' Sexes{1} 'and' Sexes{2}]).BL.DP.REM.Posthocs.(['Mut_' Sexes{1} 'Vs' Sexes{2}])  = p_vals_BL_DP_REM.Posthoc.Mut.([Sexes{1} 'vs' Sexes{2}]);

	% - SD LP -
	p_vals.(['Both' Sexes{1} 'and' Sexes{2}]).SDexcSD.LP.Wake.Posthocs.(['Mut_' Sexes{1} 'Vs' Sexes{2}]) = p_vals_SDexcSD_LP_Wake.Posthoc.Mut.([Sexes{1} 'vs' Sexes{2}]);
	p_vals.(['Both' Sexes{1} 'and' Sexes{2}]).SDexcSD.LP.NREM.Posthocs.(['Mut_' Sexes{1} 'Vs' Sexes{2}]) = p_vals_SDexcSD_LP_NREM.Posthoc.Mut.([Sexes{1} 'vs' Sexes{2}]);
	p_vals.(['Both' Sexes{1} 'and' Sexes{2}]).SDexcSD.LP.REM.Posthocs.(['Mut_' Sexes{1} 'Vs' Sexes{2}])  = p_vals_SDexcSD_LP_REM.Posthoc.Mut.([Sexes{1} 'vs' Sexes{2}]);

	% - SD DP -
	p_vals.(['Both' Sexes{1} 'and' Sexes{2}]).SD.DP.Wake.Posthocs.(['Mut_' Sexes{1} 'Vs' Sexes{2}]) = p_vals_SD_DP_Wake.Posthoc.Mut.([Sexes{1} 'vs' Sexes{2}]);
	p_vals.(['Both' Sexes{1} 'and' Sexes{2}]).SD.DP.NREM.Posthocs.(['Mut_' Sexes{1} 'Vs' Sexes{2}]) = p_vals_SD_DP_NREM.Posthoc.Mut.([Sexes{1} 'vs' Sexes{2}]);
	p_vals.(['Both' Sexes{1} 'and' Sexes{2}]).SD.DP.REM.Posthocs.(['Mut_' Sexes{1} 'Vs' Sexes{2}])  = p_vals_SD_DP_REM.Posthoc.Mut.([Sexes{1} 'vs' Sexes{2}]);


	% -- BL Male --
	sig_difference_bool.(Sexes{1}).BL.LP.Wake = sig_difference_BL_LP_Wake.(Sexes{1});
	sig_difference_bool.(Sexes{1}).BL.LP.NREM = sig_difference_BL_LP_NREM.(Sexes{1});
	sig_difference_bool.(Sexes{1}).BL.LP.REM  = sig_difference_BL_LP_REM.(Sexes{1});
	sig_difference_bool.(Sexes{1}).BL.DP.Wake = sig_difference_BL_DP_Wake.(Sexes{1});
	sig_difference_bool.(Sexes{1}).BL.DP.NREM = sig_difference_BL_DP_NREM.(Sexes{1});
	sig_difference_bool.(Sexes{1}).BL.DP.REM  = sig_difference_BL_DP_REM.(Sexes{1});

	% -- BL Female --
	sig_difference_bool.(Sexes{2}).BL.LP.Wake = sig_difference_BL_LP_Wake.(Sexes{2});
	sig_difference_bool.(Sexes{2}).BL.LP.NREM = sig_difference_BL_LP_NREM.(Sexes{2});
	sig_difference_bool.(Sexes{2}).BL.LP.REM  = sig_difference_BL_LP_REM.(Sexes{2});
	sig_difference_bool.(Sexes{2}).BL.DP.Wake = sig_difference_BL_DP_Wake.(Sexes{2});
	sig_difference_bool.(Sexes{2}).BL.DP.NREM = sig_difference_BL_DP_NREM.(Sexes{2});
	sig_difference_bool.(Sexes{2}).BL.DP.REM  = sig_difference_BL_DP_REM.(Sexes{2});

	% -- SD Male -- 
	sig_difference_bool.(Sexes{1}).SDexcSD.LP.Wake = sig_difference_SDexcSD_LP_Wake.(Sexes{1});
	sig_difference_bool.(Sexes{1}).SDexcSD.LP.NREM = sig_difference_SDexcSD_LP_NREM.(Sexes{1});
	sig_difference_bool.(Sexes{1}).SDexcSD.LP.REM  = sig_difference_SDexcSD_LP_REM.(Sexes{1});
	sig_difference_bool.(Sexes{1}).SD.DP.Wake = sig_difference_SD_DP_Wake.(Sexes{1});
	sig_difference_bool.(Sexes{1}).SD.DP.NREM = sig_difference_SD_DP_NREM.(Sexes{1});
	sig_difference_bool.(Sexes{1}).SD.DP.REM  = sig_difference_SD_DP_REM.(Sexes{1});

	% -- SD Female -- 
	sig_difference_bool.(Sexes{2}).SDexcSD.LP.Wake = sig_difference_SDexcSD_LP_Wake.(Sexes{2});
	sig_difference_bool.(Sexes{2}).SDexcSD.LP.NREM = sig_difference_SDexcSD_LP_NREM.(Sexes{2});
	sig_difference_bool.(Sexes{2}).SDexcSD.LP.REM  = sig_difference_SDexcSD_LP_REM.(Sexes{2});
	sig_difference_bool.(Sexes{2}).SD.DP.Wake = sig_difference_SD_DP_Wake.(Sexes{2});
	sig_difference_bool.(Sexes{2}).SD.DP.NREM = sig_difference_SD_DP_NREM.(Sexes{2});
	sig_difference_bool.(Sexes{2}).SD.DP.REM  = sig_difference_SD_DP_REM.(Sexes{2});

	% -- BL Male --
	if ~Female_Data_Present ANOVA_tbls.(Sexes{1}).BL.LP.Wake = ANOVA_tbl_BL_LP_Wake.(Sexes{1}); else ANOVA_tbls.(Sexes{1}).BL.LP.Wake = []; end
	if ~Female_Data_Present ANOVA_tbls.(Sexes{1}).BL.LP.NREM = ANOVA_tbl_BL_LP_NREM.(Sexes{1}); else ANOVA_tbls.(Sexes{1}).BL.LP.NREM = []; end
	if ~Female_Data_Present ANOVA_tbls.(Sexes{1}).BL.LP.REM  = ANOVA_tbl_BL_LP_REM.(Sexes{1});  else ANOVA_tbls.(Sexes{1}).BL.LP.REM  = []; end
	if ~Female_Data_Present ANOVA_tbls.(Sexes{1}).BL.DP.Wake = ANOVA_tbl_BL_DP_Wake.(Sexes{1}); else ANOVA_tbls.(Sexes{1}).BL.DP.Wake = []; end
	if ~Female_Data_Present ANOVA_tbls.(Sexes{1}).BL.DP.NREM = ANOVA_tbl_BL_DP_NREM.(Sexes{1}); else ANOVA_tbls.(Sexes{1}).BL.DP.NREM = []; end
	if ~Female_Data_Present ANOVA_tbls.(Sexes{1}).BL.DP.REM  = ANOVA_tbl_BL_DP_REM.(Sexes{1});  else ANOVA_tbls.(Sexes{1}).BL.DP.REM  = []; end

	% -- SD Male -- 
	if ~Female_Data_Present ANOVA_tbls.(Sexes{1}).SDexcSD.LP.Wake = ANOVA_tbl_SDexcSD_LP_Wake.(Sexes{1}); else ANOVA_tbls.(Sexes{1}).SDexcSD.LP.Wake = []; end
	if ~Female_Data_Present ANOVA_tbls.(Sexes{1}).SDexcSD.LP.NREM = ANOVA_tbl_SDexcSD_LP_NREM.(Sexes{1}); else ANOVA_tbls.(Sexes{1}).SDexcSD.LP.NREM = []; end
	if ~Female_Data_Present ANOVA_tbls.(Sexes{1}).SDexcSD.LP.REM  = ANOVA_tbl_SDexcSD_LP_REM.(Sexes{1});  else ANOVA_tbls.(Sexes{1}).SDexcSD.LP.REM  = []; end
	if ~Female_Data_Present ANOVA_tbls.(Sexes{1}).SD.DP.Wake      = ANOVA_tbl_SD_DP_Wake.(Sexes{1});      else ANOVA_tbls.(Sexes{1}).SD.DP.Wake      = []; end
	if ~Female_Data_Present ANOVA_tbls.(Sexes{1}).SD.DP.NREM      = ANOVA_tbl_SD_DP_NREM.(Sexes{1});      else ANOVA_tbls.(Sexes{1}).SD.DP.NREM      = []; end 
	if ~Female_Data_Present ANOVA_tbls.(Sexes{1}).SD.DP.REM       = ANOVA_tbl_SD_DP_REM.(Sexes{1});       else ANOVA_tbls.(Sexes{1}).SD.DP.REM       = []; end

	% -- BL Male and Female --
	if Female_Data_Present
		ANOVA_tbls.(['Both' Sexes{1} 'and' Sexes{2}]).BL.LP.Wake = ANOVA_tbl_BL_LP_Wake.([SexVarName 'Genotype']);  
		ANOVA_tbls.(['Both' Sexes{1} 'and' Sexes{2}]).BL.LP.NREM = ANOVA_tbl_BL_LP_NREM.([SexVarName 'Genotype']);
		ANOVA_tbls.(['Both' Sexes{1} 'and' Sexes{2}]).BL.LP.REM  = ANOVA_tbl_BL_LP_REM.([SexVarName 'Genotype']); 
		ANOVA_tbls.(['Both' Sexes{1} 'and' Sexes{2}]).BL.DP.Wake = ANOVA_tbl_BL_DP_Wake.([SexVarName 'Genotype']);  
		ANOVA_tbls.(['Both' Sexes{1} 'and' Sexes{2}]).BL.DP.NREM = ANOVA_tbl_BL_DP_NREM.([SexVarName 'Genotype']);
		ANOVA_tbls.(['Both' Sexes{1} 'and' Sexes{2}]).BL.DP.REM  = ANOVA_tbl_BL_DP_REM.([SexVarName 'Genotype']); 
	end 	

	% -- SD Male and Female --
	if Female_Data_Present
		ANOVA_tbls.(['Both' Sexes{1} 'and' Sexes{2}]).SDexcSD.LP.Wake = ANOVA_tbl_SDexcSD_LP_Wake.([SexVarName 'Genotype']);  
		ANOVA_tbls.(['Both' Sexes{1} 'and' Sexes{2}]).SDexcSD.LP.NREM = ANOVA_tbl_SDexcSD_LP_NREM.([SexVarName 'Genotype']);
		ANOVA_tbls.(['Both' Sexes{1} 'and' Sexes{2}]).SDexcSD.LP.REM  = ANOVA_tbl_SDexcSD_LP_REM.([SexVarName 'Genotype']); 
		ANOVA_tbls.(['Both' Sexes{1} 'and' Sexes{2}]).SD.DP.Wake      = ANOVA_tbl_SD_DP_Wake.([SexVarName 'Genotype']);  
		ANOVA_tbls.(['Both' Sexes{1} 'and' Sexes{2}]).SD.DP.NREM      = ANOVA_tbl_SD_DP_NREM.([SexVarName 'Genotype']);
		ANOVA_tbls.(['Both' Sexes{1} 'and' Sexes{2}]).SD.DP.REM       = ANOVA_tbl_SD_DP_REM.([SexVarName 'Genotype']); 
	end 

	% -- BL Female --
	% ANOVA_tbls.(Sexes{2}).BL.LP.Wake = ANOVA_tbl_BL_LP_Wake.(Sexes{2});
	% ANOVA_tbls.(Sexes{2}).BL.LP.NREM = ANOVA_tbl_BL_LP_NREM.(Sexes{2});
	% ANOVA_tbls.(Sexes{2}).BL.LP.REM  = ANOVA_tbl_BL_LP_REM.(Sexes{2});
	% ANOVA_tbls.(Sexes{2}).BL.DP.Wake = ANOVA_tbl_BL_DP_Wake.(Sexes{2});
	% ANOVA_tbls.(Sexes{2}).BL.DP.NREM = ANOVA_tbl_BL_DP_NREM.(Sexes{2});
	% ANOVA_tbls.(Sexes{2}).BL.DP.REM  = ANOVA_tbl_BL_DP_REM.(Sexes{2});

	% % -- SD Female -- 
	% ANOVA_tbls.(Sexes{2}).SDexcSD.LP.Wake = ANOVA_tbl_SDexcSD_LP_Wake.(Sexes{2});
	% ANOVA_tbls.(Sexes{2}).SDexcSD.LP.NREM = ANOVA_tbl_SDexcSD_LP_NREM.(Sexes{2});
	% ANOVA_tbls.(Sexes{2}).SDexcSD.LP.REM  = ANOVA_tbl_SDexcSD_LP_REM.(Sexes{2});
	% ANOVA_tbls.(Sexes{2}).SD.DP.Wake      = ANOVA_tbl_SD_DP_Wake.(Sexes{2});
	% ANOVA_tbls.(Sexes{2}).SD.DP.NREM      = ANOVA_tbl_SD_DP_NREM.(Sexes{2});
	% ANOVA_tbls.(Sexes{2}).SD.DP.REM       = ANOVA_tbl_SD_DP_REM.(Sexes{2});

% end of case where WindowLength=12
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------


elseif WindowLength==6
	
	Six_hour_avg_percentages_Male   = Twelve_hour_avg_percentages_Male;
	Six_hour_avg_percentages_Female = Twelve_hour_avg_percentages_Female;



	% If you don't have female data, set up those structs as empty
	if isempty(Twelve_hour_avg_percentages_Female) 	% I know we are doing 6-hr intervals, not 12. this is just the input. 
		Female_Data_Present = false; 
		
		% BL
		Six_hour_avg_percentages_Female.WT.BL.DPfirst6.Wake  = [];
		Six_hour_avg_percentages_Female.Mut.BL.DPfirst6.Wake = [];
		Six_hour_avg_percentages_Female.WT.BL.DPfirst6.NREM  = [];
		Six_hour_avg_percentages_Female.Mut.BL.DPfirst6.NREM = [];
		Six_hour_avg_percentages_Female.WT.BL.DPfirst6.REM   = [];
		Six_hour_avg_percentages_Female.Mut.BL.DPfirst6.REM  = [];
		Six_hour_avg_percentages_Female.WT.BL.DPlast6.Wake  = [];
		Six_hour_avg_percentages_Female.Mut.BL.DPlast6.Wake = [];
		Six_hour_avg_percentages_Female.WT.BL.DPlast6.NREM  = [];
		Six_hour_avg_percentages_Female.Mut.BL.DPlast6.NREM = [];
		Six_hour_avg_percentages_Female.WT.BL.DPlast6.REM   = [];
		Six_hour_avg_percentages_Female.Mut.BL.DPlast6.REM  = [];		

		% SD
		Six_hour_avg_percentages_Female.WT.SD.DPfirst6.Wake  = [];
		Six_hour_avg_percentages_Female.Mut.SD.DPfirst6.Wake = [];
		Six_hour_avg_percentages_Female.WT.SD.DPfirst6.NREM  = [];
		Six_hour_avg_percentages_Female.Mut.SD.DPfirst6.NREM = [];
		Six_hour_avg_percentages_Female.WT.SD.DPfirst6.REM   = [];
		Six_hour_avg_percentages_Female.Mut.SD.DPfirst6.REM  = [];
		Six_hour_avg_percentages_Female.WT.SD.DPlast6.Wake  = [];
		Six_hour_avg_percentages_Female.Mut.SD.DPlast6.Wake = [];
		Six_hour_avg_percentages_Female.WT.SD.DPlast6.NREM  = [];
		Six_hour_avg_percentages_Female.Mut.SD.DPlast6.NREM = [];
		Six_hour_avg_percentages_Female.WT.SD.DPlast6.REM   = [];
		Six_hour_avg_percentages_Female.Mut.SD.DPlast6.REM  = [];
		

		

		
	else 
		Female_Data_Present = true;
	end 

	% BL DP (first 6 hrs) Wake
	Female_WT_BL_DPfirst6_Wake  = Six_hour_avg_percentages_Female.WT.BL.DPfirst6.Wake;
	Female_Mut_BL_DPfirst6_Wake = Six_hour_avg_percentages_Female.Mut.BL.DPfirst6.Wake;
	Male_WT_BL_DPfirst6_Wake    = Six_hour_avg_percentages_Male.WT.BL.DPfirst6.Wake;
	Male_Mut_BL_DPfirst6_Wake   = Six_hour_avg_percentages_Male.Mut.BL.DPfirst6.Wake;

	% BL DP (last 6 hrs) Wake
	Female_WT_BL_DPlast6_Wake  = Six_hour_avg_percentages_Female.WT.BL.DPlast6.Wake;
	Female_Mut_BL_DPlast6_Wake = Six_hour_avg_percentages_Female.Mut.BL.DPlast6.Wake;
	Male_WT_BL_DPlast6_Wake    = Six_hour_avg_percentages_Male.WT.BL.DPlast6.Wake;
	Male_Mut_BL_DPlast6_Wake   = Six_hour_avg_percentages_Male.Mut.BL.DPlast6.Wake;

	% BL DP (first 6 hrs) NREM
	Female_WT_BL_DPfirst6_NREM  = Six_hour_avg_percentages_Female.WT.BL.DPfirst6.NREM;
	Female_Mut_BL_DPfirst6_NREM = Six_hour_avg_percentages_Female.Mut.BL.DPfirst6.NREM;
	Male_WT_BL_DPfirst6_NREM    = Six_hour_avg_percentages_Male.WT.BL.DPfirst6.NREM;
	Male_Mut_BL_DPfirst6_NREM   = Six_hour_avg_percentages_Male.Mut.BL.DPfirst6.NREM;

	% BL DP (last 6 hrs) NREM
	Female_WT_BL_DPlast6_NREM  = Six_hour_avg_percentages_Female.WT.BL.DPlast6.NREM;
	Female_Mut_BL_DPlast6_NREM = Six_hour_avg_percentages_Female.Mut.BL.DPlast6.NREM;
	Male_WT_BL_DPlast6_NREM    = Six_hour_avg_percentages_Male.WT.BL.DPlast6.NREM;
	Male_Mut_BL_DPlast6_NREM   = Six_hour_avg_percentages_Male.Mut.BL.DPlast6.NREM;

	% BL DP (first 6 hrs) REM
	Female_WT_BL_DPfirst6_REM  = Six_hour_avg_percentages_Female.WT.BL.DPfirst6.REM;
	Female_Mut_BL_DPfirst6_REM = Six_hour_avg_percentages_Female.Mut.BL.DPfirst6.REM;
	Male_WT_BL_DPfirst6_REM    = Six_hour_avg_percentages_Male.WT.BL.DPfirst6.REM;
	Male_Mut_BL_DPfirst6_REM   = Six_hour_avg_percentages_Male.Mut.BL.DPfirst6.REM;

	% BL DP (last 6 hrs) REM
	Female_WT_BL_DPlast6_REM  = Six_hour_avg_percentages_Female.WT.BL.DPlast6.REM;
	Female_Mut_BL_DPlast6_REM = Six_hour_avg_percentages_Female.Mut.BL.DPlast6.REM;
	Male_WT_BL_DPlast6_REM    = Six_hour_avg_percentages_Male.WT.BL.DPlast6.REM;
	Male_Mut_BL_DPlast6_REM   = Six_hour_avg_percentages_Male.Mut.BL.DPlast6.REM;

	% SD DP (first 6 hrs) Wake
	Female_WT_SD_DPfirst6_Wake  = Six_hour_avg_percentages_Female.WT.SD.DPfirst6.Wake;
	Female_Mut_SD_DPfirst6_Wake = Six_hour_avg_percentages_Female.Mut.SD.DPfirst6.Wake;
	Male_WT_SD_DPfirst6_Wake    = Six_hour_avg_percentages_Male.WT.SD.DPfirst6.Wake;
	Male_Mut_SD_DPfirst6_Wake   = Six_hour_avg_percentages_Male.Mut.SD.DPfirst6.Wake;

	% SD DP (last 6 hrs) Wake
	Female_WT_SD_DPlast6_Wake  = Six_hour_avg_percentages_Female.WT.SD.DPlast6.Wake;
	Female_Mut_SD_DPlast6_Wake = Six_hour_avg_percentages_Female.Mut.SD.DPlast6.Wake;
	Male_WT_SD_DPlast6_Wake    = Six_hour_avg_percentages_Male.WT.SD.DPlast6.Wake;
	Male_Mut_SD_DPlast6_Wake   = Six_hour_avg_percentages_Male.Mut.SD.DPlast6.Wake;

	% SD DP (first 6 hrs) NREM
	Female_WT_SD_DPfirst6_NREM  = Six_hour_avg_percentages_Female.WT.SD.DPfirst6.NREM;
	Female_Mut_SD_DPfirst6_NREM = Six_hour_avg_percentages_Female.Mut.SD.DPfirst6.NREM;
	Male_WT_SD_DPfirst6_NREM    = Six_hour_avg_percentages_Male.WT.SD.DPfirst6.NREM;
	Male_Mut_SD_DPfirst6_NREM   = Six_hour_avg_percentages_Male.Mut.SD.DPfirst6.NREM;

	% SD DP (last 6 hrs) NREM
	Female_WT_SD_DPlast6_NREM  = Six_hour_avg_percentages_Female.WT.SD.DPlast6.NREM;
	Female_Mut_SD_DPlast6_NREM = Six_hour_avg_percentages_Female.Mut.SD.DPlast6.NREM;
	Male_WT_SD_DPlast6_NREM    = Six_hour_avg_percentages_Male.WT.SD.DPlast6.NREM;
	Male_Mut_SD_DPlast6_NREM   = Six_hour_avg_percentages_Male.Mut.SD.DPlast6.NREM;

	% SD DP (first 6 hrs) REM
	Female_WT_SD_DPfirst6_REM  = Six_hour_avg_percentages_Female.WT.SD.DPfirst6.REM;
	Female_Mut_SD_DPfirst6_REM = Six_hour_avg_percentages_Female.Mut.SD.DPfirst6.REM;
	Male_WT_SD_DPfirst6_REM    = Six_hour_avg_percentages_Male.WT.SD.DPfirst6.REM;
	Male_Mut_SD_DPfirst6_REM   = Six_hour_avg_percentages_Male.Mut.SD.DPfirst6.REM;

	% SD DP (last 6 hrs) REM
	Female_WT_SD_DPlast6_REM  = Six_hour_avg_percentages_Female.WT.SD.DPlast6.REM;
	Female_Mut_SD_DPlast6_REM = Six_hour_avg_percentages_Female.Mut.SD.DPlast6.REM;
	Male_WT_SD_DPlast6_REM    = Six_hour_avg_percentages_Male.WT.SD.DPlast6.REM;
	Male_Mut_SD_DPlast6_REM   = Six_hour_avg_percentages_Male.Mut.SD.DPlast6.REM;
	

	% ------------------------ Baseline ---------------------------------------
	
	% ------- Dark Phase First 6 hrs --------------
	% Wake 
	[p_vals_BL_DPfirst6_Wake,sig_difference_BL_DPfirst6_Wake,ANOVA_tbl_BL_DPfirst6_Wake] = perform_2way_Anova_SexGenotype_posthoc(Female_WT_BL_DPfirst6_Wake,Female_Mut_BL_DPfirst6_Wake,Male_WT_BL_DPfirst6_Wake,Male_Mut_BL_DPfirst6_Wake,Sexes,SexVarName);
	% NREM
	[p_vals_BL_DPfirst6_NREM,sig_difference_BL_DPfirst6_NREM,ANOVA_tbl_BL_DPfirst6_NREM] = perform_2way_Anova_SexGenotype_posthoc(Female_WT_BL_DPfirst6_NREM,Female_Mut_BL_DPfirst6_NREM,Male_WT_BL_DPfirst6_NREM,Male_Mut_BL_DPfirst6_NREM,Sexes,SexVarName);
	% REM 
	[p_vals_BL_DPfirst6_REM,sig_difference_BL_DPfirst6_REM,ANOVA_tbl_BL_DPfirst6_REM]    = perform_2way_Anova_SexGenotype_posthoc(Female_WT_BL_DPfirst6_REM,Female_Mut_BL_DPfirst6_REM,Male_WT_BL_DPfirst6_REM,Male_Mut_BL_DPfirst6_REM,Sexes,SexVarName);

	% ------- Dark Phase Last 6 hrs --------------
	% Wake 
	[p_vals_BL_DPlast6_Wake,sig_difference_BL_DPlast6_Wake,ANOVA_tbl_BL_DPlast6_Wake] = perform_2way_Anova_SexGenotype_posthoc(Female_WT_BL_DPlast6_Wake,Female_Mut_BL_DPlast6_Wake,Male_WT_BL_DPlast6_Wake,Male_Mut_BL_DPlast6_Wake,Sexes,SexVarName);
	% NREM
	[p_vals_BL_DPlast6_NREM,sig_difference_BL_DPlast6_NREM,ANOVA_tbl_BL_DPlast6_NREM] = perform_2way_Anova_SexGenotype_posthoc(Female_WT_BL_DPlast6_NREM,Female_Mut_BL_DPlast6_NREM,Male_WT_BL_DPlast6_NREM,Male_Mut_BL_DPlast6_NREM,Sexes,SexVarName);
	% REM 
	[p_vals_BL_DPlast6_REM,sig_difference_BL_DPlast6_REM,ANOVA_tbl_BL_DPlast6_REM]    = perform_2way_Anova_SexGenotype_posthoc(Female_WT_BL_DPlast6_REM,Female_Mut_BL_DPlast6_REM,Male_WT_BL_DPlast6_REM,Male_Mut_BL_DPlast6_REM,Sexes,SexVarName);


	% ------------------------ Sleep Dep ---------------------------------------
	
	% ------- Dark Phase First 6 hrs --------------
	% Wake 
	[p_vals_SD_DPfirst6_Wake,sig_difference_SD_DPfirst6_Wake,ANOVA_tbl_SD_DPfirst6_Wake] = perform_2way_Anova_SexGenotype_posthoc(Female_WT_SD_DPfirst6_Wake,Female_Mut_SD_DPfirst6_Wake,Male_WT_SD_DPfirst6_Wake,Male_Mut_SD_DPfirst6_Wake,Sexes,SexVarName);
	% NREM
	[p_vals_SD_DPfirst6_NREM,sig_difference_SD_DPfirst6_NREM,ANOVA_tbl_SD_DPfirst6_NREM] = perform_2way_Anova_SexGenotype_posthoc(Female_WT_SD_DPfirst6_NREM,Female_Mut_SD_DPfirst6_NREM,Male_WT_SD_DPfirst6_NREM,Male_Mut_SD_DPfirst6_NREM,Sexes,SexVarName);
	% REM 
	[p_vals_SD_DPfirst6_REM,sig_difference_SD_DPfirst6_REM,ANOVA_tbl_SD_DPfirst6_REM]    = perform_2way_Anova_SexGenotype_posthoc(Female_WT_SD_DPfirst6_REM,Female_Mut_SD_DPfirst6_REM,Male_WT_SD_DPfirst6_REM,Male_Mut_SD_DPfirst6_REM,Sexes,SexVarName);

	% ------- Dark Phase Last 6 hrs --------------
	% Wake 
	[p_vals_SD_DPlast6_Wake,sig_difference_SD_DPlast6_Wake,ANOVA_tbl_SD_DPlast6_Wake] = perform_2way_Anova_SexGenotype_posthoc(Female_WT_SD_DPlast6_Wake,Female_Mut_SD_DPlast6_Wake,Male_WT_SD_DPlast6_Wake,Male_Mut_SD_DPlast6_Wake,Sexes,SexVarName);
	% NREM
	[p_vals_SD_DPlast6_NREM,sig_difference_SD_DPlast6_NREM,ANOVA_tbl_SD_DPlast6_NREM] = perform_2way_Anova_SexGenotype_posthoc(Female_WT_SD_DPlast6_NREM,Female_Mut_SD_DPlast6_NREM,Male_WT_SD_DPlast6_NREM,Male_Mut_SD_DPlast6_NREM,Sexes,SexVarName);
	% REM 
	[p_vals_SD_DPlast6_REM,sig_difference_SD_DPlast6_REM,ANOVA_tbl_SD_DPlast6_REM]    = perform_2way_Anova_SexGenotype_posthoc(Female_WT_SD_DPlast6_REM,Female_Mut_SD_DPlast6_REM,Male_WT_SD_DPlast6_REM,Male_Mut_SD_DPlast6_REM,Sexes,SexVarName);

	% -- BL Male --
	if Female_Data_Present p_vals.(Sexes{1}).BL.DPfirst6.Wake = p_vals_BL_DPfirst6_Wake.Posthoc.(Sexes{1}).WTvsMut; else p_vals.(Sexes{1}).BL.DPfirst6.Wake = p_vals_BL_DPfirst6_Wake.(Sexes{1}); end
	if Female_Data_Present p_vals.(Sexes{1}).BL.DPfirst6.NREM = p_vals_BL_DPfirst6_NREM.Posthoc.(Sexes{1}).WTvsMut; else p_vals.(Sexes{1}).BL.DPfirst6.NREM = p_vals_BL_DPfirst6_NREM.(Sexes{1}); end
	if Female_Data_Present p_vals.(Sexes{1}).BL.DPfirst6.REM  = p_vals_BL_DPfirst6_REM.Posthoc.(Sexes{1}).WTvsMut;  else p_vals.(Sexes{1}).BL.DPfirst6.REM  = p_vals_BL_DPfirst6_REM.(Sexes{1});  end
	if Female_Data_Present p_vals.(Sexes{1}).BL.DPlast6.Wake  = p_vals_BL_DPlast6_Wake.Posthoc.(Sexes{1}).WTvsMut;  else p_vals.(Sexes{1}).BL.DPlast6.Wake  = p_vals_BL_DPlast6_Wake.(Sexes{1});  end
	if Female_Data_Present p_vals.(Sexes{1}).BL.DPlast6.NREM  = p_vals_BL_DPlast6_NREM.Posthoc.(Sexes{1}).WTvsMut;  else p_vals.(Sexes{1}).BL.DPlast6.NREM  = p_vals_BL_DPlast6_NREM.(Sexes{1});  end
	if Female_Data_Present p_vals.(Sexes{1}).BL.DPlast6.REM   = p_vals_BL_DPlast6_REM.Posthoc.(Sexes{1}).WTvsMut;   else p_vals.(Sexes{1}).BL.DPlast6.REM   = p_vals_BL_DPlast6_REM.(Sexes{1});   end

	% -- BL Female --
	p_vals.(Sexes{2}).BL.DPfirst6.Wake = p_vals_BL_DPfirst6_Wake.Posthoc.(Sexes{2});
	p_vals.(Sexes{2}).BL.DPfirst6.NREM = p_vals_BL_DPfirst6_NREM.Posthoc.(Sexes{2});
	p_vals.(Sexes{2}).BL.DPfirst6.REM  = p_vals_BL_DPfirst6_REM.Posthoc.(Sexes{2});
	p_vals.(Sexes{2}).BL.DPlast6.Wake  = p_vals_BL_DPlast6_Wake.Posthoc.(Sexes{2});
	p_vals.(Sexes{2}).BL.DPlast6.NREM  = p_vals_BL_DPlast6_NREM.Posthoc.(Sexes{2});
	p_vals.(Sexes{2}).BL.DPlast6.REM   = p_vals_BL_DPlast6_REM.Posthoc.(Sexes{2});

	% -- BL Both --
	p_vals.(['Both' Sexes{1} 'and' Sexes{2}]).BL.DPfirst6.Wake = p_vals_BL_DPfirst6_Wake.(['Both' Sexes{1} 'and' Sexes{2}]);
	p_vals.(['Both' Sexes{1} 'and' Sexes{2}]).BL.DPfirst6.NREM = p_vals_BL_DPfirst6_NREM.(['Both' Sexes{1} 'and' Sexes{2}]);
	p_vals.(['Both' Sexes{1} 'and' Sexes{2}]).BL.DPfirst6.REM  = p_vals_BL_DPfirst6_REM.(['Both' Sexes{1} 'and' Sexes{2}]);
	p_vals.(['Both' Sexes{1} 'and' Sexes{2}]).BL.DPlast6.Wake  = p_vals_BL_DPlast6_Wake.(['Both' Sexes{1} 'and' Sexes{2}]);
	p_vals.(['Both' Sexes{1} 'and' Sexes{2}]).BL.DPlast6.NREM  = p_vals_BL_DPlast6_NREM.(['Both' Sexes{1} 'and' Sexes{2}]);
	p_vals.(['Both' Sexes{1} 'and' Sexes{2}]).BL.DPlast6.REM   = p_vals_BL_DPlast6_REM.(['Both' Sexes{1} 'and' Sexes{2}]);


	% -- SD Male -- 
	if Female_Data_Present p_vals.(Sexes{1}).SD.DPfirst6.Wake = p_vals_SD_DPfirst6_Wake.Posthoc.(Sexes{1}).WTvsMut; else p_vals.(Sexes{1}).SD.DPfirst6.Wake = p_vals_SD_DPfirst6_Wake.(Sexes{1}); end 
	if Female_Data_Present p_vals.(Sexes{1}).SD.DPfirst6.NREM = p_vals_SD_DPfirst6_NREM.Posthoc.(Sexes{1}).WTvsMut; else p_vals.(Sexes{1}).SD.DPfirst6.NREM = p_vals_SD_DPfirst6_NREM.(Sexes{1}); end
	if Female_Data_Present p_vals.(Sexes{1}).SD.DPfirst6.REM  = p_vals_SD_DPfirst6_REM.Posthoc.(Sexes{1}).WTvsMut;  else p_vals.(Sexes{1}).SD.DPfirst6.REM  = p_vals_SD_DPfirst6_REM.(Sexes{1});  end
	if Female_Data_Present p_vals.(Sexes{1}).SD.DPlast6.Wake  = p_vals_SD_DPlast6_Wake.Posthoc.(Sexes{1}).WTvsMut;  else p_vals.(Sexes{1}).SD.DPlast6.Wake  = p_vals_SD_DPlast6_Wake.(Sexes{1});  end
	if Female_Data_Present p_vals.(Sexes{1}).SD.DPlast6.NREM  = p_vals_SD_DPlast6_NREM.Posthoc.(Sexes{1}).WTvsMut;  else p_vals.(Sexes{1}).SD.DPlast6.NREM  = p_vals_SD_DPlast6_NREM.(Sexes{1});  end
	if Female_Data_Present p_vals.(Sexes{1}).SD.DPlast6.REM   = p_vals_SD_DPlast6_REM.Posthoc.(Sexes{1}).WTvsMut;   else p_vals.(Sexes{1}).SD.DPlast6.REM   = p_vals_SD_DPlast6_REM.(Sexes{1});   end

	% -- SD Female -- 
	p_vals.(Sexes{2}).SD.DPfirst6.Wake = p_vals_SD_DPfirst6_Wake.Posthoc.(Sexes{2});
	p_vals.(Sexes{2}).SD.DPfirst6.NREM = p_vals_SD_DPfirst6_NREM.Posthoc.(Sexes{2});
	p_vals.(Sexes{2}).SD.DPfirst6.REM  = p_vals_SD_DPfirst6_REM.Posthoc.(Sexes{2});
	p_vals.(Sexes{2}).SD.DPlast6.Wake  = p_vals_SD_DPlast6_Wake.Posthoc.(Sexes{2});
	p_vals.(Sexes{2}).SD.DPlast6.NREM  = p_vals_SD_DPlast6_NREM.Posthoc.(Sexes{2});
	p_vals.(Sexes{2}).SD.DPlast6.REM   = p_vals_SD_DPlast6_REM.Posthoc.(Sexes{2});

	% -- SD Both -- 
	p_vals.(['Both' Sexes{1} 'and' Sexes{2}]).SD.DPfirst6.Wake = p_vals_SD_DPfirst6_Wake.(['Both' Sexes{1} 'and' Sexes{2}]);
	p_vals.(['Both' Sexes{1} 'and' Sexes{2}]).SD.DPfirst6.NREM = p_vals_SD_DPfirst6_NREM.(['Both' Sexes{1} 'and' Sexes{2}]);
	p_vals.(['Both' Sexes{1} 'and' Sexes{2}]).SD.DPfirst6.REM  = p_vals_SD_DPfirst6_REM.(['Both' Sexes{1} 'and' Sexes{2}]);
	p_vals.(['Both' Sexes{1} 'and' Sexes{2}]).SD.DPlast6.Wake  = p_vals_SD_DPlast6_Wake.(['Both' Sexes{1} 'and' Sexes{2}]);
	p_vals.(['Both' Sexes{1} 'and' Sexes{2}]).SD.DPlast6.NREM  = p_vals_SD_DPlast6_NREM.(['Both' Sexes{1} 'and' Sexes{2}]);
	p_vals.(['Both' Sexes{1} 'and' Sexes{2}]).SD.DPlast6.REM   = p_vals_SD_DPlast6_REM.(['Both' Sexes{1} 'and' Sexes{2}]);

	% -- post-hocs for M vs F WT --
	% - BL DP first 6 hrs -
	p_vals.(['Both' Sexes{1} 'and' Sexes{2}]).BL.DPfirst6.Wake.Posthocs.(['WT_' Sexes{1} 'Vs' Sexes{2}]) = p_vals_BL_DPfirst6_Wake.Posthoc.WT.([Sexes{1} 'vs' Sexes{2}]);
	p_vals.(['Both' Sexes{1} 'and' Sexes{2}]).BL.DPfirst6.NREM.Posthocs.(['WT_' Sexes{1} 'Vs' Sexes{2}]) = p_vals_BL_DPfirst6_NREM.Posthoc.WT.([Sexes{1} 'vs' Sexes{2}]);
	p_vals.(['Both' Sexes{1} 'and' Sexes{2}]).BL.DPfirst6.REM.Posthocs.(['WT_' Sexes{1} 'Vs' Sexes{2}])  = p_vals_BL_DPfirst6_REM.Posthoc.WT.([Sexes{1} 'vs' Sexes{2}]);

	% - BL DP last 6 hrs -
	p_vals.(['Both' Sexes{1} 'and' Sexes{2}]).BL.DPlast6.Wake.Posthocs.(['WT_' Sexes{1} 'Vs' Sexes{2}]) = p_vals_BL_DPlast6_Wake.Posthoc.WT.([Sexes{1} 'vs' Sexes{2}]);
	p_vals.(['Both' Sexes{1} 'and' Sexes{2}]).BL.DPlast6.NREM.Posthocs.(['WT_' Sexes{1} 'Vs' Sexes{2}]) = p_vals_BL_DPlast6_NREM.Posthoc.WT.([Sexes{1} 'vs' Sexes{2}]);
	p_vals.(['Both' Sexes{1} 'and' Sexes{2}]).BL.DPlast6.REM.Posthocs.(['WT_' Sexes{1} 'Vs' Sexes{2}])  = p_vals_BL_DPlast6_REM.Posthoc.WT.([Sexes{1} 'vs' Sexes{2}]);

	% - SD DP first 6 hrs -
	p_vals.(['Both' Sexes{1} 'and' Sexes{2}]).SD.DPfirst6.Wake.Posthocs.(['WT_' Sexes{1} 'Vs' Sexes{2}]) = p_vals_SD_DPfirst6_Wake.Posthoc.WT.([Sexes{1} 'vs' Sexes{2}]);
	p_vals.(['Both' Sexes{1} 'and' Sexes{2}]).SD.DPfirst6.NREM.Posthocs.(['WT_' Sexes{1} 'Vs' Sexes{2}]) = p_vals_SD_DPfirst6_NREM.Posthoc.WT.([Sexes{1} 'vs' Sexes{2}]);
	p_vals.(['Both' Sexes{1} 'and' Sexes{2}]).SD.DPfirst6.REM.Posthocs.(['WT_' Sexes{1} 'Vs' Sexes{2}])  = p_vals_SD_DPfirst6_REM.Posthoc.WT.([Sexes{1} 'vs' Sexes{2}]);

	% - SD DP last 6 hrs -
	p_vals.(['Both' Sexes{1} 'and' Sexes{2}]).SD.DPlast6.Wake.Posthocs.(['WT_' Sexes{1} 'Vs' Sexes{2}]) = p_vals_SD_DPlast6_Wake.Posthoc.WT.([Sexes{1} 'vs' Sexes{2}]);
	p_vals.(['Both' Sexes{1} 'and' Sexes{2}]).SD.DPlast6.NREM.Posthocs.(['WT_' Sexes{1} 'Vs' Sexes{2}]) = p_vals_SD_DPlast6_NREM.Posthoc.WT.([Sexes{1} 'vs' Sexes{2}]);
	p_vals.(['Both' Sexes{1} 'and' Sexes{2}]).SD.DPlast6.REM.Posthocs.(['WT_' Sexes{1} 'Vs' Sexes{2}])  = p_vals_SD_DPlast6_REM.Posthoc.WT.([Sexes{1} 'vs' Sexes{2}]);

	
	% -- post-hocs for M vs F Mut --
	% - BL DP first 6 hrs -
	p_vals.(['Both' Sexes{1} 'and' Sexes{2}]).BL.DPfirst6.Wake.Posthocs.(['Mut_' Sexes{1} 'Vs' Sexes{2}]) = p_vals_BL_DPfirst6_Wake.Posthoc.Mut.([Sexes{1} 'vs' Sexes{2}]);
	p_vals.(['Both' Sexes{1} 'and' Sexes{2}]).BL.DPfirst6.NREM.Posthocs.(['Mut_' Sexes{1} 'Vs' Sexes{2}]) = p_vals_BL_DPfirst6_NREM.Posthoc.Mut.([Sexes{1} 'vs' Sexes{2}]);
	p_vals.(['Both' Sexes{1} 'and' Sexes{2}]).BL.DPfirst6.REM.Posthocs.(['Mut_' Sexes{1} 'Vs' Sexes{2}])  = p_vals_BL_DPfirst6_REM.Posthoc.Mut.([Sexes{1} 'vs' Sexes{2}]);

	% - BL DP last 6 hrs -
	p_vals.(['Both' Sexes{1} 'and' Sexes{2}]).BL.DPlast6.Wake.Posthocs.(['Mut_' Sexes{1} 'Vs' Sexes{2}]) = p_vals_BL_DPlast6_Wake.Posthoc.Mut.([Sexes{1} 'vs' Sexes{2}]);
	p_vals.(['Both' Sexes{1} 'and' Sexes{2}]).BL.DPlast6.NREM.Posthocs.(['Mut_' Sexes{1} 'Vs' Sexes{2}]) = p_vals_BL_DPlast6_NREM.Posthoc.Mut.([Sexes{1} 'vs' Sexes{2}]);
	p_vals.(['Both' Sexes{1} 'and' Sexes{2}]).BL.DPlast6.REM.Posthocs.(['Mut_' Sexes{1} 'Vs' Sexes{2}])  = p_vals_BL_DPlast6_REM.Posthoc.Mut.([Sexes{1} 'vs' Sexes{2}]);

	% - SD DP first 6 hrs -
	p_vals.(['Both' Sexes{1} 'and' Sexes{2}]).SD.DPfirst6.Wake.Posthocs.(['Mut_' Sexes{1} 'Vs' Sexes{2}]) = p_vals_SD_DPfirst6_Wake.Posthoc.Mut.([Sexes{1} 'vs' Sexes{2}]);
	p_vals.(['Both' Sexes{1} 'and' Sexes{2}]).SD.DPfirst6.NREM.Posthocs.(['Mut_' Sexes{1} 'Vs' Sexes{2}]) = p_vals_SD_DPfirst6_NREM.Posthoc.Mut.([Sexes{1} 'vs' Sexes{2}]);
	p_vals.(['Both' Sexes{1} 'and' Sexes{2}]).SD.DPfirst6.REM.Posthocs.(['Mut_' Sexes{1} 'Vs' Sexes{2}])  = p_vals_SD_DPfirst6_REM.Posthoc.Mut.([Sexes{1} 'vs' Sexes{2}]);

	% - SD DP last 6 hrs -
	p_vals.(['Both' Sexes{1} 'and' Sexes{2}]).SD.DPlast6.Wake.Posthocs.(['Mut_' Sexes{1} 'Vs' Sexes{2}]) = p_vals_SD_DPlast6_Wake.Posthoc.Mut.([Sexes{1} 'vs' Sexes{2}]);
	p_vals.(['Both' Sexes{1} 'and' Sexes{2}]).SD.DPlast6.NREM.Posthocs.(['Mut_' Sexes{1} 'Vs' Sexes{2}]) = p_vals_SD_DPlast6_NREM.Posthoc.Mut.([Sexes{1} 'vs' Sexes{2}]);
	p_vals.(['Both' Sexes{1} 'and' Sexes{2}]).SD.DPlast6.REM.Posthocs.(['Mut_' Sexes{1} 'Vs' Sexes{2}])  = p_vals_SD_DPlast6_REM.Posthoc.Mut.([Sexes{1} 'vs' Sexes{2}]);


	% -- BL Male --
	sig_difference_bool.(Sexes{1}).BL.DPfirst6.Wake = sig_difference_BL_DPfirst6_Wake.(Sexes{1});
	sig_difference_bool.(Sexes{1}).BL.DPfirst6.NREM = sig_difference_BL_DPfirst6_NREM.(Sexes{1});
	sig_difference_bool.(Sexes{1}).BL.DPfirst6.REM  = sig_difference_BL_DPfirst6_REM.(Sexes{1});
	sig_difference_bool.(Sexes{1}).BL.DPlast6.Wake  = sig_difference_BL_DPlast6_Wake.(Sexes{1});
	sig_difference_bool.(Sexes{1}).BL.DPlast6.NREM  = sig_difference_BL_DPlast6_NREM.(Sexes{1});
	sig_difference_bool.(Sexes{1}).BL.DPlast6.REM   = sig_difference_BL_DPlast6_REM.(Sexes{1});

	% -- BL Female --
	sig_difference_bool.(Sexes{2}).BL.DPfirst6.Wake = sig_difference_BL_DPfirst6_Wake.(Sexes{2});
	sig_difference_bool.(Sexes{2}).BL.DPfirst6.NREM = sig_difference_BL_DPfirst6_NREM.(Sexes{2});
	sig_difference_bool.(Sexes{2}).BL.DPfirst6.REM  = sig_difference_BL_DPfirst6_REM.(Sexes{2});
	sig_difference_bool.(Sexes{2}).BL.DPlast6.Wake  = sig_difference_BL_DPlast6_Wake.(Sexes{2});
	sig_difference_bool.(Sexes{2}).BL.DPlast6.NREM  = sig_difference_BL_DPlast6_NREM.(Sexes{2});
	sig_difference_bool.(Sexes{2}).BL.DPlast6.REM   = sig_difference_BL_DPlast6_REM.(Sexes{2});

	% -- SD Male -- 
	sig_difference_bool.(Sexes{1}).SD.DPfirst6.Wake = sig_difference_SD_DPfirst6_Wake.(Sexes{1});
	sig_difference_bool.(Sexes{1}).SD.DPfirst6.NREM = sig_difference_SD_DPfirst6_NREM.(Sexes{1});
	sig_difference_bool.(Sexes{1}).SD.DPfirst6.REM  = sig_difference_SD_DPfirst6_REM.(Sexes{1});
	sig_difference_bool.(Sexes{1}).SD.DPlast6.Wake  = sig_difference_SD_DPlast6_Wake.(Sexes{1});
	sig_difference_bool.(Sexes{1}).SD.DPlast6.NREM  = sig_difference_SD_DPlast6_NREM.(Sexes{1});
	sig_difference_bool.(Sexes{1}).SD.DPlast6.REM   = sig_difference_SD_DPlast6_REM.(Sexes{1});

	% -- SD Female -- 
	sig_difference_bool.(Sexes{2}).SD.DPfirst6.Wake = sig_difference_SD_DPfirst6_Wake.(Sexes{2});
	sig_difference_bool.(Sexes{2}).SD.DPfirst6.NREM = sig_difference_SD_DPfirst6_NREM.(Sexes{2});
	sig_difference_bool.(Sexes{2}).SD.DPfirst6.REM  = sig_difference_SD_DPfirst6_REM.(Sexes{2});
	sig_difference_bool.(Sexes{2}).SD.DPlast6.Wake = sig_difference_SD_DPlast6_Wake.(Sexes{2});
	sig_difference_bool.(Sexes{2}).SD.DPlast6.NREM = sig_difference_SD_DPlast6_NREM.(Sexes{2});
	sig_difference_bool.(Sexes{2}).SD.DPlast6.REM  = sig_difference_SD_DPlast6_REM.(Sexes{2});


	% -- BL Male --
	if ~Female_Data_Present ANOVA_tbls.(Sexes{1}).BL.DPfirst6.Wake = ANOVA_tbl_BL_DPfirst6_Wake.(Sexes{1}); else ANOVA_tbls.(Sexes{1}).BL.DPfirst6.Wake = []; end
	if ~Female_Data_Present ANOVA_tbls.(Sexes{1}).BL.DPfirst6.NREM = ANOVA_tbl_BL_DPfirst6_NREM.(Sexes{1}); else ANOVA_tbls.(Sexes{1}).BL.DPfirst6.NREM = []; end
	if ~Female_Data_Present ANOVA_tbls.(Sexes{1}).BL.DPfirst6.REM  = ANOVA_tbl_BL_DPfirst6_REM.(Sexes{1});  else ANOVA_tbls.(Sexes{1}).BL.DPfirst6.REM  = []; end
	if ~Female_Data_Present ANOVA_tbls.(Sexes{1}).BL.DPlast6.Wake  = ANOVA_tbl_BL_DPlast6_Wake.(Sexes{1});  else ANOVA_tbls.(Sexes{1}).BL.DPlast6.Wake  = []; end
	if ~Female_Data_Present ANOVA_tbls.(Sexes{1}).BL.DPlast6.NREM  = ANOVA_tbl_BL_DPlast6_NREM.(Sexes{1});  else ANOVA_tbls.(Sexes{1}).BL.DPlast6.NREM  = []; end
	if ~Female_Data_Present ANOVA_tbls.(Sexes{1}).BL.DPlast6.REM   = ANOVA_tbl_BL_DPlast6_REM.(Sexes{1});   else ANOVA_tbls.(Sexes{1}).BL.DPlast6.REM   = []; end

	% -- SD Male -- 
	if ~Female_Data_Present ANOVA_tbls.(Sexes{1}).SD.DPfirst6.Wake = ANOVA_tbl_SD_DPfirst6_Wake.(Sexes{1}); else ANOVA_tbls.(Sexes{1}).SD.DPfirst6.Wake = []; end
	if ~Female_Data_Present ANOVA_tbls.(Sexes{1}).SD.DPfirst6.NREM = ANOVA_tbl_SD_DPfirst6_NREM.(Sexes{1}); else ANOVA_tbls.(Sexes{1}).SD.DPfirst6.NREM = []; end
	if ~Female_Data_Present ANOVA_tbls.(Sexes{1}).SD.DPfirst6.REM  = ANOVA_tbl_SD_DPfirst6_REM.(Sexes{1});  else ANOVA_tbls.(Sexes{1}).SD.DPfirst6.REM  = []; end
	if ~Female_Data_Present ANOVA_tbls.(Sexes{1}).SD.DPlast6.Wake  = ANOVA_tbl_SD_DPlast6_Wake.(Sexes{1});  else ANOVA_tbls.(Sexes{1}).SD.DPlast6.Wake  = []; end
	if ~Female_Data_Present ANOVA_tbls.(Sexes{1}).SD.DPlast6.NREM  = ANOVA_tbl_SD_DPlast6_NREM.(Sexes{1});  else ANOVA_tbls.(Sexes{1}).SD.DPlast6.NREM  = []; end 
	if ~Female_Data_Present ANOVA_tbls.(Sexes{1}).SD.DPlast6.REM   = ANOVA_tbl_SD_DPlast6_REM.(Sexes{1});   else ANOVA_tbls.(Sexes{1}).SD.DPlast6.REM   = []; end

	% -- BL Male and Female --
	if Female_Data_Present
		ANOVA_tbls.(['Both' Sexes{1} 'and' Sexes{2}]).BL.DPfirst6.Wake = ANOVA_tbl_BL_DPfirst6_Wake.([SexVarName 'Genotype']);  
		ANOVA_tbls.(['Both' Sexes{1} 'and' Sexes{2}]).BL.DPfirst6.NREM = ANOVA_tbl_BL_DPfirst6_NREM.([SexVarName 'Genotype']);
		ANOVA_tbls.(['Both' Sexes{1} 'and' Sexes{2}]).BL.DPfirst6.REM  = ANOVA_tbl_BL_DPfirst6_REM.([SexVarName 'Genotype']); 
		ANOVA_tbls.(['Both' Sexes{1} 'and' Sexes{2}]).BL.DPlast6.Wake  = ANOVA_tbl_BL_DPlast6_Wake.([SexVarName 'Genotype']);  
		ANOVA_tbls.(['Both' Sexes{1} 'and' Sexes{2}]).BL.DPlast6.NREM  = ANOVA_tbl_BL_DPlast6_NREM.([SexVarName 'Genotype']);
		ANOVA_tbls.(['Both' Sexes{1} 'and' Sexes{2}]).BL.DPlast6.REM   = ANOVA_tbl_BL_DPlast6_REM.([SexVarName 'Genotype']); 
	end 	

	% -- SD Male and Female --
	if Female_Data_Present
		ANOVA_tbls.(['Both' Sexes{1} 'and' Sexes{2}]).SD.DPfirst6.Wake = ANOVA_tbl_SD_DPfirst6_Wake.([SexVarName 'Genotype']);  
		ANOVA_tbls.(['Both' Sexes{1} 'and' Sexes{2}]).SD.DPfirst6.NREM = ANOVA_tbl_SD_DPfirst6_NREM.([SexVarName 'Genotype']);
		ANOVA_tbls.(['Both' Sexes{1} 'and' Sexes{2}]).SD.DPfirst6.REM  = ANOVA_tbl_SD_DPfirst6_REM.([SexVarName 'Genotype']); 
		ANOVA_tbls.(['Both' Sexes{1} 'and' Sexes{2}]).SD.DPlast6.Wake  = ANOVA_tbl_SD_DPlast6_Wake.([SexVarName 'Genotype']);  
		ANOVA_tbls.(['Both' Sexes{1} 'and' Sexes{2}]).SD.DPlast6.NREM  = ANOVA_tbl_SD_DPlast6_NREM.([SexVarName 'Genotype']);
		ANOVA_tbls.(['Both' Sexes{1} 'and' Sexes{2}]).SD.DPlast6.REM   = ANOVA_tbl_SD_DPlast6_REM.([SexVarName 'Genotype']); 
	end 


else
	error('You entered an invalid choice for WindowLength.  It needs to be 12 or 6');
end
















































