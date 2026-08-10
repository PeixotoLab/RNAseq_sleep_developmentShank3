% little script to get all of the means and SEMs for Lizzy's resubmission Sept. 2024


%--- Figure 1 ---
% LP Wake
Fig1_LP_W_Mut_Males.Mean   = mean(Twelve_hour_avg.percentages.Male.Mut.BL.First12hrs.Wake)
Fig1_LP_W_Mut_Males.SEM    = std(Twelve_hour_avg.percentages.Male.Mut.BL.First12hrs.Wake)/sqrt(length(Twelve_hour_avg.percentages.Male.Mut.BL.First12hrs.Wake))

Fig1_LP_W_WT_Females.Mean  = mean(Twelve_hour_avg.percentages.Female.WT.BL.First12hrs.Wake)
Fig1_LP_W_WT_Females.SEM   = std(Twelve_hour_avg.percentages.Female.WT.BL.First12hrs.Wake)/sqrt(length(Twelve_hour_avg.percentages.Female.WT.BL.First12hrs.Wake))

Fig1_LP_W_Mut_Females.Mean = mean(Twelve_hour_avg.percentages.Female.Mut.BL.First12hrs.Wake)
Fig1_LP_W_Mut_Females.SEM  = std(Twelve_hour_avg.percentages.Female.Mut.BL.First12hrs.Wake)/sqrt(length(Twelve_hour_avg.percentages.Female.Mut.BL.First12hrs.Wake))

% DP Wake
Fig1_DP_W_WT_Males.Mean    = mean(Twelve_hour_avg.percentages.Male.WT.BL.Last12hrs.Wake)
Fig1_DP_W_WT_Males.SEM     = std(Twelve_hour_avg.percentages.Male.WT.BL.Last12hrs.Wake)/sqrt(length(Twelve_hour_avg.percentages.Male.WT.BL.Last12hrs.Wake))

Fig1_DP_W_Mut_Males.Mean   = mean(Twelve_hour_avg.percentages.Male.Mut.BL.Last12hrs.Wake)
Fig1_DP_W_Mut_Males.SEM    = std(Twelve_hour_avg.percentages.Male.Mut.BL.Last12hrs.Wake)/sqrt(length(Twelve_hour_avg.percentages.Male.Mut.BL.Last12hrs.Wake))

Fig1_DP_W_WT_Females.Mean  = mean(Twelve_hour_avg.percentages.Female.WT.BL.Last12hrs.Wake)
Fig1_DP_W_WT_Females.SEM   = std(Twelve_hour_avg.percentages.Female.WT.BL.Last12hrs.Wake)/sqrt(length(Twelve_hour_avg.percentages.Female.WT.BL.Last12hrs.Wake))

Fig1_DP_W_Mut_Females.Mean = mean(Twelve_hour_avg.percentages.Female.Mut.BL.Last12hrs.Wake)
Fig1_DP_W_Mut_Females.SEM  = std(Twelve_hour_avg.percentages.Female.Mut.BL.Last12hrs.Wake)/sqrt(length(Twelve_hour_avg.percentages.Female.Mut.BL.Last12hrs.Wake))


% LP NREM
Fig1_LP_N_WT_Females.Mean  = mean(Twelve_hour_avg.percentages.Female.WT.BL.First12hrs.NREM)
Fig1_LP_N_WT_Females.SEM   = std(Twelve_hour_avg.percentages.Female.WT.BL.First12hrs.NREM)/sqrt(length(Twelve_hour_avg.percentages.Female.WT.BL.First12hrs.NREM))

Fig1_LP_N_Mut_Females.Mean = mean(Twelve_hour_avg.percentages.Female.Mut.BL.First12hrs.NREM) 
Fig1_LP_N_Mut_Females.SEM  = std(Twelve_hour_avg.percentages.Female.Mut.BL.First12hrs.NREM)/sqrt(length(Twelve_hour_avg.percentages.Female.Mut.BL.First12hrs.NREM))


% DP NREM
Fig1_DP_N_WT_Males.Mean    = mean(Twelve_hour_avg.percentages.Male.WT.BL.Last12hrs.NREM)
Fig1_DP_N_WT_Males.SEM     =  std(Twelve_hour_avg.percentages.Male.WT.BL.Last12hrs.NREM)/sqrt(length(Twelve_hour_avg.percentages.Male.WT.BL.Last12hrs.NREM))

Fig1_DP_N_Mut_Males.Mean   = mean(Twelve_hour_avg.percentages.Male.Mut.BL.Last12hrs.NREM) 
Fig1_DP_N_Mut_Males.SEM    =  std(Twelve_hour_avg.percentages.Male.Mut.BL.Last12hrs.NREM)/sqrt(length(Twelve_hour_avg.percentages.Male.Mut.BL.Last12hrs.NREM))

Fig1_DP_N_WT_Females.Mean  = mean(Twelve_hour_avg.percentages.Female.WT.BL.Last12hrs.NREM)
Fig1_DP_N_WT_Females.SEM   =  std(Twelve_hour_avg.percentages.Female.WT.BL.Last12hrs.NREM)/sqrt(length(Twelve_hour_avg.percentages.Female.WT.BL.Last12hrs.NREM))

Fig1_DP_N_Mut_Females.Mean = mean(Twelve_hour_avg.percentages.Female.Mut.BL.Last12hrs.NREM)
Fig1_DP_N_Mut_Females.SEM  =  std(Twelve_hour_avg.percentages.Female.Mut.BL.Last12hrs.NREM)/sqrt(length(Twelve_hour_avg.percentages.Female.Mut.BL.Last12hrs.NREM))

% DP REM
Fig1_DP_R_Mut_Males.Mean  = mean(Twelve_hour_avg.percentages.Male.Mut.BL.Last12hrs.REM)
Fig1_DP_R_Mut_Males.SEM   =  std(Twelve_hour_avg.percentages.Male.Mut.BL.Last12hrs.REM)/sqrt(length(Twelve_hour_avg.percentages.Male.Mut.BL.Last12hrs.REM))

Fig1_DP_R_WT_Females.Mean = mean(Twelve_hour_avg.percentages.Female.WT.BL.Last12hrs.REM)
Fig1_DP_R_WT_Females.SEM  =  std(Twelve_hour_avg.percentages.Female.WT.BL.Last12hrs.REM)/sqrt(length(Twelve_hour_avg.percentages.Female.WT.BL.Last12hrs.REM))

Fig1_DP_R_Mut_Females.Mean = mean(Twelve_hour_avg.percentages.Female.Mut.BL.Last12hrs.REM)
Fig1_DP_R_Mut_Females.SEM  =  std(Twelve_hour_avg.percentages.Female.Mut.BL.Last12hrs.REM)/sqrt(length(Twelve_hour_avg.percentages.Female.Mut.BL.Last12hrs.REM))

% --------------------------------------------------------------------------------------------------------------------------------------------------------------------
% --- End of Figure 1 Means + SEMs -----------------------------------------------------------------------------------------------------------------------------------



% --- Figure 3 (sleep latency) -------------------------------
Fig3_WT_Males.Mean    = mean(latency_struct.Male.WT)
Fig3_WT_Males.SEM     =  std(latency_struct.Male.WT)/sqrt(length(latency_struct.Male.WT))

Fig3_Mut_Males.Mean   = mean(latency_struct.Male.Mut)
Fig3_Mut_Males.SEM    =  std(latency_struct.Male.Mut)/sqrt(length(latency_struct.Male.Mut))

Fig3_WT_Females.Mean  = mean(latency_struct.Female.WT) 
Fig3_WT_Females.SEM   =  std(latency_struct.Female.WT)/sqrt(length(latency_struct.Female.WT))

Fig3_Mut_Females.Mean = mean(latency_struct.Female.Mut)
Fig3_Mut_Females.SEM  =  std(latency_struct.Female.Mut)/sqrt(length(latency_struct.Female.Mut))
% ------------------------------------------------------------
% --- End of Figure 3 ----------------------------------------




% ---- Figure 5 (TIS Sleep Dep) -------------------------------
% LP Wake
Fig5_LP_W_Mut_Males.Mean   = mean(Twelve_hour_avg.percentages.Male.Mut.SDexcSD.First12hrs.Wake)
Fig5_LP_W_Mut_Males.SEM    =  std(Twelve_hour_avg.percentages.Male.Mut.SDexcSD.First12hrs.Wake)/sqrt(length(Twelve_hour_avg.percentages.Male.Mut.SD.First12hrs.Wake))

Fig5_LP_W_WT_Females.Mean  = mean(Twelve_hour_avg.percentages.Female.WT.SDexcSD.First12hrs.Wake)
Fig5_LP_W_WT_Females.SEM   =  std(Twelve_hour_avg.percentages.Female.WT.SDexcSD.First12hrs.Wake)/sqrt(length(Twelve_hour_avg.percentages.Female.WT.SD.First12hrs.Wake))

Fig5_LP_W_Mut_Females.Mean = mean(Twelve_hour_avg.percentages.Female.Mut.SDexcSD.First12hrs.Wake)
Fig5_LP_W_Mut_Females.SEM  =  std(Twelve_hour_avg.percentages.Female.Mut.SDexcSD.First12hrs.Wake)/sqrt(length(Twelve_hour_avg.percentages.Female.Mut.SD.First12hrs.Wake))

% DP Wake
Fig5_DP_W_WT_Males.Mean = mean(Twelve_hour_avg.percentages.Male.WT.SD.Last12hrs.Wake)
Fig5_DP_W_WT_Males.SEM  =  std(Twelve_hour_avg.percentages.Male.WT.SD.Last12hrs.Wake)/sqrt(length(Twelve_hour_avg.percentages.Male.WT.SD.Last12hrs.Wake))

Fig5_DP_W_Mut_Males.Mean = mean(Twelve_hour_avg.percentages.Male.Mut.SD.Last12hrs.Wake)
Fig5_DP_W_Mut_Males.SEM  =  std(Twelve_hour_avg.percentages.Male.Mut.SD.Last12hrs.Wake)/sqrt(length(Twelve_hour_avg.percentages.Male.Mut.SD.Last12hrs.Wake))

% LP NREM
Fig5_LP_N_Mut_Males.Mean = mean(Twelve_hour_avg.percentages.Male.Mut.SDexcSD.First12hrs.NREM)
Fig5_LP_N_Mut_Males.SEM  =  std(Twelve_hour_avg.percentages.Male.Mut.SDexcSD.First12hrs.NREM)/sqrt(length(Twelve_hour_avg.percentages.Male.Mut.SD.First12hrs.NREM))

Fig5_LP_N_WT_Females.Mean = mean(Twelve_hour_avg.percentages.Female.WT.SDexcSD.First12hrs.NREM)
Fig5_LP_N_WT_Females.SEM  =  std(Twelve_hour_avg.percentages.Female.WT.SDexcSD.First12hrs.NREM)/sqrt(length(Twelve_hour_avg.percentages.Female.WT.SD.First12hrs.NREM))

Fig5_LP_N_Mut_Females.Mean = mean(Twelve_hour_avg.percentages.Female.Mut.SDexcSD.First12hrs.NREM)
Fig5_LP_N_Mut_Females.SEM  =  std(Twelve_hour_avg.percentages.Female.Mut.SDexcSD.First12hrs.NREM)/sqrt(length(Twelve_hour_avg.percentages.Female.Mut.SD.First12hrs.NREM))

% DP NREM
Fig5_DP_N_WT_Males.Mean = mean(Twelve_hour_avg.percentages.Male.WT.SD.Last12hrs.NREM)
Fig5_DP_N_WT_Males.SEM  =  std(Twelve_hour_avg.percentages.Male.WT.SD.Last12hrs.NREM)/sqrt(length(Twelve_hour_avg.percentages.Male.WT.SD.Last12hrs.NREM))

Fig5_DP_N_Mut_Males.Mean = mean(Twelve_hour_avg.percentages.Male.Mut.SD.Last12hrs.NREM)
Fig5_DP_N_Mut_Males.SEM  =  std(Twelve_hour_avg.percentages.Male.Mut.SD.Last12hrs.NREM)/sqrt(length(Twelve_hour_avg.percentages.Male.Mut.SD.Last12hrs.NREM))

% LP REM
Fig5_LP_R_WT_Males.Mean = mean(Twelve_hour_avg.percentages.Male.WT.SDexcSD.First12hrs.REM)
Fig5_LP_R_WT_Males.SEM  =  std(Twelve_hour_avg.percentages.Male.WT.SDexcSD.First12hrs.REM)/sqrt(length(Twelve_hour_avg.percentages.Male.WT.SD.First12hrs.REM))

Fig5_LP_R_Mut_Males.Mean = mean(Twelve_hour_avg.percentages.Male.Mut.SDexcSD.First12hrs.REM)
Fig5_LP_R_Mut_Males.SEM  =  std(Twelve_hour_avg.percentages.Male.Mut.SDexcSD.First12hrs.REM)/sqrt(length(Twelve_hour_avg.percentages.Male.Mut.SD.First12hrs.REM))

Fig5_LP_R_WT_Females.Mean = mean(Twelve_hour_avg.percentages.Female.WT.SDexcSD.First12hrs.REM)
Fig5_LP_R_WT_Females.SEM  =  std(Twelve_hour_avg.percentages.Female.WT.SDexcSD.First12hrs.REM)/sqrt(length(Twelve_hour_avg.percentages.Female.WT.SD.First12hrs.REM))

Fig5_LP_R_Mut_Females.Mean = mean(Twelve_hour_avg.percentages.Female.Mut.SDexcSD.First12hrs.REM)
Fig5_LP_R_Mut_Females.SEM  =  std(Twelve_hour_avg.percentages.Female.Mut.SDexcSD.First12hrs.REM)/sqrt(length(Twelve_hour_avg.percentages.Female.Mut.SD.First12hrs.REM))





























