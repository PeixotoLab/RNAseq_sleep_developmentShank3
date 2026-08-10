function [Twelve_hour_avg_percentages, HourlyStatePercentages, RepMeasuresANOVA_tables_hourly_percentage, RM_model] = Time_In_State_Analysis(options)
%
% 	USAGE: out = Time_In_State_Analysis(TimeInState_struct,TimeFrameRMANOVA,inBLorSD,Sex)
%
% 	This function carries out the analysis of time in each arousal state and makes corresponding figure.
%   
% 	INPUTS: 	TimeInState_struct:  A struct. Here is an example: TimeInState.WT.Male.BL.Wake
%                                    where the first level can be WT or Mut 
%                                    second level can be Male or Female
%                                    third level can be BL or SD
%                                    fourth level can be Wake, NREM, or REM 
%
%               TimeFrameRMANOVA     A 2-element vector containing the start time and end time (in hours relative to the start of recording) to which you want to restrict analysis 
%               inBLorSD             A string representing either BL or SD or SDexcSD for Baseline or sleep dep (whole day) or sleep dep day excluding the actual sleep dep
%               Sex                  A string 'Male' or 'Female' 
% ----------------------------------------------------------------                              

arguments
    options.TimeInState
    options.SD_length_hrs    = 5
    options.TimeFrameRMANOVA = []
    options.inBLorSD         = 'Both'
    options.Sex              = 'Male'
    options.MakeFigures      = 0
    options.LegendLabels     = []
end 

TimeInState_struct = options.TimeInState;
SD_length_hrs      = options.SD_length_hrs;
TimeFrameRMANOVA   = options.TimeFrameRMANOVA;
inBLorSD           = options.inBLorSD;
Sex                = options.Sex;
MakeFigures        = options.MakeFigures;
LegendLabels       = options.LegendLabels;


% If the TimeInState struct is empty (for Female and there is no Female data), just return empty structs
if all(structfun(@isempty,TimeInState_struct.WT.(Sex).BL)) & all(structfun(@isempty,TimeInState_struct.WT.(Sex).SD)) & all(structfun(@isempty,TimeInState_struct.WT.(Sex).SDexcSD))

    Twelve_hour_avg_percentages = []; 
    HourlyStatePercentages      = [];
    RepMeasuresANOVA_tables_hourly_percentage = [];
    RM_model = [];
    return 
end 






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

% Wake24hWTBL = WT_24_hr_data_struct.Wake_24h_WT_BL;
% Wake24hWTSD = WT_24_hr_data_struct.Wake_24h_WT_SD;
% NREM24hWTBL = WT_24_hr_data_struct.NREM_24h_WT_BL;
% NREM24hWTSD = WT_24_hr_data_struct.NREM_24h_WT_SD;
% REM24hWTBL  = WT_24_hr_data_struct.REM_24h_WT_BL;
% REM24hWTSD  = WT_24_hr_data_struct.REM_24h_WT_SD;

% Wake24hMutBL = Mut_24_hr_data_struct.Wake_24h_Mut_BL;
% Wake24hMutSD = Mut_24_hr_data_struct.Wake_24h_Mut_SD;
% NREM24hMutBL = Mut_24_hr_data_struct.NREM_24h_Mut_BL;
% NREM24hMutSD = Mut_24_hr_data_struct.NREM_24h_Mut_SD; 
% REM24hMutBL  = Mut_24_hr_data_struct.REM_24h_Mut_BL;
% REM24hMutSD  = Mut_24_hr_data_struct.REM_24h_Mut_SD;



% THIS CODE MAY NOT NEED TO BE IN THIS FUNCTION.  THESE VARIABLES DONT GET USED LATER
% --- WT Baseline means and SEMS -------------------------------------
WakeWTBLmean = mean(Wake24hWTBL,2,'omitnan')';
WakeWTBLstd  = std(Wake24hWTBL','omitnan')/sqrt(size(Wake24hWTBL,2));
NREMWTBLmean = mean(NREM24hWTBL,2,'omitnan')';
NREMWTBLstd  = std(NREM24hWTBL','omitnan')/sqrt(size(NREM24hWTBL,2));
REMWTBLmean  = mean(REM24hWTBL,2,'omitnan')';
REMWTBLstd   = std(REM24hWTBL','omitnan')/sqrt(size(REM24hWTBL,2));


% ---Mutant Baseline means and SEMS ------------------------------------
WakeMutBLmean = mean(Wake24hMutBL,2,'omitnan')';
WakeMutBLstd  = std(Wake24hMutBL','omitnan')/sqrt(size(Wake24hMutBL,2));
NREMMutBLmean = mean(NREM24hMutBL,2,'omitnan')';
NREMMutBLstd  = std(NREM24hMutBL','omitnan')/sqrt(size(NREM24hMutBL,2));
REMMutBLmean  = mean(REM24hMutBL,2,'omitnan')';
REMMutBLstd   = std(REM24hMutBL','omitnan')/sqrt(size(REM24hMutBL,2));

% --- WT Sleep Dep Cnt SD,rl remove the 5h SD!! -----------------------------
WakeWTSDmean = mean(Wake24hWTSD(4:24,:),2,'omitnan')';
WakeWTSDstd  = std(Wake24hWTSD(4:24,:)','omitnan')/sqrt(size(Wake24hWTSD,2));
NREMWTSDmean = mean(NREM24hWTSD(4:24,:),2,'omitnan')';
NREMWTSDstd  = std(NREM24hWTSD(4:24,:)','omitnan')/sqrt(size(NREM24hWTSD,2));
REMWTSDmean  = mean(REM24hWTSD(4:24,:),2,'omitnan')';
REMWTSDstd   = std(REM24hWTSD(4:24,:)','omitnan')/sqrt(size(REM24hWTSD,2));


% --- Mutant Sleep Dep, remove the 5h SD!! ------------------------------------- 
WakeMutSDmean = mean(Wake24hMutSD(4:24,:),2,'omitnan')';
WakeMutSDstd  = std(Wake24hMutSD(4:24,:)','omitnan')/sqrt(size(Wake24hMutSD,2));
NREMMutSDmean = mean(NREM24hMutSD(4:24,:),2,'omitnan')';
NREMMutSDstd  = std(NREM24hMutSD(4:24,:)','omitnan')/sqrt(size(NREM24hMutSD,2));
REMMutSDmean  = mean(REM24hMutSD(4:24,:),2,'omitnan')';
REMMutSDstd   = std(REM24hMutSD(4:24,:)','omitnan')/sqrt(size(REM24hMutSD,2));


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


% --- Compute averages (% time in each state) in 12-hr chunks ---
% -First 12 hours (of baseline)
% --- WT
Wake_WT_Avg_perct_First12hrsBL_by_animal = mean(Wake24hWTBL(1:12,:))/60*100;
NREM_WT_Avg_perct_First12hrsBL_by_animal = mean(NREM24hWTBL(1:12,:))/60*100;
REM_WT_Avg_perct_First12hrsBL_by_animal  = mean(REM24hWTBL(1:12,:))/60*100;

% --- Mutant
Wake_Mut_Avg_perct_First12hrsBL_by_animal = mean(Wake24hMutBL(1:12,:))/60*100;
NREM_Mut_Avg_perct_First12hrsBL_by_animal = mean(NREM24hMutBL(1:12,:))/60*100;
REM_Mut_Avg_perct_First12hrsBL_by_animal  = mean(REM24hMutBL(1:12,:))/60*100;

% -Last (second) 12 hours (of baseline)
% --- WT
Wake_WT_Avg_perct_Last12hrsBL_by_animal = mean(Wake24hWTBL(12:24,:))/60*100;
NREM_WT_Avg_perct_Last12hrsBL_by_animal = mean(NREM24hWTBL(12:24,:))/60*100;
REM_WT_Avg_perct_Last12hrsBL_by_animal  = mean(REM24hWTBL(12:24,:))/60*100;

% --- Mutant
Wake_Mut_Avg_perct_Last12hrsBL_by_animal = mean(Wake24hMutBL(12:24,:))/60*100;
NREM_Mut_Avg_perct_Last12hrsBL_by_animal = mean(NREM24hMutBL(12:24,:))/60*100;
REM_Mut_Avg_perct_Last12hrsBL_by_animal  = mean(REM24hMutBL(12:24,:))/60*100;

% ------- Sleep Dep (entire 24 hours) -------
% -First 12 hours (of Sleep Dep, where sleep dep is the second 24 hours overall)
% --- WT
Wake_WT_Avg_perct_First12hrsSD_by_animal = mean(Wake24hWTSD(1:12,:))/60*100;
NREM_WT_Avg_perct_First12hrsSD_by_animal = mean(NREM24hWTSD(1:12,:))/60*100;
REM_WT_Avg_perct_First12hrsSD_by_animal  = mean(REM24hWTSD(1:12,:))/60*100;

% --- Mutant
Wake_Mut_Avg_perct_First12hrsSD_by_animal = mean(Wake24hMutSD(1:12,:))/60*100;
NREM_Mut_Avg_perct_First12hrsSD_by_animal = mean(NREM24hMutSD(1:12,:))/60*100;
REM_Mut_Avg_perct_First12hrsSD_by_animal  = mean(REM24hMutSD(1:12,:))/60*100;

% -Last (second) 12 hours (of Sleep Dep, where sleep dep is the second 24 hours overall)
% --- WT
Wake_WT_Avg_perct_Last12hrsSD_by_animal = mean(Wake24hWTSD(12:24,:))/60*100;
NREM_WT_Avg_perct_Last12hrsSD_by_animal = mean(NREM24hWTSD(12:24,:))/60*100;
REM_WT_Avg_perct_Last12hrsSD_by_animal  = mean(REM24hWTSD(12:24,:))/60*100;

% --- Mutant
Wake_Mut_Avg_perct_Last12hrsSD_by_animal = mean(Wake24hMutSD(12:24,:))/60*100;
NREM_Mut_Avg_perct_Last12hrsSD_by_animal = mean(NREM24hMutSD(12:24,:))/60*100;
REM_Mut_Avg_perct_Last12hrsSD_by_animal  = mean(REM24hMutSD(12:24,:))/60*100;

% ------- Sleep Dep Day (excluding the actual sleep dep time) -------
% -First 12 hours (of Sleep Dep, excluding the actual sleep dep)
% --- WT
Wake_WT_Avg_perct_First12hrsSDexcSD_by_animal = mean(Wake24hWTSD(SD_length_hrs+1:12,:))/60*100;
NREM_WT_Avg_perct_First12hrsSDexcSD_by_animal = mean(NREM24hWTSD(SD_length_hrs+1:12,:))/60*100;
REM_WT_Avg_perct_First12hrsSDexcSD_by_animal  = mean(REM24hWTSD(SD_length_hrs+1:12,:))/60*100;

% --- Mutant
Wake_Mut_Avg_perct_First12hrsSDexcSD_by_animal = mean(Wake24hMutSD(SD_length_hrs+1:12,:))/60*100;
NREM_Mut_Avg_perct_First12hrsSDexcSD_by_animal = mean(NREM24hMutSD(SD_length_hrs+1:12,:))/60*100;
REM_Mut_Avg_perct_First12hrsSDexcSD_by_animal  = mean(REM24hMutSD(SD_length_hrs+1:12,:))/60*100;




% -- Build the structs to return all the data
Twelve_hour_avg_percentages = struct();
% ---- Baseline ------------
% -- WT Baseline
Twelve_hour_avg_percentages.WT.BL.First12hrs.Wake = Wake_WT_Avg_perct_First12hrsBL_by_animal;
Twelve_hour_avg_percentages.WT.BL.First12hrs.NREM = NREM_WT_Avg_perct_First12hrsBL_by_animal;
Twelve_hour_avg_percentages.WT.BL.First12hrs.REM  = REM_WT_Avg_perct_First12hrsBL_by_animal;

Twelve_hour_avg_percentages.WT.BL.Last12hrs.Wake = Wake_WT_Avg_perct_Last12hrsBL_by_animal;
Twelve_hour_avg_percentages.WT.BL.Last12hrs.NREM = NREM_WT_Avg_perct_Last12hrsBL_by_animal;
Twelve_hour_avg_percentages.WT.BL.Last12hrs.REM  = REM_WT_Avg_perct_Last12hrsBL_by_animal;

% -- Mutant Baseline 
Twelve_hour_avg_percentages.Mut.BL.First12hrs.Wake = Wake_Mut_Avg_perct_First12hrsBL_by_animal;
Twelve_hour_avg_percentages.Mut.BL.First12hrs.NREM = NREM_Mut_Avg_perct_First12hrsBL_by_animal;
Twelve_hour_avg_percentages.Mut.BL.First12hrs.REM  = REM_Mut_Avg_perct_First12hrsBL_by_animal;

Twelve_hour_avg_percentages.Mut.BL.Last12hrs.Wake = Wake_Mut_Avg_perct_Last12hrsBL_by_animal;
Twelve_hour_avg_percentages.Mut.BL.Last12hrs.NREM = NREM_Mut_Avg_perct_Last12hrsBL_by_animal;
Twelve_hour_avg_percentages.Mut.BL.Last12hrs.REM  = REM_Mut_Avg_perct_Last12hrsBL_by_animal;


% ---- Sleep Dep ------------
% -- WT Sleep Dep
Twelve_hour_avg_percentages.WT.SD.First12hrs.Wake = Wake_WT_Avg_perct_First12hrsSD_by_animal;
Twelve_hour_avg_percentages.WT.SD.First12hrs.NREM = NREM_WT_Avg_perct_First12hrsSD_by_animal;
Twelve_hour_avg_percentages.WT.SD.First12hrs.REM  = REM_WT_Avg_perct_First12hrsSD_by_animal;

Twelve_hour_avg_percentages.WT.SDexcSD.First12hrs.Wake = Wake_WT_Avg_perct_First12hrsSDexcSD_by_animal;
Twelve_hour_avg_percentages.WT.SDexcSD.First12hrs.NREM = NREM_WT_Avg_perct_First12hrsSDexcSD_by_animal;
Twelve_hour_avg_percentages.WT.SDexcSD.First12hrs.REM  = REM_WT_Avg_perct_First12hrsSDexcSD_by_animal;

Twelve_hour_avg_percentages.WT.SD.Last12hrs.Wake = Wake_WT_Avg_perct_Last12hrsSD_by_animal;
Twelve_hour_avg_percentages.WT.SD.Last12hrs.NREM = NREM_WT_Avg_perct_Last12hrsSD_by_animal;
Twelve_hour_avg_percentages.WT.SD.Last12hrs.REM  = REM_WT_Avg_perct_Last12hrsSD_by_animal;

% -- Mutant Sleep Dep 
Twelve_hour_avg_percentages.Mut.SD.First12hrs.Wake = Wake_Mut_Avg_perct_First12hrsSD_by_animal;
Twelve_hour_avg_percentages.Mut.SD.First12hrs.NREM = NREM_Mut_Avg_perct_First12hrsSD_by_animal;
Twelve_hour_avg_percentages.Mut.SD.First12hrs.REM  = REM_Mut_Avg_perct_First12hrsSD_by_animal;

Twelve_hour_avg_percentages.Mut.SDexcSD.First12hrs.Wake = Wake_Mut_Avg_perct_First12hrsSDexcSD_by_animal;
Twelve_hour_avg_percentages.Mut.SDexcSD.First12hrs.NREM = NREM_Mut_Avg_perct_First12hrsSDexcSD_by_animal;
Twelve_hour_avg_percentages.Mut.SDexcSD.First12hrs.REM  = REM_Mut_Avg_perct_First12hrsSDexcSD_by_animal;

Twelve_hour_avg_percentages.Mut.SD.Last12hrs.Wake = Wake_Mut_Avg_perct_Last12hrsSD_by_animal;
Twelve_hour_avg_percentages.Mut.SD.Last12hrs.NREM = NREM_Mut_Avg_perct_Last12hrsSD_by_animal;
Twelve_hour_avg_percentages.Mut.SD.Last12hrs.REM  = REM_Mut_Avg_perct_Last12hrsSD_by_animal;

% ---------------------------------------------------------------------------------
% --- MAKE THE FIGURES ------------------------------------------------------------
if MakeFigures
    posA=([12 0 12 100]);
    posB=([0 0 6 100]);

    myfontsize   = 16;
    mylinewidth  = 2;
    mymarkersize = 5; 


    LPboxpos = ([0 -5 12 5]);  %[x y w h] Light period box along bottom of plot
    DPboxpos = ([12 -5 12 5]); %[x y w h] Dark period box along bottom of plot
    SDboxpos = ([0 -1 SD_length_hrs 100]);  % [x y w h] Sleep Dep hatching

    LPboxposREM = ([0 -1 12 1]);  %[x y w h] Light period box along bottom of plot REM VERSION: SHORTER
    DPboxposREM = ([12 -1 12 1]); %[x y w h] Dark period box along bottom of plot  REM VERSION: SHORTER


    FTIS_BL_RelativetoTRTa=figure;
    FTIS_BL_RelativetoTRTa.Name = 'BL Wake relative to Total Recording Time';

    % rectangle('Position',posA,'FaceColor',[0.8 0.8 0.8],'EdgeColor','none');
    % hold on

    rectangle('Position',LPboxpos,'FaceColor','w','EdgeColor','k','LineWidth',2);     % LP rectangle
    hold on
    rectangle('Position',DPboxpos,'FaceColor','k','EdgeColor','k','LineWidth',2);    % DP rectangle
    e1=errorbar(WakeWTBLRelmean,WakeWTBLRelstd,'Color','k', 'LineWidth',mylinewidth,'Marker','o','MarkerFaceColor','k','MarkerSize',mymarkersize,'CapSize',0);
    e2=errorbar(WakeMutBLRelmean,WakeMutBLRelstd,'Color','r','LineWidth',mylinewidth,'Marker','o','MarkerFaceColor','r','MarkerSize',mymarkersize,'CapSize',0);
    hold off 
    ax = FTIS_BL_RelativetoTRTa.Children;
    set(ax, 'Layer','top')
    %title('Wake, Baseline');
    ylabel({'WAKE','Time in Wake (% TRT)'});
    xlabel('Time (hours)');
    
    %ax=gca;
    l = legend([e1 e2],LegendLabels);
    %l=legend(ax,LegendLabels);
    l.Box = 'off';
    l.Location = 'None';
    l.Position = [0.7326 0.1649 0.1634 0.0798];
    l.FontSize = 12;
    set(ax, 'XTick',[1 4 7 10 13 16 19 22])
    set(ax, 'XTickLabels',{'1' '4' '7' '10' '13' '16' '19' '22'})
    xlim([0 24])
    ylim([-5 100])
    
    ax.FontSize = myfontsize;
    ax.LineWidth = mylinewidth;

    h=figure;
    h.Name = 'BL NREM relative to Total Recording Time';
    % rectangle('Position',posA,'FaceColor',[0.8 0.8 0.8],'EdgeColor','none');
    % hold on
    rectangle('Position',LPboxpos,'FaceColor','w','EdgeColor','k','LineWidth',2);     % LP rectangle
    hold on
    rectangle('Position',DPboxpos,'FaceColor','k','EdgeColor','k','LineWidth',2);    % DP rectangle
    
    errorbar(NREMMutBLRelmean,NREMMutBLRelstd,'Color','r','LineWidth',mylinewidth,'Marker','o','MarkerFaceColor','r','MarkerSize',mymarkersize,'CapSize',0);
    errorbar(NREMWTBLRelmean,NREMWTBLRelstd,'Color','k','LineWidth',mylinewidth,'Marker','o','MarkerFaceColor','k','MarkerSize',mymarkersize,'CapSize',0);
    hold off
    ax = h.Children;
    set(ax, 'Layer','top')
    %title('NREM Baseline');
    % legend('Mutant','WT');
    % legend 'boxoff';
    set(ax, 'XTick',[1 4 7 10 13 16 19 22])
    set(ax, 'XTickLabels',{'1' '4' '7' '10' '13' '16' '19' '22'})
    xlim([0 24])
    ylim([-5 100])
    ylabel({'NREM','Time in NREM (% TRT)'});
    xlabel('Time (hours)')
    
    ax.FontSize = myfontsize;
    ax.LineWidth = mylinewidth;

    H=figure;
    H.Name = 'BL REM relative to Total Recording Time';
    % rectangle('Position',posA,'FaceColor',[0.8 0.8 0.8],'EdgeColor','none');
    % hold on 
    rectangle('Position',LPboxposREM,'FaceColor','w','EdgeColor','k','LineWidth',2);     % LP rectangle
    hold on
    rectangle('Position',DPboxposREM,'FaceColor','k','EdgeColor','k','LineWidth',2);    % DP rectangle
    
    errorbar(REMMutBLRelmean,REMMutBLRelstd,'Color','r','LineWidth',mylinewidth,'Marker','o','MarkerFaceColor','r','MarkerSize',mymarkersize,'CapSize',0);
    errorbar(REMWTBLRelmean,REMWTBLRelstd,'Color','k','LineWidth',mylinewidth,'Marker','o','MarkerFaceColor','k','MarkerSize',mymarkersize,'CapSize',0);
    hold off
    ax = H.Children;
    set(ax, 'Layer','top')
    %title('REM Baseline');
    % legend('Mutant','WT');
    % legend 'boxoff';
    set(ax, 'XTick',[1 4 7 10 13 16 19 22])
    set(ax, 'XTickLabels',{'1' '4' '7' '10' '13' '16' '19' '22'})
    xlim([0 24])
    ylim ([-1 20])
    xlabel('Time (hours)')
    ylabel({'REM','Time in REM (% TRT)'});
    
    ax.FontSize = myfontsize;
    ax.LineWidth = mylinewidth;


    % ------ 3-panel figure for all of BL ---------
    H=figure;
    H.Name = 'Three-Panel BL Time in State Figure';
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
    ax = H.Children.Children;
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

    % -- BL NREM --
    nexttile
    rectangle('Position',LPboxpos,'FaceColor','w','EdgeColor','k','LineWidth',2);     % LP rectangle
    hold on
    rectangle('Position',DPboxpos,'FaceColor','k','EdgeColor','k','LineWidth',2);    % DP rectangle
    
    errorbar(NREMMutBLRelmean,NREMMutBLRelstd,'Color','r','LineWidth',mylinewidth,'Marker','o','MarkerFaceColor','r','MarkerSize',mymarkersize,'CapSize',0);
    errorbar(NREMWTBLRelmean,NREMWTBLRelstd,'Color','k','LineWidth',mylinewidth,'Marker','o','MarkerFaceColor','k','MarkerSize',mymarkersize,'CapSize',0);
    hold off
    ax = H.Children.Children(1);  
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
    % ax = gca;
    ax.FontSize = myfontsize;
    ax.LineWidth = mylinewidth;

    % -- BL REM --
    nexttile
    rectangle('Position',LPboxposREM,'FaceColor','w','EdgeColor','k','LineWidth',2);     % LP rectangle
    hold on
    rectangle('Position',DPboxposREM,'FaceColor','k','EdgeColor','k','LineWidth',2);    % DP rectangle
    
    errorbar(REMMutBLRelmean,REMMutBLRelstd,'Color','r','LineWidth',mylinewidth,'Marker','o','MarkerFaceColor','r','MarkerSize',mymarkersize,'CapSize',0);
    errorbar(REMWTBLRelmean,REMWTBLRelstd,'Color','k','LineWidth',mylinewidth,'Marker','o','MarkerFaceColor','k','MarkerSize',mymarkersize,'CapSize',0);
    hold off
    ax=H.Children.Children(1);
    set(ax, 'Layer','top')
    %title('REM Baseline');
    % legend('Mutant','WT');
    % legend 'boxoff';
    set(ax, 'XTick',[1 4 7 10 13 16 19 22])
    set(ax, 'XTickLabels',{'1' '4' '7' '10' '13' '16' '19' '22'})
    xlim([0 24])
    ylim ([-1 20])
    xlabel('Time (hours)')
    ylabel({'REM','Time in REM (% TRT)'});
    %ax = gca;
    ax.FontSize = myfontsize;
    ax.LineWidth = mylinewidth;




    % ------ SD Figures ----------------------------------------------------------------
    FRS_TIS_relativetoTRT=figure;
    FRS_TIS_relativetoTRT.Name = 'SD Wake relative to Total Recording Time';


    % rectangle('Position',posA,'FaceColor',[0.8 0.8 0.8],'EdgeColor','none');
    
    rectangle('Position',LPboxpos,'FaceColor','w','EdgeColor','k','LineWidth',2);     % LP rectangle
    hold on 
    rectangle('Position',DPboxpos,'FaceColor','k','EdgeColor','k','LineWidth',2);    % DP rectangle
    rectangle('Position',[0 -5 SD_length_hrs 5],'FaceColor',[0.5 0.5 0.5],'EdgeColor','k','LineWidth',2);  % SD rectangle
    ax=FRS_TIS_relativetoTRT.Children;
    % p2 = patch([0 0 SD_length_hrs SD_length_hrs],[-5 0 0 -5],[1 1 1]);  % SD hatching
    % p2.Parent = ax;
    % hatchfill2(p2,'single','HatchAngle',45,'HatchDensity',60,'HatchLineWidth',1.5);
   
    %rectangle('Position',posB,'FaceColor',[0.5 0.5 0.5],'EdgeColor','none');
    e1=errorbar(WakeWTSDRelmean,WakeWTSDRelstd,'Color','k','LineWidth',mylinewidth,'Marker','o','MarkerFaceColor','k','MarkerSize',mymarkersize,'CapSize',0);
    e2=errorbar(WakeMutSDRelmean,WakeMutSDRelstd,'Color','r','LineWidth',mylinewidth,'Marker','o','MarkerFaceColor','r','MarkerSize',mymarkersize,'CapSize',0);
    hold off 
    set(ax, 'Layer','top')
    %title('RS Wake')
    ylabel({'WAKE','Time in WAKE (% TRT)'});
    xlabel('Time (hours)')
    l = legend([e1 e2],LegendLabels);
    
    %l = legend(ax,LegendLabels);
    l.Box = 'off';
    l.Location = 'None';
    l.Position = [0.7326 0.1649 0.1634 0.0798];
    l.FontSize = 12;
    set(ax, 'XTick',[1 4 7 10 13 16 19 22])
    set(ax, 'XTickLabels',{'1' '4' '7' '10' '13' '16' '19' '22'})
    xlim([0 24])
    ylim([-5 100])
    
    ax.FontSize = myfontsize;
    ax.LineWidth = mylinewidth;


    % -- SD NREM ---

    H=figure;
    H.Name = 'SD NREM relative to Total Recording Time';
    % rectangle('Position',posA,'FaceColor',[0.8 0.8 0.8],'EdgeColor','none');
    % hold on
    % rectangle('Position',posB,'FaceColor',[0.5 0.5 0.5],'EdgeColor','none');
    % hold on
    rectangle('Position',LPboxpos,'FaceColor','w','EdgeColor','k','LineWidth',2);     % LP rectangle
    hold on
    rectangle('Position',DPboxpos,'FaceColor','k','EdgeColor','k','LineWidth',2);    % DP rectangle
    rectangle('Position',[0 -5 SD_length_hrs 5],'FaceColor',[0.5 0.5 0.5],'EdgeColor','k','LineWidth',2);  % SD rectangle

    errorbar(NREMMutSDRelmean,NREMMutSDRelstd,'Color','r','LineWidth',mylinewidth,'Marker','o','MarkerFaceColor','r','MarkerSize',mymarkersize,'CapSize',0);
    errorbar(NREMWTSDRelmean,NREMWTSDRelstd,'Color','k','LineWidth',mylinewidth,'Marker','o','MarkerFaceColor','k','MarkerSize',mymarkersize,'CapSize',0);
    hold off
    ax=H.Children;
    set(ax, 'Layer','top')
    %title('RS NREM');
    % legend('Mutant','WT');
    % legend 'boxoff';
    set(ax, 'XTick',[1 4 7 10 13 16 19 22])
    set(ax, 'XTickLabels',{'1' '4' '7' '10' '13' '16' '19' '22'})
    xlim([0 24])
    ylim([-5 100])
    xlabel('Time (hours)')
    ylabel({'NREM','Time in NREM (% TRT)'});
    
    ax.FontSize = myfontsize;
    ax.LineWidth = mylinewidth;


    H=figure;
    H.Name = 'SD REM relative to Total Recording Time';
    % rectangle('Position',posA,'FaceColor',[0.8 0.8 0.8],'EdgeColor','none');
    % hold on
    % rectangle('Position',posB,'FaceColor',[0.5 0.5 0.5],'EdgeColor','none');
    % hold on
    rectangle('Position',LPboxposREM,'FaceColor','w','EdgeColor','k','LineWidth',2);     % LP rectangle
    hold on
    rectangle('Position',DPboxposREM,'FaceColor','k','EdgeColor','k','LineWidth',2);    % DP rectangle
    rectangle('Position',[0 -1 SD_length_hrs 1],'FaceColor',[0.5 0.5 0.5],'EdgeColor','k','LineWidth',2);  % SD rectangle

    errorbar(REMMutSDRelmean,REMMutSDRelstd,'Color','r','LineWidth',mylinewidth,'Marker','o','MarkerFaceColor','r','MarkerSize',mymarkersize,'CapSize',0);
    errorbar(REMWTSDRelmean,REMWTSDRelstd,'Color','k','LineWidth',mylinewidth,'Marker','o','MarkerFaceColor','k','MarkerSize',mymarkersize,'CapSize',0);
    hold off 
    ax = H.Chilren;
    set(ax, 'Layer','top')
    % legend('Mutant','WT')
    % legend 'boxoff';
    %title('RS REM');
    set(ax, 'XTick',[1 4 7 10 13 16 19 22])
    set(ax, 'XTickLabels',{'1' '4' '7' '10' '13' '16' '19' '22'})
    xlim([0 24])
    ylim([-1 20])
    xlabel('Time (hours)')
    ylabel({'REM','Time in REM (% TRT)'});
    
    ax.FontSize = myfontsize;
    ax.LineWidth = mylinewidth;


    % ------ 3-panel figure for all of SD ---------
    
    H=figure;
    H.Name = 'Three-Panel SD Time in State Figure';
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

    %rectangle('Position',posB,'FaceColor',[0.5 0.5 0.5],'EdgeColor','none');
    e1=errorbar(WakeWTSDRelmean,WakeWTSDRelstd,'Color','k','LineWidth',mylinewidth,'Marker','o','MarkerFaceColor','k','MarkerSize',mymarkersize,'CapSize',0);
    e2=errorbar(WakeMutSDRelmean,WakeMutSDRelstd,'Color','r','LineWidth',mylinewidth,'Marker','o','MarkerFaceColor','r','MarkerSize',mymarkersize,'CapSize',0);
    hold off 
    ax = H.Children.Children;
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

    %-- SD NREM as part of 3-panel --
    nexttile 
    rectangle('Position',LPboxpos,'FaceColor','w','EdgeColor','k','LineWidth',2);     % LP rectangle
    hold on
    rectangle('Position',DPboxpos,'FaceColor','k','EdgeColor','k','LineWidth',2);    % DP rectangle
    rectangle('Position',[0 -5 SD_length_hrs 5],'FaceColor',[0.5 0.5 0.5],'EdgeColor','k','LineWidth',2);  % SD rectangle

    errorbar(NREMMutSDRelmean,NREMMutSDRelstd,'Color','r','LineWidth',mylinewidth,'Marker','o','MarkerFaceColor','r','MarkerSize',mymarkersize,'CapSize',0);
    errorbar(NREMWTSDRelmean,NREMWTSDRelstd,'Color','k','LineWidth',mylinewidth,'Marker','o','MarkerFaceColor','k','MarkerSize',mymarkersize,'CapSize',0);
    hold off
    ax = H.Children.Children(1);
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

    % -- SD REM as part of 3-panel --
    nexttile
    rectangle('Position',LPboxposREM,'FaceColor','w','EdgeColor','k','LineWidth',2);     % LP rectangle
    hold on
    rectangle('Position',DPboxposREM,'FaceColor','k','EdgeColor','k','LineWidth',2);    % DP rectangle
    rectangle('Position',[0 -1 SD_length_hrs 1],'FaceColor',[0.5 0.5 0.5],'EdgeColor','k','LineWidth',2);  % SD rectangle

    errorbar(REMMutSDRelmean,REMMutSDRelstd,'Color','r','LineWidth',mylinewidth,'Marker','o','MarkerFaceColor','r','MarkerSize',mymarkersize,'CapSize',0);
    errorbar(REMWTSDRelmean,REMWTSDRelstd,'Color','k','LineWidth',mylinewidth,'Marker','o','MarkerFaceColor','k','MarkerSize',mymarkersize,'CapSize',0);
    hold off 
    ax = H.Children.Children(1);
    set(ax, 'Layer','top')
    % legend('Mutant','WT')
    % legend 'boxoff';
    %title('RS REM');
    set(ax, 'XTick',[1 4 7 10 13 16 19 22])
    set(ax, 'XTickLabels',{'1' '4' '7' '10' '13' '16' '19' '22'})
    xlim([0 24])
    ylim([-1 20])
    xlabel('Time (hours)')
    ylabel({'REM','Time in REM (% TRT)'});
    
    ax.FontSize = myfontsize;
    ax.LineWidth = mylinewidth;

end 
% ---------------------------------------------------------------------------------------
% ---------------------------------------------------------------------------------------



% ---------------------------------------------------------------------------------------
% --- RUN THE STATS ---------------------------------------------------------------------

% If no timeframe is selected, do the full BL or SD day (or both)
if isempty(TimeFrameRMANOVA)
    if strcmp(inBLorSD,'BL')
        % Baseline 
        % Repeated Measures Wake Percentage by Hour (Baseline)
        WakePercByHour_WT_BL = Wake24hWTBL/60;
        WakePercByHour_Mut_BL= Wake24hMutBL/60;
        %Genotype = categorical([0 0 0 0 0 1 1 1 1 1 1]');  % Categorical is important!  otherwise ranova produces incorrect F values
        Genotype = categorical([zeros(1,size(WakePercByHour_Mut_BL,2)) ones(1,size(WakePercByHour_WT_BL,2))]');  % Categorical is important!  otherwise ranova produces incorrect F values

        Y = [WakePercByHour_Mut_BL'; WakePercByHour_WT_BL'];
        t = table(Genotype,Y(:,1),Y(:,2),Y(:,3),Y(:,4),Y(:,5),Y(:,6),Y(:,7),Y(:,8),Y(:,9),Y(:,10),Y(:,11),Y(:,12), ...
                           Y(:,13),Y(:,14),Y(:,15),Y(:,16),Y(:,17),Y(:,18),Y(:,19),Y(:,20),Y(:,21),Y(:,22),Y(:,23),Y(:,24), ...
                           'VariableNames',{'Genotype','Hour1','Hour2','Hour3','Hour4','Hour5','Hour6','Hour7','Hour8',...
                            'Hour9','Hour10','Hour11','Hour12','Hour13','Hour14','Hour15','Hour16','Hour17','Hour18','Hour19','Hour20','Hour21','Hour22','Hour23','Hour24'});

        Time = 1:24;   % Within-subjects variable
        WithinStructure = table(Time','VariableNames',{'Time'});
        WithinStructure.Time = categorical(WithinStructure.Time);

        rm_W_all_BL = fitrm(t,'Hour1-Hour24 ~ Genotype','WithinDesign',WithinStructure,'WithinModel','Time');
        RManovatbl_Wake_Percentage_Hourly_BL = ranova(rm_W_all_BL,'WithinModel','Time');

        % Repeated Measures SWS Percentage by Hour (Baseline)
        NREM_PercByHour_WT_BL  = NREM24hWTBL/60;
        NREM_PercByHour_Mut_BL = NREM24hMutBL/60;
        Y = [NREM_PercByHour_Mut_BL'; NREM_PercByHour_WT_BL'];

        t = table(Genotype,Y(:,1),Y(:,2),Y(:,3),Y(:,4),Y(:,5),Y(:,6),Y(:,7),Y(:,8),Y(:,9),Y(:,10),Y(:,11),Y(:,12), ...
                           Y(:,13),Y(:,14),Y(:,15),Y(:,16),Y(:,17),Y(:,18),Y(:,19),Y(:,20),Y(:,21),Y(:,22),Y(:,23),Y(:,24), ...
                           'VariableNames',{'Genotype','Hour1','Hour2','Hour3','Hour4','Hour5','Hour6','Hour7','Hour8',...
                            'Hour9','Hour10','Hour11','Hour12','Hour13','Hour14','Hour15','Hour16','Hour17','Hour18','Hour19','Hour20','Hour21','Hour22','Hour23','Hour24'});


        rm_N_all_BL = fitrm(t,'Hour1-Hour24 ~ Genotype','WithinDesign',WithinStructure,'WithinModel','Time');
        RManovatbl_NREM_Percentage_Hourly_BL = ranova(rm_N_all_BL,'WithinModel','Time');


        % Repeated Measures REM Percentage by Hour (Baseline)
        REM_PercByHour_WT_BL  = REM24hWTBL/60;
        REM_PercByHour_Mut_BL = REM24hMutBL/60;
        Y = [REM_PercByHour_Mut_BL'; REM_PercByHour_WT_BL'];

        t = table(Genotype,Y(:,1),Y(:,2),Y(:,3),Y(:,4),Y(:,5),Y(:,6),Y(:,7),Y(:,8),Y(:,9),Y(:,10),Y(:,11),Y(:,12), ...
                           Y(:,13),Y(:,14),Y(:,15),Y(:,16),Y(:,17),Y(:,18),Y(:,19),Y(:,20),Y(:,21),Y(:,22),Y(:,23),Y(:,24), ...
                           'VariableNames',{'Genotype','Hour1','Hour2','Hour3','Hour4','Hour5','Hour6','Hour7','Hour8',...
                            'Hour9','Hour10','Hour11','Hour12','Hour13','Hour14','Hour15','Hour16','Hour17','Hour18','Hour19','Hour20','Hour21','Hour22','Hour23','Hour24'});


        rm_R_all_BL = fitrm(t,'Hour1-Hour24 ~ Genotype','WithinDesign',WithinStructure,'WithinModel','Time');
        RManovatbl_REM_Percentage_Hourly_BL = ranova(rm_R_all_BL,'WithinModel','Time');


        % Baseline
        RepMeasuresANOVA_tables_hourly_percentage.BL.Wake = RManovatbl_Wake_Percentage_Hourly_BL;
        RepMeasuresANOVA_tables_hourly_percentage.BL.NREM = RManovatbl_NREM_Percentage_Hourly_BL;
        RepMeasuresANOVA_tables_hourly_percentage.BL.REM  = RManovatbl_REM_Percentage_Hourly_BL;

        RM_Model.Wake = rm_W_all_BL;
        RM_Model.NREM = rm_N_all_BL;
        RM_Model.REM  = rm_R_all_BL;

    elseif strcmp(inBLorSD,'SD')    
        % Sleep Dep (whole 24 hours of sleep dep day)
        % Repeated Measures Wake Percentage by Hour (Sleep Dep whole 24 hours)
        WakePercByHour_WT_SD = Wake24hWTSD/60;
        WakePercByHour_Mut_SD= Wake24hMutSD/60;
        Genotype = categorical([zeros(1,size(WakePercByHour_Mut_SD,2)) ones(1,size(WakePercByHour_WT_SD,2))]');  % Categorical is important!  otherwise ranova produces incorrect F values

        Y = [WakePercByHour_Mut_SD'; WakePercByHour_WT_SD'];
        t = table(Genotype,Y(:,1),Y(:,2),Y(:,3),Y(:,4),Y(:,5),Y(:,6),Y(:,7),Y(:,8),Y(:,9),Y(:,10),Y(:,11),Y(:,12), ...
                           Y(:,13),Y(:,14),Y(:,15),Y(:,16),Y(:,17),Y(:,18),Y(:,19),Y(:,20),Y(:,21),Y(:,22),Y(:,23),Y(:,24), ...
                           'VariableNames',{'Genotype','Hour1','Hour2','Hour3','Hour4','Hour5','Hour6','Hour7','Hour8',...
                            'Hour9','Hour10','Hour11','Hour12','Hour13','Hour14','Hour15','Hour16','Hour17','Hour18','Hour19','Hour20','Hour21','Hour22','Hour23','Hour24'});

        rm_W_all_SD = fitrm(t,'Hour1-Hour24 ~ Genotype','WithinDesign',WithinStructure,'WithinModel','Time');
        RManovatbl_Wake_Percentage_Hourly_SD = ranova(rm_W_all_SD,'WithinModel','Time');

        % Repeated Measures NREM Percentage by Hour (Sleep Dep whole 24 hours)
        NREMPercByHour_WT_SD = NREM24hWTSD/60;
        NREMPercByHour_Mut_SD= NREM24hMutSD/60;

        Y = [NREMPercByHour_Mut_SD'; NREMPercByHour_WT_SD'];
        t = table(Genotype,Y(:,1),Y(:,2),Y(:,3),Y(:,4),Y(:,5),Y(:,6),Y(:,7),Y(:,8),Y(:,9),Y(:,10),Y(:,11),Y(:,12), ...
                           Y(:,13),Y(:,14),Y(:,15),Y(:,16),Y(:,17),Y(:,18),Y(:,19),Y(:,20),Y(:,21),Y(:,22),Y(:,23),Y(:,24), ...
                           'VariableNames',{'Genotype','Hour1','Hour2','Hour3','Hour4','Hour5','Hour6','Hour7','Hour8',...
                            'Hour9','Hour10','Hour11','Hour12','Hour13','Hour14','Hour15','Hour16','Hour17','Hour18','Hour19','Hour20','Hour21','Hour22','Hour23','Hour24'});


        rm_N_all_SD = fitrm(t,'Hour1-Hour24 ~ Genotype','WithinDesign',WithinStructure,'WithinModel','Time');
        RManovatbl_NREM_Percentage_Hourly_SD = ranova(rm_N_all_SD,'WithinModel','Time');

        % Repeated Measures REM Percentage by Hour (Sleep Dep whole 24 hours)
        REMPercByHour_WT_SD = REM24hWTSD/60;
        REMPercByHour_Mut_SD= REM24hMutSD/60;

        Y = [REMPercByHour_Mut_SD'; REMPercByHour_WT_SD'];
        t = table(Genotype,Y(:,1),Y(:,2),Y(:,3),Y(:,4),Y(:,5),Y(:,6),Y(:,7),Y(:,8),Y(:,9),Y(:,10),Y(:,11),Y(:,12), ...
                           Y(:,13),Y(:,14),Y(:,15),Y(:,16),Y(:,17),Y(:,18),Y(:,19),Y(:,20),Y(:,21),Y(:,22),Y(:,23),Y(:,24), ...
                           'VariableNames',{'Genotype','Hour1','Hour2','Hour3','Hour4','Hour5','Hour6','Hour7','Hour8',...
                            'Hour9','Hour10','Hour11','Hour12','Hour13','Hour14','Hour15','Hour16','Hour17','Hour18','Hour19','Hour20','Hour21','Hour22','Hour23','Hour24'});


        rm_R_all_SD = fitrm(t,'Hour1-Hour24 ~ Genotype','WithinDesign',WithinStructure,'WithinModel','Time');
        RManovatbl_REM_Percentage_Hourly_SD = ranova(rm_R_all_SD,'WithinModel','Time');

        % Sleep Dep
        RepMeasuresANOVA_tables_hourly_percentage.SD.Wake = RManovatbl_Wake_Percentage_Hourly_SD;
        RepMeasuresANOVA_tables_hourly_percentage.SD.NREM = RManovatbl_NREM_Percentage_Hourly_SD;
        RepMeasuresANOVA_tables_hourly_percentage.SD.REM  = RManovatbl_REM_Percentage_Hourly_SD;

        RM_Model.Wake = rm_W_all_SD;
        RM_Model.NREM = rm_N_all_SD;
        RM_Model.REM  = rm_R_all_SD;


    elseif strcmp(inBLorSD,'Both')
      

      % Baseline 
        % Repeated Measures Wake Percentage by Hour (Baseline)
        WakePercByHour_WT_BL = Wake24hWTBL/60;
        WakePercByHour_Mut_BL= Wake24hMutBL/60;
        %Genotype = categorical([0 0 0 0 0 1 1 1 1 1 1]');  % Categorical is important!  otherwise ranova produces incorrect F values
        Genotype = categorical([zeros(1,size(WakePercByHour_Mut_BL,2)) ones(1,size(WakePercByHour_WT_BL,2))]');  % Categorical is important!  otherwise ranova produces incorrect F values

        Y = [WakePercByHour_Mut_BL'; WakePercByHour_WT_BL'];
        t = table(Genotype,Y(:,1),Y(:,2),Y(:,3),Y(:,4),Y(:,5),Y(:,6),Y(:,7),Y(:,8),Y(:,9),Y(:,10),Y(:,11),Y(:,12), ...
                           Y(:,13),Y(:,14),Y(:,15),Y(:,16),Y(:,17),Y(:,18),Y(:,19),Y(:,20),Y(:,21),Y(:,22),Y(:,23),Y(:,24), ...
                           'VariableNames',{'Genotype','Hour1','Hour2','Hour3','Hour4','Hour5','Hour6','Hour7','Hour8',...
                            'Hour9','Hour10','Hour11','Hour12','Hour13','Hour14','Hour15','Hour16','Hour17','Hour18','Hour19','Hour20','Hour21','Hour22','Hour23','Hour24'});

        Time = 1:24;   % Within-subjects variable
        WithinStructure = table(Time','VariableNames',{'Time'});
        WithinStructure.Time = categorical(WithinStructure.Time);

        rm_W_all_BL = fitrm(t,'Hour1-Hour24 ~ Genotype','WithinDesign',WithinStructure,'WithinModel','Time');
        RManovatbl_Wake_Percentage_Hourly_BL = ranova(rm_W_all_BL,'WithinModel','Time');

        % Repeated Measures SWS Percentage by Hour (Baseline)
        NREM_PercByHour_WT_BL  = NREM24hWTBL/60;
        NREM_PercByHour_Mut_BL = NREM24hMutBL/60;
        Y = [NREM_PercByHour_Mut_BL'; NREM_PercByHour_WT_BL'];

        t = table(Genotype,Y(:,1),Y(:,2),Y(:,3),Y(:,4),Y(:,5),Y(:,6),Y(:,7),Y(:,8),Y(:,9),Y(:,10),Y(:,11),Y(:,12), ...
                           Y(:,13),Y(:,14),Y(:,15),Y(:,16),Y(:,17),Y(:,18),Y(:,19),Y(:,20),Y(:,21),Y(:,22),Y(:,23),Y(:,24), ...
                           'VariableNames',{'Genotype','Hour1','Hour2','Hour3','Hour4','Hour5','Hour6','Hour7','Hour8',...
                            'Hour9','Hour10','Hour11','Hour12','Hour13','Hour14','Hour15','Hour16','Hour17','Hour18','Hour19','Hour20','Hour21','Hour22','Hour23','Hour24'});


        rm_N_all_BL = fitrm(t,'Hour1-Hour24 ~ Genotype','WithinDesign',WithinStructure,'WithinModel','Time');
        RManovatbl_NREM_Percentage_Hourly_BL = ranova(rm_N_all_BL,'WithinModel','Time');


        % Repeated Measures REM Percentage by Hour (Baseline)
        REM_PercByHour_WT_BL  = REM24hWTBL/60;
        REM_PercByHour_Mut_BL = REM24hMutBL/60;
        Y = [REM_PercByHour_Mut_BL'; REM_PercByHour_WT_BL'];

        t = table(Genotype,Y(:,1),Y(:,2),Y(:,3),Y(:,4),Y(:,5),Y(:,6),Y(:,7),Y(:,8),Y(:,9),Y(:,10),Y(:,11),Y(:,12), ...
                           Y(:,13),Y(:,14),Y(:,15),Y(:,16),Y(:,17),Y(:,18),Y(:,19),Y(:,20),Y(:,21),Y(:,22),Y(:,23),Y(:,24), ...
                           'VariableNames',{'Genotype','Hour1','Hour2','Hour3','Hour4','Hour5','Hour6','Hour7','Hour8',...
                            'Hour9','Hour10','Hour11','Hour12','Hour13','Hour14','Hour15','Hour16','Hour17','Hour18','Hour19','Hour20','Hour21','Hour22','Hour23','Hour24'});


        rm_R_all_BL = fitrm(t,'Hour1-Hour24 ~ Genotype','WithinDesign',WithinStructure,'WithinModel','Time');
        RManovatbl_REM_Percentage_Hourly_BL = ranova(rm_R_all_BL,'WithinModel','Time');


    % Sleep Dep (whole 24 hours of sleep dep day)
        % Repeated Measures Wake Percentage by Hour (Sleep Dep whole 24 hours)
        WakePercByHour_WT_SD = Wake24hWTSD/60;
        WakePercByHour_Mut_SD= Wake24hMutSD/60;
        Genotype = categorical([zeros(1,size(WakePercByHour_Mut_SD,2)) ones(1,size(WakePercByHour_WT_SD,2))]');  % Categorical is important!  otherwise ranova produces incorrect F values

        Y = [WakePercByHour_Mut_SD'; WakePercByHour_WT_SD'];
        t = table(Genotype,Y(:,1),Y(:,2),Y(:,3),Y(:,4),Y(:,5),Y(:,6),Y(:,7),Y(:,8),Y(:,9),Y(:,10),Y(:,11),Y(:,12), ...
                           Y(:,13),Y(:,14),Y(:,15),Y(:,16),Y(:,17),Y(:,18),Y(:,19),Y(:,20),Y(:,21),Y(:,22),Y(:,23),Y(:,24), ...
                           'VariableNames',{'Genotype','Hour1','Hour2','Hour3','Hour4','Hour5','Hour6','Hour7','Hour8',...
                            'Hour9','Hour10','Hour11','Hour12','Hour13','Hour14','Hour15','Hour16','Hour17','Hour18','Hour19','Hour20','Hour21','Hour22','Hour23','Hour24'});

        rm_W_all_SD = fitrm(t,'Hour1-Hour24 ~ Genotype','WithinDesign',WithinStructure,'WithinModel','Time');
        RManovatbl_Wake_Percentage_Hourly_SD = ranova(rm_W_all_SD,'WithinModel','Time');

        % Repeated Measures NREM Percentage by Hour (Sleep Dep whole 24 hours)
        NREMPercByHour_WT_SD = NREM24hWTSD/60;
        NREMPercByHour_Mut_SD= NREM24hMutSD/60;

        Y = [NREMPercByHour_Mut_SD'; NREMPercByHour_WT_SD'];
        t = table(Genotype,Y(:,1),Y(:,2),Y(:,3),Y(:,4),Y(:,5),Y(:,6),Y(:,7),Y(:,8),Y(:,9),Y(:,10),Y(:,11),Y(:,12), ...
                           Y(:,13),Y(:,14),Y(:,15),Y(:,16),Y(:,17),Y(:,18),Y(:,19),Y(:,20),Y(:,21),Y(:,22),Y(:,23),Y(:,24), ...
                           'VariableNames',{'Genotype','Hour1','Hour2','Hour3','Hour4','Hour5','Hour6','Hour7','Hour8',...
                            'Hour9','Hour10','Hour11','Hour12','Hour13','Hour14','Hour15','Hour16','Hour17','Hour18','Hour19','Hour20','Hour21','Hour22','Hour23','Hour24'});


        rm_N_all_SD = fitrm(t,'Hour1-Hour24 ~ Genotype','WithinDesign',WithinStructure,'WithinModel','Time');
        RManovatbl_NREM_Percentage_Hourly_SD = ranova(rm_N_all_SD,'WithinModel','Time');

        % Repeated Measures REM Percentage by Hour (Sleep Dep whole 24 hours)
        REMPercByHour_WT_SD = REM24hWTSD/60;
        REMPercByHour_Mut_SD= REM24hMutSD/60;

        Y = [REMPercByHour_Mut_SD'; REMPercByHour_WT_SD'];
        t = table(Genotype,Y(:,1),Y(:,2),Y(:,3),Y(:,4),Y(:,5),Y(:,6),Y(:,7),Y(:,8),Y(:,9),Y(:,10),Y(:,11),Y(:,12), ...
                           Y(:,13),Y(:,14),Y(:,15),Y(:,16),Y(:,17),Y(:,18),Y(:,19),Y(:,20),Y(:,21),Y(:,22),Y(:,23),Y(:,24), ...
                           'VariableNames',{'Genotype','Hour1','Hour2','Hour3','Hour4','Hour5','Hour6','Hour7','Hour8',...
                            'Hour9','Hour10','Hour11','Hour12','Hour13','Hour14','Hour15','Hour16','Hour17','Hour18','Hour19','Hour20','Hour21','Hour22','Hour23','Hour24'});


        rm_R_all_SD = fitrm(t,'Hour1-Hour24 ~ Genotype','WithinDesign',WithinStructure,'WithinModel','Time');
        RManovatbl_REM_Percentage_Hourly_SD = ranova(rm_R_all_SD,'WithinModel','Time');


        % Baseline
        RepMeasuresANOVA_tables_hourly_percentage.BL.Wake = RManovatbl_Wake_Percentage_Hourly_BL;
        RepMeasuresANOVA_tables_hourly_percentage.BL.NREM = RManovatbl_NREM_Percentage_Hourly_BL;
        RepMeasuresANOVA_tables_hourly_percentage.BL.REM  = RManovatbl_REM_Percentage_Hourly_BL;

        % Sleep Dep
        RepMeasuresANOVA_tables_hourly_percentage.SD.Wake = RManovatbl_Wake_Percentage_Hourly_SD;
        RepMeasuresANOVA_tables_hourly_percentage.SD.NREM = RManovatbl_NREM_Percentage_Hourly_SD;
        RepMeasuresANOVA_tables_hourly_percentage.SD.REM  = RManovatbl_REM_Percentage_Hourly_SD;

        RM_model.BL.Wake = rm_W_all_BL;
        RM_model.BL.NREM = rm_N_all_BL;
        RM_model.BL.REM  = rm_R_all_BL;
        RM_model.SD.Wake = rm_W_all_SD;
        RM_model.SD.NREM = rm_N_all_SD;
        RM_model.SD.REM  = rm_R_all_SD;

    end % end of choosing BL or SD or Both
end     % end of if isempty(TimeFrameRMANOVA)
% ------------------------------------------------------------------------------------------------------


% Compare BL to Sleep Dep?  






% --- If requested, run a repeated measures ANOVA on a custom timeframe, given as hours from start of recording. [1 12] means first 12 hours. don't start with 0
% if nargin==3
% 	error('You called Time_In_State_Analysis.m with 3 arguments. Either 2 or 4 are acceptable.')
% end 
if ~isempty(TimeFrameRMANOVA)
    if strcmp(inBLorSD,'BL')
        % Initialize variables that will be used for all 3 RMANOVA calls.  
		% Repeated Measures Wake Percentage by Hour (Baseline)
        WakePercByHour_WT_BL = Wake24hWTBL/60;
        WakePercByHour_Mut_BL= Wake24hMutBL/60;
        %Genotype = categorical([0 0 0 0 0 1 1 1 1 1 1]');  % Categorical is important!  otherwise ranova produces incorrect F values
		Genotype = categorical([zeros(1,size(WakePercByHour_Mut_BL,2)) ones(1,size(WakePercByHour_WT_BL,2))]');  % Categorical is important!  otherwise ranova produces incorrect F values

		Time = TimeFrameRMANOVA(1):TimeFrameRMANOVA(2);   % Within-subjects variable
		WithinStructure = table(Time','VariableNames',{'Time'});
		WithinStructure.Time = categorical(WithinStructure.Time);

		rm_model_string = ['Hour',num2str(TimeFrameRMANOVA(1)),'-Hour',num2str(TimeFrameRMANOVA(2)),' ~ Genotype']; 


        Y = [WakePercByHour_Mut_BL'; WakePercByHour_WT_BL'];
    	t_W_custom_BL = table(Genotype,Y(:,1),Y(:,2),Y(:,3),Y(:,4),Y(:,5),Y(:,6),Y(:,7),Y(:,8),Y(:,9),Y(:,10),Y(:,11),Y(:,12), ...
                   	   				   Y(:,13),Y(:,14),Y(:,15),Y(:,16),Y(:,17),Y(:,18),Y(:,19),Y(:,20),Y(:,21),Y(:,22),Y(:,23),Y(:,24), ...
                   	   				     'VariableNames',{'Genotype','Hour1','Hour2','Hour3','Hour4','Hour5','Hour6','Hour7','Hour8',...
                       				     'Hour9','Hour10','Hour11','Hour12','Hour13','Hour14','Hour15','Hour16','Hour17','Hour18','Hour19','Hour20','Hour21','Hour22','Hour23','Hour24'});

		
		rm_W_custom_BL = fitrm(t_W_custom_BL,rm_model_string,'WithinDesign',WithinStructure,'WithinModel','Time');
		RManovatbl_Wake_Percentage_Hourly_BL_custom = ranova(rm_W_custom_BL,'WithinModel','Time');

		% Repeated Measures SWS Percentage by Hour (Baseline)
		NREM_PercByHour_WT_BL  = NREM24hWTBL/60;
		NREM_PercByHour_Mut_BL = NREM24hMutBL/60;
		Y = [NREM_PercByHour_Mut_BL'; NREM_PercByHour_WT_BL'];

		t_N_custom_BL = table(Genotype,Y(:,1),Y(:,2),Y(:,3),Y(:,4),Y(:,5),Y(:,6),Y(:,7),Y(:,8),Y(:,9),Y(:,10),Y(:,11),Y(:,12), ...
                   	   					Y(:,13),Y(:,14),Y(:,15),Y(:,16),Y(:,17),Y(:,18),Y(:,19),Y(:,20),Y(:,21),Y(:,22),Y(:,23),Y(:,24), ...
                       					'VariableNames',{'Genotype','Hour1','Hour2','Hour3','Hour4','Hour5','Hour6','Hour7','Hour8',...
                       					'Hour9','Hour10','Hour11','Hour12','Hour13','Hour14','Hour15','Hour16','Hour17','Hour18','Hour19','Hour20','Hour21','Hour22','Hour23','Hour24'});


		rm_N_custom_BL = fitrm(t_N_custom_BL,rm_model_string,'WithinDesign',WithinStructure,'WithinModel','Time');
		RManovatbl_NREM_Percentage_Hourly_BL_custom = ranova(rm_N_custom_BL,'WithinModel','Time');


		% Repeated Measures REM Percentage by Hour (Baseline)
		REM_PercByHour_WT_BL  = REM24hWTBL/60;
		REM_PercByHour_Mut_BL = REM24hMutBL/60;
		Y = [REM_PercByHour_Mut_BL'; REM_PercByHour_WT_BL'];

		t_R_custom_BL = table(Genotype,Y(:,1),Y(:,2),Y(:,3),Y(:,4),Y(:,5),Y(:,6),Y(:,7),Y(:,8),Y(:,9),Y(:,10),Y(:,11),Y(:,12), ...
                   	   	   			   Y(:,13),Y(:,14),Y(:,15),Y(:,16),Y(:,17),Y(:,18),Y(:,19),Y(:,20),Y(:,21),Y(:,22),Y(:,23),Y(:,24), ...
                           			   'VariableNames',{'Genotype','Hour1','Hour2','Hour3','Hour4','Hour5','Hour6','Hour7','Hour8',...
                           			   'Hour9','Hour10','Hour11','Hour12','Hour13','Hour14','Hour15','Hour16','Hour17','Hour18','Hour19','Hour20','Hour21','Hour22','Hour23','Hour24'});


		rm_R_custom_BL = fitrm(t_R_custom_BL,rm_model_string,'WithinDesign',WithinStructure,'WithinModel','Time');
		RManovatbl_REM_Percentage_Hourly_BL_custom = ranova(rm_R_custom_BL,'WithinModel','Time');

		% Baseline
		RepMeasuresANOVA_tables_hourly_percentage.BL_custom.Wake = RManovatbl_Wake_Percentage_Hourly_BL_custom;
		RepMeasuresANOVA_tables_hourly_percentage.BL_custom.NREM = RManovatbl_NREM_Percentage_Hourly_BL_custom;
		RepMeasuresANOVA_tables_hourly_percentage.BL_custom.REM  = RManovatbl_REM_Percentage_Hourly_BL_custom;

        RM_model.Wake = rm_W_custom_BL;
        RM_model.NREM = rm_N_custom_BL;
        RM_model.REM  = rm_R_custom_BL; 

	elseif strcmp(inBLorSD,'SD')
        % --- Sleep Dep (whole 24 hours, but using hours TimeFrameRMANOVA(1) to TimeFrameRMANOVA(2)) ----
        % Wake 
        WakePercByHour_WT_SD  = Wake24hWTSD/60;
        WakePercByHour_Mut_SD = Wake24hMutSD/60;
		% Sleep Dep
		% Repeated Measures Wake Percentage by Hour (Sleep Dep)
		Genotype = categorical([zeros(1,size(WakePercByHour_Mut_SD,2)) ones(1,size(WakePercByHour_WT_SD,2))]');  % Categorical is important!  otherwise ranova produces incorrect F values
		
		Time = TimeFrameRMANOVA(1):TimeFrameRMANOVA(2);   % Within-subjects variable
		
		WithinStructure = table(Time','VariableNames',{'Time'});
		WithinStructure.Time = categorical(WithinStructure.Time);

		rm_model_string = ['Hour',num2str(TimeFrameRMANOVA(1)),'-Hour',num2str(TimeFrameRMANOVA(2)),' ~ Genotype']; 

		
		
		
		Y = [WakePercByHour_Mut_SD'; WakePercByHour_WT_SD'];
		t_W_custom_SD = table(Genotype,Y(:,1),Y(:,2),Y(:,3),Y(:,4),Y(:,5),Y(:,6),Y(:,7),Y(:,8),Y(:,9),Y(:,10),Y(:,11),Y(:,12), ...
                   		   			   Y(:,13),Y(:,14),Y(:,15),Y(:,16),Y(:,17),Y(:,18),Y(:,19),Y(:,20),Y(:,21),Y(:,22),Y(:,23),Y(:,24), ...
                   		   			   'VariableNames',{'Genotype','Hour1','Hour2','Hour3','Hour4','Hour5','Hour6','Hour7','Hour8',...
                    	   			   'Hour9','Hour10','Hour11','Hour12','Hour13','Hour14','Hour15','Hour16','Hour17','Hour18','Hour19','Hour20','Hour21','Hour22','Hour23','Hour24'});

		rm_W_custom_SD = fitrm(t_W_custom_SD,rm_model_string,'WithinDesign',WithinStructure,'WithinModel','Time');
		RManovatbl_Wake_Percentage_Hourly_SD_custom = ranova(rm_W_custom_SD,'WithinModel','Time');

		% Repeated Measures NREM Percentage by Hour (Sleep Dep)
		NREMPercByHour_WT_SD = NREM24hWTSD/60;
		NREMPercByHour_Mut_SD= NREM24hMutSD/60;

		Y = [NREMPercByHour_Mut_SD'; NREMPercByHour_WT_SD'];
		t_N_custom_SD = table(Genotype,Y(:,1),Y(:,2),Y(:,3),Y(:,4),Y(:,5),Y(:,6),Y(:,7),Y(:,8),Y(:,9),Y(:,10),Y(:,11),Y(:,12), ...
                   		   			   Y(:,13),Y(:,14),Y(:,15),Y(:,16),Y(:,17),Y(:,18),Y(:,19),Y(:,20),Y(:,21),Y(:,22),Y(:,23),Y(:,24), ...
                   		   			   'VariableNames',{'Genotype','Hour1','Hour2','Hour3','Hour4','Hour5','Hour6','Hour7','Hour8',...
                    	   			   'Hour9','Hour10','Hour11','Hour12','Hour13','Hour14','Hour15','Hour16','Hour17','Hour18','Hour19','Hour20','Hour21','Hour22','Hour23','Hour24'});


		rm_N_custom_SD = fitrm(t_N_custom_SD,rm_model_string,'WithinDesign',WithinStructure,'WithinModel','Time');
		RManovatbl_NREM_Percentage_Hourly_SD_custom = ranova(rm_N_custom_SD,'WithinModel','Time');

		% Repeated Measures REM Percentage by Hour (Sleep Dep)
		REMPercByHour_WT_SD = REM24hWTSD/60;
		REMPercByHour_Mut_SD= REM24hMutSD/60;

		Y = [REMPercByHour_Mut_SD'; REMPercByHour_WT_SD'];
		t_R_custom_SD = table(Genotype,Y(:,1),Y(:,2),Y(:,3),Y(:,4),Y(:,5),Y(:,6),Y(:,7),Y(:,8),Y(:,9),Y(:,10),Y(:,11),Y(:,12), ...
                   		   			   Y(:,13),Y(:,14),Y(:,15),Y(:,16),Y(:,17),Y(:,18),Y(:,19),Y(:,20),Y(:,21),Y(:,22),Y(:,23),Y(:,24), ...
                   		   			   'VariableNames',{'Genotype','Hour1','Hour2','Hour3','Hour4','Hour5','Hour6','Hour7','Hour8',...
                    	   			   'Hour9','Hour10','Hour11','Hour12','Hour13','Hour14','Hour15','Hour16','Hour17','Hour18','Hour19','Hour20','Hour21','Hour22','Hour23','Hour24'});


		rm_R_custom_SD = fitrm(t_R_custom_SD,rm_model_string,'WithinDesign',WithinStructure,'WithinModel','Time');
		RManovatbl_REM_Percentage_Hourly_SD_custom = ranova(rm_R_custom_SD,'WithinModel','Time');
		% ------------------------------------------------------------------------------------------------------

		
		RepMeasuresANOVA_tables_hourly_percentage.SD_custom.Wake = RManovatbl_Wake_Percentage_Hourly_SD_custom;
		RepMeasuresANOVA_tables_hourly_percentage.SD_custom.NREM = RManovatbl_NREM_Percentage_Hourly_SD_custom;
		RepMeasuresANOVA_tables_hourly_percentage.SD_custom.REM  = RManovatbl_REM_Percentage_Hourly_SD_custom;

        RM_model.Wake = rm_W_custom_SD;
        RM_model.NREM = rm_N_custom_SD;
        RM_model.REM  = rm_R_custom_SD;
    

	elseif strcmp(inBLorSD,'Both')
		% --- Baseline --- 
		
		% Initialize variables that will be used for all 3 RMANOVA calls.  
		%Genotype = categorical([0 0 0 0 0 1 1 1 1 1 1]');  % Categorical is important!  otherwise ranova produces incorrect F values
		Genotype = categorical([zeros(1,size(WakePercByHour_Mut_BL,2)) ones(1,size(WakePercByHour_WT_BL,2))]');  % Categorical is important!  otherwise ranova produces incorrect F values

		Time = TimeFrameRMANOVA(1):TimeFrameRMANOVA(2);   % Within-subjects variable
		WithinStructure = table(Time','VariableNames',{'Time'});
		WithinStructure.Time = categorical(WithinStructure.Time);

		rm_model_string = ['Hour',num2str(TimeFrameRMANOVA(1)),'-Hour',num2str(TimeFrameRMANOVA(2)),' ~ Genotype']; 




		% Repeated Measures Wake Percentage by Hour (Baseline)
		WakePercByHour_WT_BL = Wake24hWTBL/60;
		WakePercByHour_Mut_BL= Wake24hMutBL/60;
		Y = [WakePercByHour_Mut_BL'; WakePercByHour_WT_BL'];
		t_W_custom_BL = table(Genotype,Y(:,1),Y(:,2),Y(:,3),Y(:,4),Y(:,5),Y(:,6),Y(:,7),Y(:,8),Y(:,9),Y(:,10),Y(:,11),Y(:,12), ...
                   	   				   Y(:,13),Y(:,14),Y(:,15),Y(:,16),Y(:,17),Y(:,18),Y(:,19),Y(:,20),Y(:,21),Y(:,22),Y(:,23),Y(:,24), ...
                   	   				   'VariableNames',{'Genotype','Hour1','Hour2','Hour3','Hour4','Hour5','Hour6','Hour7','Hour8',...
                       				   'Hour9','Hour10','Hour11','Hour12','Hour13','Hour14','Hour15','Hour16','Hour17','Hour18','Hour19','Hour20','Hour21','Hour22','Hour23','Hour24'});

		
		rm_W_custom_BL = fitrm(t_W_custom_BL,rm_model_string,'WithinDesign',WithinStructure,'WithinModel','Time');
		RManovatbl_Wake_Percentage_Hourly_BL_custom = ranova(rm_W_custom_BL,'WithinModel','Time');

		% Repeated Measures SWS Percentage by Hour (Baseline)
		NREM_PercByHour_WT_BL  = NREM24hWTBL/60;
		NREM_PercByHour_Mut_BL = NREM24hMutBL/60;
		Y = [NREM_PercByHour_Mut_BL'; NREM_PercByHour_WT_BL'];

		t_N_custom_BL = table(Genotype,Y(:,1),Y(:,2),Y(:,3),Y(:,4),Y(:,5),Y(:,6),Y(:,7),Y(:,8),Y(:,9),Y(:,10),Y(:,11),Y(:,12), ...
                   	   					Y(:,13),Y(:,14),Y(:,15),Y(:,16),Y(:,17),Y(:,18),Y(:,19),Y(:,20),Y(:,21),Y(:,22),Y(:,23),Y(:,24), ...
                       					'VariableNames',{'Genotype','Hour1','Hour2','Hour3','Hour4','Hour5','Hour6','Hour7','Hour8',...
                       					'Hour9','Hour10','Hour11','Hour12','Hour13','Hour14','Hour15','Hour16','Hour17','Hour18','Hour19','Hour20','Hour21','Hour22','Hour23','Hour24'});


		rm_N_custom_BL = fitrm(t_N_custom_BL,rm_model_string,'WithinDesign',WithinStructure,'WithinModel','Time');
		RManovatbl_NREM_Percentage_Hourly_BL_custom = ranova(rm_N_custom_BL,'WithinModel','Time');


		% Repeated Measures REM Percentage by Hour (Baseline)
		REM_PercByHour_WT_BL  = REM24hWTBL/60;
		REM_PercByHour_Mut_BL = REM24hMutBL/60;
		Y = [REM_PercByHour_Mut_BL'; REM_PercByHour_WT_BL'];

		t_R_custom_BL = table(Genotype,Y(:,1),Y(:,2),Y(:,3),Y(:,4),Y(:,5),Y(:,6),Y(:,7),Y(:,8),Y(:,9),Y(:,10),Y(:,11),Y(:,12), ...
                   	   	   			   Y(:,13),Y(:,14),Y(:,15),Y(:,16),Y(:,17),Y(:,18),Y(:,19),Y(:,20),Y(:,21),Y(:,22),Y(:,23),Y(:,24), ...
                           			   'VariableNames',{'Genotype','Hour1','Hour2','Hour3','Hour4','Hour5','Hour6','Hour7','Hour8',...
                           			   'Hour9','Hour10','Hour11','Hour12','Hour13','Hour14','Hour15','Hour16','Hour17','Hour18','Hour19','Hour20','Hour21','Hour22','Hour23','Hour24'});


		rm_R_custom_BL = fitrm(t_R_custom_BL,rm_model_string,'WithinDesign',WithinStructure,'WithinModel','Time');
		RManovatbl_REM_Percentage_Hourly_BL_custom = ranova(rm_R_custom_BL,'WithinModel','Time');

		% --- Sleep Dep ----

		% Sleep Dep
		% Repeated Measures Wake Percentage by Hour (Sleep Dep)
		
		Genotype = categorical([zeros(1,size(WakePercByHour_Mut_SD,2)) ones(1,size(WakePercByHour_WT_SD,2))]');  % Categorical is important!  otherwise ranova produces incorrect F values

		% Wake 
		WakePercByHour_WT_SD  = Wake24hWTSD/60;
		WakePercByHour_Mut_SD = Wake24hMutSD/60;
		
		Y = [WakePercByHour_Mut_SD'; WakePercByHour_WT_SD'];
		t_W_custom_SD = table(Genotype,Y(:,1),Y(:,2),Y(:,3),Y(:,4),Y(:,5),Y(:,6),Y(:,7),Y(:,8),Y(:,9),Y(:,10),Y(:,11),Y(:,12), ...
                   		   			   Y(:,13),Y(:,14),Y(:,15),Y(:,16),Y(:,17),Y(:,18),Y(:,19),Y(:,20),Y(:,21),Y(:,22),Y(:,23),Y(:,24), ...
                   		   			   'VariableNames',{'Genotype','Hour1','Hour2','Hour3','Hour4','Hour5','Hour6','Hour7','Hour8',...
                    	   			   'Hour9','Hour10','Hour11','Hour12','Hour13','Hour14','Hour15','Hour16','Hour17','Hour18','Hour19','Hour20','Hour21','Hour22','Hour23','Hour24'});

		rm_W_custom_SD = fitrm(t_W_custom_SD,rm_model_string,'WithinDesign',WithinStructure,'WithinModel','Time');
		RManovatbl_Wake_Percentage_Hourly_SD_custom = ranova(rm_W_custom_SD,'WithinModel','Time');

		% Repeated Measures NREM Percentage by Hour (Sleep Dep)
		NREMPercByHour_WT_SD = NREM24hWTSD/60;
		NREMPercByHour_Mut_SD= NREM24hMutSD/60;

		Y = [NREMPercByHour_Mut_SD'; NREMPercByHour_WT_SD'];
		t_N_custom_SD = table(Genotype,Y(:,1),Y(:,2),Y(:,3),Y(:,4),Y(:,5),Y(:,6),Y(:,7),Y(:,8),Y(:,9),Y(:,10),Y(:,11),Y(:,12), ...
                   		   			   Y(:,13),Y(:,14),Y(:,15),Y(:,16),Y(:,17),Y(:,18),Y(:,19),Y(:,20),Y(:,21),Y(:,22),Y(:,23),Y(:,24), ...
                   		   			   'VariableNames',{'Genotype','Hour1','Hour2','Hour3','Hour4','Hour5','Hour6','Hour7','Hour8',...
                    	   			   'Hour9','Hour10','Hour11','Hour12','Hour13','Hour14','Hour15','Hour16','Hour17','Hour18','Hour19','Hour20','Hour21','Hour22','Hour23','Hour24'});


		rm_N_custom_SD = fitrm(t_N_custom_SD,rm_model_string,'WithinDesign',WithinStructure,'WithinModel','Time');
		RManovatbl_NREM_Percentage_Hourly_SD_custom = ranova(rm_N_custom_SD,'WithinModel','Time');

		% Repeated Measures REM Percentage by Hour (Sleep Dep)
		REMPercByHour_WT_SD = REM24hWTSD/60;
		REMPercByHour_Mut_SD= REM24hMutSD/60;

		Y = [REMPercByHour_Mut_SD'; REMPercByHour_WT_SD'];
		t_R_custom_SD = table(Genotype,Y(:,1),Y(:,2),Y(:,3),Y(:,4),Y(:,5),Y(:,6),Y(:,7),Y(:,8),Y(:,9),Y(:,10),Y(:,11),Y(:,12), ...
                   		   			   Y(:,13),Y(:,14),Y(:,15),Y(:,16),Y(:,17),Y(:,18),Y(:,19),Y(:,20),Y(:,21),Y(:,22),Y(:,23),Y(:,24), ...
                   		   			   'VariableNames',{'Genotype','Hour1','Hour2','Hour3','Hour4','Hour5','Hour6','Hour7','Hour8',...
                    	   			   'Hour9','Hour10','Hour11','Hour12','Hour13','Hour14','Hour15','Hour16','Hour17','Hour18','Hour19','Hour20','Hour21','Hour22','Hour23','Hour24'});


		rm_R_custom_SD = fitrm(t_R_custom_SD,rm_model_string,'WithinDesign',WithinStructure,'WithinModel','Time');
		RManovatbl_REM_Percentage_Hourly_SD_custom = ranova(rm_R_custom_SD,'WithinModel','Time');
		% ------------------------------------------------------------------------------------------------------

		% Baseline
		RepMeasuresANOVA_tables_hourly_percentage.BL_custom.Wake = RManovatbl_Wake_Percentage_Hourly_BL_custom;
		RepMeasuresANOVA_tables_hourly_percentage.BL_custom.NREM = RManovatbl_NREM_Percentage_Hourly_BL_custom;
		RepMeasuresANOVA_tables_hourly_percentage.BL_custom.REM  = RManovatbl_REM_Percentage_Hourly_BL_custom;

		% SD
		RepMeasuresANOVA_tables_hourly_percentage.SD_custom.Wake = RManovatbl_Wake_Percentage_Hourly_SD_custom;
		RepMeasuresANOVA_tables_hourly_percentage.SD_custom.NREM = RManovatbl_NREM_Percentage_Hourly_SD_custom;
		RepMeasuresANOVA_tables_hourly_percentage.SD_custom.REM  = RManovatbl_REM_Percentage_Hourly_SD_custom;

        RM_Model.BL.Wake = rm_W_custom_BL;
        RM_Model.BL.NREM = rm_N_custom_BL;
        RM_Model.BL.REM  = rm_R_custom_BL;
        RM_Model.SD.Wake = rm_W_custom_SD;
        RM_Model.SD.NREM = rm_N_custom_SD;
        RM_Model.SD.REM  = rm_R_custom_SD; 

	end 
end % end of if ~isempty(TimeFrameRMANOVA)



% --------- Set up HourlyStatePercentages struct ----------------
% -- WT BL All24 hrs --
HourlyStatePercentages.Means.WT.BL.All24hrs.Wake = WakeWTBLRelmean;
HourlyStatePercentages.Means.WT.BL.All24hrs.NREM = NREMWTBLRelmean;
HourlyStatePercentages.Means.WT.BL.All24hrs.REM  = REMWTBLRelmean;

% -- WT BL First 12 hrs -- 
HourlyStatePercentages.Means.WT.BL.First12hrs.Wake = WakeWTBLRelmean(1:12);
HourlyStatePercentages.Means.WT.BL.First12hrs.NREM = NREMWTBLRelmean(1:12);
HourlyStatePercentages.Means.WT.BL.First12hrs.REM  = REMWTBLRelmean(1:12);

% -- WT BL Last 12 hrs --
HourlyStatePercentages.Means.WT.BL.Last12hrs.Wake = WakeWTBLRelmean(13:24);
HourlyStatePercentages.Means.WT.BL.Last12hrs.NREM = NREMWTBLRelmean(13:24);
HourlyStatePercentages.Means.WT.BL.Last12hrs.REM  = REMWTBLRelmean(13:24);



% -- WT SD All24 hrs --
HourlyStatePercentages.Means.WT.SD.All24hrs.Wake = WakeWTSDRelmean;
HourlyStatePercentages.Means.WT.SD.All24hrs.NREM = NREMWTSDRelmean;
HourlyStatePercentages.Means.WT.SD.All24hrs.REM  = REMWTSDRelmean;

% -- WT SD All24 hrs --
HourlyStatePercentages.Means.WT.SD.All24hrexcSD.Wake = WakeWTSDRelmean(SD_length_hrs+1:24);
HourlyStatePercentages.Means.WT.SD.All24hrexcSD.NREM = NREMWTSDRelmean(SD_length_hrs+1:24);
HourlyStatePercentages.Means.WT.SD.All24hrexcSD.REM  = REMWTSDRelmean(SD_length_hrs+1:24);

% -- WT SD First 12 hrs --
HourlyStatePercentages.Means.WT.SD.First12hrs.Wake = WakeWTSDRelmean(1:12);
HourlyStatePercentages.Means.WT.SD.First12hrs.NREM = NREMWTSDRelmean(1:12);
HourlyStatePercentages.Means.WT.SD.First12hrs.REM  = REMWTSDRelmean(1:12);
% -- WT SD Last 12 hrs --
HourlyStatePercentages.Means.WT.SD.Last12hrs.Wake = WakeWTSDRelmean(13:24);
HourlyStatePercentages.Means.WT.SD.Last12hrs.NREM = NREMWTSDRelmean(13:24);
HourlyStatePercentages.Means.WT.SD.Last12hrs.REM  = REMWTSDRelmean(13:24);



% -- Mut BL All 24 hrs -- 
HourlyStatePercentages.Means.Mut.BL.All24hrs.Wake = WakeMutBLRelmean;
HourlyStatePercentages.Means.Mut.BL.All24hrs.NREM = NREMMutBLRelmean;
HourlyStatePercentages.Means.Mut.BL.All24hrs.REM  = REMMutBLRelmean;

% -- Mut BL First 12 hrs --
HourlyStatePercentages.Means.Mut.BL.First12hrs.Wake = WakeMutBLRelmean(1:12);
HourlyStatePercentages.Means.Mut.BL.First12hrs.NREM = NREMMutBLRelmean(1:12);
HourlyStatePercentages.Means.Mut.BL.First12hrs.REM  = REMMutBLRelmean(1:12);
% -- Mut BL Last 12 hrs --
HourlyStatePercentages.Means.Mut.BL.Last12hrs.Wake = WakeMutBLRelmean(13:24);
HourlyStatePercentages.Means.Mut.BL.Last12hrs.NREM = NREMMutBLRelmean(13:24);
HourlyStatePercentages.Means.Mut.BL.Last12hrs.REM  = REMMutBLRelmean(13:24);


% -- Mut SD All 24 hrs -- 
HourlyStatePercentages.Means.Mut.SD.All24hrs.Wake = WakeMutSDRelmean;
HourlyStatePercentages.Means.Mut.SD.All24hrs.NREM = NREMMutSDRelmean;
HourlyStatePercentages.Means.Mut.SD.All24hrs.REM  = REMMutSDRelmean;

% -- Mut SD All 24 hrs (excluding actual SD)-- 
HourlyStatePercentages.Means.Mut.SD.All24hrexcSD.Wake = WakeMutSDRelmean(SD_length_hrs+1:24);
HourlyStatePercentages.Means.Mut.SD.All24hrexcSD.NREM = NREMMutSDRelmean(SD_length_hrs+1:24);
HourlyStatePercentages.Means.Mut.SD.All24hrexcSD.REM  = REMMutSDRelmean(SD_length_hrs+1:24);

% -- Mut SD First 12 hrs --
HourlyStatePercentages.Means.Mut.SD.First12hrs.Wake = WakeMutSDRelmean(1:12);
HourlyStatePercentages.Means.Mut.SD.First12hrs.NREM = NREMMutSDRelmean(1:12);
HourlyStatePercentages.Means.Mut.SD.First12hrs.REM  = REMMutSDRelmean(1:12);
% -- Mut SD Last 12 hrs -- 
HourlyStatePercentages.Means.Mut.SD.Last12hrs.Wake = WakeMutSDRelmean(13:24);
HourlyStatePercentages.Means.Mut.SD.Last12hrs.NREM = NREMMutSDRelmean(13:24);
HourlyStatePercentages.Means.Mut.SD.Last12hrs.REM  = REMMutSDRelmean(13:24);
% ---------------------------------------------------------------------------------
% ---------------------------------------------------------------------------------
% --- End of Means ----------------------------------------------------------------




% --------- Set up HourlyStatePercentages SEMs ----------------
% -- WT BL All24 hrs --
HourlyStatePercentages.SEMs.WT.BL.All24hrs.Wake = WakeWTBLRelstd;
HourlyStatePercentages.SEMs.WT.BL.All24hrs.NREM = NREMWTBLRelstd;
HourlyStatePercentages.SEMs.WT.BL.All24hrs.REM  = REMWTBLRelstd;

% -- WT BL First 12 hrs -- 
HourlyStatePercentages.SEMs.WT.BL.First12hrs.Wake = WakeWTBLRelstd(1:12);
HourlyStatePercentages.SEMs.WT.BL.First12hrs.NREM = NREMWTBLRelstd(1:12);
HourlyStatePercentages.SEMs.WT.BL.First12hrs.REM  = REMWTBLRelstd(1:12);

% -- WT BL Last 12 hrs --
HourlyStatePercentages.SEMs.WT.BL.Last12hrs.Wake = WakeWTBLRelstd(13:24);
HourlyStatePercentages.SEMs.WT.BL.Last12hrs.NREM = NREMWTBLRelstd(13:24);
HourlyStatePercentages.SEMs.WT.BL.Last12hrs.REM  = REMWTBLRelstd(13:24);



% -- WT SD All24 hrs --
HourlyStatePercentages.SEMs.WT.SD.All24hrs.Wake = WakeWTSDRelstd;
HourlyStatePercentages.SEMs.WT.SD.All24hrs.NREM = NREMWTSDRelstd;
HourlyStatePercentages.SEMs.WT.SD.All24hrs.REM  = REMWTSDRelstd;

% -- WT SD All24 hrs (excluding actual SD) --
HourlyStatePercentages.SEMs.WT.SD.All24hrexcSD.Wake = WakeWTSDRelstd(SD_length_hrs+1:24);
HourlyStatePercentages.SEMs.WT.SD.All24hrexcSD.NREM = NREMWTSDRelstd(SD_length_hrs+1:24);
HourlyStatePercentages.SEMs.WT.SD.All24hrexcSD.REM  = REMWTSDRelstd(SD_length_hrs+1:24);


% -- WT SD First 12 hrs --
HourlyStatePercentages.SEMs.WT.SD.First12hrs.Wake = WakeWTSDRelstd(1:12);
HourlyStatePercentages.SEMs.WT.SD.First12hrs.NREM = NREMWTSDRelstd(1:12);
HourlyStatePercentages.SEMs.WT.SD.First12hrs.REM  = REMWTSDRelstd(1:12);
% -- WT SD Last 12 hrs --
HourlyStatePercentages.SEMs.WT.SD.Last12hrs.Wake = WakeWTSDRelstd(13:24);
HourlyStatePercentages.SEMs.WT.SD.Last12hrs.NREM = NREMWTSDRelstd(13:24);
HourlyStatePercentages.SEMs.WT.SD.Last12hrs.REM  = REMWTSDRelstd(13:24);



% -- Mut BL All 24 hrs -- 
HourlyStatePercentages.SEMs.Mut.BL.All24hrs.Wake = WakeMutBLRelstd;
HourlyStatePercentages.SEMs.Mut.BL.All24hrs.NREM = NREMMutBLRelstd;
HourlyStatePercentages.SEMs.Mut.BL.All24hrs.REM  = REMMutBLRelstd;

% -- Mut BL First 12 hrs --
HourlyStatePercentages.SEMs.Mut.BL.First12hrs.Wake = WakeMutBLRelstd(1:12);
HourlyStatePercentages.SEMs.Mut.BL.First12hrs.NREM = NREMMutBLRelstd(1:12);
HourlyStatePercentages.SEMs.Mut.BL.First12hrs.REM  = REMMutBLRelstd(1:12);
% -- Mut BL Last 12 hrs --
HourlyStatePercentages.SEMs.Mut.BL.Last12hrs.Wake = WakeMutBLRelstd(13:24);
HourlyStatePercentages.SEMs.Mut.BL.Last12hrs.NREM = NREMMutBLRelstd(13:24);
HourlyStatePercentages.SEMs.Mut.BL.Last12hrs.REM  = REMMutBLRelstd(13:24);


% -- Mut SD All 24 hrs -- 
HourlyStatePercentages.SEMs.Mut.SD.All24hrs.Wake = WakeMutSDRelstd;
HourlyStatePercentages.SEMs.Mut.SD.All24hrs.NREM = NREMMutSDRelstd;
HourlyStatePercentages.SEMs.Mut.SD.All24hrs.REM  = REMMutSDRelstd;

% -- Mut SD All 24 hrs (without actual SD)-- 
HourlyStatePercentages.SEMs.Mut.SD.All24hrexcSD.Wake = WakeMutSDRelstd(SD_length_hrs+1:24);
HourlyStatePercentages.SEMs.Mut.SD.All24hrexcSD.NREM = NREMMutSDRelstd(SD_length_hrs+1:24);
HourlyStatePercentages.SEMs.Mut.SD.All24hrexcSD.REM  = REMMutSDRelstd(SD_length_hrs+1:24);

% -- Mut SD First 12 hrs --
HourlyStatePercentages.SEMs.Mut.SD.First12hrs.Wake = WakeMutSDRelstd(1:12);
HourlyStatePercentages.SEMs.Mut.SD.First12hrs.NREM = NREMMutSDRelstd(1:12);
HourlyStatePercentages.SEMs.Mut.SD.First12hrs.REM  = REMMutSDRelstd(1:12);
% -- Mut SD Last 12 hrs -- 
HourlyStatePercentages.SEMs.Mut.SD.Last12hrs.Wake = WakeMutSDRelstd(13:24);
HourlyStatePercentages.SEMs.Mut.SD.Last12hrs.NREM = NREMMutSDRelstd(13:24);
HourlyStatePercentages.SEMs.Mut.SD.Last12hrs.REM  = REMMutSDRelstd(13:24);
% ---------------------------------------------------------------------------------
% ---------------------------------------------------------------------------------
% --- End of SDs ----------------------------------------------------------------



% ---- Put all of the RM fitrm outputs in a struct to return ---
% RM_models.BL.Wake = rm_W_custom_BL;
% RM_models.BL.NREM = rm_N_custom_BL;
% RM_models.BL.REM  = rm_R_custom_BL;
% RM_models.SD.Wake = rm_W_custom_SD;
% RM_models.SD.NREM = rm_N_custom_SD;
% RM_models.SD.REM  = rm_R_custom_SD;






















































