function fig_handle = plot_time_in_state_M_vs_F_WT_vs_Mut(Twelve_hour_avg_percentages_Male,Twelve_hour_avg_percentages_Female,p_vals_struct,BLorSD,First12hrsOrLast12hrs,SleepStage)
%
% USAGE: fig_handle = plot_time_in_state_M_vs_F_WT_vs_Mut(Twelve_hour_avg_percentages_Male,Twelve_hour_avg_percentages_Female,BLorSD,First12hrsOrLast12hrs,SleepStage)
%
% This function makes a plot with dots at 4 horizontal locations: Male+WT, Male+Mutant, Female+WT, Female+Mutant, 
% vertical axis is  percentage of time in the particular arousal state (comparing means across Sex+Genotype combos)
% either for the light phase (first12hrs) or dark phase (Last12hrs) and BL or SD.
%
% INPUTS: SleepStage: options are 'Wake', 'NREM' or 'REM'


if strcmp(SleepStage,'WAKE')
	SleepStage_for_ylabel = upper(SleepStage);  % convert Wake to WAKE for y label boldface part
else
	SleepStage_for_ylabel = SleepStage; 
end 

% Extract the relevant p-values from the p-vals struct
if strcmp(First12hrsOrLast12hrs,'First12hrs')
  	LPorDP = 'LP';
  elseif strcmp(First12hrsOrLast12hrs,'Last12hrs')
  	LPorDP = 'DP';
  else
  	error('Invalid entry for First12hrsOrLast12hrs in plot_time_in_state_M_vs_F_WT_vs_Mut.')
end 
if isfield(p_vals_struct.(BLorSD),LPorDP) % The sleep dep day does not have p-values for the light phase, so set all p-vals to empty in this case.
	p_val.Male.WTvsMut     = p_vals_struct.(BLorSD).(LPorDP).(SleepStage).Posthoc.Male.WTvsMut;
	p_val.Female.WTvsMut   = p_vals_struct.(BLorSD).(LPorDP).(SleepStage).Posthoc.Female.WTvsMut; 
	p_val.WT.MaleVsFemale  = p_vals_struct.(BLorSD).(LPorDP).(SleepStage).Posthoc.WT.MalevsFemale;
	p_val.Mut.MaleVsFemale = p_vals_struct.(BLorSD).(LPorDP).(SleepStage).Posthoc.Mut.MalevsFemale; 
else 
	p_val.Male.WTvsMut     = [];
	p_val.Female.WTvsMut   = []; 
	p_val.WT.MaleVsFemale  = [];
	p_val.Mut.MaleVsFemale = []; 
end 
% --- Counts per group -----
N_male_WT   = length(Twelve_hour_avg_percentages_Male.WT.(BLorSD).(First12hrsOrLast12hrs).(SleepStage));
N_female_WT = length(Twelve_hour_avg_percentages_Female.WT.(BLorSD).(First12hrsOrLast12hrs).(SleepStage)); 

N_male_Mut   = length(Twelve_hour_avg_percentages_Male.Mut.(BLorSD).(First12hrsOrLast12hrs).(SleepStage));
N_female_Mut = length(Twelve_hour_avg_percentages_Female.Mut.(BLorSD).(First12hrsOrLast12hrs).(SleepStage)); 

% --- Means ------------------------------
Mean_M_WT  = mean(Twelve_hour_avg_percentages_Male.WT.(BLorSD).(First12hrsOrLast12hrs).(SleepStage));
Mean_M_Mut = mean(Twelve_hour_avg_percentages_Male.Mut.(BLorSD).(First12hrsOrLast12hrs).(SleepStage));
Mean_F_WT  = mean(Twelve_hour_avg_percentages_Female.WT.(BLorSD).(First12hrsOrLast12hrs).(SleepStage));
Mean_F_Mut = mean(Twelve_hour_avg_percentages_Female.Mut.(BLorSD).(First12hrsOrLast12hrs).(SleepStage));


% -- First 12 hours --
% Mean_M_WT_first12hr_Wake  = mean(Twelve_hour_avg_percentages_Male.WT.BL.First12hrs.Wake);
% Mean_M_Mut_first12hr_Wake = mean(Twelve_hour_avg_percentages_Male.Mut.BL.First12hrs.Wake);
% Mean_F_WT_first12hr_Wake  = mean(Twelve_hour_avg_percentages_Female.WT.BL.First12hrs.Wake);
% Mean_F_Mut_first12hr_Wake = mean(Twelve_hour_avg_percentages_Female.Mut.BL.First12hrs.Wake);

% Mean_M_WT_first12hr_NREM  = mean(Twelve_hour_avg_percentages_Male.WT.BL.First12hrs.NREM);
% Mean_M_Mut_first12hr_NREM = mean(Twelve_hour_avg_percentages_Male.Mut.BL.First12hrs.NREM);
% Mean_F_WT_first12hr_NREM  = mean(Twelve_hour_avg_percentages_Female.WT.BL.First12hrs.NREM);
% Mean_F_Mut_first12hr_NREM = mean(Twelve_hour_avg_percentages_Female.Mut.BL.First12hrs.NREM);


% Mean_M_WT_first12hr_REM  = mean(Twelve_hour_avg_percentages_Male.WT.BL.First12hrs.REM);
% Mean_M_Mut_first12hr_REM = mean(Twelve_hour_avg_percentages_Male.Mut.BL.First12hrs.REM);
% Mean_F_WT_first12hr_REM  = mean(Twelve_hour_avg_percentages_Female.WT.BL.First12hrs.REM);
% Mean_F_Mut_first12hr_REM = mean(Twelve_hour_avg_percentages_Female.Mut.BL.First12hrs.REM);

% % -- Second 12 hours --
% Mean_M_WT_last12hr_Wake  = mean(Twelve_hour_avg_percentages_Male.WT.BL.Last12hrs.Wake);
% Mean_M_Mut_last12hr_Wake = mean(Twelve_hour_avg_percentages_Male.Mut.BL.Last12hrs.Wake);
% Mean_F_WT_last12hr_Wake  = mean(Twelve_hour_avg_percentages_Female.WT.BL.Last12hrs.Wake);
% Mean_F_Mut_last12hr_Wake = mean(Twelve_hour_avg_percentages_Female.Mut.BL.Last12hrs.Wake);

% Mean_M_WT_last12hr_NREM  = mean(Twelve_hour_avg_percentages_Male.WT.BL.Last12hrs.NREM);
% Mean_M_Mut_last12hr_NREM = mean(Twelve_hour_avg_percentages_Male.Mut.BL.Last12hrs.NREM);
% Mean_F_WT_last12hr_NREM  = mean(Twelve_hour_avg_percentages_Female.WT.BL.Last12hrs.NREM);
% Mean_F_Mut_last12hr_NREM = mean(Twelve_hour_avg_percentages_Female.Mut.BL.Last12hrs.NREM);

% Mean_M_WT_last12hr_REM  = mean(Twelve_hour_avg_percentages_Male.WT.BL.Last12hrs.REM);
% Mean_M_Mut_last12hr_REM = mean(Twelve_hour_avg_percentages_Male.Mut.BL.Last12hrs.REM);
% Mean_F_WT_last12hr_REM  = mean(Twelve_hour_avg_percentages_Female.WT.BL.Last12hrs.REM);
% Mean_F_Mut_last12hr_REM = mean(Twelve_hour_avg_percentages_Female.Mut.BL.Last12hrs.REM);

% -----------------------------------------------------------
figure
f=gcf;
f.Renderer = 'painters';
p=plot(1*ones(1,N_male_WT),Twelve_hour_avg_percentages_Male.WT.(BLorSD).(First12hrsOrLast12hrs).(SleepStage),  'ko',2*ones(1,N_male_Mut),Twelve_hour_avg_percentages_Male.Mut.(BLorSD).(First12hrsOrLast12hrs).(SleepStage),  'ro', ...
	 3*ones(1,N_female_WT),Twelve_hour_avg_percentages_Female.WT.(BLorSD).(First12hrsOrLast12hrs).(SleepStage),'ko',4*ones(1,N_female_Mut),Twelve_hour_avg_percentages_Female.Mut.(BLorSD).(First12hrsOrLast12hrs).(SleepStage),'ro');
ax=gca();
ax.XLim = [0 5];
ax.YLim = [0 100];
ax.YTick = 0:25:100;
ax.XTick = [];
ax.FontSize = 14;
ax.Box = 'off';
ax.Color = 'none';
ax.LineWidth = 2;
ax.TickDir = 'out';
ax.YLabel.String = {strcat('{\bf ',SleepStage_for_ylabel,'}');['Time in ', SleepStage, ' (% TRT)']};
ax.XTick = [1.5 3.5];
ax.XTickLabel = {'Males','Females'};
ax.XAxis.TickLength = [0 0];
if strcmp(First12hrsOrLast12hrs,'Last12hrs')
	ax.Color = [0.9 0.9 0.9];  % Use a light gray background if plotting dark phase
end  
for i=1:4
	p(i).LineWidth = 3;
	p(i).MarkerSize = 12;
end 
p(1).MarkerFaceColor = 'k';
p(2).MarkerFaceColor = 'r';
hold on
line1 = line([0.85 1.15],[Mean_M_WT Mean_M_WT]);
line1.Color =[0.8 0.8 0.8];
line1.LineWidth = 6;
line2 = line([1.85 2.15],[Mean_M_Mut Mean_M_Mut]);
line2.Color =[0.8 0.8 0.8];
line2.LineWidth = 6;
line3 = line([2.85 3.15],[Mean_F_WT Mean_F_WT]);
line3.Color =[0.8 0.8 0.8];
line3.LineWidth = 6;
line4 = line([3.85 4.15],[Mean_F_Mut Mean_F_Mut]);
line4.Color =[0.8 0.8 0.8];
line4.LineWidth = 6;


% --- add * and # if significant interactions ---
vert_fudge_factor  = 5;
asterisk_FontSize  = 24;
sig_line_thickness = 1.5;

% Male WTvsMut
if p_val.Male.WTvsMut < 0.05    & p_val.Male.WTvsMut >=0.01
	symbol.MaleWTVsMut = '*';
elseif p_val.Male.WTvsMut <0.01 & p_val.Male.WTvsMut >=0.001
	symbol.MaleWTVsMut = '**';
elseif p_val.Male.WTvsMut < 0.001
	symbol.MaleWTVsMut = '***';
end 

if p_val.Male.WTvsMut < 0.05
	y_loc = max([Twelve_hour_avg_percentages_Male.WT.(BLorSD).(First12hrsOrLast12hrs).(SleepStage) ...
	 			 Twelve_hour_avg_percentages_Male.Mut.(BLorSD).(First12hrsOrLast12hrs).(SleepStage)])+vert_fudge_factor;
	
	text(ax,1.5,y_loc,symbol.MaleWTVsMut,'FontSize',asterisk_FontSize,'HorizontalAlignment','center');
	plot(ax,[1,2],[y_loc-2,y_loc-2],'k','LineWidth',sig_line_thickness);
end 


% Female WTvsMut
if p_val.Female.WTvsMut < 0.05    & p_val.Female.WTvsMut >=0.01
	symbol.FemaleWTVsMut = '*';
elseif p_val.Female.WTvsMut <0.01 & p_val.Female.WTvsMut >=0.001
	symbol.FemaleWTVsMut = '**';
elseif p_val.Female.WTvsMut < 0.001
	symbol.FemaleWTVsMut = '***';
end 

if p_val.Female.WTvsMut < 0.05
	y_loc = max([Twelve_hour_avg_percentages_Female.WT.(BLorSD).(First12hrsOrLast12hrs).(SleepStage) ...
	 			 Twelve_hour_avg_percentages_Female.Mut.(BLorSD).(First12hrsOrLast12hrs).(SleepStage)])+vert_fudge_factor;
	
	text(ax,3.5,y_loc,symbol.FemaleWTVsMut,'FontSize',asterisk_FontSize,'HorizontalAlignment','center');
	plot(ax,[3,4],[y_loc-2,y_loc-2],'k','LineWidth',sig_line_thickness);
end 

% WT Male vs Female
if p_val.WT.MaleVsFemale < 0.05    & p_val.WT.MaleVsFemale >=0.01
	symbol.WTMaleVsFemale = '#';
elseif p_val.WT.MaleVsFemale <0.01 & p_val.WT.MaleVsFemale >=0.001
	symbol.WTMaleVsFemale = '##';
elseif p_val.WT.MaleVsFemale < 0.001
	symbol.WTMaleVsFemale = '###';
end 

if p_val.WT.MaleVsFemale < 0.05
	y_loc = max([Twelve_hour_avg_percentages_Male.WT.(BLorSD).(First12hrsOrLast12hrs).(SleepStage) ...
	 			 Twelve_hour_avg_percentages_Female.WT.(BLorSD).(First12hrsOrLast12hrs).(SleepStage)])+vert_fudge_factor;
	
	text(ax,2,y_loc,symbol.WTMaleVsFemale,'FontSize',asterisk_FontSize,'Color','red','HorizontalAlignment','center');
	plot(ax,[1,3],[y_loc-2,y_loc-2],'r','LineWidth',sig_line_thickness);
end

% Mut Male vs Female
if p_val.Mut.MaleVsFemale < 0.05    & p_val.Mut.MaleVsFemale >=0.01
	symbol.MutMaleVsFemale = '#';
elseif p_val.Mut.MaleVsFemale <0.01 & p_val.Mut.MaleVsFemale >=0.001
	symbol.MutMaleVsFemale = '##';
elseif p_val.Mut.MaleVsFemale < 0.001
	symbol.MutMaleVsFemale = '###';
end 

if p_val.Mut.MaleVsFemale < 0.05
	y_loc = max([Twelve_hour_avg_percentages_Male.Mut.(BLorSD).(First12hrsOrLast12hrs).(SleepStage) ...
	 			 Twelve_hour_avg_percentages_Female.Mut.(BLorSD).(First12hrsOrLast12hrs).(SleepStage)])+vert_fudge_factor;
	
	text(ax,3,y_loc,symbol.MutMaleVsFemale,'FontSize',asterisk_FontSize,'Color','red','HorizontalAlignment','center');
	plot(ax,[2,4],[y_loc-2,y_loc-2],'r','LineWidth',sig_line_thickness);
end

hold off
if strcmp(SleepStage,'Wake')  % Only add the legend if Wake (since it's the top graph)
	l=legend([p(1) p(2)],'Wildtype','Shank3^{\DeltaC}');
	l.Box = 'off';
	l.Location = 'northwest';
	l.FontSize = 14;
end 
ax.Title.String = strcat(SleepStage,{' '},BLorSD,{' '},First12hrsOrLast12hrs); % Add an informative title (feel free to delete this line)




% % put all of the figures as subplots in a bigger figure
% % Get a list of all of the open figures
% % figlist=get(groot,'Children');
 
% % newfig=figure;
% % tcl=tiledlayout(newfig,'flow')
 
% % for i = 1:numel(figlist)
% %     figure(figlist(i));
% %     ax=gca;
% %     ax.Parent=tcl;
% %     ax.Layout.Tile=i;
% % end

fig_handle = f;  % So you return the figure handle
