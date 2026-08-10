function fig_handle = plot_output_var_M_vs_F_WT_vs_Mut(Twelve_hour_avg_values_Male,Twelve_hour_avg_values_Female,outcome_var,p_vals_struct,BLorSD,First12hrsOrLast12hrs,SleepStage,Sexes,LegendLabels)
%
% USAGE: fig_handle = plot_output_var_M_vs_F_WT_vs_Mut(Twelve_hour_avg_values_Male,Twelve_hour_avg_values_Female,outcome_var,p_vals_struct,BLorSD,First12hrsOrLast12hrs,SleepStage)
%
% This function makes a plot with dots at 4 horizontal locations: Male+WT, Male+Mutant, Female+WT, Female+Mutant, 
% vertical axis is  one of the following (based on outcome_var): 
% 							percentage of time in the particular arousal state (comparing means across Sex+Genotype combos)
% 							bout counts    
% 							bout durations 

% either for the light phase (first12hrs) or dark phase (Last12hrs) and BL or SD.
%
% INPUTS: SleepStage: 		options are 'Wake', 'NREM' or 'REM'
%					LegendLabels: 	a 2-element cell array like {'WT','Mut'} for instance


% Handle the case where no female data is present
if isempty(Twelve_hour_avg_values_Female)
	NoFemaleData = true; 
	Twelve_hour_avg_values_Female.WT.(BLorSD).(First12hrsOrLast12hrs).(SleepStage) =[];
	Twelve_hour_avg_values_Female.Mut.(BLorSD).(First12hrsOrLast12hrs).(SleepStage) =[];
else 
	NoFemaleData = false;
end 

if strcmp(SleepStage,'Wake')
	SleepStage_for_ylabel = upper(SleepStage);  % convert Wake to WAKE for y label boldface part
else
	SleepStage_for_ylabel = SleepStage; 
end 


% Set up labels for the vertical axis, based on the output variable and the arousal state
if strcmp(outcome_var,'Percentages')
	y_label_string = {SleepStage_for_ylabel;['(% TRT)']};
elseif strcmp(outcome_var,'Bout_Counts')
	y_label_string = {SleepStage_for_ylabel;['Avg Bouts/hr']};
elseif strcmp(outcome_var,'Bout_Durations')
	y_label_string = {SleepStage_for_ylabel;['Bout Duration (min)']};
else
	error('In plot_output_var_M_vs_F_WT_vs_Mut.m: you entered an outcome_var that is invalid.  Options are ''Percentages'', ''Bout_Counts'', or ''Bout_Durations'' ')
end 



% Extract the relevant p-values from the p-vals struct
if strcmp(First12hrsOrLast12hrs,'First12hrs')
  	LPorDP = 'LP';
  elseif strcmp(First12hrsOrLast12hrs,'Last12hrs')
  	LPorDP = 'DP';
  elseif strcmp(First12hrsOrLast12hrs,'DPfirst6')
  	LPorDP = 'DPfirst6';
  elseif strcmp(First12hrsOrLast12hrs,'DPlast6')
  	LPorDP = 'DPlast6';
  else
  	error('Invalid entry for First12hrsOrLast12hrs in plot_time_in_state_M_vs_F_WT_vs_Mut. Options are First12hrs,Last12hrs,DPfirst6hrs,DPlast6hrs.')
end 


% p_val.Male.WTvsMut     = p_vals_struct.(BLorSD).(LPorDP).(SleepStage).Posthoc.Male.WTvsMut;
% p_val.Female.WTvsMut   = p_vals_struct.(BLorSD).(LPorDP).(SleepStage).Posthoc.Female.WTvsMut; 
% p_val.WT.MaleVsFemale  = p_vals_struct.(BLorSD).(LPorDP).(SleepStage).Posthoc.WT.MalevsFemale;
% p_val.Mut.MaleVsFemale = p_vals_struct.(BLorSD).(LPorDP).(SleepStage).Posthoc.Mut.MalevsFemale; 

if ~strcmp(outcome_var,'Percentages')  % If not doing Percentages
	if NoFemaleData
		p_val.(Sexes{1}).WTvsMut.Adjusted   = p_vals_struct.(Sexes{1}).(BLorSD).(LPorDP).(SleepStage).ANOVA.Genotype; % plotting code needs adj or un, but if only one sex, these are the same
		p_val.(Sexes{1}).WTvsMut.Unadjusted = p_vals_struct.(Sexes{1}).(BLorSD).(LPorDP).(SleepStage).ANOVA.Genotype;
	else 
		p_val.(Sexes{1}).WTvsMut = p_vals_struct.(Sexes{1}).(BLorSD).(LPorDP).(SleepStage);   % if Female data present, we did a 2-way ANOVA first, then only post-hocs if p-vals were small
	end  

	p_val.(Sexes{2}).WTvsMut = p_vals_struct.(Sexes{2}).(BLorSD).(LPorDP).(SleepStage).WTvsMut; 
	p_val.WT.([Sexes{1} 'Vs' Sexes{2}])  = p_vals_struct.(['Both' Sexes{1} 'and' Sexes{2}]).(BLorSD).(LPorDP).(SleepStage).Posthocs.(['WT_'  Sexes{1} 'Vs' Sexes{2}]);
	p_val.Mut.([Sexes{1} 'Vs' Sexes{2}]) = p_vals_struct.(['Both' Sexes{1} 'and' Sexes{2}]).(BLorSD).(LPorDP).(SleepStage).Posthocs.(['Mut_' Sexes{1} 'Vs' Sexes{2}]); 

else 	% for the case where Percentages is output var
 	if NoFemaleData
		%use p vals from repeated measures ANOVAs (need adjusted and un because of plotting code later. Same if no female data present)
		p_val.(Sexes{1}).WTvsMut.Adjusted   = p_vals_struct.(Sexes{1}).(BLorSD).(LPorDP).(SleepStage).Posthocs.WTvsMut;
		p_val.(Sexes{1}).WTvsMut.Unadjusted = p_vals_struct.(Sexes{1}).(BLorSD).(LPorDP).(SleepStage).Posthocs.WTvsMut;
		p_val.(Sexes{2}).WTvsMut   = [];
		p_val.WT.([Sexes{1} 'Vs' Sexes{2}])  = [];
		p_val.Mut.([Sexes{1} 'Vs' Sexes{2}]) = [];

	else   % Percentages and both sexes present
		 p_val.WT.([Sexes{1}  'Vs' Sexes{2}]) = p_vals_struct.(BLorSD).(LPorDP).(SleepStage).Posthocs.(['WT_'  Sexes{1} 'vs' Sexes{2}]);
		 p_val.Mut.([Sexes{1} 'Vs' Sexes{2}]) = p_vals_struct.(BLorSD).(LPorDP).(SleepStage).Posthocs.(['Mut_' Sexes{1} 'vs' Sexes{2}]);
		 p_val.(Sexes{2}).WTvsMut = p_vals_struct.(BLorSD).(LPorDP).(SleepStage).Posthocs.([Sexes{2} '_WTvsMut']);
		 p_val.(Sexes{1}).WTvsMut = p_vals_struct.(BLorSD).(LPorDP).(SleepStage).Posthocs.([Sexes{1} '_WTvsMut']);
	end
end


% --- Counts per group -----
N_male_WT   = length(Twelve_hour_avg_values_Male.WT.(BLorSD).(First12hrsOrLast12hrs).(SleepStage));
N_female_WT = length(Twelve_hour_avg_values_Female.WT.(BLorSD).(First12hrsOrLast12hrs).(SleepStage)); 

N_male_Mut   = length(Twelve_hour_avg_values_Male.Mut.(BLorSD).(First12hrsOrLast12hrs).(SleepStage));
N_female_Mut = length(Twelve_hour_avg_values_Female.Mut.(BLorSD).(First12hrsOrLast12hrs).(SleepStage)); 

% --- Means ------------------------------
Mean_M_WT  = mean(Twelve_hour_avg_values_Male.WT.(BLorSD).(First12hrsOrLast12hrs).(SleepStage),'omitmissing');
Mean_M_Mut = mean(Twelve_hour_avg_values_Male.Mut.(BLorSD).(First12hrsOrLast12hrs).(SleepStage),'omitmissing');
Mean_F_WT  = mean(Twelve_hour_avg_values_Female.WT.(BLorSD).(First12hrsOrLast12hrs).(SleepStage),'omitmissing');
Mean_F_Mut = mean(Twelve_hour_avg_values_Female.Mut.(BLorSD).(First12hrsOrLast12hrs).(SleepStage),'omitmissing');


% -----------------------------------------------------------
%figure

axis_linewidth = 2; 
myfontsize = 16;
mymarkersize   = 120;
myjitterwidth  = 0.2;

f=figure;
f.Renderer = 'painters';

% if NoFemaleData
% 	f.Position = [300 100 800 700];
% end 



% p=plot(1*ones(1,N_male_WT),Twelve_hour_avg_values_Male.WT.(BLorSD).(First12hrsOrLast12hrs).(SleepStage),  'ko',2*ones(1,N_male_Mut),Twelve_hour_avg_values_Male.Mut.(BLorSD).(First12hrsOrLast12hrs).(SleepStage),  'ro', ...
% 	 3*ones(1,N_female_WT),Twelve_hour_avg_values_Female.WT.(BLorSD).(First12hrsOrLast12hrs).(SleepStage),'ko',4*ones(1,N_female_Mut),Twelve_hour_avg_values_Female.Mut.(BLorSD).(First12hrsOrLast12hrs).(SleepStage),'ro');
% jitter plot code
sw1 = swarmchart(1*ones(1,N_male_WT),Twelve_hour_avg_values_Male.WT.(BLorSD).(First12hrsOrLast12hrs).(SleepStage),mymarkersize);
hold on 
sw2 = swarmchart(2*ones(1,N_male_Mut),Twelve_hour_avg_values_Male.Mut.(BLorSD).(First12hrsOrLast12hrs).(SleepStage),mymarkersize);
sw3 = swarmchart(3*ones(1,N_female_WT),Twelve_hour_avg_values_Female.WT.(BLorSD).(First12hrsOrLast12hrs).(SleepStage),0.6*mymarkersize);
sw4 = swarmchart(4*ones(1,N_female_Mut),Twelve_hour_avg_values_Female.Mut.(BLorSD).(First12hrsOrLast12hrs).(SleepStage),0.6*mymarkersize);
sw1.MarkerFaceColor = 'k';  % black for WT
sw2.MarkerFaceColor = 'r';  % red for Mut
sw3.MarkerFaceColor = 'none';  
sw4.MarkerFaceColor = 'none';  
sw1.XJitterWidth = myjitterwidth;
sw2.XJitterWidth = myjitterwidth;
sw3.XJitterWidth = myjitterwidth;
sw4.XJitterWidth = myjitterwidth;
sw1.XJitter = 'none';    						% Turn off jitter
sw2.XJitter = 'none';    						% Turn off jitter
sw3.XJitter = 'none';    						% Turn off jitter
sw4.XJitter = 'none';    						% Turn off jitter

sw1.MarkerEdgeColor = 'none';  
sw2.MarkerEdgeColor = 'none';  
sw3.MarkerEdgeColor = 'k';  
sw4.MarkerEdgeColor = 'r';
sw3.LineWidth = 3;
sw4.LineWidth = 3;

ax=f.Children;
ax.XLim = [0 5];
%ax.YLim = [0 100];
all_data_in_plot = [Twelve_hour_avg_values_Male.WT.(BLorSD).(First12hrsOrLast12hrs).(SleepStage),  Twelve_hour_avg_values_Female.WT.(BLorSD).(First12hrsOrLast12hrs).(SleepStage), ...
										Twelve_hour_avg_values_Male.Mut.(BLorSD).(First12hrsOrLast12hrs).(SleepStage), Twelve_hour_avg_values_Female.Mut.(BLorSD).(First12hrsOrLast12hrs).(SleepStage)];
ax.YLim  = [0 1.2*max(all_data_in_plot)];
%ax.YTick = 0:25:100;
ax.XTick = [];
ax.FontSize = myfontsize;
ax.FontName = 'DejaVu Sans';
ax.Box = 'off';
ax.Color = 'none';
ax.LineWidth = axis_linewidth;
ax.TickDir = 'out';
%ax.YLabel.String = {strcat('{\bf ',SleepStage_for_ylabel,'}');['Time in ', SleepStage, ' (% TRT)']};
ax.YLabel.String = y_label_string;

ax.XTick = [1.5 3.5];
ax.XTickLabel = {Sexes{1},Sexes{2}};
ax.XAxis.TickLength = [0 0];
if strcmp(First12hrsOrLast12hrs,'Last12hrs')
	ax.Color = [0.9 0.9 0.9];  % Use a light gray background if plotting dark phase
end  
% for i=1:length(p)
% 	p(i).LineWidth = 3;
% 	p(i).MarkerSize = 12;
% end 
% p(1).MarkerFaceColor = 'k';
% p(2).MarkerFaceColor = 'r';
%hold on
line1 = line([0.85 1.15],[Mean_M_WT Mean_M_WT]);
line1.Color =[0.3 0.3 0.3];
line1.LineWidth = 6;
line2 = line([1.85 2.15],[Mean_M_Mut Mean_M_Mut]);
line2.Color =[0.3 0.3 0.3];
line2.LineWidth = 6;
line3 = line([2.85 3.15],[Mean_F_WT Mean_F_WT]);
line3.Color =[0.3 0.3 0.3];
line3.LineWidth = 6;
line4 = line([3.85 4.15],[Mean_F_Mut Mean_F_Mut]);
line4.Color =[0.3 0.3 0.3];
line4.LineWidth = 6;


% --- add * and # if significant interactions ---
vert_fudge_factor1 = 5;
vert_fudge_factor2 = 8;
asterisk_FontSize  = 28;
hashtag_FontSize   = 16;
sig_line_thickness = 1.5;
diff_geno_color    = [0.3 0.3 0.3];
diff_sex_color     = [0.3 0.3 0.3];  

% Male WTvsMut
if p_val.(Sexes{1}).WTvsMut.Adjusted < 0.05    & p_val.(Sexes{1}).WTvsMut.Adjusted >=0.01
	symbol.MaleWTVsMut = '*';
elseif p_val.(Sexes{1}).WTvsMut.Adjusted <0.01 & p_val.(Sexes{1}).WTvsMut.Adjusted >=0.001
	symbol.MaleWTVsMut = '**';
elseif p_val.(Sexes{1}).WTvsMut.Adjusted < 0.001
	symbol.MaleWTVsMut = '***';
end 

y_loc = 1.1*max(all_data_in_plot);

% add the p-value for males WT vs Mut (if you only have Male data)
if NoFemaleData
	if p_val.(Sexes{1}).WTvsMut.Adjusted < 0.01 p_val_displayed = '< 0.01'; else p_val_displayed = ['= ',num2str(p_val.(Sexes{1}).WTvsMut.Adjusted,'%.2f')]; end  
	text(ax,1.5,ax.YLim(2),['p value ',p_val_displayed],'FontSize',14,'HorizontalAlignment','center');
end

if ~NoFemaleData  % If female data are present
	
	% add the p-value for females WT vs Mut
	% if p_val.Female.WTvsMut < 0.01 p_val_displayed = '< 0.01'; else p_val_displayed = ['= ',num2str(p_val.Female.WTvsMut,'%.2f')]; end 
	% 	text(ax,3.5,ax.YLim(2),['p value ',p_val_displayed],'FontSize',14,'HorizontalAlignment','center');
	% end 

	% Set up the symbols

	% Female WTvsMut
	if p_val.(Sexes{2}).WTvsMut.Adjusted < 0.05    & p_val.(Sexes{2}).WTvsMut.Adjusted >=0.01
		symbol.FemaleWTVsMut = '*';
	elseif p_val.(Sexes{2}).WTvsMut.Adjusted <0.01 & p_val.(Sexes{2}).WTvsMut.Adjusted >=0.001
		symbol.FemaleWTVsMut = '**';
	elseif p_val.(Sexes{2}).WTvsMut.Adjusted < 0.001
		symbol.FemaleWTVsMut = '***';
	end 

	% WT Male vs Female
	if p_val.WT.([Sexes{1} 'Vs' Sexes{2}]).Adjusted < 0.05    & p_val.WT.([Sexes{1} 'Vs' Sexes{2}]).Adjusted >=0.01
		symbol.WTMaleVsFemale = '#';
	elseif p_val.WT.([Sexes{1} 'Vs' Sexes{2}]).Adjusted <0.01 & p_val.WT.([Sexes{1} 'Vs' Sexes{2}]).Adjusted >=0.001
		symbol.WTMaleVsFemale = '##';
	elseif p_val.WT.([Sexes{1} 'Vs' Sexes{2}]).Adjusted < 0.001
		symbol.WTMaleVsFemale = '###';
	end 

	% Mut Male vs Female
	if p_val.Mut.([Sexes{1} 'Vs' Sexes{2}]).Adjusted < 0.05    & p_val.Mut.([Sexes{1} 'Vs' Sexes{2}]).Adjusted >=0.01
		symbol.MutMaleVsFemale = '#';
	elseif p_val.Mut.([Sexes{1} 'Vs' Sexes{2}]).Adjusted <0.01 & p_val.Mut.([Sexes{1} 'Vs' Sexes{2}]).Adjusted >=0.001
		symbol.MutMaleVsFemale = '##';
	elseif p_val.Mut.([Sexes{1} 'Vs' Sexes{2}]).Adjusted < 0.001
		symbol.MutMaleVsFemale = '###';
	end 

	% add black asterisks for genotype differences within Sex
	% Male asterisks
	if p_val.(Sexes{1}).WTvsMut.Adjusted < 0.05
		y_loc = max([Twelve_hour_avg_values_Male.WT.(BLorSD).(First12hrsOrLast12hrs).(SleepStage) ...
		 			 Twelve_hour_avg_values_Male.Mut.(BLorSD).(First12hrsOrLast12hrs).(SleepStage)])+vert_fudge_factor1;
			
		y_loc = 1.1*max(all_data_in_plot);
		if y_loc>ax.YLim(2) y_loc=ax.YLim(2); end 
		text(ax,1.5,y_loc,symbol.MaleWTVsMut,'FontSize',asterisk_FontSize,'HorizontalAlignment','center','Color',diff_geno_color);
		plot(ax,[1,2],[y_loc,y_loc],'Color',diff_geno_color,'LineWidth',sig_line_thickness);
	end 

	% Female asterisks
	if p_val.(Sexes{2}).WTvsMut.Adjusted < 0.05
		y_loc = max([Twelve_hour_avg_values_Female.WT.(BLorSD).(First12hrsOrLast12hrs).(SleepStage) ...
		 			 Twelve_hour_avg_values_Female.Mut.(BLorSD).(First12hrsOrLast12hrs).(SleepStage)])+vert_fudge_factor1;
		y_loc = 1.1*max(all_data_in_plot);
		if y_loc>ax.YLim(2) y_loc=ax.YLim(2); end 
		text(ax,3.5,y_loc,symbol.FemaleWTVsMut,'FontSize',asterisk_FontSize,'HorizontalAlignment','center','Color',diff_geno_color);
		plot(ax,[3,4],[y_loc,y_loc],'Color',diff_geno_color,'LineWidth',sig_line_thickness);
	end 

	% WT # M vs F
	if p_val.WT.([Sexes{1} 'Vs' Sexes{2}]).Adjusted < 0.05
		y_loc = max([Twelve_hour_avg_values_Male.WT.(BLorSD).(First12hrsOrLast12hrs).(SleepStage) ...
		 			 Twelve_hour_avg_values_Female.WT.(BLorSD).(First12hrsOrLast12hrs).(SleepStage)])+vert_fudge_factor2;
		if y_loc>ax.YLim(2) y_loc=ax.YLim(2); end 
		text(ax,2,y_loc,symbol.WTMaleVsFemale,'FontSize',hashtag_FontSize,'Color','red','HorizontalAlignment','center','VerticalAlignment','bottom','Color',diff_sex_color);
		plot(ax,[1,3],[y_loc,y_loc],'Color',diff_sex_color,'LineWidth',sig_line_thickness);
	end

	% Mut # M vs F 
	if p_val.Mut.([Sexes{1} 'Vs' Sexes{2}]).Adjusted < 0.05
		y_loc = max([Twelve_hour_avg_values_Male.Mut.(BLorSD).(First12hrsOrLast12hrs).(SleepStage) ...
		 			 Twelve_hour_avg_values_Female.Mut.(BLorSD).(First12hrsOrLast12hrs).(SleepStage)])+vert_fudge_factor2;
		if y_loc>ax.YLim(2) y_loc=ax.YLim(2); end 
		text(ax,3,y_loc,symbol.MutMaleVsFemale,'FontSize',hashtag_FontSize,'Color','red','HorizontalAlignment','center','VerticalAlignment','bottom','Color',diff_sex_color);
		plot(ax,[2,4],[y_loc,y_loc],'Color',diff_sex_color,'LineWidth',sig_line_thickness);
	end
end % end of if female data are present 


hold off
if strcmp(SleepStage,'Wake')  % Only add the legend if Wake (since it's the top graph)
	%l=legend([p(1) p(2)],LegendLabels{1},LegendLabels{2});
	[l,objh] = legend([sw1 sw2],LegendLabels{1},LegendLabels{2},'AutoUpdate','off');
	objhl    = findobj(objh, 'type', 'patch'); % objects of legend of type patch
	set(objhl, 'Markersize', 11);  % set the size of circles in legend.  you kind of need to guess
	l.Box = 'off';
	l.Location = 'northwest';
	l.FontSize = myfontsize;
else 
	dummy_legend = legend('none');  		% This is to keep combine_time_in_state_panels_into_one_big_figure.m to not complain
	set(dummy_legend,'visible','off'); 
end 
%ax.Title.String = strcat(SleepStage,{' '},BLorSD,{' '},First12hrsOrLast12hrs); % Add an informative title (feel free to delete this line)


% Handle the case where there is no female data change the x-axis limits so we don't see the "Female" empty spot
if NoFemaleData
	ax.XLim = [0 3];
	ax.XAxis.TickValues = [1 2];  % change the label on the x-axis from "Males" to "WT" and "Mut"
	ax.XAxis.TickLabels = LegendLabels; %{'WT';'Mut'}; 
end  



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