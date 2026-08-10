function [latency_struct,latency_fig_handles,latency_p_values,latency_sig_difference_bool,ANOVA_table_latency] = sleep_latency_analysis(Scores_table,sleep_dep_length_hours,epoch_duration_secs,firstNREM_episode_duration_epochs,NREM_char,LegendLabels,Sexes,SexVarName)
%
% USAGE: [latency_fig_handles,latency_p_values,latency_sig_difference_bool] = sleep_latency_analysis(Scores_struct,sleep_dep_length_hours,epoch_duration_secs,firstNREM_episode_duration_epochs)
%
%
% This function reads in the state scores and computes the latency to the first NREM episode 
% that is at least firstNREM_episode_duration_epochs long. 
% 
% INPUTS:  	Scores_struct:  					a table containing the following columns: Scores, Sex, WTorMut, BLorSD 
% 			sleep_dep_length_hours:   			assumed to start at the beginning of the data passed into this function.  a scalar
% 			epoch_duration_secs: 				the epoch duration in seconds (used to score the data)
% 			firstNREM_episode_duration_epochs	the number of contiguous NREM epochs that are counted as the first NREM episode.
%
% OUTPUTS:  latency_fig_handles: 			a cell array of handles to the figures
%			latency_p_values:				a vector (or struct)? of the p-values comparing latency between WT and Mut for instance 


ind_first_epoch_after_sleep_dep = sleep_dep_length_hours*60*60/epoch_duration_secs + 1;  



% Pull apart the WT and Mutant data (SD only)
ind.WTSD  = (Scores_table.WTorMut=='WT'  & Scores_table.BLorSD=='SD');
ind.MutSD = (Scores_table.WTorMut=='Mut' & Scores_table.BLorSD=='SD');

ind.Male.WTSD    = (Scores_table.Sex==Sexes{1} & Scores_table.WTorMut=='WT'  & Scores_table.BLorSD=='SD');
ind.Male.MutSD   = (Scores_table.Sex==Sexes{1} & Scores_table.WTorMut=='Mut'  & Scores_table.BLorSD=='SD');

ind.Female.WTSD  = (Scores_table.Sex==Sexes{2} & Scores_table.WTorMut=='WT'  & Scores_table.BLorSD=='SD');
ind.Female.MutSD = (Scores_table.Sex==Sexes{2} & Scores_table.WTorMut=='Mut' & Scores_table.BLorSD=='SD');


Scores_WTSD = Scores_table.Scores(ind.WTSD);

Scores_Male_WT  = Scores_table.Scores(ind.Male.WTSD);
Scores_Male_Mut = Scores_table.Scores(ind.Male.MutSD);

Scores_Female_WT  = Scores_table.Scores(ind.Female.WTSD);
Scores_Female_Mut = Scores_table.Scores(ind.Female.MutSD);


% Do you have female data?
if isempty(Scores_Female_WT) & isempty(Scores_Female_Mut)
	FemaleDataPresent = 0;
else 
	FemaleDataPresent = 1;
end  




% -- WT --------------------
% Male
for i=1:length(Scores_Male_WT)  % loop over animals
	% pull out the scores 
	% calculate the time to the first NREM episode that is at least firstNREM_episode_duration_epochs long 

	scores_after_SD = Scores_Male_WT{i}(ind_first_epoch_after_sleep_dep:end);
	N_locs = ismember(scores_after_SD,NREM_char); % boolean for locations that are any of NREM_char
	Nruns_start_stop = contiguous(N_locs,1);  %  
	Nruns_start_stop = Nruns_start_stop{1,2}; % this is a N by 2 matrix where N is the number of runs found
	Nruns_length = Nruns_start_stop(:,2)-Nruns_start_stop(:,1) + 1; % length of this is the number of runs found
	ind_first_Nrun = find(Nruns_length>=firstNREM_episode_duration_epochs);
	latency.Male.WT(i) = Nruns_start_stop(ind_first_Nrun(1),1);
	latency.Male.WT(i) = latency.Male.WT(i)*epoch_duration_secs/60;  % Convert to minutes

end 

% Female
if FemaleDataPresent
	for i=1:length(Scores_Female_WT)
		
		scores_after_SD = Scores_Female_WT{i}(ind_first_epoch_after_sleep_dep:end);
		N_locs = ismember(scores_after_SD,NREM_char); % boolean for locations that are any of NREM_char
		Nruns_start_stop = contiguous(N_locs,1);  %  
		Nruns_start_stop = Nruns_start_stop{1,2}; % this is a N by 2 matrix where N is the number of runs found
		Nruns_length = Nruns_start_stop(:,2)-Nruns_start_stop(:,1) + 1; % length of this is the number of runs found
		ind_first_Nrun = find(Nruns_length>=firstNREM_episode_duration_epochs);
		latency.Female.WT(i) = Nruns_start_stop(ind_first_Nrun(1),1);
		latency.Female.WT(i) = latency.Female.WT(i)*epoch_duration_secs/60;  % Convert to minutes
	end 
else 
	latency.Female.WT = [];
end  
% -------------------------

% -- Mutant -------------
% Male
for i=1:length(Scores_Male_Mut)  % loop over animals
	% pull out the scores 
	% calculate the time to the first NREM episode that is at least firstNREM_episode_duration_epochs long 

	scores_after_SD = Scores_Male_Mut{i}(ind_first_epoch_after_sleep_dep:end);
	N_locs = ismember(scores_after_SD,NREM_char); % boolean for locations that are any of NREM_char
	Nruns_start_stop = contiguous(N_locs,1);  %  
	Nruns_start_stop = Nruns_start_stop{1,2}; % this is a N by 2 matrix where N is the number of runs found
	Nruns_length = Nruns_start_stop(:,2)-Nruns_start_stop(:,1) + 1; % length of this is the number of runs found
	ind_first_Nrun = find(Nruns_length>=firstNREM_episode_duration_epochs);
	latency.Male.Mut(i) = Nruns_start_stop(ind_first_Nrun(1),1);
	latency.Male.Mut(i) = latency.Male.Mut(i)*epoch_duration_secs/60;  % Convert to minutes

end 

% Female
if FemaleDataPresent
	for i=1:length(Scores_Female_Mut)
		scores_after_SD = Scores_Female_Mut{i}(ind_first_epoch_after_sleep_dep:end);
		N_locs = ismember(scores_after_SD,NREM_char); % boolean for locations that are any of NREM_char
		Nruns_start_stop = contiguous(N_locs,1);  %  
		Nruns_start_stop = Nruns_start_stop{1,2}; % this is a N by 2 matrix where N is the number of runs found
		Nruns_length = Nruns_start_stop(:,2)-Nruns_start_stop(:,1) + 1; % length of this is the number of runs found
		ind_first_Nrun = find(Nruns_length>=firstNREM_episode_duration_epochs);
		latency.Female.Mut(i) = Nruns_start_stop(ind_first_Nrun(1),1);
		latency.Female.Mut(i) = latency.Female.Mut(i)*epoch_duration_secs/60;  % Convert to minutes
	end 
else 
	latency.Female.Mut = [];
end 


% compute averages and sems
N_male_WT    = length(latency.Male.WT);
N_male_Mut   = length(latency.Male.Mut);
N_female_WT  = length(latency.Female.WT);
N_female_Mut = length(latency.Female.Mut);

Mean_M_WT  = mean(latency.Male.WT);
Mean_M_Mut = mean(latency.Male.Mut);
Mean_F_WT  = mean(latency.Female.WT);
Mean_F_Mut = mean(latency.Female.Mut); 


% do stats and compute p-values (so we can add annotations to the figure if needed)
[p_vals_latency,sig_difference_latency,ANOVA_table_latency] = perform_2way_Anova_SexGenotype_posthoc(latency.Female.WT,latency.Female.Mut,latency.Male.WT,latency.Male.Mut,Sexes,SexVarName);



latency_p_values=p_vals_latency;
latency_sig_difference_bool = sig_difference_latency;

% Since code below looks for adjusted p-values (from case with both sexes), if you have only males
% call the p-values adjusted and unadjusted
if ~FemaleDataPresent  % If no female data
	% Male WTvs Mut
	WTvsMut = latency_p_values.Posthoc.Male.WTvsMut;
	rmfield(latency_p_values.Posthoc.Male,'WTvsMut');  	% remove the field
	latency_p_values.Posthoc.Male.WTvsMut = struct;		% create the field
	latency_p_values.Posthoc.Male.WTvsMut.Adjusted   = WTvsMut;
	latency_p_values.Posthoc.Male.WTvsMut.Unadjusted = WTvsMut;

	% Female, WT vs Mut (set to empty)
	rmfield(latency_p_values.Posthoc.Female,'WTvsMut');	% remove the field
	latency_p_values.Posthoc.Female.WTvsMut = struct('Adjusted',[],'Unadjusted',[]); 	% create the field
	

	%WT, Male Vs Female (set to empty)
	rmfield(latency_p_values.Posthoc.WT,'MalevsFemale'); 	% remove field
	latency_p_values.Posthoc.WT.MalevsFemale = struct('Adjusted',[],'Unadjusted',[]);		% create field
	

	%Mut, Male Vs Female (set to empty)
	rmfield(latency_p_values.Posthoc.Mut,'MalevsFemale'); 	% remove field
	latency_p_values.Posthoc.Mut.MalevsFemale = struct('Adjusted',[],'Unadjusted',[]);		% create field
	

end



% make plot
figure
axis_linewidth = 2;
max_y = max([latency.Male.WT latency.Male.Mut latency.Female.WT latency.Female.Mut]);

f=gcf;
f.Name = 'Sleep Latency Analysis';
f.Renderer = 'painters';
p=plot(1*ones(1,N_male_WT),latency.Male.WT,'ko',2*ones(1,N_male_Mut),latency.Male.Mut,'ro', ...
	   4*ones(1,N_female_WT),latency.Female.WT,'ko',5*ones(1,N_female_Mut),latency.Female.Mut,'ro');
ax=f.Children;
ax.XLim = [0 6];
%ax.YLim = [0 90];
ax.YLim = [0 1.2*max_y];
%ax.YTick = 0:30:90;
ax.XTick = [];
ax.FontSize = 14;
ax.Box = 'off';
ax.Color = 'none';
ax.LineWidth = axis_linewidth;
ax.TickDir = 'out';
ax.YLabel.String = 'Latency to NREMS (min)';
ax.XTick = [1.5 4.5];
ax.XTickLabel = {Sexes{1},Sexes{2}};
ax.XAxis.TickLength = [0 0];

for i=1:length(p)
	p(i).LineWidth = 3;
	p(i).MarkerSize = 12;
end 
p(1).MarkerFaceColor = 'k';
p(2).MarkerFaceColor = 'r';
hold on
line1 = line([0.85 1.15],[Mean_M_WT Mean_M_WT]);
line1.Color =[0.4 0.4 0.4];
line1.LineWidth = 6;
line2 = line([1.85 2.15],[Mean_M_Mut Mean_M_Mut]);
line2.Color =[0.4 0.4 0.4];
line2.LineWidth = 6;
line3 = line([3.85 4.15],[Mean_F_WT Mean_F_WT]);
line3.Color =[0.4 0.4 0.4];
line3.LineWidth = 6;
line4 = line([4.85 5.15],[Mean_F_Mut Mean_F_Mut]);
line4.Color =[0.4 0.4 0.4];
line4.LineWidth = 6;
%hold off


ax.Title.String = 'Sleep Latency'; % Add an informative title (feel free to delete this line)

% add * and # if significant interactions
vert_fudge_factor = 5;
asterisk_FontSize = 24;
sig_line_thickness = 1.5;
% Male WTvsMut
if latency_p_values.Posthoc.(Sexes{1}).WTvsMut.Adjusted < 0.05 & latency_p_values.Posthoc.(Sexes{1}).WTvsMut.Adjusted >=0.01
	symbol.MaleWTVsMut = '*';
elseif latency_p_values.Posthoc.(Sexes{1}).WTvsMut.Adjusted <0.01 & latency_p_values.Posthoc.(Sexes{1}).WTvsMut.Adjusted >=0.001
	symbol.MaleWTVsMut = '**';
elseif latency_p_values.Posthoc.(Sexes{1}).WTvsMut.Adjusted < 0.001
	symbol.MaleWTVsMut = '***';
end 

if latency_p_values.Posthoc.(Sexes{1}).WTvsMut.Adjusted < 0.05
	y_loc = max([latency.Male.WT latency.(Sexes{1}).Mut])+vert_fudge_factor;
	text(ax,1.5,y_loc,symbol.MaleWTVsMut,'FontSize',asterisk_FontSize,'HorizontalAlignment','center');
	plot(ax,[1,2],[y_loc-2,y_loc-2],'k','LineWidth',sig_line_thickness);
end 

% Female WTvsMut
if latency_p_values.Posthoc.(Sexes{2}).WTvsMut.Adjusted < 0.05 & latency_p_values.Posthoc.(Sexes{2}).WTvsMut.Adjusted >=0.01
	symbol.FemaleWTVsMut = '*';
elseif latency_p_values.Posthoc.(Sexes{2}).WTvsMut.Adjusted <0.01 & latency_p_values.Posthoc.(Sexes{2}).WTvsMut.Adjusted >=0.001
	symbol.FemaleWTVsMut = '**';
elseif latency_p_values.Posthoc.(Sexes{2}).WTvsMut.Adjusted < 0.001
	symbol.FemaleWTVsMut = '***';
end 
if latency_p_values.Posthoc.(Sexes{2}).WTvsMut.Adjusted < 0.05
	y_loc = max([latency.Female.WT latency.Female.Mut])+vert_fudge_factor;
	text(ax,4.5,y_loc,symbol.FemaleWTVsMut,'FontSize',asterisk_FontSize,'HorizontalAlignment','center');
	plot(ax,[4,5],[y_loc-2,y_loc-2],'k','LineWidth',sig_line_thickness);
end 

% WT Male Vs Female
if latency_p_values.Posthoc.WT.([Sexes{1} 'vs' Sexes{2}]).Adjusted < 0.05 & latency_p_values.Posthoc.WT.([Sexes{1} 'vs' Sexes{2}]).Adjusted >=0.01
	symbol.WTMaleVsFemale = '#';
elseif latency_p_values.Posthoc.WT.([Sexes{1} 'vs' Sexes{2}]).Adjusted <0.01 & latency_p_values.Posthoc.WT.([Sexes{1} 'vs' Sexes{2}]).Adjusted >=0.001
	symbol.WTMaleVsFemale = '##';
elseif latency_p_values.Posthoc.WT.([Sexes{1} 'vs' Sexes{2}]).Adjusted < 0.001
	symbol.WTMaleVsFemale = '###';
end 
if latency_p_values.Posthoc.WT.([Sexes{1} 'vs' Sexes{2}]).Adjusted < 0.05
	y_loc = max([latency.Female.WT latency.Male.WT])+vert_fudge_factor;
	text(ax,2.5,y_loc,symbol.WTMaleVsFemale,'FontSize',asterisk_FontSize,'Color','red','HorizontalAlignment','center');
	plot(ax,[1,4],[y_loc-2,y_loc-2],'r','LineWidth',sig_line_thickness);
end 


% Mut Male Vs Female
if latency_p_values.Posthoc.Mut.([Sexes{1} 'vs' Sexes{2}]).Adjusted < 0.05 & latency_p_values.Posthoc.Mut.([Sexes{1} 'vs' Sexes{2}]).Adjusted >=0.01
	symbol.MutMaleVsFemale = '#';
elseif latency_p_values.Posthoc.Mut.([Sexes{1} 'vs' Sexes{2}]).Adjusted <0.01 & latency_p_values.Posthoc.Mut.([Sexes{1} 'vs' Sexes{2}]).Adjusted >=0.001
	symbol.MutMaleVsFemale = '##';
elseif latency_p_values.Posthoc.Mut.([Sexes{1} 'vs' Sexes{2}]).Adjusted < 0.001
	symbol.MutMaleVsFemale = '###';
end 
if latency_p_values.Posthoc.Mut.([Sexes{1} 'vs' Sexes{2}]).Adjusted < 0.05
	y_loc = max([latency.Female.Mut latency.Male.Mut])+vert_fudge_factor;
	text(ax,3.5,y_loc,symbol.MutMaleVsFemale,'FontSize',asterisk_FontSize,'Color','red','HorizontalAlignment','center');
	plot(ax,[2,5],[y_loc-2,y_loc-2],'r','LineWidth',sig_line_thickness);
end 

l=legend([p(1) p(2)],LegendLabels{1},LegendLabels{2});
l.Box = 'off';
l.Location = 'northwest';
l.FontSize = 14;



% Handle the case where there is no female data change the x-axis limits so we don't see the "Female" empty spot
if FemaleDataPresent==0
	ax.XLim = [0 3];
	ax.XAxis.TickValues = [1 2];  % change the label on the x-axis from "Males" to "WT" and "Mut"
	ax.XAxis.TickLabels = {'WT';'Mut'}; 
end  


hold off 
% -----------------------------------------------------------------------------
% --- End of code for calculating and plotting latency to NREMS in minutes ----









latency_fig_handles=f;






latency_struct = latency;  
