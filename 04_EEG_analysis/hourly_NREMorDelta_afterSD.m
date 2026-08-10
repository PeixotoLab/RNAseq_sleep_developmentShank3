function out = hourly_NREM_Delta_or_Sigma_afterSD(options)
%
% This function reads in Hourly values of either Delta or sigma during either NREM or 



arguments 
	options.Hourly_power_inSTATE
	options.DeltaOrTheta
	options.WorNREM
end 

Hourly_power_inSTATE = options.Hourly_power_inSTATE;
DeltaOrTheta         = options.DeltaOrTheta;
WorNREM              = options.WorNREM; 



power_inSTATE_WT_Male_normalization    = mean(Hourly_power_inSTATE.WT.Male.BL ,2,'omitnan');
power_inSTATE_Mut_Male_normalization   = mean(Hourly_power_inSTATE.Mut.Male.BL,2,'omitnan');
power_inSTATE_WT_Female_normalization  = mean(Hourly_power_inSTATE.WT.Female.BL,2,'omitnan');
power_inSTATE_Mut_Female_normalization = mean(Hourly_power_inSTATE.Mut.Female.BL,2,'omitnan');  
	
Normalized_power.WT.Male    = (Hourly_power_inSTATE.WT.Male.SD./power_inSTATE_WT_Male_normalization)*100;
Normalized_power.Mut.Male   = (Hourly_power_inSTATE.Mut.Male.SD./power_inSTATE_Mut_Male_normalization)*100;
Normalized_power.WT.Female  = (Hourly_power_inSTATE.WT.Female.SD./power_inSTATE_WT_Female_normalization)*100;
Normalized_power.Mut.Female = (Hourly_power_inSTATE.Mut.Female.SD./power_inSTATE_Mut_Female_normalization)*100;  

% Next compute the means and SEs to make the plots  (use 6:12) since it's only hours 6-12 of the SD day.  (end of light period)
Means_norm_power.WT.Male    = mean(Normalized_power.WT.Male,1,'omitnan');
Means_norm_power.Mut.Male   = mean(Normalized_power.Mut.Male,1,'omitnan');
Means_norm_power.WT.Female  = mean(Normalized_power.WT.Female,1,'omitnan');
Means_norm_power.Mut.Female = mean(Normalized_power.Mut.Female,1,'omitnan');

SEMs_norm_power.WT.Male     = std(Normalized_power.WT.Male,0,1,'omitnan')./sqrt(sum(~isnan(Normalized_power.WT.Male),1));
SEMs_norm_power.Mut.Male    = std(Normalized_power.Mut.Male,0,1,'omitnan')./sqrt(sum(~isnan(Normalized_power.Mut.Male),1)); 
SEMs_norm_power.WT.Female   = std(Normalized_power.WT.Female,0,1,'omitnan')./sqrt(sum(~isnan(Normalized_power.WT.Female),1));
SEMs_norm_power.Mut.Female  = std(Normalized_power.Mut.Female,0,1,'omitnan')./sqrt(sum(~isnan(Normalized_power.Mut.Female),1));


% Make the figure
line_thick = 2;
filled_marker_size = 18;
open_marker_size   = 6;
figure('Renderer', 'painters', 'Position', [600 400 650 320])
f=gcf;
f.Renderer = 'painters';
subplot(1,2,1)    				% Male panel 
pm=plot(6:12,Means_norm_power.WT.Male(6:12),'k.',6:12,Means_norm_power.Mut.Male(6:12),'r.');
pm(1).LineWidth=line_thick;
pm(2).LineWidth=line_thick;
pm(1).MarkerSize = filled_marker_size;
pm(2).MarkerSize = filled_marker_size;
hold on 
eb1=errorbar(6:12,Means_norm_power.WT.Male(6:12),SEMs_norm_power.WT.Male(6:12),'k');
eb2=errorbar(6:12,Means_norm_power.Mut.Male(6:12),SEMs_norm_power.Mut.Male(6:12),'r');
eb1.LineWidth = line_thick;
eb2.LineWidth = line_thick;
eb1.CapSize = 0;
eb2.CapSize = 0;
yline(100,'--','LineWidth',3);
r1 = rectangle('Position',[0 75 5 2]);
r1.FaceColor = 'red';
r1.LineWidth = 1;
r2 = rectangle('Position',[5 75 7 2]);
r2.LineWidth = 1;
hold off
ax=gca;
ax.XLim = [0 12];
ax.YLim = [75 175];
ax.Box= 'off';
ax.LineWidth=1.0;
ax.Title.String = 'Males';
ax.XTick = [ 0 6 12];
ax.YTick = [75:25:175];
ax.XLabel.String = 'Hour';
ax.YLabel.String = 'Normalized NREM Delta';
ax.YLabel.FontWeight = 'bold';
ax.YLabel.FontSize   = 14;
ax.FontSize = 16;

subplot(1,2,2) 				% Female panel
pm=plot(6:12,Means_norm_power.WT.Female(6:12),'ko',6:12,Means_norm_power.Mut.Female(6:12),'ro');
pm(1).LineWidth=line_thick;
pm(2).LineWidth=line_thick;
pm(1).MarkerSize = open_marker_size;
pm(2).MarkerSize = open_marker_size;

hold on 
eb1=errorbar(6:12,Means_norm_power.WT.Female(6:12),SEMs_norm_power.WT.Female(6:12),'o-','MarkerFaceColor','w','MarkerEdgeColor','k','Color','k');
eb2=errorbar(6:12,Means_norm_power.Mut.Female(6:12),SEMs_norm_power.Mut.Female(6:12),'o-','MarkerFaceColor','w','MarkerEdgeColor','r','Color','r');
eb1.LineWidth = line_thick;
eb2.LineWidth = line_thick;
eb1.CapSize = 0;
eb2.CapSize = 0;
yline(100,'--','LineWidth',3);
r1 = rectangle('Position',[0 75 5 2]);
r1.FaceColor = 'red';
r1.LineWidth = 1;
r2 = rectangle('Position',[5 75 7 2]);
r2.LineWidth = 1;
hold off
ax=gca;
ax.XLim = [0 12];
ax.YLim = [75 175];
ax.Box= 'off';
ax.LineWidth=1.0;
ax.Title.String = 'Females';
ax.XTick = [ 0 6 12];
ax.YTick = [75:25:175];
ax.XLabel.String = 'Hour';
ax.FontSize = 16;
%ax.YLabel.String = 'Normalized NREM Delta';
% ax.YLabel.FontWeight = 'bold';
% ax.YLabel.FontSize   = 14;

% --- Carry out the repeated measures ANOVA to compare the two curves -------------
% ---------------------------------------------------------------------------------
% --- Males, WT vs Mut ---
Mut_Male_6_12hr = Normalized_power.Mut.Male(:,6:12);
WT_Male_6_12hr  = Normalized_power.WT.Male(:,6:12);
%Genotype = categorical([0 0 0 0 0 1 1 1 1 1 1]');  % Categorical is important!  otherwise ranova produces incorrect F values
Genotype = categorical([zeros(1,size(Mut_Male_6_12hr,1)) ones(1,size(WT_Male_6_12hr,1))]');  % Categorical is important!  otherwise ranova produces incorrect F values

Y = [Mut_Male_6_12hr; WT_Male_6_12hr];
t = table(Genotype,Y(:,1),Y(:,2),Y(:,3),Y(:,4),Y(:,5),Y(:,6),Y(:,7), ...
                   'VariableNames',{'Genotype','Hour6','Hour7','Hour8','Hour9','Hour10','Hour11','Hour12'});

Time = 6:12;   % Within-subjects variable
WithinStructure = table(Time','VariableNames',{'Time'});
WithinStructure.Time = categorical(WithinStructure.Time);

rm = fitrm(t,'Hour6-Hour12 ~ Genotype','WithinDesign',WithinStructure,'WithinModel','Time');
RManovatbl_SleepPress_Male_WTVsMut = ranova(rm,'WithinModel','Time')

% follow up with post-hoc pairwise comparisons (and Sidak correction)
if RManovatbl_SleepPress_Male_WTVsMut.pValue(2) < 0.05
	disp('Do post-hoc comparisons here.')
end 


% --- Females, WT vs Mut ---
Mut_Female_6_12hr = Normalized_power.Mut.Female(:,6:12);
WT_Female_6_12hr  = Normalized_power.WT.Female(:,6:12);
%Genotype = categorical([0 0 0 0 0 1 1 1 1 1 1]');  % Categorical is important!  otherwise ranova produces incorrect F values
Genotype = categorical([zeros(1,size(Mut_Female_6_12hr,1)) ones(1,size(WT_Female_6_12hr,1))]');  % Categorical is important!  otherwise ranova produces incorrect F values

Y = [Mut_Female_6_12hr; WT_Female_6_12hr];
t = table(Genotype,Y(:,1),Y(:,2),Y(:,3),Y(:,4),Y(:,5),Y(:,6),Y(:,7), ...
                   'VariableNames',{'Genotype','Hour6','Hour7','Hour8','Hour9','Hour10','Hour11','Hour12'});

Time = 6:12;   % Within-subjects variable
WithinStructure = table(Time','VariableNames',{'Time'});
WithinStructure.Time = categorical(WithinStructure.Time);

rm = fitrm(t,'Hour6-Hour12 ~ Genotype','WithinDesign',WithinStructure,'WithinModel','Time');
RManovatbl_SleepPress_Feale_WTVsMut = ranova(rm,'WithinModel','Time')
if RM





















