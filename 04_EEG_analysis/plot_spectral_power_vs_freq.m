function fig_handle = plot_spectral_power_vs_freq(options)
%
% USAGE: fig_handle = plot_spectral_power_vs_freq(Means_struct,SEMS_struct,sex,BLorSD,state,LegendLabels,EEG_bin_edges,EEGLowerLimitHz,Normalization) 
% INPUTS:  



arguments
    options.Means_struct     = []
    options.SEMS_struct      = []
    options.sex              = []
    options.BLorSD           = []
    options.state            = []
    options.LegendLabels     = []
    options.EEG_bin_edges 	 = []
    options.EEGLowerLimitHz  = []
    options.Normalization    = []
end 

Means_struct    = options.Means_struct;
SEMS_struct     = options.SEMS_struct;
sex       	    = options.sex;
BLorSD          = options.BLorSD;
state 		    = options.state;
LegendLabels    = options.LegendLabels;
EEG_bin_edges   = options.EEG_bin_edges;
EEGLowerLimitHz = options.EEGLowerLimitHz;
Normalization   = options.Normalization;





yWT  = Means_struct.(sex).WT.(BLorSD).(state);
yMut = Means_struct.(sex).Mut.(BLorSD).(state);

SEMS_WT  = SEMS_struct.(sex).WT.(BLorSD).(state);
SEMS_Mut = SEMS_struct.(sex).Mut.(BLorSD).(state); 

% if strcmp(Normalization,'AreaUnderCurve')
% 	ylabel_string = {'Power'; 'Normalized to Area Under Curve'};
% elseif strcmp(Normalization,'MeanPowerAllStates')
% 	ylabel_string = {'Power Normalized to Mean'; 'Power in BL (%)'};
% end
ylabel_string ={'Normalized';'Power'};  % MJR 1.25:  simplify so it looks nicer



if strcmp(state,'Wake')
	title_string = 'WAKE';
elseif strcmp(state,'NREM')
	title_string = 'NREM';
elseif strcmp(state,'REM')
	title_string = 'REM';
else
	error('You entered an invalid state value in plot_spectral_power_vs_freq.m ')
end 

% if strcmp(BLorSD,'BL')
% 	title_string = strcat(title_string,' Baseline');
% elseif strcmp(BLorSD,'SD')
% 	title_string = strcat(title_string,' Sleep Dep');
% elseif strcmp(BLorSD,'BLFirst12hr')
% 	title_string = strcat(title_string,' Baseline First 12 hr');
% elseif strcmp(BLorSD,'BLLast12hr')
% 	title_string = strcat(title_string,' Baseline Last 12 hr');
% elseif strcmp(BLorSD,'SDFirst12hr')
% 	title_string = strcat(title_string,' Sleep Dep First 12 hr');
% elseif strcmp(BLorSD,'SDLast12hr')
% 	title_string = strcat(title_string,' Sleep Dep Last 12 hr');
% else
% 	error('You entered an invalid Baseline or Sleep Dep value in plot_spectral_power_vs_freq.m. Options are ''BL'',''SD'', ''BLFirst12hr'',''BLLast12hr'',''SDFirst12hr'', or ''SDLast12hr''. ')
% end 

% if strcmp(sex,'Male')
% 	title_string = strcat(title_string,' Male');
% 	mylinestyle = '-';
% elseif strcmp(sex,'Female')
% 	title_string = strcat(title_string,' Female');
% 	mylinestyle = '--';
% else
% 	error('You entered an invalid value for sex in plot_spectral_power_vs_freq.m. Options are ''Male'' or ''Female''. ')
% end 

title_string = [title_string,' ',sex];
mylinestyle  = '-';


H=figure;

mylinewidth_ax    = 2;
mylinewidth_curve = 1.5;
myfontsize  	  = 16;

H.Name = strcat('Spectral Power ', state);
[~, indexOfNearestLowerLim] = min(abs(EEG_bin_edges-EEGLowerLimitHz));
xvals = EEG_bin_edges(indexOfNearestLowerLim:end);

plot(xvals,yWT,'Color','k','LineWidth',mylinewidth_curve,'LineStyle',mylinestyle)  
hold on 
plot(xvals,yMut,'Color','r','LineWidth',mylinewidth_curve,'LineStyle',mylinestyle)
ts_W_WT  = tinv([0.025 0.975],length(yWT)-1);
ts_W_Mut = tinv([0.025 0.975],length(yMut)-1);
CI_WT    = ts_W_WT(2)*SEMS_WT;
CI_Mut   = ts_W_Mut(2)*SEMS_Mut;
patch([xvals fliplr(xvals)],[yWT-CI_WT   fliplr(yWT + CI_WT)],  'k','FaceAlpha',0.2,'EdgeColor','none')
patch([xvals fliplr(xvals)],[yMut-CI_Mut fliplr(yMut + CI_Mut)],'r','FaceAlpha',0.2,'EdgeColor','none')    
hold off
ax=gca;
ax.FontSize = myfontsize;
ax.FontName = 'DejaVu Sans';
ax.XLim = [0 20];
ax.XLabel.String = 'Frequency (Hz)';
%ax.XTick=0:5:max(xvals);
ax.XTick=0:5:20;
ax.YLabel.String = ylabel_string; %'Normalized Power';
%ax.YLabel.FontSize=14;
ax.Box = 'off';
ax.LineWidth = mylinewidth_ax;
ax.Title.String     = title_string;
ax.Title.FontWeight = 'normal';

% shade the background if dark period
if contains(BLorSD,'Last12hr')
	ax.Color = [0.9 0.9 0.9];  % Use a light gray background if plotting dark phase
end 
% Make the y axis ticks reflect the normalization (ie if normalized to mean power in BL, use percentages)
% if strcmp(Normalization,'MeanPowerAllStates')
% 	disp('ax.YTickLabel before modifying: ')
% 	disp(ax.YTickLabel)
% 	ax.YTickLabel = strcat(cellstr(string(str2double(ax.YTickLabel))),'%');
% end


if strcmp(state,'Wake') & strcmp(sex,'Male')
	legend(LegendLabels)			% Was legend('Wildtype','Shank3\Delta C')
	legend boxoff
end 




fig_handle = gcf;