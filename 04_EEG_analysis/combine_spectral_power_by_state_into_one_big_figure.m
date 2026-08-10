function combine_spectral_power_by_state_into_one_big_figure(FemaleDataPresent,Normalization,FigNames)
% This little script makes two figures out of 12 existing figures. 
% specifically, it makes a figure with 6 panels (W light phase, N light phase, R light phase, then same for dark phase)
% by combining the 6 figures (of each individual panel) just made into one large figure
% then it does the same thing for the the SD. 
% resulting in two new figures  
% This is specifically for spectral power figures.  Top row is Male.  Bottom row is Female.  Columns from left to right
% are Wake, NREM, REM 
% For each figure, we are combining light and dark phases.  





% Combine these 6 BL figures into one large figure with 6 panels
figlist   = get(groot,'Children');
if FemaleDataPresent
    figlistBL = figlist(7:12); % Only take the last 6 (they are ordered backwards!)
else 
    figlistBL = figlist(4:6);  % If no female data, take the last 3  
end 
figlistBL = flipud(figlistBL);  % flipped to get the ordering correct
%figlistBL = [figlistBL(1) figlistBL(4) figlistBL(2) figlistBL(5) figlistBL(3) figlistBL(6)];  % the ordering got weird, this fixes it


% Combine these 6 SD figures into one large figure with 6 panels
if FemaleDataPresent
    figlistSD = figlist(1:6);  % Only take the first 6 (they are ordered backwards!)
else 
    figlistSD = figlist(1:3);  % If no female data, take the first 3.
end 

figlistSD = flipud(figlistSD);  % flipped to get the ordering correct
%figlistSD = [figlistSD(1) figlistSD(4) figlistSD(2) figlistSD(5) figlistSD(3) figlistSD(6)];  % the ordering got weird, this fixes it

newfigBL = figure;
newfigBL.Name = FigNames{1};

figBL = gcf;
figBL.Position = [860 460 880 380];
newfigBL.Renderer = 'painters';
% tcl    = tiledlayout(newfig,'flow');
if FemaleDataPresent
    tclBL = tiledlayout(newfigBL,2,3);
else 
    tclBL = tiledlayout(newfigBL,1,3);
end 

newfigSD = figure;
newfigSD.Name = FigNames{2};
figSD = gcf;
figSD.Position = [860 460 880 380];
newfigSD.Renderer = 'painters';
% tcl    = tiledlayout(newfig,'flow');
if FemaleDataPresent
    tclSD = tiledlayout(newfigSD,2,3);
else 
    tclSD = tiledlayout(newfigSD,1,3); 
end 


% Determine vertical axis limits to make them the same between BL and SD (and between first12 and last12) by state
% Case with only Male data 
if ~FemaleDataPresent
    max_BL = max([figlistBL(1).CurrentAxes.YLim,figlistBL(2).CurrentAxes.YLim,figlistBL(3).CurrentAxes.YLim]);
    max_SD = max([figlistSD(1).CurrentAxes.YLim,figlistSD(2).CurrentAxes.YLim,figlistSD(3).CurrentAxes.YLim]);
else 
    max_BL = max([figlistBL(1).CurrentAxes.YLim,figlistBL(2).CurrentAxes.YLim,figlistBL(3).CurrentAxes.YLim,figlistBL(4).CurrentAxes.YLim,figlistBL(5).CurrentAxes.YLim,figlistBL(6).CurrentAxes.YLim]);
    max_SD = max([figlistSD(1).CurrentAxes.YLim,figlistSD(2).CurrentAxes.YLim,figlistSD(3).CurrentAxes.YLim,figlistSD(4).CurrentAxes.YLim,figlistSD(5).CurrentAxes.YLim,figlistSD(6).CurrentAxes.YLim]);   
end



for i = 1:numel(figlistBL)
    figure(figlistBL(i));
    L = findobj(figlistBL(i),'type','line');  % get the line so we can make it dashed if Female data
    ax=gca;
    ax.Parent=tclBL;
    ax.Layout.Tile=i;
    ax.YLim = [0 max_BL];
    if i==1 | i==2 | i==3               % remove the x axis label for top 3 plots
        ax.XLabel.String = [];
    end 
    if i==2 | i==3 | i==5 | i==6        % remove the y axis label for 4 rightmost plots 
        ax.YLabel.String = [];
    end 
    if i==4 | i==5 | i==6               % dashed lines for bottom 3 plots (female?)
        L(1).LineStyle = '--';
        L(2).LineStyle = '--';
    end
    if strcmp(Normalization,'MeanPowerAllStates')
        %ax.YTickLabel = strcat(cellstr(string(str2double(ax.YTickLabel))),'%'); % add % here and not earlier
    end 
    % if i==1
    %     ax.Title.String = 'Baseline';
    % end 
end 

% add title above the subplots
if contains(FigNames{1},'LP')
    title(tclBL,'Light Period','FontSize',18)
elseif contains(FigNames{1},'DP')
    title(tclBL,'Dark Period','FontSize',18)
end 
% ---- End of BL  --------------------------------------

% --- Start of SD --------------------------------------
for i = 1:numel(figlistSD)
    figure(figlistSD(i));
    L = findobj(figlistSD(i),'type','line');  % get the line so we can make it dashed if Female data
    ax=gca;
    ax.Parent=tclSD;
    ax.Layout.Tile=i;
    ax.YLim = [0 max_SD];
    if i==1 | i==2 | i==3
        ax.XLabel.String = [];
    end
    if i==2 | i==3 | i==5 | i==6       % remove the y axis label for 4 rightmost plots 
        ax.YLabel.String = [];
    end 
    if i==4 | i==5 | i==6               % dashed lines for bottom 3 plots (female?)
        L(1).LineStyle = '--';
        L(2).LineStyle = '--';
    end
    if strcmp(Normalization,'MeanPowerAllStates')
        %ax.YTickLabel = strcat(cellstr(string(str2double(ax.YTickLabel))),'%'); % Add % here and not earlier
    end
    % if i==1
    %     ax.Title.String = 'Sleep Dep';
    % end 
end

% add title above the subplots
if contains(FigNames{1},'LP')
    title(tclSD,'Light Period','FontSize',18)
elseif contains(FigNames{1},'DP')
    title(tclSD,'Dark Period','FontSize',18)
end 

% close all of the empty figures that were made empty when they were combined into larger figures
close(figlistBL)
close(figlistSD)


