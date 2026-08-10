function combine_time_in_state_panels_into_one_big_figure(metric)
% This little script makes two figures out of 12 existing figures. 
% specifically, it makes a figure with 6 panels (W light phase, N light phase, R light phase, then same for dark phase)
% by combining the 6 figures (of each individual panel) just made into one large figure
% then it does the same thing for the the SD. 
% resulting in two new figures  





% Combine these 6 BL figures into one large figure with 6 panels
figlist   = get(groot,'Children');
figlistBL = figlist(7:12);  % Only take the last 6 (they are ordered backwards!)
figlistBL = flipud(figlistBL);  % flipped to get the ordering correct
figlistBL = [figlistBL(1) figlistBL(4) figlistBL(2) figlistBL(5) figlistBL(3) figlistBL(6)];  % the ordering got weird, this fixes it


% Combine these 6 SD figures into one large figure with 6 panels
figlistSD = figlist(1:6);  % Only take the last 6 (they are ordered backwards!)
figlistSD = flipud(figlistSD);  % flipped to get the ordering correct
figlistSD = [figlistSD(1) figlistSD(4) figlistSD(2) figlistSD(5) figlistSD(3) figlistSD(6)];  % the ordering got weird, this fixes it

newfigBL = figure;
newfigBL.Name = strcat(metric,' BL Six Panels');
newfigBL.Renderer = 'painters';
%newfigBL.Position = [350 267 1035 549];
newfigBL.Position = [300 100 800 700];
% tcl    = tiledlayout(newfig,'flow');
tclBL    = tiledlayout(newfigBL,3,2);

newfigSD = figure;
newfigSD.Name = strcat(metric,' SD Six Panels');
newfigSD.Renderer = 'painters';
%newfigSD.Position = [340 260 1035 549];
newfigSD.Position = [300 100 800 700];
% tcl    = tiledlayout(newfig,'flow');
tclSD    = tiledlayout(newfigSD,3,2);


% Determine vertical axis limits to make them the same between BL and SD (and between first12 and last12) by state
%Wake
wake_max = max([figlistBL(1).CurrentAxes.YLim,figlistBL(2).CurrentAxes.YLim,figlistSD(1).CurrentAxes.YLim,figlistSD(2).CurrentAxes.YLim]);
nrem_max = max([figlistBL(3).CurrentAxes.YLim,figlistBL(4).CurrentAxes.YLim,figlistSD(3).CurrentAxes.YLim,figlistSD(4).CurrentAxes.YLim]);
rem_max  = max([figlistBL(5).CurrentAxes.YLim,figlistBL(6).CurrentAxes.YLim,figlistSD(5).CurrentAxes.YLim,figlistSD(6).CurrentAxes.YLim]);



% -- BL --
for i = 1:numel(figlistBL)
    disp(['i = ', num2str(i)])
    figure(figlistBL(i));
    ax=gca;
    ax.Parent=tclBL;
    ax.Layout.Tile=i;
    if i==1
        ax.Title.String = 'Light Period';
        ax.Title.FontWeight = 'normal';
        ax.Title.FontName   = 'DejaVu Sans';   % python matplotlib like fonts
        ax.Title.FontSize = 13;
        ax.YLim = [0 wake_max];
    end 
    if i==2
        ax.Title.String = 'Dark Period';
        ax.Title.FontWeight = 'normal';
        ax.Title.FontName   = 'DejaVu Sans';   % python matplotlib like fonts
        ax.Title.FontSize = 13;
        ax.YLim = [0 wake_max];
    end
    % determine state and set ax.YLim
    if i==2
        ax.YLim = [0 wake_max];
    end 
    if i==3 | i==4
        ax.YLim = [0 nrem_max];
    end
    if i==5 | i==6
        ax.YLim = [0 rem_max];
    end 
end 


% -- SD ---
for i = 1:numel(figlistSD)
    figure(figlistSD(i));
    ax=gca;
    ax.Parent=tclSD;
    ax.Layout.Tile=i;
    if i==1
        ax.Title.String = 'Light Period';
        ax.Title.FontWeight = 'normal';
        ax.Title.FontName   = 'DejaVu Sans';   % python matplotlib like fonts
        ax.Title.FontSize   = 13;
        ax.YLim = [ 0 wake_max];
    end 
    if i==2
        ax.Title.String = 'Dark Period';
        ax.Title.FontWeight = 'normal';
        ax.Title.FontName   = 'DejaVu Sans';   % python matplotlib like fonts
        ax.Title.FontSize   = 13;
        ax.YLim = [0 wake_max];
    end
    % determine state and set ax.YLim
    if i==2
        ax.YLim = [0 wake_max];
    end 
    if i==3 | i==4
        ax.YLim = [0 nrem_max];
    end
    if i==5 | i==6
        ax.YLim = [0 rem_max];
    end 
end


% close all of the empty figures that were made empty when they were combined into larger figures
close(figlistBL)
close(figlistSD)


