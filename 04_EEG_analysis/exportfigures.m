FolderName = '/Users/caitlinottaway/Library/CloudStorage/OneDrive-WashingtonStateUniversity(email.wsu.edu)/Pexioto/EEG_Analysis/P24_P30_P45_P60/Figures/P45_figures for only spectra with animals used according to lizzy spreadsheet/' %destination path
FigList = findobj(allchild(0), 'flat', 'Type', 'figure');

for iFig = 1:length(FigList)
    FigHandle = FigList(iFig);
    ax = findobj(FigHandle, 'type', 'axes');
    
    % Initialize figure name
    FigName = sprintf('Figure_%d', iFig);  % Default fallback name
    
    savefig(FigHandle, fullfile(FolderName, [FigName '.fig']));
end
