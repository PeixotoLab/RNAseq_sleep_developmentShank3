function [artefact_percentages_Table, overall_av_perc,exit_code] = calc_artifact_percentages(Scores,Sexes)




% go through each field of Scores and calculate the percentage of epochs 
% that were scored as artifact for each recording.  


% Male WT BL
for i=1:length(Scores.(Sexes{1}).WT.BL)
	artefact_percentages.(Sexes{1}).WT.BL(i) = sum(contains(Scores.(Sexes{1}).WT.BL{i}, '*'))/length(Scores.(Sexes{1}).WT.BL{i});
end

% Male WT SD
for i=1:length(Scores.(Sexes{1}).WT.SD)
	artefact_percentages.(Sexes{1}).WT.SD(i) = sum(contains(Scores.(Sexes{1}).WT.SD{i}, '*'))/length(Scores.(Sexes{1}).WT.SD{i});
end

% Male Mut BL
for i=1:length(Scores.(Sexes{1}).Mut.BL)
	artefact_percentages.(Sexes{1}).Mut.BL(i) = sum(contains(Scores.(Sexes{1}).Mut.BL{i}, '*'))/length(Scores.(Sexes{1}).Mut.BL{i});
end

% Male Mut SD
for i=1:length(Scores.(Sexes{1}).Mut.SD)
	artefact_percentages.(Sexes{1}).Mut.SD(i) = sum(contains(Scores.(Sexes{1}).Mut.SD{i}, '*'))/length(Scores.(Sexes{1}).Mut.SD{i});
end

% if Female data
if ~isempty(Scores.(Sexes{2}).WT.BL)
	% Female WT BL
	for i=1:length(Scores.(Sexes{2}).WT.BL)
		artefact_percentages.(Sexes{2}).WT.BL(i) = sum(contains(Scores.(Sexes{2}).WT.BL{i}, '*'))/length(Scores.(Sexes{2}).WT.BL{i});
	end
	% Female WT SD
	for i=1:length(Scores.(Sexes{2}).WT.SD)
		artefact_percentages.(Sexes{2}).WT.SD(i) = sum(contains(Scores.(Sexes{2}).WT.SD{i}, '*'))/length(Scores.(Sexes{2}).WT.SD{i});
	end
	% Female Mut BL
	for i=1:length(Scores.(Sexes{2}).Mut.BL)
		artefact_percentages.(Sexes{2}).Mut.BL(i) = sum(contains(Scores.(Sexes{2}).Mut.BL{i}, '*'))/length(Scores.(Sexes{2}).Mut.BL{i});
	end
	% Female Mut SD
	for i=1:length(Scores.(Sexes{2}).Mut.SD)
		artefact_percentages.(Sexes{2}).Mut.SD(i) = sum(contains(Scores.(Sexes{2}).Mut.SD{i}, '*'))/length(Scores.(Sexes{2}).Mut.SD{i});
	end
else
	artefact_percentages.(Sexes{2}).WT.BL = [];
	artefact_percentages.(Sexes{2}).WT.SD = [];
	artefact_percentages.(Sexes{2}).Mut.BL= [];
	artefact_percentages.(Sexes{2}).Mut.SD= [];
end



% Calculate the overall average of the percentage of epochs that were scored as artifact

aggregated = [artefact_percentages.(Sexes{1}).WT.BL   artefact_percentages.(Sexes{1}).WT.SD ...
			  artefact_percentages.(Sexes{1}).Mut.BL  artefact_percentages.(Sexes{1}).Mut.SD...
			  artefact_percentages.(Sexes{2}).WT.BL artefact_percentages.(Sexes{2}).WT.SD ...
			  artefact_percentages.(Sexes{2}).Mut.BL artefact_percentages.(Sexes{2}).Mut.SD];

artefact_percentages.Average = mean(aggregated);

% Make a table to put all of the artefact percentages in and return
artefact_percentages_Table = table;
artefact_percentages_Table.Sex           = [repmat('M',length(artefact_percentages.(Sexes{1}).WT.BL) + length(artefact_percentages.(Sexes{1}).WT.SD) + length(artefact_percentages.(Sexes{1}).Mut.BL) + length(artefact_percentages.(Sexes{1}).Mut.SD),1); repmat('F',length(artefact_percentages.(Sexes{2}).WT.BL) + length(artefact_percentages.(Sexes{2}).WT.SD) + length(artefact_percentages.(Sexes{2}).Mut.BL) + length(artefact_percentages.(Sexes{2}).Mut.SD),1)]; 
artefact_percentages_Table.Genotype      = [repmat('WT ',length(artefact_percentages.(Sexes{1}).WT.BL)+length(artefact_percentages.(Sexes{1}).WT.SD),1); repmat('Mut',length(artefact_percentages.(Sexes{1}).Mut.BL)+length(artefact_percentages.(Sexes{1}).Mut.SD),1); repmat('WT ',length(artefact_percentages.(Sexes{2}).WT.BL)+length(artefact_percentages.(Sexes{2}).WT.SD),1); repmat('Mut',length(artefact_percentages.(Sexes{2}).Mut.BL)+length(artefact_percentages.(Sexes{2}).Mut.SD),1)]; 
artefact_percentages_Table.BLorSD        = [repmat('BL',length(artefact_percentages.(Sexes{1}).WT.BL),1);repmat('SD',length(artefact_percentages.(Sexes{1}).WT.SD),1); repmat('BL',length(artefact_percentages.(Sexes{1}).Mut.BL),1); repmat('SD',length(artefact_percentages.(Sexes{1}).Mut.SD),1);repmat('BL',length(artefact_percentages.(Sexes{2}).WT.BL),1);repmat('SD',length(artefact_percentages.(Sexes{2}).WT.SD),1); repmat('BL',length(artefact_percentages.(Sexes{2}).Mut.BL),1); repmat('SD',length(artefact_percentages.(Sexes{2}).Mut.SD),1)];
artefact_percentages_Table.Artefact_percentage = [artefact_percentages.(Sexes{1}).WT.BL'; artefact_percentages.(Sexes{1}).WT.SD'; artefact_percentages.(Sexes{1}).Mut.BL'; artefact_percentages.(Sexes{1}).Mut.SD'; artefact_percentages.(Sexes{2}).WT.BL'; artefact_percentages.(Sexes{2}).WT.SD'; artefact_percentages.(Sexes{2}).Mut.BL'; artefact_percentages.(Sexes{2}).Mut.SD']; 

artefact_percentages_Table.Artefact_percentage = strcat(num2str(artefact_percentages_Table.Artefact_percentage*100, '%.1f'), '%');

overall_av_perc = artefact_percentages.Average;



if overall_av_perc > 0.20 | sum(aggregated>0.3)>0
    answer = questdlg('Across all recordings, the average percentage of artifacts is > 20 percent or there is at least one recording with > 30 percent artefact. Would you like me to stop and tell you which file(s) have too many artifacts?','Artifacts Issue','STOP. And tell me.','Keep Going.  It''s fine.','Keep Going.  It''s fine.');
    drawnow; pause(0.05);
    if strcmp(answer,'STOP. And tell me.')
        disp('Here are the files and each files'' percentage of artifacts:')
        disp(artefact_percentages_Table)
        exit_code = 1;
    else
    	exit_code = 0;
    end
else 
	exit_code = 0;
end


 