function efficiency = calc_sleep_dep_efficiency(Scores,Sexes,SD_length_hrs,epoch_duration_secs)
%
% USAGE: calc_sleep_dep_efficiency(Scores,SD_length_hrs)
%
% This function calculates the efficiency of the Sleep Dep
% calculated as the percentage of the sleep dep that the animal was actually awake. 
%
%
% INPUTS: 	Scores 				A struct with the following format:  Scores.Male.WT.SD 
% 								each of these fields contains N entries that are Mx1 where
%								N is the number of animals, and M is the # of epochs in each recording.
% 			SD_length_hrs 		self-explanatory
%			epoch_duration_secs	self-explanatory
%
% OUTPUT: 	The function throws an error if either of the following happens: 1) a SD recording
% 			has sleep dep efficiency < 90%, or a BL recording has more than 90% Wake during the 
% 			the first SD_length_hrs of the recording.  


cut_off = .90;  

num_epochs_in_SD = SD_length_hrs*60*60/epoch_duration_secs; 
% --- Make sure all SD days have at least 90% wakefulness during SD -------------------------------
% -- Male WT --
for i=1:length(Scores.(Sexes{1}).WT.SD)
	scores_during_SD = Scores.(Sexes{1}).WT.SD{i}(1:num_epochs_in_SD);
	num_non_artefact_epochs_in_SD = length(scores_during_SD) - length(find(strcmp(scores_during_SD, 'M'))); 
	efficiency.M.WT(i) = sum(strcmp(scores_during_SD, 'W') | strcmp(scores_during_SD,'W*') | strcmp(scores_during_SD, 'WA') | strcmp(scores_during_SD,'WA*'))/num_non_artefact_epochs_in_SD;
end 

% -- Male Mut --
for i=1:length(Scores.(Sexes{1}).Mut.SD)
	scores_during_SD = Scores.(Sexes{1}).Mut.SD{i}(1:num_epochs_in_SD);
	num_non_artefact_epochs_in_SD = length(scores_during_SD) - length(find(strcmp(scores_during_SD, 'M'))); 
	efficiency.M.Mut(i) = sum(strcmp(scores_during_SD, 'W') | strcmp(scores_during_SD,'W*') | strcmp(scores_during_SD, 'WA') | strcmp(scores_during_SD,'WA*'))/num_non_artefact_epochs_in_SD;
end 

% -- Female WT --
if ~isempty(Scores.(Sexes{2}).WT.SD)
	for i=1:length(Scores.(Sexes{2}).WT.SD)
		scores_during_SD = Scores.(Sexes{2}).WT.SD{i}(1:num_epochs_in_SD);
		num_non_artefact_epochs_in_SD = length(scores_during_SD) - length(find(strcmp(scores_during_SD, 'M'))); 
		efficiency.F.WT(i) = sum(strcmp(scores_during_SD, 'W') | strcmp(scores_during_SD,'W*') | strcmp(scores_during_SD, 'WA') | strcmp(scores_during_SD,'WA*'))/num_non_artefact_epochs_in_SD;
	end 
else
	efficiency.F.WT = [];
end 

% -- Female Mut --
if ~isempty(Scores.(Sexes{2}).Mut.SD)
	for i=1:length(Scores.(Sexes{2}).Mut.SD)
		scores_during_SD = Scores.(Sexes{2}).Mut.SD{i}(1:num_epochs_in_SD);
		num_non_artefact_epochs_in_SD = length(scores_during_SD) - length(find(strcmp(scores_during_SD, 'M'))); 
		efficiency.F.Mut(i) = sum(strcmp(scores_during_SD, 'W') | strcmp(scores_during_SD,'W*') | strcmp(scores_during_SD, 'WA') | strcmp(scores_during_SD,'WA*'))/num_non_artefact_epochs_in_SD;
	end 
else
	efficiency.F.Mut = [];
end 


err_mess1 = '';
err_mess2 = '';
err_mess3 = '';
err_mess4 = '';

if sum(efficiency.M.WT < cut_off) > 0 
	%error(['You have ',num2str(sum(efficiency.M.WT < cut_off)),' Male WT SD recording(s) with sleep dep efficiency less than ',num2str(100*cut_off),' %.'])
	err_mess1 = append(Sexes{1}, ' WT had ',num2str(sum(efficiency.M.WT < cut_off)),' recording(s) with sleep dep efficiency less than ',num2str(100*cut_off),' percent.\n');
end 

if sum(efficiency.M.Mut < cut_off) > 0 
	%error(['You have ',num2str(sum(efficiency.M.Mut < cut_off)),' Male Mut SD recording(s) with sleep dep efficiency less than ',num2str(100*cut_off),' %.'])
	err_mess2 = append(Sexes{1}, ' Mut had ',num2str(sum(efficiency.M.Mut < cut_off)),' recording(s) with sleep dep efficiency less than ',num2str(100*cut_off),' percent.\n');
end

if ~isempty(efficiency.F.WT) & sum(efficiency.F.WT < cut_off) > 0 
	%error(['You have ',num2str(sum(efficiency.F.WT < cut_off)),' Female WT SD recording(s) with sleep dep efficiency less than ',num2str(100*cut_off),' %.'])
	err_mess3 = append(Sexes{2}, ' WT had ',num2str(sum(efficiency.F.WT < cut_off)),' recording(s) with sleep dep efficiency less than ',num2str(100*cut_off),' percent.\n');
end 

if ~isempty(efficiency.F.Mut) & sum(efficiency.F.Mut < cut_off) > 0 
	%error(['You have ',num2str(sum(efficiency.F.Mut < cut_off)),' Female Mut SD recording(s) with sleep dep efficiency less than ',num2str(100*cut_off),' %.'])
	err_mess4 = append(Sexes{2}, ' Mut had ',num2str(sum(efficiency.F.Mut < cut_off)),' recording(s) with sleep dep efficiency less than ',num2str(100*cut_off),' percent.\n');
end

error_message = [err_mess1 , ...
				 err_mess2 , ...
				 err_mess3 , ...
				 err_mess4];

% This is now done in the main script so you can tell the user which recording(s) have the low efficiency.  
% if ~isempty(error_message) 
% 	answer = questdlg('At least one of the recordings had a sleep dep efficiency less than 90 percent. What should I do?','Sleep Dep Efficiency Issue','STOP','Keep Going.  It''s fine.','Keep Going.  It''s fine.');
% 	drawnow; pause(0.05);  % this prevents questdlg from hanging
% 	if strcmp(answer,'STOP')
% 		error('myerror:with:colon',error_message)
% 	end 
% end 
% --- End of Checking SD day ----------------------------------------------------------------------
% -------------------------------------------------------------------------------------------------

% --- Make sure all BL days have less than 90% wakefulness during the first SD_length_hrs hours ---
% -- Male WT --
for i=1:length(Scores.(Sexes{1}).WT.BL)
	scores_during_BL = Scores.(Sexes{1}).WT.BL{i}(1:num_epochs_in_SD);
	PercentW_BL.M.WT(i) = sum(strcmp(scores_during_BL, 'W') | strcmp(scores_during_BL,'W*') | strcmp(scores_during_BL, 'WA') | strcmp(scores_during_BL,'WA*'))/num_epochs_in_SD;
end 

% -- Male Mut --
for i=1:length(Scores.(Sexes{1}).Mut.BL)
	scores_during_BL = Scores.(Sexes{1}).Mut.BL{i}(1:num_epochs_in_SD);
	PercentW_BL.M.Mut(i) = sum(strcmp(scores_during_BL, 'W') | strcmp(scores_during_BL,'W*') | strcmp(scores_during_BL, 'WA') | strcmp(scores_during_BL,'WA*'))/num_epochs_in_SD;
end 

% -- Female WT --
if ~isempty(Scores.(Sexes{2}).WT.BL)
	for i=1:length(Scores.(Sexes{2}).WT.BL)
		scores_during_BL = Scores.(Sexes{2}).WT.BL{i}(1:num_epochs_in_SD);
		PercentW_BL.F.WT(i) = sum(strcmp(scores_during_BL, 'W') | strcmp(scores_during_BL,'W*') | strcmp(scores_during_BL, 'WA') | strcmp(scores_during_BL,'WA*'))/num_epochs_in_SD;
	end 
else
	PercentW_BL.F.WT = [];
end 

% -- Female Mut --
if ~isempty(Scores.(Sexes{2}).Mut.BL)
	for i=1:length(Scores.(Sexes{2}).Mut.BL)
		scores_during_BL = Scores.(Sexes{2}).Mut.BL{i}(1:num_epochs_in_SD);
		PercentW_BL.F.Mut(i) = sum(strcmp(scores_during_BL, 'W') | strcmp(scores_during_BL,'W*') | strcmp(scores_during_BL, 'WA') | strcmp(scores_during_BL,'WA*'))/num_epochs_in_SD;
	end 
else
	PercentW_BL.F.Mut = [];
end 


err_messBL1 = '';
err_messBL2 = '';
err_messBL3 = '';
err_messBL4 = '';

if sum(PercentW_BL.M.WT > cut_off) > 0 
	err_messBL1 = append('Male WT had ',num2str(sum(PercentW_BL.M.WT > cut_off)),' recording(s) with Baseline Wakefulness more than ',num2str(100*cut_off),' percent. Perhaps a BL and SD file got mixed up.\n');
end 

if sum(PercentW_BL.M.Mut > cut_off) > 0 
	err_messBL2 = append('Male Mut had ',num2str(sum(PercentW_BL.M.Mut > cut_off)),' recording(s) with Baseline Wakefulness more than ',num2str(100*cut_off),' percent. Perhaps a BL and SD file got mixed up.\n');
end

if ~isempty(PercentW_BL.F.WT) & sum(PercentW_BL.F.WT > cut_off) > 0 
	err_messBL3 = append('Female WT had ',num2str(sum(PercentW_BL.F.WT > cut_off)),' recording(s) with Baseline Wakefulness more than ',num2str(100*cut_off),' percent. Perhaps a BL and SD file got mixed up.\n');
end 

if ~isempty(PercentW_BL.F.Mut) & sum(PercentW_BL.F.Mut > cut_off) > 0 
	err_messBL4 = append('Female Mut had ',num2str(sum(PercentW_BL.F.Mut > cut_off)),' recording(s) with Baseline Wakefulness more than ',num2str(100*cut_off),' percent. Perhaps a BL and SD file got mixed up.\n');
end

error_messageBL = [err_messBL1 , ...
				   err_messBL2 , ...
				   err_messBL3 , ...
				   err_messBL4];
if ~isempty(error_messageBL)
	error('myerror:with:colonBL',error_messageBL)
end 


