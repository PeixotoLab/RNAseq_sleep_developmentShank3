function Avg_bout_counts = calculate_avg_bout_counts(options)
%
% USAGE: Avg_bout_counts = calculate_avg_bout_counts(Scores,Sex,WTorMut,BLorSD,bout_mins,TimeFrame)
%
%
% INPUTS:	Scores: struct created in ExtractEEGData. Scores.Male.WT.BL is a 1xN cell array where N
%					is the number of animals, and each cell is an nx1 cell where n is the number of epochs in a day,
%					typically 21600 for 4-second epochs.  
%
%			Sex:  		'M' or 'F'
% 			WTorMut:  	Wildtype or Mutant
% 			BLorSD:  	baseline or sleep dep
% 			bout_mins:	bout minimums (what is considered a W bout?) set in Bout_Minimums in the main script
% 			TimeFrame:  a 2-element vector specifying the hours from the beginning of that day (BL or SD)
% 						examples:  [1 12] first 12 hours, [1 6] first six hours, [13 24] last 12 hours

arguments
    options.Scores
    options.Sex              = 'Male';
    options.WTorMut 		 
    options.BLorSD   
    options.BoutMinimums
    options.TimeFrame  
end 

Scores   	= options.Scores;
Sex       	= options.Sex;
WTorMut		= options.WTorMut;
BLorSD   	= options.BLorSD;
bout_mins 	= options.BoutMinimums;
TimeFrame	= options.TimeFrame;




L = length(Scores{1});  % should be 21600 if 4-second epochs and a 24-hr recording

if mod(L,24) ~= 0 
	error('Scores does not have the right number of epochs.  It needs to be divisible by 24.')
end 
num_epochs_in_an_hr = L/24;
a = TimeFrame(1);
b = TimeFrame(2);
time_ind = [(a-1)*num_epochs_in_an_hr+1:num_epochs_in_an_hr*b];
%time_ind = [(a-1)*(L/24)+1:L/(24/b)];


% % set up the struct field name based on the timeframe chosen
% if isequal(TimeFrame,[1 12])
% 	TimeFrametime_field_str = "First12hrs";
% elseif isequal(TimeFrame,[1 6])
% 	TimeFrametime_field_str = "Hrs1_6";
% elseif isequal(TimeFrame,[7 12])
% 	TimeFrametime_field_str = "Hrs7_12";
% elseif isequal(TimeFrame,[13 24])
% 	TimeFrametime_field_str = "Last12hrs";
% elseif isequal(TimeFrame,[13 18])
% 	TimeFrametime_field_str = "Hrs13_18";
% elseif isequal(TimeFrame,[19 24])
% 	TimeFrametime_field_str = "Hrs19_24";
% else
% 	error("In calculate_avg_bout_counts.m: you chose a TimeFrame that I was not expecting.  Try again.")
% end 



% First confirm that all entries in Scores have the same length.  If not, throw an error
if ~all(cellfun(@(e) isequal(size(Scores{1}), size(e)) , Scores(2:end)))
	error('In calculate_avg_bout_counts.m: the scores cells are not the same size')
end  

for i=1:length(Scores)

    scores_this_segment = Scores{i}(time_ind);

    WAKEEpoch = (strcmp(scores_this_segment,'WA'))| (strcmp(scores_this_segment,'W')) | (strcmp(scores_this_segment,'W*'));
    NREMEpoch = (strcmp(scores_this_segment,'N')) | (strcmp(scores_this_segment,'NR'))| (strcmp(scores_this_segment,'N*'));
    REMEpoch  = (strcmp(scores_this_segment,'R')) | (strcmp(scores_this_segment,'RR'))| (strcmp(scores_this_segment,'R*'));

    scores_num = WAKEEpoch + 2*NREMEpoch + 3*REMEpoch;  % Convert W, N, and R to 1, 2, and 3

    runs = contiguous(scores_num);
    W_idx = find([runs{:,1}]==1);
    N_idx = find([runs{:,1}]==2);
    R_idx = find([runs{:,1}]==3);

    % Wake
    if ~isempty(W_idx) 					% If there are some W runs
	    W_bout_startstops = runs{W_idx,2};
	    W_bout_durations  = W_bout_startstops(:,2)-W_bout_startstops(:,1) + 1;
	    W_bout_startstops = W_bout_startstops(W_bout_durations>=bout_mins.W,:);  % Keep only those bout startstops that are longer than the mininum bout length
	else 
		W_bout_startstops = [];   		% handle the case where no W in this window
	end  
    
    % NREM
    if ~isempty(N_idx)
	    N_bout_startstops = runs{N_idx,2};
	    N_bout_durations  = N_bout_startstops(:,2)-N_bout_startstops(:,1) + 1;
	    N_bout_startstops = N_bout_startstops(N_bout_durations>=bout_mins.N,:);  % Keep only those bout startstops that are longer than the mininum bout length
	else 
		N_bout_startstops = [];  		% handle the case where no N in this window
	end 

	% REM
	if ~isempty(R_idx)
	    R_bout_startstops = runs{R_idx,2};
	    R_bout_durations  = R_bout_startstops(:,2)-R_bout_startstops(:,1) + 1;
	    R_bout_startstops = R_bout_startstops(R_bout_durations>=bout_mins.R,:);  % Keep only those bout startstops that are longer than the mininum bout length
	else 
		R_bout_startstops = [];  		% handle the case where no R in this window
	end 

    Avg_bout_counts.Wake(i) = size(W_bout_startstops,1)/(TimeFrame(2)-TimeFrame(1)+1); % the number of wake bouts for this animal (divide by hours to get bouts/hour)
    Avg_bout_counts.NREM(i) = size(N_bout_startstops,1)/(TimeFrame(2)-TimeFrame(1)+1); % the number of NREM bouts for this animal (divide by hours to get bouts/hour)
    Avg_bout_counts.REM(i)  = size(R_bout_startstops,1)/(TimeFrame(2)-TimeFrame(1)+1); % the number of REM bouts for this animal (divide by hours to get bouts/hour)


    Avg_bout_durations.Wake(i) = mean(W_bout_durations)*epoch_duration/60; % the average wake bout duration (in minutes) for this animal 
    Avg_bout_durations.NREM(i) = mean(N_bout_durations)*epoch_duration/60; % the average NREM bout duration (in minutes) for this animal 
    Avg_bout_durations.REM(i)  = mean(R_bout_durations)*epoch_duration/60; % the average REM bout duration (in minutes) for this animal 

	% TESTING:  return raw counts, not normalized by the number of hours
	% Avg_bout_counts.Wake(i) = size(W_bout_startstops,1); % the number of wake bouts for this animal 
    % Avg_bout_counts.NREM(i) = size(N_bout_startstops,1); % the number of NREM bouts for this animal 
    % Avg_bout_counts.REM(i)  = size(R_bout_startstops,1); % the number of REM bouts for this animal 

end 