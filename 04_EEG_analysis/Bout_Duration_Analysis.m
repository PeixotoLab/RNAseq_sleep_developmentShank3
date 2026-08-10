function [Twelve_hour_avg_bout_durations] = Bout_Duration_Analysis(options)
%
% 	USAGE: [Twelve_hour_avg_bout_durations, TwoWayANOVA_tables_bout_durations] = Bout_Duration_Analysis(TimeInState_struct,TimeFrameRMANOVA,inBLorSD,Sex)
%
% 	This function carries out the analysis of bout each arousal state and makes corresponding figure.
%   
% 	INPUTS: 	Scores:             A struct. Here is an example: Scores.Male.WT.BL
%                                    where the first level can be Male or Female 
%                                    second level can be WT or Mut  
%                                    third level can be BL or SD.  This should be for a 24-hour period
%                                    
%
%               Sex:                A string 'Male' or 'Female' 
%               sleep_dep_length:   duration in hours of the sleep dep
%               epoch_duration:     in seconds
%               BoutMinimums:       a struct with 3 fields: W, N, and R. Values are the minimum number of consecutive epochs of that state to count as a bout.  
% 
%
%  OUTPUTS:     Twelve_hour_avg_bout_counts:        A struct with the following form: 
%                                                   Twelve_hour_avg_bout_counts.Female.WT.BL.First12hrs.Wake 
%                                                   and this contains 
% ----------------------------------------------------------------                              

arguments
    options.Scores
    options.Sex              = 'Male'
    options.sleep_dep_length = 5;
    options.epoch_duration   = 4;
    options.BoutMinimums
    options.WindowLengthHrs  = 12;
end 

Scores           = options.Scores;
Sex              = options.Sex;
sleep_dep_length = options.sleep_dep_length;
epoch_duration   = options.epoch_duration;
bout_mins        = options.BoutMinimums;
WindowLengthHrs  = options.WindowLengthHrs;


% W_bout_min = 0;
% N_bout_min = 0;
% R_bout_min = 0;  


% If the Scores struct is empty (for Female and there is no Female data), just return empty structs
if isempty(Scores.(Sex).WT.BL) & isempty(Scores.(Sex).WT.SD) 

    Twelve_hour_avg_bout_durations    = [];
    TwoWayANOVA_tables_bout_durations = [];

    return 
end 


% Sanity check: Do you have the right number of epochs for a 24-hour recording?  
for i=1:length(Scores.(Sex).WT.BL)
    if length(Scores.(Sex).WT.BL{i}) ~= (24*60*60)/epoch_duration
        error('In funcion Bout_Duration_Analysis, a  WT BL recording has the wrong number of bouts.') 
    end 
end 
for i=1:length(Scores.(Sex).WT.SD)
    if length(Scores.(Sex).WT.SD{i}) ~= (24*60*60)/epoch_duration
        error('In funcion Bout_Duration_Analysis, a  WT SD recording has the wrong number of bouts.') 
    end 
end 

% Build the struct to return the output
Twelve_hour_avg_bout_durations = struct();
% ---- Baseline ------------


% --- WT -------------------------------
% ------------------------------------------------------------------------------
% ---- WT BL First12Hrs ----  (do all 3 states here)
for i=1:length(Scores.(Sex).WT.BL)

    scores_this_segment = Scores.(Sex).WT.BL{i}(1:length(Scores.(Sex).WT.BL{i})/2);

    WAKEEpoch = (strcmp(scores_this_segment,'WA'))| (strcmp(scores_this_segment,'W')) | (strcmp(scores_this_segment,'W*'));
    NREMEpoch = (strcmp(scores_this_segment,'N')) | (strcmp(scores_this_segment,'NR'))| (strcmp(scores_this_segment,'N*'));
    REMEpoch  = (strcmp(scores_this_segment,'R')) | (strcmp(scores_this_segment,'RR'))| (strcmp(scores_this_segment,'R*'));

    scores_num = WAKEEpoch + 2*NREMEpoch + 3*REMEpoch;  % Convert W, N, and R to 1, 2, and 3

    runs = contiguous(scores_num);
    W_idx = find([runs{:,1}]==1);
    N_idx = find([runs{:,1}]==2);
    R_idx = find([runs{:,1}]==3);

    W_bout_startstops = runs{W_idx,2};
    W_bout_durations  = W_bout_startstops(:,2)-W_bout_startstops(:,1) + 1;
	W_bout_durations  = W_bout_durations(W_bout_durations>=bout_mins.W);

	N_bout_startstops = runs{N_idx,2};
    N_bout_durations  = N_bout_startstops(:,2)-N_bout_startstops(:,1) + 1;
	N_bout_durations  = N_bout_durations(N_bout_durations>=bout_mins.N);

	R_bout_startstops = runs{R_idx,2};
    R_bout_durations  = R_bout_startstops(:,2)-R_bout_startstops(:,1) + 1;
	R_bout_durations  = R_bout_durations(R_bout_durations>=bout_mins.R);



    Twelve_hour_avg_bout_durations.WT.BL.First12hrs.Wake(i) = mean(W_bout_durations)*epoch_duration/60; % the average wake bout duration (in minutes) for this animal 
    Twelve_hour_avg_bout_durations.WT.BL.First12hrs.NREM(i) = mean(N_bout_durations)*epoch_duration/60; % the average NREM bout duration (in minutes) for this animal 
    Twelve_hour_avg_bout_durations.WT.BL.First12hrs.REM(i)  = mean(R_bout_durations)*epoch_duration/60; % the average REM bout duration (in minutes) for this animal 

end 
% --- WT BL Last12Hrs ----    (do all 3 states here)
for i=1:length(Scores.(Sex).WT.BL)

    scores_this_segment = Scores.(Sex).WT.BL{i}(length(Scores.(Sex).WT.BL{i})/2+1:length(Scores.(Sex).WT.BL{i}));

    WAKEEpoch = (strcmp(scores_this_segment,'WA'))| (strcmp(scores_this_segment,'W')) | (strcmp(scores_this_segment,'W*')); % these return boolean vectors
    NREMEpoch = (strcmp(scores_this_segment,'N')) | (strcmp(scores_this_segment,'NR'))| (strcmp(scores_this_segment,'N*'));
    REMEpoch  = (strcmp(scores_this_segment,'R')) | (strcmp(scores_this_segment,'RR'))| (strcmp(scores_this_segment,'R*'));

    scores_num = WAKEEpoch + 2*NREMEpoch + 3*REMEpoch;  % Convert W, N, and R to 1, 2, and 3

    runs = contiguous(scores_num);  % help contiguous to see outputs
    W_idx = find([runs{:,1}]==1);
    N_idx = find([runs{:,1}]==2);
    R_idx = find([runs{:,1}]==3);

	W_bout_startstops = runs{W_idx,2};
    W_bout_durations  = W_bout_startstops(:,2)-W_bout_startstops(:,1) + 1;
	W_bout_durations  = W_bout_durations(W_bout_durations>=bout_mins.W);

	N_bout_startstops = runs{N_idx,2};
    N_bout_durations  = N_bout_startstops(:,2)-N_bout_startstops(:,1) + 1;
	N_bout_durations  = N_bout_durations(N_bout_durations>=bout_mins.N);

	R_bout_startstops = runs{R_idx,2};
    R_bout_durations  = R_bout_startstops(:,2)-R_bout_startstops(:,1) + 1;
	R_bout_durations  = R_bout_durations(R_bout_durations>=bout_mins.R);


    Twelve_hour_avg_bout_durations.WT.BL.Last12hrs.Wake(i) = mean(W_bout_durations)*epoch_duration/60; % the average wake bout duration (in minutes)  for this animal 
    Twelve_hour_avg_bout_durations.WT.BL.Last12hrs.NREM(i) = mean(N_bout_durations)*epoch_duration/60; % the average NREM bout duration (in minutes)  for this animal 
    Twelve_hour_avg_bout_durations.WT.BL.Last12hrs.REM(i)  = mean(R_bout_durations)*epoch_duration/60; % the average REM bout duration  (in minutes) for this animal

end 

% --- WT SD First12Hrs ---    (do all 3 states here)
for i=1:length(Scores.(Sex).WT.SD)

    scores_this_segment = Scores.(Sex).WT.SD{i}(1:length(Scores.(Sex).WT.SD{i})/2);

    WAKEEpoch = (strcmp(scores_this_segment,'WA'))| (strcmp(scores_this_segment,'W')) | (strcmp(scores_this_segment,'W*'));
    NREMEpoch = (strcmp(scores_this_segment,'N')) | (strcmp(scores_this_segment,'NR'))| (strcmp(scores_this_segment,'N*'));
    REMEpoch  = (strcmp(scores_this_segment,'R')) | (strcmp(scores_this_segment,'RR'))| (strcmp(scores_this_segment,'R*'));

    scores_num = WAKEEpoch + 2*NREMEpoch + 3*REMEpoch;  % Convert W, N, and R to 1, 2, and 3

    runs = contiguous(scores_num);
    W_idx = find([runs{:,1}]==1);
    N_idx = find([runs{:,1}]==2);
    R_idx = find([runs{:,1}]==3);

	W_bout_startstops = runs{W_idx,2};
    W_bout_durations  = W_bout_startstops(:,2)-W_bout_startstops(:,1) + 1;
	W_bout_durations  = W_bout_durations(W_bout_durations>=bout_mins.W);

	N_bout_startstops = runs{N_idx,2};
    N_bout_durations  = N_bout_startstops(:,2)-N_bout_startstops(:,1) + 1;
	N_bout_durations  = N_bout_durations(N_bout_durations>=bout_mins.N);

	R_bout_startstops = runs{R_idx,2};
    R_bout_durations  = R_bout_startstops(:,2)-R_bout_startstops(:,1) + 1;
	R_bout_durations  = R_bout_durations(R_bout_durations>=bout_mins.R);

    Twelve_hour_avg_bout_durations.WT.SD.First12hrs.Wake(i) = mean(W_bout_durations)*epoch_duration/60; % the average wake bout duration (in minutes)  for this animal 
    Twelve_hour_avg_bout_durations.WT.SD.First12hrs.NREM(i) = mean(N_bout_durations)*epoch_duration/60; % the average NREM bout duration (in minutes)  for this animal 
    Twelve_hour_avg_bout_durations.WT.SD.First12hrs.REM(i)  = mean(R_bout_durations)*epoch_duration/60; % the average REM bout duration  (in minutes) for this animal

end 
% --- WT SD Last12Hrs ----    (do all 3 states here)
for i=1:length(Scores.(Sex).WT.SD)

    scores_this_segment = Scores.(Sex).WT.SD{i}(length(Scores.(Sex).WT.SD{i})/2+1:length(Scores.(Sex).WT.SD{i}));

    WAKEEpoch = (strcmp(scores_this_segment,'WA'))| (strcmp(scores_this_segment,'W')) | (strcmp(scores_this_segment,'W*')); % these return boolean vectors
    NREMEpoch = (strcmp(scores_this_segment,'N')) | (strcmp(scores_this_segment,'NR'))| (strcmp(scores_this_segment,'N*'));
    REMEpoch  = (strcmp(scores_this_segment,'R')) | (strcmp(scores_this_segment,'RR'))| (strcmp(scores_this_segment,'R*'));

    scores_num = WAKEEpoch + 2*NREMEpoch + 3*REMEpoch;  % Convert W, N, and R to 1, 2, and 3

    runs = contiguous(scores_num);  % help contiguous to see outputs
    W_idx = find([runs{:,1}]==1);
    N_idx = find([runs{:,1}]==2);
    R_idx = find([runs{:,1}]==3);

	W_bout_startstops = runs{W_idx,2};
    W_bout_durations  = W_bout_startstops(:,2)-W_bout_startstops(:,1) + 1;
	W_bout_durations  = W_bout_durations(W_bout_durations>=bout_mins.W);

	N_bout_startstops = runs{N_idx,2};
    N_bout_durations  = N_bout_startstops(:,2)-N_bout_startstops(:,1) + 1;
	N_bout_durations  = N_bout_durations(N_bout_durations>=bout_mins.N);

	R_bout_startstops = runs{R_idx,2};
    R_bout_durations  = R_bout_startstops(:,2)-R_bout_startstops(:,1) + 1;
	R_bout_durations  = R_bout_durations(R_bout_durations>=bout_mins.R);

    Twelve_hour_avg_bout_durations.WT.SD.Last12hrs.Wake(i) = mean(W_bout_durations)*epoch_duration/60; % the average wake bout duration (in minutes)  for this animal 
    Twelve_hour_avg_bout_durations.WT.SD.Last12hrs.NREM(i) = mean(N_bout_durations)*epoch_duration/60; % the average NREM bout duration (in minutes)  for this animal 
    Twelve_hour_avg_bout_durations.WT.SD.Last12hrs.REM(i)  = mean(R_bout_durations)*epoch_duration/60; % the average REM bout duration  (in minutes) for this animal 

end 

% --- WT SD First12Hrs (excluding actual SD) ---    (do all 3 states here)
for i=1:length(Scores.(Sex).WT.SD)

    num_epochs_in_SD = sleep_dep_length*60*60/epoch_duration;
    scores_this_segment = Scores.(Sex).WT.SD{i}(num_epochs_in_SD+1:length(Scores.(Sex).WT.SD{i})/2);

    WAKEEpoch = (strcmp(scores_this_segment,'WA'))| (strcmp(scores_this_segment,'W')) | (strcmp(scores_this_segment,'W*'));
    NREMEpoch = (strcmp(scores_this_segment,'N')) | (strcmp(scores_this_segment,'NR'))| (strcmp(scores_this_segment,'N*'));
    REMEpoch  = (strcmp(scores_this_segment,'R')) | (strcmp(scores_this_segment,'RR'))| (strcmp(scores_this_segment,'R*'));

    scores_num = WAKEEpoch + 2*NREMEpoch + 3*REMEpoch;  % Convert W, N, and R to 1, 2, and 3

    runs = contiguous(scores_num);
    W_idx = find([runs{:,1}]==1);
    N_idx = find([runs{:,1}]==2);
    R_idx = find([runs{:,1}]==3);

    W_bout_startstops = runs{W_idx,2};
    W_bout_durations  = W_bout_startstops(:,2)-W_bout_startstops(:,1) + 1;
	W_bout_durations  = W_bout_durations(W_bout_durations>=bout_mins.W);

	N_bout_startstops = runs{N_idx,2};
    N_bout_durations  = N_bout_startstops(:,2)-N_bout_startstops(:,1) + 1;
	N_bout_durations  = N_bout_durations(N_bout_durations>=bout_mins.N);

	R_bout_startstops = runs{R_idx,2};
    R_bout_durations  = R_bout_startstops(:,2)-R_bout_startstops(:,1) + 1;
	R_bout_durations  = R_bout_durations(R_bout_durations>=bout_mins.R);


    Twelve_hour_avg_bout_durations.WT.SDexcSD.First12hrs.Wake(i) = mean(W_bout_durations)*epoch_duration/60; % the average wake bout duration (in minutes)  for this animal 
    Twelve_hour_avg_bout_durations.WT.SDexcSD.First12hrs.NREM(i) = mean(N_bout_durations)*epoch_duration/60; % the average NREM bout duration (in minutes)  for this animal 
    Twelve_hour_avg_bout_durations.WT.SDexcSD.First12hrs.REM(i)  = mean(R_bout_durations)*epoch_duration/60; % the average REM bout duration  (in minutes) for this animal

end 


% --- Mutant -------------------------------
% ------------------------------------------------------------------------------
% ---- Mut BL First12Hrs ----  (do all 3 states here)
for i=1:length(Scores.(Sex).Mut.BL)

    scores_this_segment = Scores.(Sex).Mut.BL{i}(1:length(Scores.(Sex).Mut.BL{i})/2);

    WAKEEpoch = (strcmp(scores_this_segment,'WA'))| (strcmp(scores_this_segment,'W')) | (strcmp(scores_this_segment,'W*'));
    NREMEpoch = (strcmp(scores_this_segment,'N')) | (strcmp(scores_this_segment,'NR'))| (strcmp(scores_this_segment,'N*'));
    REMEpoch  = (strcmp(scores_this_segment,'R')) | (strcmp(scores_this_segment,'RR'))| (strcmp(scores_this_segment,'R*'));

    scores_num = WAKEEpoch + 2*NREMEpoch + 3*REMEpoch;  % Convert W, N, and R to 1, 2, and 3

    runs = contiguous(scores_num);
    W_idx = find([runs{:,1}]==1);
    N_idx = find([runs{:,1}]==2);
    R_idx = find([runs{:,1}]==3);

	W_bout_startstops = runs{W_idx,2};
    W_bout_durations  = W_bout_startstops(:,2)-W_bout_startstops(:,1) + 1;
	W_bout_durations  = W_bout_durations(W_bout_durations>=bout_mins.W);

	N_bout_startstops = runs{N_idx,2};
    N_bout_durations  = N_bout_startstops(:,2)-N_bout_startstops(:,1) + 1;
	N_bout_durations  = N_bout_durations(N_bout_durations>=bout_mins.N);

	R_bout_startstops = runs{R_idx,2};
    R_bout_durations  = R_bout_startstops(:,2)-R_bout_startstops(:,1) + 1;
	R_bout_durations  = R_bout_durations(R_bout_durations>=bout_mins.R);

	Twelve_hour_avg_bout_durations.Mut.BL.First12hrs.Wake(i) = mean(W_bout_durations)*epoch_duration/60; % the average wake bout duration (in minutes)  for this animal 
    Twelve_hour_avg_bout_durations.Mut.BL.First12hrs.NREM(i) = mean(N_bout_durations)*epoch_duration/60; % the average NREM bout duration (in minutes)  for this animal 
    Twelve_hour_avg_bout_durations.Mut.BL.First12hrs.REM(i)  = mean(R_bout_durations)*epoch_duration/60; % the average REM bout duration  (in minutes) for this animal

end 
% --- Mut BL Last12Hrs ----    (do all 3 states here)
for i=1:length(Scores.(Sex).Mut.BL)

    scores_this_segment = Scores.(Sex).Mut.BL{i}(length(Scores.(Sex).Mut.BL{i})/2+1:length(Scores.(Sex).Mut.BL{i}));

    WAKEEpoch = (strcmp(scores_this_segment,'WA'))| (strcmp(scores_this_segment,'W')) | (strcmp(scores_this_segment,'W*')); % these return boolean vectors
    NREMEpoch = (strcmp(scores_this_segment,'N')) | (strcmp(scores_this_segment,'NR'))| (strcmp(scores_this_segment,'N*'));
    REMEpoch  = (strcmp(scores_this_segment,'R')) | (strcmp(scores_this_segment,'RR'))| (strcmp(scores_this_segment,'R*'));

    scores_num = WAKEEpoch + 2*NREMEpoch + 3*REMEpoch;  % Convert W, N, and R to 1, 2, and 3

    runs = contiguous(scores_num);  % help contiguous to see outputs
    W_idx = find([runs{:,1}]==1);
    N_idx = find([runs{:,1}]==2);
    R_idx = find([runs{:,1}]==3);

	W_bout_startstops = runs{W_idx,2};
    W_bout_durations  = W_bout_startstops(:,2)-W_bout_startstops(:,1) + 1;
	W_bout_durations  = W_bout_durations(W_bout_durations>=bout_mins.W);

	N_bout_startstops = runs{N_idx,2};
    N_bout_durations  = N_bout_startstops(:,2)-N_bout_startstops(:,1) + 1;
	N_bout_durations  = N_bout_durations(N_bout_durations>=bout_mins.N);

	R_bout_startstops = runs{R_idx,2};
    R_bout_durations  = R_bout_startstops(:,2)-R_bout_startstops(:,1) + 1;
	R_bout_durations  = R_bout_durations(R_bout_durations>=bout_mins.R);

    Twelve_hour_avg_bout_durations.Mut.BL.Last12hrs.Wake(i) = mean(W_bout_durations)*epoch_duration/60; % the average wake bout duration (in minutes)  for this animal 
    Twelve_hour_avg_bout_durations.Mut.BL.Last12hrs.NREM(i) = mean(N_bout_durations)*epoch_duration/60; % the average NREM bout duration (in minutes)  for this animal 
    Twelve_hour_avg_bout_durations.Mut.BL.Last12hrs.REM(i)  = mean(R_bout_durations)*epoch_duration/60; % the average REM bout duration  (in minutes) for this animal

end 

% --- Mut SD First12Hrs ---    (do all 3 states here)
for i=1:length(Scores.(Sex).Mut.SD)

    scores_this_segment = Scores.(Sex).Mut.SD{i}(1:length(Scores.(Sex).Mut.SD{i})/2);

    WAKEEpoch = (strcmp(scores_this_segment,'WA'))| (strcmp(scores_this_segment,'W')) | (strcmp(scores_this_segment,'W*'));
    NREMEpoch = (strcmp(scores_this_segment,'N')) | (strcmp(scores_this_segment,'NR'))| (strcmp(scores_this_segment,'N*'));
    REMEpoch  = (strcmp(scores_this_segment,'R')) | (strcmp(scores_this_segment,'RR'))| (strcmp(scores_this_segment,'R*'));

    scores_num = WAKEEpoch + 2*NREMEpoch + 3*REMEpoch;  % Convert W, N, and R to 1, 2, and 3

    runs = contiguous(scores_num);
    W_idx = find([runs{:,1}]==1);
    N_idx = find([runs{:,1}]==2);
    R_idx = find([runs{:,1}]==3);

    W_bout_startstops = runs{W_idx,2};
    W_bout_durations  = W_bout_startstops(:,2)-W_bout_startstops(:,1) + 1;
	W_bout_durations  = W_bout_durations(W_bout_durations>=bout_mins.W);

	N_bout_startstops = runs{N_idx,2};
    N_bout_durations  = N_bout_startstops(:,2)-N_bout_startstops(:,1) + 1;
	N_bout_durations  = N_bout_durations(N_bout_durations>=bout_mins.N);

	R_bout_startstops = runs{R_idx,2};
    R_bout_durations  = R_bout_startstops(:,2)-R_bout_startstops(:,1) + 1;
	R_bout_durations  = R_bout_durations(R_bout_durations>=bout_mins.R);


    Twelve_hour_avg_bout_durations.Mut.SD.First12hrs.Wake(i) = mean(W_bout_durations)*epoch_duration/60; % the average wake bout duration (in minutes)  for this animal 
    Twelve_hour_avg_bout_durations.Mut.SD.First12hrs.NREM(i) = mean(N_bout_durations)*epoch_duration/60; % the average NREM bout duration (in minutes)  for this animal 
    Twelve_hour_avg_bout_durations.Mut.SD.First12hrs.REM(i)  = mean(R_bout_durations)*epoch_duration/60; % the average REM bout duration  (in minutes) for this animal

end 
% --- Mut SD Last12Hrs ----    (do all 3 states here)
for i=1:length(Scores.(Sex).Mut.SD)

    scores_this_segment = Scores.(Sex).Mut.SD{i}(length(Scores.(Sex).Mut.SD{i})/2+1:length(Scores.(Sex).Mut.SD{i}));

    WAKEEpoch = (strcmp(scores_this_segment,'WA'))| (strcmp(scores_this_segment,'W')) | (strcmp(scores_this_segment,'W*')); % these return boolean vectors
    NREMEpoch = (strcmp(scores_this_segment,'N')) | (strcmp(scores_this_segment,'NR'))| (strcmp(scores_this_segment,'N*'));
    REMEpoch  = (strcmp(scores_this_segment,'R')) | (strcmp(scores_this_segment,'RR'))| (strcmp(scores_this_segment,'R*'));

    scores_num = WAKEEpoch + 2*NREMEpoch + 3*REMEpoch;  % Convert W, N, and R to 1, 2, and 3

    runs = contiguous(scores_num);  % help contiguous to see outputs
    W_idx = find([runs{:,1}]==1);
    N_idx = find([runs{:,1}]==2);
    R_idx = find([runs{:,1}]==3);

    W_bout_startstops = runs{W_idx,2};
    W_bout_durations  = W_bout_startstops(:,2)-W_bout_startstops(:,1) + 1;
	W_bout_durations  = W_bout_durations(W_bout_durations>=bout_mins.W);

	N_bout_startstops = runs{N_idx,2};
    N_bout_durations  = N_bout_startstops(:,2)-N_bout_startstops(:,1) + 1;
	N_bout_durations  = N_bout_durations(N_bout_durations>=bout_mins.N);

	R_bout_startstops = runs{R_idx,2};
    R_bout_durations  = R_bout_startstops(:,2)-R_bout_startstops(:,1) + 1;
	R_bout_durations  = R_bout_durations(R_bout_durations>=bout_mins.R);


    Twelve_hour_avg_bout_durations.Mut.SD.Last12hrs.Wake(i) = mean(W_bout_durations)*epoch_duration/60; % the average wake bout duration (in minutes)  for this animal 
    Twelve_hour_avg_bout_durations.Mut.SD.Last12hrs.NREM(i) = mean(N_bout_durations)*epoch_duration/60; % the average NREM bout duration (in minutes)  for this animal 
    Twelve_hour_avg_bout_durations.Mut.SD.Last12hrs.REM(i)  = mean(R_bout_durations)*epoch_duration/60; % the average REM bout duration  (in minutes) for this animal

end 

% --- Mut SD (excluding the actual SD) First12Hrs ---    (do all 3 states here)
for i=1:length(Scores.(Sex).Mut.SD)

    num_epochs_in_SD = sleep_dep_length*60*60/epoch_duration;
    scores_this_segment = Scores.(Sex).Mut.SD{i}(num_epochs_in_SD+1:length(Scores.(Sex).Mut.SD{i})/2);

    WAKEEpoch = (strcmp(scores_this_segment,'WA'))| (strcmp(scores_this_segment,'W')) | (strcmp(scores_this_segment,'W*'));
    NREMEpoch = (strcmp(scores_this_segment,'N')) | (strcmp(scores_this_segment,'NR'))| (strcmp(scores_this_segment,'N*'));
    REMEpoch  = (strcmp(scores_this_segment,'R')) | (strcmp(scores_this_segment,'RR'))| (strcmp(scores_this_segment,'R*'));

    scores_num = WAKEEpoch + 2*NREMEpoch + 3*REMEpoch;  % Convert W, N, and R to 1, 2, and 3

    runs = contiguous(scores_num);
    W_idx = find([runs{:,1}]==1);
    N_idx = find([runs{:,1}]==2);
    R_idx = find([runs{:,1}]==3);

	W_bout_startstops = runs{W_idx,2};
    W_bout_durations  = W_bout_startstops(:,2)-W_bout_startstops(:,1) + 1;
	W_bout_durations  = W_bout_durations(W_bout_durations>=bout_mins.W);
	
	N_bout_startstops = runs{N_idx,2};
    N_bout_durations  = N_bout_startstops(:,2)-N_bout_startstops(:,1) + 1;
	N_bout_durations  = N_bout_durations(N_bout_durations>=bout_mins.N);

	R_bout_startstops = runs{R_idx,2};
    R_bout_durations  = R_bout_startstops(:,2)-R_bout_startstops(:,1) + 1;
	R_bout_durations  = R_bout_durations(R_bout_durations>=bout_mins.R);


    Twelve_hour_avg_bout_durations.Mut.SDexcSD.First12hrs.Wake(i) = mean(W_bout_durations)*epoch_duration/60; % the average wake bout duration (in minutes)  for this animal 
    Twelve_hour_avg_bout_durations.Mut.SDexcSD.First12hrs.NREM(i) = mean(N_bout_durations)*epoch_duration/60; % the average NREM bout duration (in minutes)  for this animal 
    Twelve_hour_avg_bout_durations.Mut.SDexcSD.First12hrs.REM(i)  = mean(R_bout_durations)*epoch_duration/60; % the average REM bout duration  (in minutes) for this animal


end 
% ----------------------------------------------------------------------



