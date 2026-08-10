function [Avg_bout_counts] = Bout_Number_Analysis(options)
%
% 	USAGE: [Twelve_hour_avg_bouts, TwoWayANOVA_tables_bout_duration] = Bout_Number_Analysis(TimeInState_struct,TimeFrameRMANOVA,inBLorSD,Sex)
%
% 	This function carries out the analysis of bout each arousal state and makes corresponding figure.
%   
% 	INPUTS: 	Scores:  A struct. Here is an example: Scores.Male.WT.BL
%                                    where the first level can be Male or Female 
%                                    second level can be WT or Mut  
%                                    third level can be BL or SD.  This should be for a 24-hour period
%                                    
%
%               Sex:     A string 'Male' or 'Female'  
%               sleep_dep_length:   duration in hours of the sleep dep
%               epoch_duration:     in seconds
%               BoutMinimums:       a struct with 3 fields: W, N, and R. Values are the minimum number of consecutive epochs of that state to count as a bout.  
% 
% 
%
%  OUTPUTS:     Avg_bout_counts:        A struct with the following form: 
%                                                   Avg_bout_counts.Female.WT.BL.First12hrs.Wake 
%                                                   and this contains bout counts normalized by hour 
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



% If the Scores struct is empty (for Female and there is no Female data), just return empty structs
if isempty(Scores.(Sex).WT.BL) & isempty(Scores.(Sex).WT.SD) 

    Avg_bout_counts    = [];
    TwoWayANOVA_tables_bout_counts = [];

    return 
end 


% Sanity check: Do you have the right number of epochs for a 24-hour recording?  
for i=1:length(Scores.(Sex).WT.BL)
    if length(Scores.(Sex).WT.BL{i}) ~= (24*60*60)/epoch_duration
        error('In funcion Bout_Number_Analysis, a  WT BL recording has the wrong number of bouts.') 
    end 
end 
for i=1:length(Scores.(Sex).WT.SD)
    if length(Scores.(Sex).WT.SD{i}) ~= (24*60*60)/epoch_duration
        error('In funcion Bout_Number_Analysis, a  WT SD recording has the wrong number of bouts.') 
    end 
end 

% Build the struct to return the output
Avg_bout_counts = struct();
% ---- Baseline ------------

if WindowLengthHrs==12
    % --- WT -------------------------------
    % ------------------------------------------------------------------------------
    % ---- WT BL First12Hrs ----  
    [Avg_bout_counts.WT.BL.First12hrs, ~] = calculate_bout_counts_duration(Scores=Scores.(Sex).WT.BL,Sex=Sex,WTorMut='WT',BLorSD='BL',BoutMinimums=bout_mins,epoch_duration=epoch_duration,TimeFrame=[1 12]);

    % --- WT BL Last12Hrs   ----    
    [Avg_bout_counts.WT.BL.Last12hrs,~] = calculate_bout_counts_duration(Scores=Scores.(Sex).WT.BL,Sex=Sex,WTorMut='WT',BLorSD='BL',BoutMinimums=bout_mins,epoch_duration=epoch_duration,TimeFrame=[13 24]);


    % --- WT SD First12Hrs  ----    
    [Avg_bout_counts.WT.SD.First12hrs,~] = calculate_bout_counts_duration(Scores=Scores.(Sex).WT.SD,Sex=Sex,WTorMut='WT',BLorSD='SD',BoutMinimums=bout_mins,epoch_duration=epoch_duration,TimeFrame=[1 12]);

    % --- WT SD Last12Hrs   ----    
    [Avg_bout_counts.WT.SD.Last12hrs,~] = calculate_bout_counts_duration(Scores=Scores.(Sex).WT.SD,Sex=Sex,WTorMut='WT',BLorSD='SD',BoutMinimums=bout_mins,epoch_duration=epoch_duration,TimeFrame=[13 24]);

    % --- WT SD First12Hrs (excluding actual SD) ---    
    [Avg_bout_counts.WT.SDexcSD.First12hrs,~] = calculate_bout_counts_duration(Scores=Scores.(Sex).WT.SD,Sex=Sex,WTorMut='WT',BLorSD='SD',BoutMinimums=bout_mins,epoch_duration=epoch_duration,TimeFrame=[sleep_dep_length+1 12]);



    % --- Mutant -------------------------------
    % ------------------------------------------------------------------------------
    % ---- Mut BL First12Hrs ----  
    [Avg_bout_counts.Mut.BL.First12hrs,~] = calculate_avg_bout_counts(Scores=Scores.(Sex).Mut.BL,Sex=Sex,WTorMut='Mut',BLorSD='BL',BoutMinimums=bout_mins,epoch_duration=epoch_duration,TimeFrame=[1 12]);

    % --- Mut BL Last12Hrs ----    
    [Avg_bout_counts.Mut.BL.Last12hrs,~] = calculate_avg_bout_counts(Scores=Scores.(Sex).Mut.BL,Sex=Sex,WTorMut='Mut',BLorSD='BL',BoutMinimums=bout_mins,epoch_duration=epoch_duration,TimeFrame=[13 24]);

    
    % --- Mut SD First12Hrs ---    
    [Avg_bout_counts.Mut.SD.First12hrs,~] = calculate_avg_bout_counts(Scores=Scores.(Sex).Mut.SD,Sex=Sex,WTorMut='Mut',BLorSD='SD',BoutMinimums=bout_mins,epoch_duration=epoch_duration,TimeFrame=[1 12]);

    % --- Mut SD Last12Hrs ----    
    [Avg_bout_counts.Mut.SD.Last12hrs,~] = calculate_avg_bout_counts(Scores=Scores.(Sex).Mut.SD,Sex=Sex,WTorMut='Mut',BLorSD='SD',BoutMinimums=bout_mins,epoch_duration=epoch_duration,TimeFrame=[13 24]);

    % --- Mut SD (excluding the actual SD) First12Hrs ---    
    [Avg_bout_counts.Mut.SDexcSD.First12hrs,~] = calculate_avg_bout_counts(Scores=Scores.(Sex).Mut.SD,Sex=Sex,WTorMut='Mut',BLorSD='SD',BoutMinimums=bout_mins,epoch_duration=epoch_duration,TimeFrame=[sleep_dep_length+1 12]);
end 

% ----------------------------------------------------------------------
% ----------------------------------------------------------------------


%  -- If requested, analyze the dark period in 6-hr segments
if WindowLengthHrs==6
   % --- WT -------------------------------
    % ------------------------------------------------------------------------------
    % ---- WT BL First12Hrs ----  
    % Avg_bout_counts.WT.BL.First12hrs = calculate_avg_bout_counts(Scores=Scores.(Sex).WT.BL,Sex='Male',WTorMut='WT',BLorSD='BL',BoutMinimums=bout_mins,TimeFrame=[1 12]);

    % --- WT BL Dark Period First 6Hrs   ----    
    Avg_bout_counts.WT.BL.DPfirst6 = calculate_avg_bout_counts(Scores=Scores.(Sex).WT.BL,Sex=Sex,WTorMut='WT',BLorSD='BL',BoutMinimums=bout_mins,TimeFrame=[13 18]);

    % --- WT BL Dark Period Last 6Hrs   ----    
    Avg_bout_counts.WT.BL.DPlast6 = calculate_avg_bout_counts(Scores=Scores.(Sex).WT.BL,Sex=Sex,WTorMut='WT',BLorSD='BL',BoutMinimums=bout_mins,TimeFrame=[19 24]);


    % --- WT SD First12Hrs  ----    
    %Avg_bout_counts.WT.SD.First12hrs = calculate_avg_bout_counts(Scores=Scores.(Sex).WT.SD,Sex=Sex,WTorMut='WT',BLorSD='SD',BoutMinimums=bout_mins,TimeFrame=[1 12]);

    % --- WT SD Dark Period First 6Hrs   ----    
    Avg_bout_counts.WT.SD.DPfirst6 = calculate_avg_bout_counts(Scores=Scores.(Sex).WT.SD,Sex=Sex,WTorMut='WT',BLorSD='SD',BoutMinimums=bout_mins,TimeFrame=[13 18]);

    % --- WT SD Dark Period Last 6Hrs   ----    
    Avg_bout_counts.WT.SD.DPlast6 = calculate_avg_bout_counts(Scores=Scores.(Sex).WT.SD,Sex=Sex,WTorMut='WT',BLorSD='SD',BoutMinimums=bout_mins,TimeFrame=[19 24]);

    


    % --- Mutant -------------------------------
    % ------------------------------------------------------------------------------
    % ---- Mut BL First12Hrs ----  
    %Avg_bout_counts.Mut.BL.First12hrs = calculate_avg_bout_counts(Scores=Scores.(Sex).Mut.BL,Sex=Sex,WTorMut='WT',BLorSD='BL',BoutMinimums=bout_mins,TimeFrame=[1 12]);

    % --- Mut BL Dark Period First 6Hrs ----    
    Avg_bout_counts.Mut.BL.DPfirst6 = calculate_avg_bout_counts(Scores=Scores.(Sex).Mut.BL,Sex=Sex,WTorMut='Mut',BLorSD='BL',BoutMinimums=bout_mins,TimeFrame=[13 18]);

    % --- Mut BL Dark Period Last 6Hrs ----    
    Avg_bout_counts.Mut.BL.DPlast6 = calculate_avg_bout_counts(Scores=Scores.(Sex).Mut.BL,Sex=Sex,WTorMut='Mut',BLorSD='BL',BoutMinimums=bout_mins,TimeFrame=[19 24]);

    
    % --- Mut SD First12Hrs ---    
    %Avg_bout_counts.Mut.SD.First12hrs = calculate_avg_bout_counts(Scores=Scores.(Sex).Mut.SD,Sex=Sex,WTorMut='Mut',BLorSD='BL',BoutMinimums=bout_mins,TimeFrame=[1 12]);

    % --- Mut SD Dark Period First 6Hrs ----    
    Avg_bout_counts.Mut.SD.DPfirst6 = calculate_avg_bout_counts(Scores=Scores.(Sex).Mut.SD,Sex=Sex,WTorMut='Mut',BLorSD='SD',BoutMinimums=bout_mins,TimeFrame=[13 18]);

    % --- Mut SD Dark Period Last 6Hrs  ---    
    Avg_bout_counts.Mut.SD.DPlast6 = calculate_avg_bout_counts(Scores=Scores.(Sex).Mut.SD,Sex=Sex,WTorMut='Mut',BLorSD='SD',BoutMinimums=bout_mins,TimeFrame=[19 24]);
end 
% -- End of 6-hour window case --  




