function all_scores_table = build_scores_table(Scores_struct,Sexes)
%
% USAGE: all_scores_table = build_scores_table(Scores_struct)
%
% This function reads in the struct Scores_struct that was created by calls to Extract_EEG_data.m
% and changes the data from a multi-level struct into a table 
% with the following columns:  scores 	Sex		WTorMut		BLorSD
%
% INPUT: 		Scores_struct 		with the following structure: Scores_struct.Male.WT.SD or Scores_struct.Female.Mut.BL  
% 
% This way it will be easy to filter for the different combinations later on. 


% Initialize  the table
all_scores_table = table('Size',[1,4],...
    'VariableNames',{'Scores','Sex','WTorMut','BLorSD'},...
    'VariableTypes',{'cell','cellstr','cellstr','cellstr'});

% Male 
% Male WT BL
for i=1:length(Scores_struct.(Sexes{1}).WT.BL)
	scores_cell = {Scores_struct.(Sexes{1}).WT.BL{i},Sexes{1},'WT','BL'};
	all_scores_table = [all_scores_table; scores_cell];  
end 

% Male WT SD
for i=1:length(Scores_struct.(Sexes{1}).WT.SD)
	scores_cell = {Scores_struct.(Sexes{1}).WT.SD{i},Sexes{1},'WT','SD'};
	all_scores_table = [all_scores_table; scores_cell];  
end 

% Male Mut BL
for i=1:length(Scores_struct.(Sexes{1}).Mut.BL)
	scores_cell = {Scores_struct.(Sexes{1}).Mut.BL{i},Sexes{1},'Mut','BL'};
	all_scores_table = [all_scores_table; scores_cell];
end   


% Male Mut SD
for i=1:length(Scores_struct.(Sexes{1}).Mut.SD)
	scores_cell = {Scores_struct.(Sexes{1}).Mut.SD{i},Sexes{1},'Mut','SD'};
	all_scores_table = [all_scores_table; scores_cell];
end   



% -- Female --  
% Female WT BL
for i=1:length(Scores_struct.(Sexes{2}).WT.BL)
	scores_cell = {Scores_struct.(Sexes{2}).WT.BL{i},Sexes{2},'WT','BL'};
	all_scores_table = [all_scores_table; scores_cell];  
end 

% Female WT SD
for i=1:length(Scores_struct.(Sexes{2}).WT.SD)
	scores_cell = {Scores_struct.(Sexes{2}).WT.SD{i},Sexes{2},'WT','SD'};
	all_scores_table = [all_scores_table; scores_cell];  
end 

% Female Mut BL
for i=1:length(Scores_struct.(Sexes{2}).Mut.BL)
	scores_cell = {Scores_struct.(Sexes{2}).Mut.BL{i},Sexes{2},'Mut','BL'};
	all_scores_table = [all_scores_table; scores_cell];  
end 

% Female Mut SD
for i=1:length(Scores_struct.(Sexes{2}).Mut.SD)
	scores_cell = {Scores_struct.(Sexes{2}).Mut.SD{i},Sexes{2},'Mut','SD'};
	all_scores_table = [all_scores_table; scores_cell];  
end 

all_scores_table(1,:) = [];  % remove first row since it is empty
all_scores_table.BLorSD  = categorical(all_scores_table.BLorSD);
all_scores_table.Sex     = categorical(all_scores_table.Sex);
all_scores_table.WTorMut = categorical(all_scores_table.WTorMut);











