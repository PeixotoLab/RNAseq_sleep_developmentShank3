% ConvertSleepSignCSVtoMatlabMAT.m
% 
% This script reads in a directory where your csv files live and then reads in each 
% .csv file and converts it to a MATLAB table and saves it as EEG in a matlab .mat file. 
%
% CHANGE the lines that begin with path_csv and path_mat to be what you need. 
%
%
% You will need to call this script for each of the following directories where you have data:
%
% Baseline WT Male
% Baseline Mut Male
% SD WT Male
% SD Mut Male
% Baseline WT Female
% Baseline Mut Female
% SD WT Female
% SD Mut Female 

tic;
%%
close all;
clc;
clear all;


%% path
%path_csv = uigetdir([], 'Please choose a directory containing the .csv files you wish to convert.')
path_csv ='/Users/caitlinottaway/Library/CloudStorage/OneDrive-WashingtonStateUniversity(email.wsu.edu)/Pexioto/EEG_Analysis/P24_P30_P45_P60/P45_CSV/todo'
%path_csv = '/Users/caitlinottaway/Library/CloudStorage/OneDrive-WashingtonStateUniversity(email.wsu.edu)/Pexioto/EEG_Analysis/P24_P30_Reanalysis/P30_CSVs/todo/';  % directory where the .csv files live
path_mat = '/Users/caitlinottaway/Library/CloudStorage/OneDrive-WashingtonStateUniversity(email.wsu.edu)/Pexioto/EEG_Analysis/P24_P30_P45_P60/P45_MAT';  % directory where you want to save the .mat files (it can be the same!)

myfiles = dir(strcat(path_csv,'/*.csv'))
file = {myfiles.name};

numfile = length(file);


%% Read in each file and convert to a .mat file with the _FFT.mat format
for i= 1:numfile
    clear filename
    filename=fullfile([path_csv+"/"+file{i}])
    data = mylocalread(filename);
    data_varnames = data{1}(22,1:end-1);
    data=data{1,1}(23:end,1:end-1);
    %% convert from cell array to table, give column names based on variable

    EEG = cell2table(data,'VariableNames',data_varnames);
    %% extract the file name to make saving easier

    fname=extractBefore(file{i},'FFT');
    filename=fullfile([path_mat+"/"+fname+'_FFT'+'.mat']) % for the long files
    save(filename,'EEG')
    clear EEG data_varnames data filename
end
toc;




% ------------------------------------------------------------------------------------------------------------
% ------------------------------------------------------------------------------------------------------------
% ------------------------------------------------------------------------------------------------------------
%% local functions, DOH

% this is my new go-to for general purpose reading of text files. it's a little slow, but it gets the job done
% it puts all data separated by certain delimeters from the text files into an array of cell matrixes
%%
function out = mylocalread(f)
    n=length(f);
    out=cell(n,1);
    for i=1:n
        fid=fopen(f{i}); %open the file
        j=1;
        while ~feof(fid) %keep reading lines until the file end of file (feof)
            line=fgetl(fid); %get next line
            line=strsplit(line,{', ',',','\t'}); % edit these delimeters when necessary ,'_'
            if ~iscell(line), line={line}; end %make sure its a cell
            m=length(line);
            for k=1:m, line(k)=mylocalfunc(line{k}); end %found out that this is faster than cellfun
            out{i}(j,1:m)=line; 
            j=j+1;
        end
        fclose(fid);
%         out{i}=cellfun(@mylocalfunc,out{i});
    end
end

% determine if data from is numeric. if so, convert to numeric
function B = mylocalfunc(A)
    if all(ismember(A, '0123456789+-.eEdD')) && ~isempty(A)
        B={str2double(A)};
    elseif strcmp(A,'None')
        B={NaN};
    else
        B={A};
    end
end

% ------------------------------------------------------------------------------------------------------------
% ------------------------------------------------------------------------------------------------------------
%% End of local functions, DOH

