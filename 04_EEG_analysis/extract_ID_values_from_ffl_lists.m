function IDs = extract_ID_values_from_ffl_lists(ffl_list,delimiter)
%
% this little function reads in the file list and finds the IDs by
% looking at the text between underscores.  That text between 
% underscores that is different between files is the ID.  Every
% other text between underscores is the same.  

for i=1:length(ffl_list)
	commonparts(i,:) = strsplit(ffl_list{i},delimiter);
end 
for j=1:size(commonparts,2) % loop over columns to find same
	if length(unique(commonparts(:,j))) > 1
		IDs = commonparts(:,j);  	% cell array of ID values
	end
end


