function Tnew = change_ANOVA_table_names(T,SexVarName)
%
% usage:  Tnew = change_ANOVA_table_names(T,SexVarName)
% 
% This function changes the row lables of an ANOVA table to change 
% the variable 'Sex' to whatever you want it to be (ie 'Age' if you are treating age like sex)



% now find the text 'Sex' in each of these and replace it with SexVarName
for i = 1:numel(T.Properties.RowNames)
    if ischar(T.Properties.RowNames{i}) || isstring(T.Properties.RowNames{i})
        T.Properties.RowNames{i} = strrep(T.Properties.RowNames{i}, 'Sex', SexVarName);
    end
end

Tnew = T;  