function [p_vals,sig_difference,ANOVA_table] = perform_2way_Anova_SexGenotype_posthoc(Female_WT,Female_Mut,Male_WT,Male_Mut,Sexes,SexVarName)
% 
% This function reads in the data in structs like Twelve_hour_avg_percentages_Female.WT.BL.First12hrs.Wake, etc
% and runs a 2-way ANOVA for sex and genotype.  If any effects are significant, it follows up with a post-hoc 
% test 
% This can be used to add asterisks or hashtags to plots 



% --- Male data ------------------------------
% --- Perform 1-way ANOVA (or a simple t-test)
if isempty(Female_WT) & isempty(Female_Mut)
	Output_variable = [Male_WT';Male_Mut'];
	Genotype = [repmat({'WT'},length(Male_WT),1); ...
	   			repmat({'Mut'},length(Male_Mut),1)];

	% I know I don't need to use anovan for this, but it was easier to copy from below
	[p,tbl,stats,terms] = anovan(Output_variable,{Genotype},'varnames','Genotype','display','off');

	[h,p_ttest_Male_WTVsMut] = ttest2(Male_WT,Male_Mut);  % 1-way ANOVA is the same as a t-test

	% Sanity check: 1-way ANOVA and ttest2 should give you the same p-value
	if abs(p-p_ttest_Male_WTVsMut)>1e-5 
		error('In perform_2way_Anova_SexGenotype_posthoc: one-way ANOVA and ttest2 gave different p-values.')
	end 


	p_vals.(Sexes{1}).ANOVA.(SexVarName) = zeros(0);
	p_vals.(Sexes{1}).ANOVA.Genotype = p(1);
	p_vals.(Sexes{1}).ANOVA.(['GenotypeX' SexVarName]) = zeros(0);

	p_vals.Posthoc.(Sexes{1}).WTvsMut = p(1);
	p_vals.Posthoc.(Sexes{2}).WTvsMut = zeros(0);
	p_vals.Posthoc.WT.([Sexes{1} 'vs' Sexes{2}])  = zeros(0);
	p_vals.Posthoc.Mut.([Sexes{1} 'vs' Sexes{2}]) = zeros(0);

	sig_difference.(Sexes{1}).WTvsMut = p_ttest_Male_WTVsMut     < 0.05;
	sig_difference.(Sexes{2}).WTvsMut = [];
	sig_difference.WT.([Sexes{1} 'vs' Sexes{2}])  = [];
	sig_difference.Mut.([Sexes{1} 'vs' Sexes{2}]) = [];

	ANOVA_table.(Sexes{1}).Genotype = tbl;

	p_vals.(Sexes{2}).ANOVA.Genotype = zeros(0);
	p_vals.(['Both' Sexes{1} 'and' Sexes{2}]).ANOVA.Genotype = zeros(0); 
	ANOVA_table.(Sexes{2}).Genotype 	= []; 

% --- The Case Where Both Sexes are present -----------------------------------
% create a vector of all the concatenated output data
elseif ~isempty(Female_WT) & ~isempty(Female_Mut)
	Output_variable = [Female_WT';...
				   Female_Mut';...
				   Male_WT';...
				   Male_Mut'];

	sex = [repmat({'Female'},length(Female_WT),1); ...
	   repmat({'Female'},length(Female_Mut),1); ...
	   repmat({'Male'},length(Male_WT),1); ...
	   repmat({'Male'},length(Male_Mut),1)];


	% Set up Genotype vector to be 'WT','WT','WT'... 'Mut','Mut','Mut'
	Genotype = [repmat({'WT'},length(Female_WT),1); ...
		repmat({'Mut'},length(Female_Mut),1); ...
		repmat({'WT'},length(Male_WT),1); ...
		repmat({'Mut'},length(Male_Mut),1)];

	sex = categorical(sex);
	Genotype = categorical(Genotype);




	[p,tbl,stats,terms] = anovan(Output_variable,{sex Genotype},'model',2,'varnames',{'Sex','Genotype'},'display','off');

	% Label the just-created ANOVA table with the arousal state and perhaps the metric (time in state, etc)
	% TODO
	% headingObj = findall(0,'Type','uicontrol','Tag','Heading');
	% headingObj(1).String = 'New String for ANOVA';    			% (1) means the most recent ANOVA table

	% adjust for multiple comparisons using BH correction
	%p = mafdr(p,'BHFDR','true');    % q contains the adjusted p values


	% Set up the post-hoc values as empty
	p_ttest_WT_MaleVsFemale  = zeros(0);
	p_ttest_Mut_MaleVsFemale = zeros(0);
	p_ttest_Male_WTVsMut     = zeros(0);
	p_ttest_Female_WTVsMut   = zeros(0);

	if p(1) < 0.05   % if sex is a significant factor, do t-tests between sexes for each genotype
		[h,p_ttest_WT_MaleVsFemale]  = ttest2(Female_WT,Male_WT);
		[h,p_ttest_Mut_MaleVsFemale] = ttest2(Female_Mut,Male_Mut);
	end 

	% alpha_sid = 1-(1-0.05)^(1/2);  % this is the Sidak (or Dunn-Sidak) correction to an alpha of 0.05

	% if p_ttest_WT_MaleVsFemale < alpha_sid
	% 	disp('Sig difference between sexes WT')
	% end 

	% if p_ttest_Mut_MaleVsFemale < alpha_sid
	% 	disp('Sig difference between sexes Mut')
	% end 

	if p(2) < 0.05 % if genotype is a significant factor, do t-tests between genotypes for each sex
		[h,p_ttest_Male_WTVsMut]   = ttest2(Male_WT,Male_Mut);
		[h,p_ttest_Female_WTVsMut] = ttest2(Female_WT,Female_Mut);
	end 

	if p(3) < 0.05 % if the interaction term is significant, do t-tests between genotypes for each sex and between sex for each genotype
		[h,p_ttest_WT_MaleVsFemale]  = ttest2(Female_WT,Male_WT);
		[h,p_ttest_Mut_MaleVsFemale] = ttest2(Female_Mut,Male_Mut);
		[h,p_ttest_Male_WTVsMut]     = ttest2(Male_WT,Male_Mut);
		[h,p_ttest_Female_WTVsMut]   = ttest2(Female_WT,Female_Mut);
	end 

	% Adjust the p-values after doing the post-hocs, based on how many tests you did


	% Separately, do an ANOVA just on Female data (with Genotype as a factor), just like for Males
	% Output_variableF = [Female_WT';Female_Mut'];
	% GenotypeF = [repmat({'WT'},length(Female_WT),1); ...
	%    			repmat({'Mut'},length(Female_Mut),1)];

	% % I know I don't need to use anovan for this, but it was easier to copy from below
	% [pF,tblF,stats,terms] = anovan(Output_variableF,{GenotypeF},'varnames','Genotype','display','off');

	% [h,p_ttest_Female_WTVsMut] = ttest2(Female_WT,Female_Mut);  % 1-way ANOVA is the same as a t-test

	% % Sanity check: 1-way ANOVA and ttest2 should give you the same p-value
	% if abs(pF-p_ttest_Female_WTVsMut)>1e-5 
	% 	error('In perform_2way_Anova_SexGenotype_posthoc: one-way ANOVA and ttest2 for female data gave different p-values.')
	% end 

	% ANOVA_table.(Sexes{2}).Genotype = tblF;

	% p_vals.(Sexes{2}).ANOVA.Sex 	 	 = [];
	% p_vals.(Sexes{2}).ANOVA.Genotype 	 = pF(1);
	% p_vals.(Sexes{2}).ANOVA.GenotypeXSex = [];




	p_vals.(['Both' Sexes{1} 'and' Sexes{2}]).ANOVA.(SexVarName) = p(1);
	p_vals.(['Both' Sexes{1} 'and' Sexes{2}]).ANOVA.Genotype     = p(2);
	p_vals.(['Both' Sexes{1} 'and' Sexes{2}]).ANOVA.(['GenotypeX' SexVarName]) = p(3);


	% set up a stencil for which p_ttest values are non-empty
	% T.nonempty = isempty([p_ttest_WT_MaleVsFemale p_ttest_Mut_MaleVsFemale ...
	% 	p_ttest_Male_WTVsMut p_ttest_Female_WTVsMut]);
	% concatenate all p_ttest values into one vector
	T = table;
	T.pvals = cell(4,1);
	T.pvals{1} = p_ttest_WT_MaleVsFemale;
	T.pvals{2} = p_ttest_Mut_MaleVsFemale; 
	T.pvals{3} = p_ttest_Male_WTVsMut;
	T.pvals{4} = p_ttest_Female_WTVsMut;
	T.name{1}  = 'WT_MaleVsFemale';
	T.name{2}  = 'Mut_MaleVsFemale';
	T.name{3}  = 'Male_WTVsMut';
	T.name{4}  = 'Female_WTVsMut';
	T.name = categorical(T.name);
	T.isempty = cellfun('isempty',T.pvals);

	% run the vector through mafdr
	if sum(T.isempty)<4   % if at least one of the post-hocs were computed, adjust with mafdr
		adj_pvals = mafdr(cell2mat(T.pvals),'BHFDR','true');
		T.pvals(T.isempty==0) = num2cell(adj_pvals);
	end 

	% return both the unadjusted and adjusted p-values
	% Set up structures first
	% Male WT vs Mut
	M_WTvsMut_cell = cell(2);
	M_WTvsMut_cell{1,1} = 'Unadjusted';
	M_WTvsMut_cell{2,1} = p_ttest_Male_WTVsMut;  
	M_WTvsMut_cell{1,2} = 'Adjusted';
	M_WTvsMut_cell{2,2} = T.pvals{T.name=='Male_WTVsMut'};

	% Female WT vs Mut
	F_WTvsMut_cell = cell(2);
	F_WTvsMut_cell{1,1} = 'Unadjusted';
	F_WTvsMut_cell{2,1} = p_ttest_Female_WTVsMut;  
	F_WTvsMut_cell{1,2} = 'Adjusted';
	F_WTvsMut_cell{2,2} = T.pvals{T.name=='Female_WTVsMut'};

	% WT Male vs Female
	WT_MvsF_cell = cell(2);
	WT_MvsF_cell{1,1} = 'Unadjusted';
	WT_MvsF_cell{2,1} = p_ttest_WT_MaleVsFemale;
	WT_MvsF_cell{1,2} = 'Adjusted';
	WT_MvsF_cell{2,2} = T.pvals{T.name=='WT_MaleVsFemale'};

	% Mut Male vs Female
	Mut_MvsF_cell = cell(2);
	Mut_MvsF_cell{1,1} = 'Unadjusted';
	Mut_MvsF_cell{2,1} = p_ttest_Mut_MaleVsFemale;
	Mut_MvsF_cell{1,2} = 'Adjusted';
	Mut_MvsF_cell{2,2} = T.pvals{T.name=='Mut_MaleVsFemale'};

	p_vals.Posthoc.(Sexes{1}).WTvsMut     = cell2table(M_WTvsMut_cell(2,:),'VariableNames',M_WTvsMut_cell(1,:));
	p_vals.Posthoc.(Sexes{2}).WTvsMut   = cell2table(F_WTvsMut_cell(2,:),'VariableNames',F_WTvsMut_cell(1,:));
	p_vals.Posthoc.WT.([Sexes{1} 'vs' Sexes{2}])  = cell2table(WT_MvsF_cell(2,:),'VariableNames',WT_MvsF_cell(1,:));
	p_vals.Posthoc.Mut.([Sexes{1} 'vs' Sexes{2}]) = cell2table(Mut_MvsF_cell(2,:),'VariableNames',Mut_MvsF_cell(1,:));

	% if either value is empty, make it NaN not a 0x0 cell because of if statements later. Matlab complains
	if iscell(p_vals.Posthoc.(Sexes{1}).WTvsMut.Unadjusted) 
		p_vals.Posthoc.(Sexes{1}).WTvsMut.Unadjusted = NaN;
	end 
	if iscell(p_vals.Posthoc.(Sexes{1}).WTvsMut.Adjusted) 
		p_vals.Posthoc.(Sexes{1}).WTvsMut.Adjusted = NaN;
	end 
	if iscell(p_vals.Posthoc.(Sexes{2}).WTvsMut.Unadjusted) 
		p_vals.Posthoc.(Sexes{2}).WTvsMut.Unadjusted = NaN;
	end 
	if iscell(p_vals.Posthoc.(Sexes{2}).WTvsMut.Adjusted) 
		p_vals.Posthoc.(Sexes{2}).WTvsMut.Adjusted = NaN;
	end 
	if iscell(p_vals.Posthoc.WT.([Sexes{1} 'vs' Sexes{2}]).Unadjusted)
		p_vals.Posthoc.WT.([Sexes{1} 'vs' Sexes{2}]).Unadjusted = NaN;
	end
	if iscell(p_vals.Posthoc.WT.([Sexes{1} 'vs' Sexes{2}]).Adjusted)
		p_vals.Posthoc.WT.([Sexes{1} 'vs' Sexes{2}]).Adjusted = NaN;
	end
	if iscell(p_vals.Posthoc.Mut.([Sexes{1} 'vs' Sexes{2}]).Unadjusted)
		p_vals.Posthoc.Mut.([Sexes{1} 'vs' Sexes{2}]).Unadjusted = NaN;
	end
	if iscell(p_vals.Posthoc.Mut.([Sexes{1} 'vs' Sexes{2}]).Adjusted)
		p_vals.Posthoc.Mut.([Sexes{1} 'vs' Sexes{2}]).Adjusted = NaN;
	end

	% p_vals.Posthoc.(Sexes{1}).WTvsMut     = T.pvals{T.name=='Male_WTVsMut'};
	% p_vals.Posthoc.(Sexes{2}).WTvsMut   = T.pvals{T.name=='Female_WTVsMut'};
	% p_vals.Posthoc.WT.([Sexes{1} 'vs' Sexes{2}])  = T.pvals{T.name=='WT_MaleVsFemale'};
	% p_vals.Posthoc.Mut.([Sexes{1} 'vs' Sexes{2}]) = T.pvals{T.name=='Mut_MaleVsFemale'};

	% p_vals.Posthoc.(Sexes{1}).WTvsMut     = p_ttest_Male_WTVsMut;
	% p_vals.Posthoc.(Sexes{2}).WTvsMut   = p_ttest_Female_WTVsMut;
	% p_vals.Posthoc.WT.([Sexes{1} 'vs' Sexes{2}])  = p_ttest_WT_MaleVsFemale;
	% p_vals.Posthoc.Mut.([Sexes{1} 'vs' Sexes{2}]) = p_ttest_Mut_MaleVsFemale;

	sig_difference.(Sexes{1}).WTvsMut     = p_ttest_Male_WTVsMut     < 0.05;
	sig_difference.(Sexes{2}).WTvsMut   = p_ttest_Female_WTVsMut   < 0.05;
	sig_difference.WT.([Sexes{1} 'vs' Sexes{2}])  = p_ttest_WT_MaleVsFemale  < 0.05;
	sig_difference.Mut.([Sexes{1} 'vs' Sexes{2}]) = p_ttest_Mut_MaleVsFemale < 0.05;

	ANOVA_table.([SexVarName 'Genotype']) = tbl;

	% Fix the ANOVA table to use the correct variable name for Sex (Age, etc)
	if(~strcmp(SexVarName,'Sex'))
		ANOVA_table.([SexVarName 'Genotype'])(:,1) = strrep(ANOVA_table.([SexVarName 'Genotype'])(:,1),'Sex',SexVarName);		
	end


else 
	error('Something went wrong in perform_2way_Anova_SexGenotype_posthoc.m: You did not have Male+Female data or only Male.')
end 
