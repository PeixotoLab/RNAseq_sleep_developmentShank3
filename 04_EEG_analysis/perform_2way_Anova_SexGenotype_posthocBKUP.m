% [p_vals] = perform_2way_Anova_SexGenotype_posthoc(Female_WT,Female_Mut,Male_WT,Male_Mut)
% 
% This function reads in the data in structs like Twelve_hour_avg_percentages_Female.WT.BL.First12hrs.Wake, etc
% and runs a 2-way ANOVA for sex and genotype.  If any effects are significant, it follows up with a post-hoc 
% test 
% This can be used to add asterisks or hashtags to plots 





% create a vector of all the concatenated output data

WakePercentage = [Twelve_hour_avg_percentages_Female.WT.BL.First12hrs.Wake';...
				  Twelve_hour_avg_percentages_Female.Mut.BL.First12hrs.Wake';...
				  Twelve_hour_avg_percentages_Male.WT.BL.First12hrs.Wake';...
				  Twelve_hour_avg_percentages_Male.Mut.BL.First12hrs.Wake'];

sex = [repmat({'Female'},length(Twelve_hour_avg_percentages_Female.WT.BL.First12hrs.Wake),1); ...
	   repmat({'Female'},length(Twelve_hour_avg_percentages_Female.Mut.BL.First12hrs.Wake),1); ...
	   repmat({'Male'},length(Twelve_hour_avg_percentages_Male.WT.BL.First12hrs.Wake),1); ...
	   repmat({'Male'},length(Twelve_hour_avg_percentages_Male.Mut.BL.First12hrs.Wake),1)];


% Set up Genotype vector to be 'WT','WT','WT'... 'Mut','Mut','Mut'
Genotype = [repmat({'WT'},length(Twelve_hour_avg_percentages_Female.WT.BL.First12hrs.Wake),1); ...
	   repmat({'Mut'},length(Twelve_hour_avg_percentages_Female.Mut.BL.First12hrs.Wake),1); ...
	   repmat({'WT'},length(Twelve_hour_avg_percentages_Male.WT.BL.First12hrs.Wake),1); ...
	   repmat({'Mut'},length(Twelve_hour_avg_percentages_Male.Mut.BL.First12hrs.Wake),1)];
sex = categorical(sex);
Genotype = categorical(Genotype);




[p,tbl,stats,terms] = anovan(WakePercentage,{sex Genotype},'model',2,'varnames',{'Sex','Genotype'})

if p(1) < 0.05  % if sex is a significant factor, do t-tests between sexes for each genotype
	[h,p_sex_WT]=ttest2(Twelve_hour_avg_percentages_Female.WT.BL.First12hrs.Wake,Twelve_hour_avg_percentages_Male.WT.BL.First12hrs.Wake);
	[h,p_sex_Mut]=ttest2(Twelve_hour_avg_percentages_Female.Mut.BL.First12hrs.Wake,Twelve_hour_avg_percentages_Male.Mut.BL.First12hrs.Wake);
	end 
alpha_sid = 1-(1-0.05)^(1/2);
if p_sex_WT < alpha_sid
	disp('Sig difference between sexes WT')
end 

if p_sex_Mut < alpha_sid
	disp('Sig difference between sexes Mut')
end 

if p(2) < 0.05 % if genotype is a significant factor, do t-tests between genotypes for each sex
	[h,p_geno_Male]=ttest2(Twelve_hour_avg_percentages_Male.WT.BL.First12hrs.Wake,Twelve_hour_avg_percentages_Male.Mut.BL.First12hrs.Wake);
	[h,p_geno_Female]=ttest2(Twelve_hour_avg_percentages_Female.WT.BL.First12hrs.Wake,Twelve_hour_avg_percentages_Female.Mut.BL.First12hrs.Wake);
	end 

p_vals.Male.WTvsMut     = p_geno_Male;
p_vals.Female.WTvsMut   = p_geno_Female;
p_vals.WT.MalevsFemale  = p_sex_WT;
p_vals.Mut.MalevsFemale = p_sex_Mut;



