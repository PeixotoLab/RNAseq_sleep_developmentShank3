function [RM_ANOVA_tables,PostHocs_Recovery] = hourly_NREM_Delta_or_Sigma_afterSD(options)
%
% This function reads in Hourly values of either Delta or sigma during either NREM or 
% 
% NOTE:  This assumes the filenames have underscores.  This is how I find the ID values


arguments 
	options.Hourly_power_inSTATE
	options.DeltaOrTheta
	options.WorNREM
	options.ffl_list_WT_BL_Males
	options.ffl_list_Mut_BL_Males
	options.ffl_list_WT_BL_Females
	options.ffl_list_Mut_BL_Females
	options.ffl_list_WT_SD_Males
	options.ffl_list_Mut_SD_Males
	options.ffl_list_WT_SD_Females
	options.ffl_list_Mut_SD_Females
	options.Plot_BL_too = 0
	options.BinSizeHrs  = 1
	options.SD_length_hrs
	options.Sexes
end 

Hourly_power_inSTATE    = options.Hourly_power_inSTATE;
DeltaOrTheta            = options.DeltaOrTheta;
Sexes   			    = options.Sexes;
WorNREM                 = options.WorNREM; 
ffl_list_WT_BL_Males    = options.ffl_list_WT_BL_Males;
ffl_list_Mut_BL_Males   = options.ffl_list_Mut_BL_Males;
ffl_list_WT_BL_Females  = options.ffl_list_WT_BL_Females;
ffl_list_Mut_BL_Females = options.ffl_list_Mut_BL_Females;
ffl_list_WT_SD_Males    = options.ffl_list_WT_SD_Males;
ffl_list_Mut_SD_Males   = options.ffl_list_Mut_SD_Males;
ffl_list_WT_SD_Females  = options.ffl_list_WT_SD_Females;
ffl_list_Mut_SD_Females = options.ffl_list_Mut_SD_Females;
Plot_BL_too 			= options.Plot_BL_too;
BinSizeHrs     			= options.BinSizeHrs; 
SD_length_hrs 			= options.SD_length_hrs;


% Is there female data present? 
if isempty(Hourly_power_inSTATE.WT.(Sexes{2}).BL)
	FemaleDataPresent = 0;
else 
	FemaleDataPresent = 1;
end 

% --- First, determine which recordings are present in BL and SD ---
delimiter = '_';

% -- Males --
IDs_WT_BL_Males  = extract_ID_values_from_ffl_lists(ffl_list_WT_BL_Males,delimiter);
IDs_WT_SD_Males  = extract_ID_values_from_ffl_lists(ffl_list_WT_SD_Males,delimiter);
IDs_Mut_BL_Males = extract_ID_values_from_ffl_lists(ffl_list_Mut_BL_Males,delimiter);
IDs_Mut_SD_Males = extract_ID_values_from_ffl_lists(ffl_list_Mut_SD_Males,delimiter);

if FemaleDataPresent
	% -- Females --
	IDs_WT_BL_Females  = extract_ID_values_from_ffl_lists(ffl_list_WT_BL_Females,delimiter);
	IDs_WT_SD_Females  = extract_ID_values_from_ffl_lists(ffl_list_WT_SD_Females,delimiter);
	IDs_Mut_BL_Females = extract_ID_values_from_ffl_lists(ffl_list_Mut_BL_Females,delimiter);
	IDs_Mut_SD_Females = extract_ID_values_from_ffl_lists(ffl_list_Mut_SD_Females,delimiter);
end 

% -- Get list of IDs present in both BL and SD
IDs_BothBLandSD_WT_Male  = intersect(IDs_WT_BL_Males,IDs_WT_SD_Males);
IDs_BothBLandSD_Mut_Male = intersect(IDs_Mut_BL_Males,IDs_Mut_SD_Males);

if FemaleDataPresent
 	IDs_BothBLandSD_WT_Female  = intersect(IDs_WT_BL_Females,IDs_WT_SD_Females);
	IDs_BothBLandSD_Mut_Female = intersect(IDs_Mut_BL_Females,IDs_Mut_SD_Females);
end 


% Normalize to the last 4 hours of the LP in Baseline (hours 9-12)?
ind_last4hrs_BL = 9:12;
%ind_last4hrs_BL = 1:24


% Male WT
for i=1:length(IDs_BothBLandSD_WT_Male)
	indBL = find(strcmp(IDs_WT_BL_Males, IDs_BothBLandSD_WT_Male{i}));
	indSD = find(strcmp(IDs_WT_SD_Males, IDs_BothBLandSD_WT_Male{i}));
	power_inSTATE_WT_Male_normalization(i) = mean(Hourly_power_inSTATE.WT.(Sexes{1}).BL(indBL,ind_last4hrs_BL),'omitnan');
	Normalized_power.SD.WT.Male(i,:) = (Hourly_power_inSTATE.WT.(Sexes{1}).SD(indSD,:)./power_inSTATE_WT_Male_normalization(i))*100;
	Normalized_power.BL.WT.Male(i,:) = (Hourly_power_inSTATE.WT.(Sexes{1}).BL(indBL,:)./power_inSTATE_WT_Male_normalization(i))*100;
end 

% Male Mut
for i=1:length(IDs_BothBLandSD_Mut_Male)
	indBL = find(strcmp(IDs_Mut_BL_Males, IDs_BothBLandSD_Mut_Male{i}));
	indSD = find(strcmp(IDs_Mut_SD_Males, IDs_BothBLandSD_Mut_Male{i}));
	power_inSTATE_Mut_Male_normalization(i) = mean(Hourly_power_inSTATE.Mut.(Sexes{1}).BL(indBL,ind_last4hrs_BL),'omitnan');
	Normalized_power.SD.Mut.Male(i,:) = (Hourly_power_inSTATE.Mut.(Sexes{1}).SD(indSD,:)./power_inSTATE_Mut_Male_normalization(i))*100;
	Normalized_power.BL.Mut.Male(i,:) = (Hourly_power_inSTATE.Mut.(Sexes{1}).BL(indBL,:)./power_inSTATE_Mut_Male_normalization(i))*100;
end 

if FemaleDataPresent
	% Female WT
	for i=1:length(IDs_BothBLandSD_WT_Female)
		indBL = find(strcmp(IDs_WT_BL_Females, IDs_BothBLandSD_WT_Female{i}));
		indSD = find(strcmp(IDs_WT_SD_Females, IDs_BothBLandSD_WT_Female{i}));
		power_inSTATE_WT_Female_normalization(i) = mean(Hourly_power_inSTATE.WT.(Sexes{2}).BL(indBL,ind_last4hrs_BL),'omitnan');
		Normalized_power.SD.WT.Female(i,:) = (Hourly_power_inSTATE.WT.(Sexes{2}).SD(indSD,:)./power_inSTATE_WT_Female_normalization(i))*100;
		Normalized_power.BL.WT.Female(i,:) = (Hourly_power_inSTATE.WT.(Sexes{2}).BL(indBL,:)./power_inSTATE_WT_Female_normalization(i))*100;
	end 

	% Female Mut
	for i=1:length(IDs_BothBLandSD_Mut_Female)
		indBL = find(strcmp(IDs_Mut_BL_Females, IDs_BothBLandSD_Mut_Female{i}));
		indSD = find(strcmp(IDs_Mut_SD_Females, IDs_BothBLandSD_Mut_Female{i}));
		power_inSTATE_Mut_Female_normalization(i) = mean(Hourly_power_inSTATE.Mut.(Sexes{2}).BL(indBL,ind_last4hrs_BL),'omitnan');
		Normalized_power.SD.Mut.Female(i,:) = (Hourly_power_inSTATE.Mut.(Sexes{2}).SD(indSD,:)./power_inSTATE_Mut_Female_normalization(i))*100;
		Normalized_power.BL.Mut.Female(i,:) = (Hourly_power_inSTATE.Mut.(Sexes{2}).BL(indBL,:)./power_inSTATE_Mut_Female_normalization(i))*100;
	end 
else 
	Normalized_power.SD.WT.Female  = [];
	Normalized_power.SD.Mut.Female = [];
	Normalized_power.BL.WT.Female  = [];
	Normalized_power.BL.Mut.Female = [];
end 

% power_inSTATE_WT_Male_normalization    = mean(Hourly_power_inSTATE.WT.Male.BL ,2,'omitnan');
% power_inSTATE_Mut_Male_normalization   = mean(Hourly_power_inSTATE.Mut.Male.BL,2,'omitnan');
% power_inSTATE_WT_Female_normalization  = mean(Hourly_power_inSTATE.WT.Female.BL,2,'omitnan');
% power_inSTATE_Mut_Female_normalization = mean(Hourly_power_inSTATE.Mut.Female.BL,2,'omitnan');  
	
% Normalized_power.WT.Male    = (Hourly_power_inSTATE.WT.Male.SD./power_inSTATE_WT_Male_normalization)*100;
% Normalized_power.Mut.Male   = (Hourly_power_inSTATE.Mut.Male.SD./power_inSTATE_Mut_Male_normalization)*100;
% Normalized_power.WT.Female  = (Hourly_power_inSTATE.WT.Female.SD./power_inSTATE_WT_Female_normalization)*100;
% Normalized_power.Mut.Female = (Hourly_power_inSTATE.Mut.Female.SD./power_inSTATE_Mut_Female_normalization)*100;  

% Next compute the means and SEs to make the plots  
% - SD -
Means_norm_power.SD.WT.Male    = mean(Normalized_power.SD.WT.Male,1,'omitnan');
Means_norm_power.SD.Mut.Male   = mean(Normalized_power.SD.Mut.Male,1,'omitnan');
Means_norm_power.SD.WT.Female  = mean(Normalized_power.SD.WT.Female,1,'omitnan');
Means_norm_power.SD.Mut.Female = mean(Normalized_power.SD.Mut.Female,1,'omitnan');

% - BL -
Means_norm_power.BL.WT.Male    = mean(Normalized_power.BL.WT.Male,1,'omitnan');
Means_norm_power.BL.Mut.Male   = mean(Normalized_power.BL.Mut.Male,1,'omitnan');
Means_norm_power.BL.WT.Female  = mean(Normalized_power.BL.WT.Female,1,'omitnan');
Means_norm_power.BL.Mut.Female = mean(Normalized_power.BL.Mut.Female,1,'omitnan');

% - SD -
SEMs_norm_power.SD.WT.Male     = std(Normalized_power.SD.WT.Male,0,1,'omitnan')./sqrt(sum(~isnan(Normalized_power.SD.WT.Male),1));
SEMs_norm_power.SD.Mut.Male    = std(Normalized_power.SD.Mut.Male,0,1,'omitnan')./sqrt(sum(~isnan(Normalized_power.SD.Mut.Male),1)); 
SEMs_norm_power.SD.WT.Female   = std(Normalized_power.SD.WT.Female,0,1,'omitnan')./sqrt(sum(~isnan(Normalized_power.SD.WT.Female),1));
SEMs_norm_power.SD.Mut.Female  = std(Normalized_power.SD.Mut.Female,0,1,'omitnan')./sqrt(sum(~isnan(Normalized_power.SD.Mut.Female),1));

% - BL -
SEMs_norm_power.BL.WT.Male     = std(Normalized_power.BL.WT.Male,0,1,'omitnan')./sqrt(sum(~isnan(Normalized_power.BL.WT.Male),1));
SEMs_norm_power.BL.Mut.Male    = std(Normalized_power.BL.Mut.Male,0,1,'omitnan')./sqrt(sum(~isnan(Normalized_power.BL.Mut.Male),1)); 
SEMs_norm_power.BL.WT.Female   = std(Normalized_power.BL.WT.Female,0,1,'omitnan')./sqrt(sum(~isnan(Normalized_power.BL.WT.Female),1));
SEMs_norm_power.BL.Mut.Female  = std(Normalized_power.BL.Mut.Female,0,1,'omitnan')./sqrt(sum(~isnan(Normalized_power.BL.Mut.Female),1));


% --- Carry out the repeated measures ANOVA to compare the two curves -------------
% ---------------------------------------------------------------------------------
if BinSizeHrs == 1
	% --- Males, WT vs Mut ---
	
	Mut_Male_afterSD = Normalized_power.SD.Mut.Male(:,SD_length_hrs+1:12);  % Uncomment these lines to make more general. not assuming 5 hr SD
	WT_Male_afterSD  = Normalized_power.SD.WT.Male(:,SD_length_hrs+1:12);
	Genotype = categorical([zeros(1,size(Mut_Male_afterSD,1)) ones(1,size(WT_Male_afterSD,1))]');  % Categorical is important!  otherwise ranova produces incorrect F values

	Y = [Mut_Male_afterSD; WT_Male_afterSD];
	t = array2table(Y);
	t = addvars(t,Genotype,'Before',1);
	

	Time = SD_length_hrs+1:12;   % Within-subjects variable
	
	WithinStructure = table(Time','VariableNames',{'Time'});
	WithinStructure.Time = categorical(WithinStructure.Time);

	mystr = strcat('Y1-Y',num2str(12-SD_length_hrs));  %Y1-Y7 for instance
	rm = fitrm(t,strcat(mystr,' ~ Genotype'),'WithinDesign',WithinStructure,'WithinModel','Time');
	RManovatbl_SleepPress_Male_WTVsMut = ranova(rm,'WithinModel','Time');

	Male_WTvsMut_logical = zeros(1,12-SD_length_hrs);

	% If the RM anova gives sig result for genotype, follow up with post-hoc pairwise comparisons (and BH correction)
	if RManovatbl_SleepPress_Male_WTVsMut.pValue('Genotype') < 0.05 | RManovatbl_SleepPress_Male_WTVsMut.pValue('Genotype:Time') < 0.05
		postHoc_SleepPress_Male_WTVsMut = multcompare(rm,'Genotype','By','Time','ComparisonType','lsd');  % Post-hocs 
		postHoc_SleepPress_Male_WTVsMut.pValue = mafdr(postHoc_SleepPress_Male_WTVsMut.pValue,'BHFDR', true); % BH adjustment to p-values
		[c,ia,ib] = unique(postHoc_SleepPress_Male_WTVsMut(:,1),'rows');
		postHoc_SleepPress_Male_WTVsMut = postHoc_SleepPress_Male_WTVsMut(ia,:);   	% These two rows remove repeated rows in first column since Time is repeated 2X.

		Male_WTvsMut_logical = postHoc_SleepPress_Male_WTVsMut.pValue < 0.05;  		% logical vector 1 if pVal<0.05, 0 otherwise
	end 


	% --- Females, WT vs Mut ---
	if FemaleDataPresent
		Mut_Female_afterSD = Normalized_power.SD.Mut.Female(:,SD_length_hrs+1:12);
		WT_Female_afterSD  = Normalized_power.SD.WT.Female(:,SD_length_hrs+1:12);
		Genotype = categorical([zeros(1,size(Mut_Female_afterSD,1)) ones(1,size(WT_Female_afterSD,1))]');  % Categorical is important!  otherwise ranova produces incorrect F values

		Y = [Mut_Female_afterSD; WT_Female_afterSD];
		t = array2table(Y);
		t = addvars(t,Genotype,'Before',1);
		
		Time = SD_length_hrs+1:12;   % Within-subjects variable
		WithinStructure = table(Time','VariableNames',{'Time'});
		WithinStructure.Time = categorical(WithinStructure.Time);

		mystr = strcat('Y1-Y',num2str(12-SD_length_hrs));  % Y1-Y7 for instance
		rm = fitrm(t,strcat(mystr,' ~ Genotype'),'WithinDesign',WithinStructure,'WithinModel','Time');
		RManovatbl_SleepPress_Female_WTVsMut = ranova(rm,'WithinModel','Time');

		Female_WTvsMut_logical = zeros(1,12-SD_length_hrs);
		% If the RM anova gives sig result for genotype, follow up with post-hoc pairwise comparisons (and Sidak correction)
		if RManovatbl_SleepPress_Female_WTVsMut.pValue('Genotype') < 0.05 | RManovatbl_SleepPress_Female_WTVsMut.pValue('Genotype:Time') < 0.05
			% Post-hoc comparisons, Female, WT vs Mut
			postHoc_SleepPress_Female_WTVsMut = multcompare(rm,'Genotype','By','Time','ComparisonType','lsd');  % Using no correction 
			postHoc_SleepPress_Female_WTVsMut.pValue = mafdr(postHoc_SleepPress_Female_WTVsMut.pValue,'BHFDR', true); % BH adjustment to p-values
			[c,ia,ib] = unique(postHoc_SleepPress_Female_WTVsMut(:,1),'rows');
			postHoc_SleepPress_Female_WTVsMut = postHoc_SleepPress_Female_WTVsMut(ia,:);   	% These two rows remove repeated rows in first column since Time is repeated 2X.

			Female_WTvsMut_logical = postHoc_SleepPress_Female_WTVsMut.pValue < 0.05;  		% logical vector 1 if pVal<0.05, 0 otherwise
		end  
	else 
		RManovatbl_SleepPress_Female_WTVsMut = []; 
		RManovatbl_SleepPress_WT_MVsF        = [];
		RManovatbl_SleepPress_Mut_MVsF 		 = []; 	
	end % if female data present

	% --- Male and Female together, WT vs Mut --- 
	% Mut_Female_6_12hr = Normalized_power.Mut.Female(:,6:12);
	% WT_Female_6_12hr  = Normalized_power.WT.Female(:,6:12);
	% Mut_Male_6_12hr = Normalized_power.Mut.Male(:,6:12);
	% WT_Male_6_12hr  = Normalized_power.WT.Male(:,6:12);
	%Genotype = categorical([0 0 0 0 0 1 1 1 1 1 1]');  % Categorical is important!  otherwise ranova produces incorrect F values
	if FemaleDataPresent
		Genotype = categorical([zeros(1,size(Mut_Female_afterSD,1)) ones(1,size(WT_Female_afterSD,1)) ...
								zeros(1,size(Mut_Male_afterSD,1))   ones(1,size(WT_Male_afterSD,1))]');  % Categorical is important!  otherwise ranova produces incorrect F values
		Sex = categorical([zeros(1,size(Mut_Female_afterSD,1)+size(WT_Female_afterSD,1)) ones(1,size(Mut_Male_afterSD,1)+size(WT_Male_afterSD,1))]');


		Y = [Mut_Female_afterSD; WT_Female_afterSD; Mut_Male_afterSD; WT_Male_afterSD];
		t = array2table(Y);
		t = addvars(t,Sex,'Before',1);
		t = addvars(t,Genotype,'Before',1);

		
		Time = SD_length_hrs+1:12;   % Within-subjects variable
		WithinStructure = table(Time','VariableNames',{'Time'});
		WithinStructure.Time = categorical(WithinStructure.Time);

		mystr = strcat('Y1-Y',num2str(12-SD_length_hrs));  % Y1-Y7 for instance
		rm = fitrm(t,strcat(mystr,' ~ Genotype*Sex'),'WithinDesign',WithinStructure,'WithinModel','Time');
		RManovatbl_SleepPress_MaleVsFemale_WTVsMut = ranova(rm,'WithinModel','Time');
		%mauchly(rm)
		% Follow up with post-hocs if Genotype is a significant factor
		% Female data
		% Female_WTvsMut_logical = zeros(1,7);
		% % If the RM anova gives sig result for genotype, follow up with post-hoc pairwise comparisons (and Sidak correction)
		% if RManovatbl_SleepPress_MaleVsFemale_WTVsMut.pValue('Genotype') < 0.05 | RManovatbl_SleepPress_MaleVsFemale_WTVsMut.pValue('Genotype:Sex') < 0.05
		% 	% Post-hoc comparisons, Female, WT vs Mut
		% 	postHoc_SleepPress_Female_WTVsMut = multcompare(rm,'Genotype','By','Time','ComparisonType','dunn-sidak');  % Using Sidak correction. 
		% 	[c,ia,ib] = unique(postHoc_SleepPress_Female_WTVsMut(:,1),'rows');
		% 	postHoc_SleepPress_Female_WTVsMut = postHoc_SleepPress_Female_WTVsMut(ia,:);   	% These two rows remove repeated rows in first column since Time is repeated 2X.

		% 	Female_WTvsMut_logical = postHoc_SleepPress_Female_WTVsMut.pValue < 0.05;  		% logical vector 1 if pVal<0.05, 0 otherwise
		% end 

		% --- Compare Sex within Genotype (one for WT, one for Mutant) ---
		% - WT -
		Sex = categorical([zeros(1,size(WT_Female_afterSD,1)) ones(1,size(WT_Male_afterSD,1))]');

		Y = [WT_Female_afterSD; WT_Male_afterSD];
		t = array2table(Y);
		t = addvars(t,Sex,'Before',1);

		% t = table(Sex,Y(:,1),Y(:,2),Y(:,3),Y(:,4),Y(:,5),Y(:,6),Y(:,7), ...
		%                    'VariableNames',{'Sex','Hour6','Hour7','Hour8','Hour9','Hour10','Hour11','Hour12'});
		Time = SD_length_hrs+1:12;   % Within-subjects variable
		WithinStructure = table(Time','VariableNames',{'Time'});
		WithinStructure.Time = categorical(WithinStructure.Time);

		mystr = strcat('Y1-Y',num2str(12-SD_length_hrs));  % Y1-Y7 for instance
		rm = fitrm(t,strcat(mystr,' ~ Sex'),'WithinDesign',WithinStructure,'WithinModel','Time');
		RManovatbl_SleepPress_WT_MVsF = ranova(rm,'WithinModel','Time');

		% - Mutant -
		Sex = categorical([zeros(1,size(Mut_Female_afterSD,1)) ones(1,size(Mut_Male_afterSD,1))]');

		Y = [Mut_Female_afterSD; Mut_Male_afterSD];
		t = array2table(Y);
		t = addvars(t,Sex,'Before',1);

		% t = table(Sex,Y(:,1),Y(:,2),Y(:,3),Y(:,4),Y(:,5),Y(:,6),Y(:,7), ...
		%                    'VariableNames',{'Sex','Hour6','Hour7','Hour8','Hour9','Hour10','Hour11','Hour12'});
		Time = SD_length_hrs+1:12;   % Within-subjects variable
		WithinStructure = table(Time','VariableNames',{'Time'});
		WithinStructure.Time = categorical(WithinStructure.Time);

		mystr = strcat('Y1-Y',num2str(12-SD_length_hrs));  % Y1-Y7 for instance
		rm = fitrm(t,strcat(mystr,' ~ Sex'),'WithinDesign',WithinStructure,'WithinModel','Time');
		RManovatbl_SleepPress_Mut_MVsF = ranova(rm,'WithinModel','Time');


		% -- Adjust the p-values for the 4 post-hocs done:  WT_MvsF, Mut_MvsF, M_WTvsMut, F_WTvsMut
		% - First make vectors for each group of 4 p-values
		Intercept_p_vals = [RManovatbl_SleepPress_WT_MVsF{'(Intercept)','pValue'}      ...
							RManovatbl_SleepPress_Mut_MVsF{'(Intercept)','pValue'}     ...
							RManovatbl_SleepPress_Male_WTVsMut{'(Intercept)','pValue'} ...
							RManovatbl_SleepPress_Female_WTVsMut{'(Intercept)','pValue'}];

		Group_p_vals = [RManovatbl_SleepPress_WT_MVsF{'Sex','pValue'}  ...
						RManovatbl_SleepPress_Mut_MVsF{'Sex','pValue'} ...
						RManovatbl_SleepPress_Male_WTVsMut{'Genotype','pValue'} ...
						RManovatbl_SleepPress_Female_WTVsMut{'Genotype','pValue'}];

		InterXTime_p_vals = [RManovatbl_SleepPress_WT_MVsF{'(Intercept):Time','pValue'}      ...
							RManovatbl_SleepPress_Mut_MVsF{'(Intercept):Time','pValue'}     ...
							RManovatbl_SleepPress_Male_WTVsMut{'(Intercept):Time','pValue'} ...
							RManovatbl_SleepPress_Female_WTVsMut{'(Intercept):Time','pValue'}];

		GroupXTime_p_vals = [RManovatbl_SleepPress_WT_MVsF{'Sex:Time','pValue'}  ...
							 RManovatbl_SleepPress_Mut_MVsF{'Sex:Time','pValue'} ...
							 RManovatbl_SleepPress_Male_WTVsMut{'Genotype:Time','pValue'} ...
							 RManovatbl_SleepPress_Female_WTVsMut{'Genotype:Time','pValue'}];

		% Adjust each group of p-values using mafdr
		Intercept_p_valsADJ = mafdr(Intercept_p_vals,'BHFDR',true);
		Group_p_valsADJ 	= mafdr(Group_p_vals,'BHFDR',true);
		InterXTime_pvalsADJ = mafdr(InterXTime_p_vals,'BHFDR',true);
		GroupXTime_p_valsADJ= mafdr(GroupXTime_p_vals,'BHFDR',true);

		% Add adjusted p-values (padded with NaNs) to RM anova tables
		RManovatbl_SleepPress_WT_MVsF  	   	 = addvars(RManovatbl_SleepPress_WT_MVsF, 		[Intercept_p_valsADJ(1) Group_p_valsADJ(1) NaN InterXTime_pvalsADJ(1) GroupXTime_p_valsADJ(1) NaN]', 'After', 'pValue', 'NewVariableNames', 'pValueBH_adjusted');
		RManovatbl_SleepPress_Mut_MVsF 	   	 = addvars(RManovatbl_SleepPress_Mut_MVsF,		[Intercept_p_valsADJ(2) Group_p_valsADJ(2) NaN InterXTime_pvalsADJ(2) GroupXTime_p_valsADJ(2) NaN]', 'After', 'pValue', 'NewVariableNames', 'pValueBH_adjusted');
		RManovatbl_SleepPress_Male_WTVsMut 	 = addvars(RManovatbl_SleepPress_Male_WTVsMut,  [Intercept_p_valsADJ(3) Group_p_valsADJ(3) NaN InterXTime_pvalsADJ(3) GroupXTime_p_valsADJ(3) NaN]', 'After', 'pValue', 'NewVariableNames', 'pValueBH_adjusted');
		RManovatbl_SleepPress_Female_WTVsMut = addvars(RManovatbl_SleepPress_Female_WTVsMut,[Intercept_p_valsADJ(4) Group_p_valsADJ(4) NaN InterXTime_pvalsADJ(4) GroupXTime_p_valsADJ(4) NaN]', 'After', 'pValue', 'NewVariableNames', 'pValueBH_adjusted'); 


	else 
		RManovatbl_SleepPress_MaleVsFemale_WTVsMut = [];
		RManovatbl_SleepPress_WT_MVsF        = [];
		RManovatbl_SleepPress_Mut_MVsF 		 = [];
	end  
	% -- end of BinSizeHrs = 1 case 


% -- BinSizeHrs = 2 case ---------------------------------
elseif BinSizeHrs==2
	% --- Males, WT vs Mut ---
	Mut_Male_6_12hr = Normalized_power.SD.Mut.Male(:,6:12);
	WT_Male_6_12hr  = Normalized_power.SD.WT.Male(:,6:12);
	
	Mut_Male_2hrbins = [mean(Mut_Male_6_12hr(:,1:2),2) mean(Mut_Male_6_12hr(:,3:4),2) mean(Mut_Male_6_12hr(:,5:6),2)];
	WT_Male_2hrbins  = [mean(WT_Male_6_12hr(:,1:2),2)  mean(WT_Male_6_12hr(:,3:4),2)  mean(WT_Male_6_12hr(:,5:6),2)];

	Genotype = categorical([zeros(1,size(Mut_Male_2hrbins,1)) ones(1,size(WT_Male_2hrbins,1))]');  % Categorical is important!  otherwise ranova produces incorrect F values

	Y = [Mut_Male_2hrbins; WT_Male_2hrbins];
	t = table(Genotype,Y(:,1),Y(:,2),Y(:,3), ...
	                   'VariableNames',{'Genotype','Hour67','Hour89','Hour1011'});

	TimeBins = 1:3;   % Within-subjects variable
	WithinStructure = table(TimeBins','VariableNames',{'Time'});
	WithinStructure.Time = categorical(WithinStructure.Time);

	rm = fitrm(t,'Hour67-Hour1011 ~ Genotype','WithinDesign',WithinStructure,'WithinModel','Time');
	RManovatbl_SleepPress_2hrBins_Male_WTVsMut = ranova(rm,'WithinModel','Time');

	Male_WTvsMut_2hrBins_logical = zeros(1,3);

	% If the RM anova gives sig result for genotype, follow up with post-hoc pairwise comparisons (and BH correction)
	if RManovatbl_SleepPress_2hrBins_Male_WTVsMut.pValue('Genotype') < 0.05 | RManovatbl_SleepPress_2hrBins_Male_WTVsMut.pValue('Genotype:Time') < 0.05 
		postHoc_SleepPress_2hrBins_Male_WTVsMut = multcompare(rm,'Genotype','By','Time','ComparisonType','lsd');  % Post-hocs 
		postHoc_SleepPress_2hrBins_Male_WTVsMut.pValue = mafdr(postHoc_SleepPress_2hrBins_Male_WTVsMut.pValue,'BHFDR', true); % BH adjustment to p-values
		[c,ia,ib] = unique(postHoc_SleepPress_2hrBins_Male_WTVsMut(:,1),'rows');
		postHoc_SleepPress_2hrBins_Male_WTVsMut = postHoc_SleepPress_2hrBins_Male_WTVsMut(ia,:);   	% These two rows remove repeated rows in first column since Time is repeated 2X.

		Male_WTvsMut_2hrBins_logical = postHoc_SleepPress_2hrBins_Male_WTVsMut.pValue < 0.05;  		% logical vector 1 if pVal<0.05, 0 otherwise
	end 


	% --- Females, WT vs Mut ---
	if FemaleDataPresent
		Mut_Female_6_12hr = Normalized_power.SD.Mut.Female(:,6:12);
		WT_Female_6_12hr  = Normalized_power.SD.WT.Female(:,6:12);
		
		Mut_Female_2hrbins = [mean(Mut_Female_6_12hr(:,1:2),2) mean(Mut_Female_6_12hr(:,3:4),2) mean(Mut_Female_6_12hr(:,5:6),2)];
		WT_Female_2hrbins  = [mean(WT_Female_6_12hr(:,1:2),2)  mean(WT_Female_6_12hr(:,3:4),2)  mean(WT_Female_6_12hr(:,5:6),2)];


		%Genotype = categorical([0 0 0 0 0 1 1 1 1 1 1]');  % Categorical is important!  otherwise ranova produces incorrect F values
		Genotype = categorical([zeros(1,size(Mut_Female_2hrbins,1)) ones(1,size(WT_Female_2hrbins,1))]');  % Categorical is important!  otherwise ranova produces incorrect F values

		Y = [Mut_Female_2hrbins; WT_Female_2hrbins];
		t = table(Genotype,Y(:,1),Y(:,2),Y(:,3), ...
		                   'VariableNames',{'Genotype','Hour67','Hour89','Hour1011'});

		TimeBins = 1:3;   % Within-subjects variable
		WithinStructure = table(TimeBins','VariableNames',{'Time'});
		WithinStructure.Time = categorical(WithinStructure.Time);

		rm = fitrm(t,'Hour67-Hour1011 ~ Genotype','WithinDesign',WithinStructure,'WithinModel','Time');
		RManovatbl_SleepPress_2hrBins_Female_WTVsMut = ranova(rm,'WithinModel','Time');

		Female_WTvsMut_2hrBins_logical = zeros(1,3);
		% If the RM anova gives sig result for genotype, follow up with post-hoc pairwise comparisons (and Sidak correction)
		if RManovatbl_SleepPress_2hrBins_Female_WTVsMut.pValue('Genotype') < 0.05 | RManovatbl_SleepPress_2hrBins_Female_WTVsMut.pValue('Genotype:Time') < 0.05
			% Post-hoc comparisons, Female, WT vs Mut
			postHoc_SleepPress_2hrBins_Female_WTVsMut = multcompare(rm,'Genotype','By','Time','ComparisonType','lsd');  % Using no correction 
			postHoc_SleepPress_2hrBins_Female_WTVsMut.pValue = mafdr(postHoc_SleepPress_2hrBins_Female_WTVsMut.pValue,'BHFDR', true); % BH adjustment to p-values
			[c,ia,ib] = unique(postHoc_SleepPress_2hrBins_Female_WTVsMut(:,1),'rows');
			postHoc_SleepPress_2hrBins_Female_WTVsMut = postHoc_SleepPress_2hrBins_Female_WTVsMut(ia,:);   	% These two rows remove repeated rows in first column since Time is repeated 2X.

			Female_WTvsMut_2hrBins_logical = postHoc_SleepPress_2hrBins_Female_WTVsMut.pValue < 0.05;  		% logical vector 1 if pVal<0.05, 0 otherwise
		end  
	else 
		RManovatbl_SleepPress_Female_WTVsMut = [];
		RManovatbl_SleepPress_WT_MVsF        = [];
		RManovatbl_SleepPress_Mut_MVsF 		 = [];

	end % if female data present

	% --- Male and Female together, WT vs Mut --- 
	% Mut_Female_6_12hr = Normalized_power.Mut.Female(:,6:12);
	% WT_Female_6_12hr  = Normalized_power.WT.Female(:,6:12);
	% Mut_Male_6_12hr = Normalized_power.Mut.Male(:,6:12);
	% WT_Male_6_12hr  = Normalized_power.WT.Male(:,6:12);
	%Genotype = categorical([0 0 0 0 0 1 1 1 1 1 1]');  % Categorical is important!  otherwise ranova produces incorrect F values
	if FemaleDataPresent
		Genotype = categorical([zeros(1,size(Mut_Female_2hrbins,1)) ones(1,size(WT_Female_2hrbins,1)) ...
								zeros(1,size(Mut_Male_2hrbins,1))   ones(1,size(WT_Male_2hrbins,1))]');  % Categorical is important!  otherwise ranova produces incorrect F values
		Sex = categorical([zeros(1,size(Mut_Female_2hrbins,1)+size(WT_Female_2hrbins,1)) ones(1,size(Mut_Male_2hrbins,1)+size(WT_Male_2hrbins,1))]');


		Y = [Mut_Female_2hrbins; WT_Female_2hrbins; Mut_Male_2hrbins; WT_Male_2hrbins];
		t = table(Genotype,Sex,Y(:,1),Y(:,2),Y(:,3), ...
		                   'VariableNames',{'Genotype','Sex','Hour67','Hour89','Hour1011'});

		TimeBins = 1:3;   % Within-subjects variable
		WithinStructure = table(TimeBins','VariableNames',{'Time'});
		WithinStructure.Time = categorical(WithinStructure.Time);

		rm = fitrm(t,'Hour67-Hour1011 ~ Genotype*Sex','WithinDesign',WithinStructure,'WithinModel','Time');
		RManovatbl_SleepPress_2hrBins_MaleVsFemale_WTVsMut = ranova(rm,'WithinModel','Time');
		%mauchly(rm)
		% Follow up with post-hocs if Genotype is a significant factor
		% Female data
		% Female_WTvsMut_logical = zeros(1,7);
		% % If the RM anova gives sig result for genotype, follow up with post-hoc pairwise comparisons (and Sidak correction)
		% if RManovatbl_SleepPress_MaleVsFemale_WTVsMut.pValue('Genotype') < 0.05 | RManovatbl_SleepPress_MaleVsFemale_WTVsMut.pValue('Genotype:Sex') < 0.05
		% 	% Post-hoc comparisons, Female, WT vs Mut
		% 	postHoc_SleepPress_Female_WTVsMut = multcompare(rm,'Genotype','By','Time','ComparisonType','dunn-sidak');  % Using Sidak correction. 
		% 	[c,ia,ib] = unique(postHoc_SleepPress_Female_WTVsMut(:,1),'rows');
		% 	postHoc_SleepPress_Female_WTVsMut = postHoc_SleepPress_Female_WTVsMut(ia,:);   	% These two rows remove repeated rows in first column since Time is repeated 2X.

		% 	Female_WTvsMut_logical = postHoc_SleepPress_Female_WTVsMut.pValue < 0.05;  		% logical vector 1 if pVal<0.05, 0 otherwise
		% end 

		% --- Compare Sex within Genotype (one for WT, one for Mutant) ---
		% - WT -
		Sex = categorical([zeros(1,size(WT_Female_2hrbins,1)) ones(1,size(WT_Male_2hrbins,1))]');

		Y = [WT_Female_2hrbins; WT_Male_2hrbins];
		t = table(Sex,Y(:,1),Y(:,2),Y(:,3), ...
		                   'VariableNames',{'Sex','Hour67','Hour89','Hour1011'});
		TimeBins = 1:3;   % Within-subjects variable
		WithinStructure = table(TimeBins','VariableNames',{'Time'});
		WithinStructure.Time = categorical(WithinStructure.Time);

		rm = fitrm(t,'Hour67-Hour1011 ~ Sex','WithinDesign',WithinStructure,'WithinModel','Time');
		RManovatbl_SleepPress_2hrBins_WT_MVsF = ranova(rm,'WithinModel','Time');

		% - Mutant -
		Sex = categorical([zeros(1,size(Mut_Female_2hrbins,1)) ones(1,size(Mut_Male_2hrbins,1))]');

		Y = [Mut_Female_2hrbins; Mut_Male_2hrbins];
		t = table(Sex,Y(:,1),Y(:,2),Y(:,3), ...
		                   'VariableNames',{'Sex','Hour67','Hour89','Hour1011'});
		TimeBins = 1:3;   % Within-subjects variable
		WithinStructure = table(TimeBins','VariableNames',{'Time'});
		WithinStructure.Time = categorical(WithinStructure.Time);

		rm = fitrm(t,'Hour67-Hour1011 ~ Sex','WithinDesign',WithinStructure,'WithinModel','Time');
		RManovatbl_SleepPress_2hrBins_Mut_MVsF = ranova(rm,'WithinModel','Time');


		% -- Adjust the p-values for the 4 post-hocs done:  WT_MvsF, Mut_MvsF, M_WTvsMut, F_WTvsMut
		% - First make vectors for each group of 4 p-values
		Intercept_p_vals = [RManovatbl_SleepPress_2hrBins_WT_MVsF{'(Intercept)','pValue'}      ...
							RManovatbl_SleepPress_2hrBins_Mut_MVsF{'(Intercept)','pValue'}     ...
							RManovatbl_SleepPress_2hrBins_Male_WTVsMut{'(Intercept)','pValue'} ...
							RManovatbl_SleepPress_2hrBins_Female_WTVsMut{'(Intercept)','pValue'}];

		Group_p_vals = [RManovatbl_SleepPress_2hrBins_WT_MVsF{'Sex','pValue'}  ...
						RManovatbl_SleepPress_2hrBins_Mut_MVsF{'Sex','pValue'} ...
						RManovatbl_SleepPress_2hrBins_Male_WTVsMut{'Genotype','pValue'} ...
						RManovatbl_SleepPress_2hrBins_Female_WTVsMut{'Genotype','pValue'}];

		InterXTime_p_vals = [RManovatbl_SleepPress_2hrBins_WT_MVsF{'(Intercept):Time','pValue'}      ...
							 RManovatbl_SleepPress_2hrBins_Mut_MVsF{'(Intercept):Time','pValue'}     ...
							 RManovatbl_SleepPress_2hrBins_Male_WTVsMut{'(Intercept):Time','pValue'} ...
							 RManovatbl_SleepPress_2hrBins_Female_WTVsMut{'(Intercept):Time','pValue'}];

		GroupXTime_p_vals = [RManovatbl_SleepPress_2hrBins_WT_MVsF{'Sex:Time','pValue'}  ...
							 RManovatbl_SleepPress_2hrBins_Mut_MVsF{'Sex:Time','pValue'} ...
							 RManovatbl_SleepPress_2hrBins_Male_WTVsMut{'Genotype:Time','pValue'} ...
							 RManovatbl_SleepPress_2hrBins_Female_WTVsMut{'Genotype:Time','pValue'}];

		% Adjust each group of p-values using mafdr
		Intercept_p_valsADJ = mafdr(Intercept_p_vals,'BHFDR',true);
		Group_p_valsADJ 	= mafdr(Group_p_vals,'BHFDR',true);
		InterXTime_pvalsADJ = mafdr(InterXTime_p_vals,'BHFDR',true);
		GroupXTime_p_valsADJ= mafdr(GroupXTime_p_vals,'BHFDR',true);

		% Add adjusted p-values (padded with NaNs) to RM anova tables
		RManovatbl_SleepPress_2hrBins_WT_MVsF  	   	 = addvars(RManovatbl_SleepPress_2hrBins_WT_MVsF, 		[Intercept_p_valsADJ(1) Group_p_valsADJ(1) NaN InterXTime_pvalsADJ(1) GroupXTime_p_valsADJ(1) NaN]', 'After', 'pValue', 'NewVariableNames', 'pValueBH_adjusted');
		RManovatbl_SleepPress_2hrBins_Mut_MVsF 	   	 = addvars(RManovatbl_SleepPress_2hrBins_Mut_MVsF,		[Intercept_p_valsADJ(2) Group_p_valsADJ(2) NaN InterXTime_pvalsADJ(2) GroupXTime_p_valsADJ(2) NaN]', 'After', 'pValue', 'NewVariableNames', 'pValueBH_adjusted');
		RManovatbl_SleepPress_2hrBins_Male_WTVsMut 	 = addvars(RManovatbl_SleepPress_2hrBins_Male_WTVsMut,  [Intercept_p_valsADJ(3) Group_p_valsADJ(3) NaN InterXTime_pvalsADJ(3) GroupXTime_p_valsADJ(3) NaN]', 'After', 'pValue', 'NewVariableNames', 'pValueBH_adjusted');
		RManovatbl_SleepPress_2hrBins_Female_WTVsMut = addvars(RManovatbl_SleepPress_2hrBins_Female_WTVsMut,[Intercept_p_valsADJ(4) Group_p_valsADJ(4) NaN InterXTime_pvalsADJ(4) GroupXTime_p_valsADJ(4) NaN]', 'After', 'pValue', 'NewVariableNames', 'pValueBH_adjusted'); 


	else 
		RManovatbl_SleepPress_2hrBins_MaleVsFemale_WTVsMut = [];
	end
else
	error('You entered a BinSizeHrs other than 1 or 2.  We only have code for 1 or 2.')
end

% ---------------------------------------------------------------------------------------------------------------------
% --- End of repeated measures ANOVA ----------------------------------------------------------------------------------


% --- Make the SD figure -------------------------------------------------------------------------------------------------
% ---------------------------------------------------------------------------------------------------------------------
% -- 1-hr bin case --
if BinSizeHrs ==1
	line_thick = 2;
	filled_marker_size = 18;
	open_marker_size   = 6;
	ymaxM = max([max(Means_norm_power.SD.WT.Male(SD_length_hrs+1:12) + SEMs_norm_power.SD.WT.Male(SD_length_hrs+1:12)) max(Means_norm_power.SD.Mut.Male(SD_length_hrs+1:12) + SEMs_norm_power.SD.Mut.Male(SD_length_hrs+1:12))]) + 5;

	if FemaleDataPresent
		H=figure('Renderer', 'painters', 'Position', [600 400 650 320]);
		H.Name = 'NREM Power After SD';
		H.Renderer = 'painters';
	else
		H=figure('Renderer', 'painters', 'Position', [600 400 325 320]);
		H.Name = 'NREM Power After SD';
		H.Renderer = 'painters';
	end 

	
	if FemaleDataPresent
		subplot(1,2,1)    				% Male panel
	end  
	pm=plot(SD_length_hrs+1:12,Means_norm_power.SD.WT.Male(SD_length_hrs+1:12),'k.',SD_length_hrs+1:12,Means_norm_power.SD.Mut.Male(SD_length_hrs+1:12),'r.');
	pm(1).LineWidth=line_thick;
	pm(2).LineWidth=line_thick;
	pm(1).MarkerSize = filled_marker_size;
	pm(2).MarkerSize = filled_marker_size;
	hold on 
	eb1=errorbar(SD_length_hrs+1:12,Means_norm_power.SD.WT.Male(SD_length_hrs+1:12),SEMs_norm_power.SD.WT.Male(SD_length_hrs+1:12),'k');
	eb2=errorbar(SD_length_hrs+1:12,Means_norm_power.SD.Mut.Male(SD_length_hrs+1:12),SEMs_norm_power.SD.Mut.Male(SD_length_hrs+1:12),'r');
	eb1.LineWidth = line_thick;
	eb2.LineWidth = line_thick;
	eb1.CapSize = 0;
	eb2.CapSize = 0;
	yline(100,'--','LineWidth',3);
	% r1 = rectangle('Position',[0 75 5 2]);
	% r1.FaceColor = 'red';
	% r1.LineWidth = 1;
	p1 = patch([0 0 SD_length_hrs SD_length_hrs],[75 77 77 75],[1 1 1]);
	%p1.FaceColor = [1 1 1];  % white face color
	hatchfill2(p1,'single','HatchAngle',45,'HatchDensity',65,'HatchLineWidth',1.5);
	r2 = rectangle('Position',[SD_length_hrs 75 12-SD_length_hrs 2]);
	r2.LineWidth = line_thick;
	% if Male_WTvsMut_logical is 1, add an asterisk above the data
	for i=SD_length_hrs+1:12
		if Male_WTvsMut_logical(i-SD_length_hrs)
			text(i,max([Means_norm_power.SD.WT.Male(i)+SEMs_norm_power.SD.WT.Male(i),Means_norm_power.SD.Mut.Male(i)+SEMs_norm_power.SD.Mut.Male(i)])+3,'*','FontSize',30,'HorizontalAlignment','center');
		end 
	end 

	hold off
	%axM=gca;
	axM = H.Children;
	axM.XLim = [0 12];
	axM.YLim = [75 ymaxM];
	axM.Box= 'off';
	axM.LineWidth=line_thick;
	axM.Title.String = Sexes{1};
	axM.Title.FontWeight = 'normal';
	axM.XTick = [ 0 6 12];
	axM.YTick = [75:25:ymaxM];
	axM.XLabel.String = 'Hour';
	axM.YLabel.String = 'Normalized NREM Delta';
	%axM.YLabel.FontWeight = 'bold';
	axM.YLabel.FontSize   = 14;
	axM.FontSize = 16;
	l=legend([eb1 eb2],[Sexes{1} ' WT'],[Sexes{1} ' Mut']);
	l.Location = 'northwest';
	l.Box = 'off';
	l.FontSize = 14;

	 
	if FemaleDataPresent
		subplot(1,2,2) 				% Female panel
		
		ymaxF = max([max(Means_norm_power.SD.WT.Female(SD_length_hrs+1:12) + SEMs_norm_power.SD.WT.Female(SD_length_hrs+1:12)) max(Means_norm_power.SD.Mut.Female(SD_length_hrs+1:12) + SEMs_norm_power.SD.Mut.Female(SD_length_hrs+1:12))]) + 5;
		ymaxBOTH = max([ymaxM ymaxF]);

		
		pm=plot(SD_length_hrs+1:12,Means_norm_power.SD.WT.Female(SD_length_hrs+1:12),'ko',SD_length_hrs+1:12,Means_norm_power.SD.Mut.Female(SD_length_hrs+1:12),'ro');
		pm(1).LineWidth=line_thick;
		pm(2).LineWidth=line_thick;
		pm(1).MarkerSize = open_marker_size;
		pm(2).MarkerSize = open_marker_size;

		hold on 
		eb1=errorbar(SD_length_hrs+1:12,Means_norm_power.SD.WT.Female(SD_length_hrs+1:12),SEMs_norm_power.SD.WT.Female(SD_length_hrs+1:12),'o--','MarkerFaceColor','w','MarkerEdgeColor','k','Color','k');
		eb2=errorbar(SD_length_hrs+1:12,Means_norm_power.SD.Mut.Female(SD_length_hrs+1:12),SEMs_norm_power.SD.Mut.Female(SD_length_hrs+1:12),'o--','MarkerFaceColor','w','MarkerEdgeColor','r','Color','r');
		eb1.LineWidth = line_thick;
		eb2.LineWidth = line_thick;
		eb1.CapSize = 0;
		eb2.CapSize = 0;
		yline(100,'--','LineWidth',3);
		% r1 = rectangle('Position',[0 75 5 2]);
		% r1.FaceColor = 'red';
		% r1.LineWidth = 1;
		p1 = patch([0 0 SD_length_hrs SD_length_hrs],[75 77 77 75],[1 1 1]);
		hatchfill2(p1,'single','HatchAngle',45,'HatchDensity',65,'HatchLineWidth',1.5);
		r2 = rectangle('Position',[SD_length_hrs 75 12-SD_length_hrs 2]);
		r2.LineWidth = line_thick;
		
		% if Female_WTvsMut_logical is 1, add an asterisk above the data
		for i=SD_length_hrs+1:12
			if Female_WTvsMut_logical(i-SD_length_hrs)
				text(i,max([Means_norm_power.SD.WT.Female(i)+SEMs_norm_power.SD.WT.Female(i),Means_norm_power.SD.Mut.Female(i)+SEMs_norm_power.SD.Mut.Female(i)])+3,'*','FontSize',30,'HorizontalAlignment','center');
			end 
		end 

		%axF=gca;
		axF = H.Children(1);    % was H.Children(2), but this does not work if you add a legend above. Most recent item in H.Children is first. 
		axF.XLim = [0 12];
		axF.YLim = [75 ymaxBOTH];
		axF.Box= 'off';
		axF.LineWidth=line_thick;
		axF.Title.String = Sexes{2};
		axF.Title.FontWeight = 'normal';
		axF.XTick = [ 0 6 12];
		axF.YTick = [75:25:ymaxBOTH];
		axF.XLabel.String = 'Hour';
		axF.YLabel.String = 'Normalized NREM Delta';
		%axF.YLabel.FontWeight = 'bold';
		axF.YLabel.FontSize   = 14;
		axF.FontSize = 16;

		% axM.YLim = [75 ymaxBOTH];  % Set both M and F vertical axis to same value
		% axM.YTick = [75:25:ymaxBOTH];

		hold off
		l=legend([eb1 eb2],[Sexes{2} ' WT'],[Sexes{2} ' Mut']);
		l.Location = 'northwest';
		l.Box = 'off';
		l.FontSize = 14;


		subplot(1,2,1)  % go back to the first subplot to set the axis limits
		axM.YLim = [75 ymaxBOTH];  % Set both M and F vertical axis to same value
		axM.YTick = [75:25:ymaxBOTH];

		%ax.YLabel.String = 'Normalized NREM Delta';
		% ax.YLabel.FontWeight = 'bold';
		% ax.YLabel.FontSize   = 14;

		%  ---- Make one plot that has all 4 curves on it.  ---
		HB=figure('Renderer', 'painters');
		HB.Name = 'NREM Power after SD (both sexes)';
		% Male
		pm=plot(SD_length_hrs+1:12,Means_norm_power.SD.WT.Male(SD_length_hrs+1:12),'k.',SD_length_hrs+1:12,Means_norm_power.SD.Mut.Male(SD_length_hrs+1:12),'r.');
		pm(1).LineWidth=line_thick;
		pm(2).LineWidth=line_thick;
		pm(1).MarkerSize = filled_marker_size;
		pm(2).MarkerSize = filled_marker_size;
		hold on 
		% Female
		pf=plot(SD_length_hrs+1:12,Means_norm_power.SD.WT.Female(SD_length_hrs+1:12),'ko',SD_length_hrs+1:12,Means_norm_power.SD.Mut.Female(SD_length_hrs+1:12),'ro');
		pf(1).LineWidth=line_thick;
		pf(2).LineWidth=line_thick;
		pf(1).MarkerSize = open_marker_size;
		pf(2).MarkerSize = open_marker_size;

		% Male
		eb1m=errorbar(SD_length_hrs+1:12,Means_norm_power.SD.WT.Male(SD_length_hrs+1:12),SEMs_norm_power.SD.WT.Male(SD_length_hrs+1:12),'k');
		eb2m=errorbar(SD_length_hrs+1:12,Means_norm_power.SD.Mut.Male(SD_length_hrs+1:12),SEMs_norm_power.SD.Mut.Male(SD_length_hrs+1:12),'r');
		eb1m.LineWidth = line_thick;
		eb2m.LineWidth = line_thick;
		eb1m.CapSize = 0;
		eb2m.CapSize = 0;

		% Female 
		eb1f=errorbar(SD_length_hrs+1:12,Means_norm_power.SD.WT.Female(SD_length_hrs+1:12),SEMs_norm_power.SD.WT.Female(SD_length_hrs+1:12),'o--','MarkerFaceColor','w','MarkerEdgeColor','k','Color','k');
		eb2f=errorbar(SD_length_hrs+1:12,Means_norm_power.SD.Mut.Female(SD_length_hrs+1:12),SEMs_norm_power.SD.Mut.Female(SD_length_hrs+1:12),'o--','MarkerFaceColor','w','MarkerEdgeColor','r','Color','r');
		eb1f.LineWidth = line_thick;
		eb2f.LineWidth = line_thick;
		eb1f.CapSize = 0;
		eb2f.CapSize = 0;

		yline(100,'--','LineWidth',3);
		p1 = patch([0 0 SD_length_hrs SD_length_hrs],[75 77 77 75],[1 1 1]);
		%p1.FaceColor = [1 1 1];  % white face color
		hatchfill2(p1,'single','HatchAngle',45,'HatchDensity',65,'HatchLineWidth',1.5);
		r2 = rectangle('Position',[SD_length_hrs 75 12-SD_length_hrs 2]);
		r2.LineWidth = line_thick;
		hold off 
		
		%axB=gca;
		axB = HB.Children;
		axB.XLim = [0 12];
		axB.YLim = [75 axB.YLim(2)];
		axB.Box= 'off';
		axB.LineWidth=line_thick;
		axB.Title.String = [Sexes{1} ' and ' Sexes{2}];
		axB.Title.FontWeight = 'normal';
		axB.XTick = [ 0 6 12];
		axB.YTick = [75:25:ymaxBOTH];
		axB.XLabel.String = 'Hour';
		axB.YLabel.String = 'Normalized NREM Delta';
		%axF.YLabel.FontWeight = 'bold';
		axB.FontSize = 16;
		l=legend([eb1m eb2m eb1f eb2f],[Sexes{1} ' WT'],[Sexes{1} ' Mut'],[Sexes{2} ' WT'],[Sexes{2} ' Mut']);
		l.Location = 'northwest';
		l.Box = 'off';
		l.FontSize = 14;


		% -------  Make a plot with 2 panels, each with 2 curves: WT M vs F and Mut M vs F -------
		HB2=figure('Renderer', 'painters');
		HB2.Name = 'NREM Power after SD (separated by sex)';
		
		% - WT M vs F -
		subplot(1,2,1)
		ymaxBOTH = max([max(Means_norm_power.SD.WT.Female(SD_length_hrs+1:12) + SEMs_norm_power.SD.WT.Female(SD_length_hrs+1:12)) max(Means_norm_power.SD.WT.Male(SD_length_hrs+1:12) + SEMs_norm_power.SD.WT.Male(SD_length_hrs+1:12))]) + 5;
		

		% may not need these 5 lines
		% pm=plot(6:12,Means_norm_power.SD.WT.Female(6:12),'ko',6:12,Means_norm_power.SD.WT.Male(6:12),'k');
		% pm(1).LineWidth=line_thick;
		% pm(2).LineWidth=line_thick;
		% pm(1).MarkerSize = open_marker_size;
		% pm(2).MarkerSize = open_marker_size;

		hold on 
		eb1=errorbar(SD_length_hrs+1:12,Means_norm_power.SD.WT.Female(SD_length_hrs+1:12),SEMs_norm_power.SD.WT.Female(SD_length_hrs+1:12),'o--','MarkerFaceColor','w','MarkerEdgeColor','k','Color','k');
		eb2=errorbar(SD_length_hrs+1:12,Means_norm_power.SD.WT.Male(SD_length_hrs+1:12),  SEMs_norm_power.SD.WT.Male(SD_length_hrs+1:12),'o-','MarkerFaceColor','k','MarkerEdgeColor','k','Color','k');
		eb1.LineWidth = line_thick;
		eb2.LineWidth = line_thick;
		eb1.CapSize = 0;
		eb2.CapSize = 0;
		yline(100,'--','LineWidth',3);
		% r1 = rectangle('Position',[0 75 5 2]);
		% r1.FaceColor = 'red';
		% r1.LineWidth = 1;
		p1 = patch([0 0 SD_length_hrs SD_length_hrs],[75 77 77 75],[1 1 1]);
		hatchfill2(p1,'single','HatchAngle',45,'HatchDensity',65,'HatchLineWidth',1.5);
		r2 = rectangle('Position',[SD_length_hrs 75 12-SD_length_hrs 2]);
		r2.LineWidth = line_thick;
		
		% if WT_MvsF_logical is 1, add an asterisk above the data
		% for i=6:12
		% 	if WT_MvsF_logical(i-5)
		% 		text(i,max([Means_norm_power.SD.WT.Female(i)+SEMs_norm_power.SD.WT.Female(i),Means_norm_power.SD.Mut.Female(i)+SEMs_norm_power.SD.Mut.Female(i)])+3,'*','FontSize',30,'HorizontalAlignment','center');
		% 	end 
		% end 

		%axWTMvF=gca;
		axWTMvF = HB2.Children;
		axWTMvF.XLim = [0 12];
		axWTMvF.YLim = [75 ymaxBOTH];
		axWTMvF.Box= 'off';
		axWTMvF.LineWidth=line_thick;
		axWTMvF.Title.String = ['WT ' Sexes{1} ' V ' Sexes{2}];
		axWTMvF.Title.FontWeight = 'normal';
		axWTMvF.XTick = [ 0 6 12];
		axWTMvF.YTick = [75:25:ymaxBOTH];
		axWTMvF.XLabel.String = 'Hour';
		axWTMvF.YLabel.String = 'Normalized NREM Delta';
		%axF.YLabel.FontWeight = 'bold';
		axWTMvF.FontSize = 16;

		% axM.YLim = [75 ymaxBOTH];  % Set both M and F vertical axis to same value
		% axM.YTick = [75:25:ymaxBOTH];

		hold off
		l=legend([eb1 eb2],Sexes{2},Sexes{1});
		l.Location = 'northwest';
		l.Box = 'off';
		l.FontSize = 14;
		
		% - Mut M vs F -
		subplot(1,2,2)
		ymaxBOTH = max([max(Means_norm_power.SD.Mut.Female(SD_length_hrs+1:12) + SEMs_norm_power.SD.Mut.Female(SD_length_hrs+1:12)) max(Means_norm_power.SD.Mut.Male(SD_length_hrs+1:12) + SEMs_norm_power.SD.Mut.Male(SD_length_hrs+1:12))]) + 5;
		

		% may not need these 5 lines
		% pm=plot(6:12,Means_norm_power.SD.WT.Female(6:12),'ko',6:12,Means_norm_power.SD.WT.Male(6:12),'k');
		% pm(1).LineWidth=line_thick;
		% pm(2).LineWidth=line_thick;
		% pm(1).MarkerSize = open_marker_size;
		% pm(2).MarkerSize = open_marker_size;

		hold on 
		eb3=errorbar(SD_length_hrs+1:12,Means_norm_power.SD.Mut.Female(SD_length_hrs+1:12),SEMs_norm_power.SD.Mut.Female(SD_length_hrs+1:12),'o--','MarkerFaceColor','w','MarkerEdgeColor','r','Color','r');
		eb4=errorbar(SD_length_hrs+1:12,Means_norm_power.SD.Mut.Male(SD_length_hrs+1:12),  SEMs_norm_power.SD.Mut.Male(SD_length_hrs+1:12),'o-','MarkerFaceColor','r','MarkerEdgeColor','r','Color','r');
		eb3.LineWidth = line_thick;
		eb4.LineWidth = line_thick;
		eb3.CapSize = 0;
		eb4.CapSize = 0;
		yline(100,'--','LineWidth',3);
		% r1 = rectangle('Position',[0 75 5 2]);
		% r1.FaceColor = 'red';
		% r1.LineWidth = 1;
		p1 = patch([0 0 SD_length_hrs SD_length_hrs],[75 77 77 75],[1 1 1]);
		hatchfill2(p1,'single','HatchAngle',45,'HatchDensity',65,'HatchLineWidth',1.5);
		r2 = rectangle('Position',[SD_length_hrs 75 12-SD_length_hrs 2]);
		r2.LineWidth = line_thick;
		
		% if WT_MvsF_logical is 1, add an asterisk above the data
		% for i=6:12
		% 	if WT_MvsF_logical(i-5)
		% 		text(i,max([Means_norm_power.SD.WT.Female(i)+SEMs_norm_power.SD.WT.Female(i),Means_norm_power.SD.Mut.Female(i)+SEMs_norm_power.SD.Mut.Female(i)])+3,'*','FontSize',30,'HorizontalAlignment','center');
		% 	end 
		% end 

		%axMutMvF=gca;
		
		axMutMvF = HB2.Children(1);
		axMutMvF.XLim = [0 12];
		axMutMvF.YLim = [75 ymaxBOTH];
		axMutMvF.Box= 'off';
		axMutMvF.LineWidth=line_thick;
		axMutMvF.Title.String = ['Mut ' Sexes{1} ' V ' Sexes{2}];
		axMutMvF.Title.FontWeight = 'normal';
		axMutMvF.XTick = [ 0 6 12];
		axMutMvF.YTick = [75:25:ymaxBOTH];
		axMutMvF.XLabel.String = 'Hour';
		axMutMvF.YLabel.String = 'Normalized NREM Delta';
		%axF.YLabel.FontWeight = 'bold';
		axMutMvF.FontSize = 16;
		hold off

		l=legend([eb3 eb4],Sexes{2},Sexes{1});
		l.Location = 'northwest';
		l.Box = 'off';
		l.FontSize = 14;



	end % end if female data present 

	% ---- 2 Hours Bins ---------------------
elseif BinSizeHrs ==2 
	line_thick = 2;
	filled_marker_size = 18;
	open_marker_size   = 6;
	
	% Average to get 2-hr bins, not 1-hr
	Means_norm_power_2hrbins.WT.Male = [mean(Means_norm_power.SD.WT.Male(6:7)) mean(Means_norm_power.SD.WT.Male(8:9)) mean(Means_norm_power.SD.WT.Male(10:11))];
	SEMs_norm_power_2hrbins.WT.Male  = [mean(SEMs_norm_power.SD.WT.Male(6:7))  mean(SEMs_norm_power.SD.WT.Male(8:9))  mean(SEMs_norm_power.SD.WT.Male(10:11))];

	Means_norm_power_2hrbins.Mut.Male = [mean(Means_norm_power.SD.Mut.Male(6:7)) mean(Means_norm_power.SD.Mut.Male(8:9)) mean(Means_norm_power.SD.Mut.Male(10:11))];
	SEMs_norm_power_2hrbins.Mut.Male  = [mean(SEMs_norm_power.SD.Mut.Male(6:7))  mean(SEMs_norm_power.SD.Mut.Male(8:9))  mean(SEMs_norm_power.SD.Mut.Male(10:11))];

	ymaxM = max([max(Means_norm_power_2hrbins.WT.Male + SEMs_norm_power_2hrbins.WT.Male) max(Means_norm_power_2hrbins.Mut.Male + SEMs_norm_power_2hrbins.Mut.Male)]) + 5;

	if FemaleDataPresent
		H=figure('Renderer', 'painters', 'Position', [600 400 650 320]);
		H.Name = 'NREM Power After SD (2-hr bins)';
		H.Renderer = 'painters';
	else
		H=figure('Renderer', 'painters', 'Position', [600 400 325 320]);
		H.Name = 'NREM Power After SD (2-hr bins)';
		H.Renderer = 'painters';
	end 

	
	if FemaleDataPresent
		subplot(1,2,1)    				% Male panel
	end 
	xvals = [6.5 8.5 10.5];
	pm=plot(xvals,Means_norm_power_2hrbins.WT.Male,'k.',xvals,Means_norm_power_2hrbins.Mut.Male,'r.');
	pm(1).LineWidth=line_thick;
	pm(2).LineWidth=line_thick;
	pm(1).MarkerSize = filled_marker_size;
	pm(2).MarkerSize = filled_marker_size;
	hold on 
	eb1=errorbar(xvals,Means_norm_power_2hrbins.WT.Male,SEMs_norm_power_2hrbins.WT.Male,'k');
	eb2=errorbar(xvals,Means_norm_power_2hrbins.Mut.Male,SEMs_norm_power_2hrbins.Mut.Male,'r');
	eb1.LineWidth = line_thick;
	eb2.LineWidth = line_thick;
	eb1.CapSize = 0;
	eb2.CapSize = 0;
	yline(100,'--','LineWidth',3);
	% r1 = rectangle('Position',[0 75 5 2]);
	% r1.FaceColor = 'red';
	% r1.LineWidth = 1;
	p1 = patch([0 0 SD_length_hrs SD_length_hrs],[75 77 77 75],[1 1 1]);
	%p1.FaceColor = [1 1 1];  % white face color
	hatchfill2(p1,'single','HatchAngle',45,'HatchDensity',65,'HatchLineWidth',1.5);
	r2 = rectangle('Position',[SD_length_hrs 75 12-SD_length_hrs 2]);
	r2.LineWidth = line_thick;
	% if Male_WTvsMut_logical is 1, add an asterisk above the data
	for i=1:3
		if Male_WTvsMut_2hrBins_logical(i)
			text(xvals(i),max([Means_norm_power_2hrbins.WT.Male(i)+SEMs_norm_power_2hrbins.WT.Male(i),Means_norm_power_2hrbins.Mut.Male(i)+SEMs_norm_power_2hrbins.Mut.Male(i)])+3,'*','FontSize',30,'HorizontalAlignment','center');
		end 
	end 

	hold off
	%axM=gca;
	axM = H.Children;
	axM.XLim = [0 12];
	axM.YLim = [75 ymaxM];
	axM.Box= 'off';
	axM.LineWidth=line_thick;
	axM.Title.String = [Sexes{1} ' (2-hr bins)'];
	axM.Title.FontWeight = 'normal';
	axM.XTick = [ 0 6 12];
	axM.YTick = [75:25:ymaxM];
	axM.XLabel.String = 'Hour';
	axM.YLabel.String = 'Normalized NREM Delta';
	%axM.YLabel.FontWeight = 'bold';
	axM.YLabel.FontSize   = 14;
	axM.FontSize = 16;
	l=legend([eb1 eb2],[Sexes{1} ' WT'],[Sexes{1} ' Mut']);
	l.Location = 'northwest';
	l.Box = 'off';
	l.FontSize = 14;

	 
	if FemaleDataPresent
		subplot(1,2,2) 				% Female panel
		
		% Average to get 2-hr bins, not 1-hr
		Means_norm_power_2hrbins.WT.Female = [mean(Means_norm_power.SD.WT.Female(6:7)) mean(Means_norm_power.SD.WT.Female(8:9)) mean(Means_norm_power.SD.WT.Female(10:11))];
		SEMs_norm_power_2hrbins.WT.Female  = [mean(SEMs_norm_power.SD.WT.Female(6:7))  mean(SEMs_norm_power.SD.WT.Female(8:9))  mean(SEMs_norm_power.SD.WT.Female(10:11))];

		Means_norm_power_2hrbins.Mut.Female = [mean(Means_norm_power.SD.Mut.Female(6:7)) mean(Means_norm_power.SD.Mut.Female(8:9)) mean(Means_norm_power.SD.Mut.Female(10:11))];
		SEMs_norm_power_2hrbins.Mut.Female  = [mean(SEMs_norm_power.SD.Mut.Female(6:7))  mean(SEMs_norm_power.SD.Mut.Female(8:9))  mean(SEMs_norm_power.SD.Mut.Female(10:11))];



		ymaxF = max([max(Means_norm_power_2hrbins.WT.Female + SEMs_norm_power_2hrbins.WT.Female) max(Means_norm_power_2hrbins.Mut.Female + SEMs_norm_power_2hrbins.Mut.Female)]) + 5;
		ymaxBOTH = max([ymaxM ymaxF]);

		
		pm=plot(xvals,Means_norm_power_2hrbins.WT.Female,'ko',xvals,Means_norm_power_2hrbins.Mut.Female,'ro');
		pm(1).LineWidth=line_thick;
		pm(2).LineWidth=line_thick;
		pm(1).MarkerSize = open_marker_size;
		pm(2).MarkerSize = open_marker_size;

		hold on 
		eb1=errorbar(xvals,Means_norm_power_2hrbins.WT.Female,SEMs_norm_power_2hrbins.WT.Female,'o--','MarkerFaceColor','w','MarkerEdgeColor','k','Color','k');
		eb2=errorbar(xvals,Means_norm_power_2hrbins.Mut.Female,SEMs_norm_power_2hrbins.Mut.Female,'o--','MarkerFaceColor','w','MarkerEdgeColor','r','Color','r');
		eb1.LineWidth = line_thick;
		eb2.LineWidth = line_thick;
		eb1.CapSize = 0;
		eb2.CapSize = 0;
		yline(100,'--','LineWidth',3);
		% r1 = rectangle('Position',[0 75 5 2]);
		% r1.FaceColor = 'red';
		% r1.LineWidth = 1;
		p1 = patch([0 0 SD_length_hrs SD_length_hrs],[75 77 77 75],[1 1 1]);
		hatchfill2(p1,'single','HatchAngle',45,'HatchDensity',65,'HatchLineWidth',1.5);
		r2 = rectangle('Position',[SD_length_hrs 75 12-SD_length_hrs 2]);
		r2.LineWidth = line_thick;
		
		% if Female_WTvsMut_logical is 1, add an asterisk above the data
		for i=1:3
			if Female_WTvsMut_2hrBins_logical(i)
				text(xvals(i),max([Means_norm_power_2hrbins.WT.Female(i)+SEMs_norm_power_2hrbins.WT.Female(i),Means_norm_power_2hrbins.Mut.Female(i)+SEMs_norm_power_2hrbins.Mut.Female(i)])+3,'*','FontSize',30,'HorizontalAlignment','center');
			end 
		end 

		%axF=gca;
		axF = H.Children(1);
		axF.XLim = [0 12];
		axF.YLim = [75 ymaxBOTH];
		axF.Box= 'off';
		axF.LineWidth=line_thick;
		axF.Title.String = [Sexes{2} ' (2-hr bins)'];
		axF.Title.FontWeight = 'normal';
		axF.XTick = [ 0 6 12];
		axF.YTick = [75:25:ymaxBOTH];
		axF.XLabel.String = 'Hour';
		axF.YLabel.String = 'Normalized NREM Delta';
		%axF.YLabel.FontWeight = 'bold';
		axF.YLabel.FontSize   = 14;
		axF.FontSize = 16;

		% axM.YLim = [75 ymaxBOTH];  % Set both M and F vertical axis to same value
		% axM.YTick = [75:25:ymaxBOTH];

		hold off
		l=legend([eb1 eb2],[Sexes{2} ' WT'],[Sexes{2} ' Mut']);
		l.Location = 'northwest';
		l.Box = 'off';
		l.FontSize = 14;


		subplot(1,2,1)  % go back to the first subplot to set the axis limits
		axM.YLim = [75 ymaxBOTH];  % Set both M and F vertical axis to same value
		axM.YTick = [75:25:ymaxBOTH];

		%ax.YLabel.String = 'Normalized NREM Delta';
		% ax.YLabel.FontWeight = 'bold';
		% ax.YLabel.FontSize   = 14;

		%  ---- Make one plot that has all 4 curves on it.  ---
		HB=figure('Renderer', 'painters');
		HB.Name = 'NREM Power after SD (both sexes) 2-hr bins';
		% Male
		pm=plot(xvals,Means_norm_power_2hrbins.WT.Male,'k.',xvals,Means_norm_power_2hrbins.Mut.Male,'r.');
		pm(1).LineWidth=line_thick;
		pm(2).LineWidth=line_thick;
		pm(1).MarkerSize = filled_marker_size;
		pm(2).MarkerSize = filled_marker_size;
		hold on 
		% Female
		pf=plot(xvals,Means_norm_power_2hrbins.WT.Female,'ko',xvals,Means_norm_power_2hrbins.Mut.Female,'ro');
		pf(1).LineWidth=line_thick;
		pf(2).LineWidth=line_thick;
		pf(1).MarkerSize = open_marker_size;
		pf(2).MarkerSize = open_marker_size;

		% Male
		eb1m=errorbar(xvals,Means_norm_power_2hrbins.WT.Male,SEMs_norm_power_2hrbins.WT.Male,'k');
		eb2m=errorbar(xvals,Means_norm_power_2hrbins.Mut.Male,SEMs_norm_power_2hrbins.Mut.Male,'r');
		eb1m.LineWidth = line_thick;
		eb2m.LineWidth = line_thick;
		eb1m.CapSize = 0;
		eb2m.CapSize = 0;

		% Female 
		eb1f=errorbar(xvals,Means_norm_power_2hrbins.WT.Female,SEMs_norm_power_2hrbins.WT.Female,'o--','MarkerFaceColor','w','MarkerEdgeColor','k','Color','k');
		eb2f=errorbar(xvals,Means_norm_power_2hrbins.Mut.Female,SEMs_norm_power_2hrbins.Mut.Female,'o--','MarkerFaceColor','w','MarkerEdgeColor','r','Color','r');
		eb1f.LineWidth = line_thick;
		eb2f.LineWidth = line_thick;
		eb1f.CapSize = 0;
		eb2f.CapSize = 0;

		yline(100,'--','LineWidth',3);
		p1 = patch([0 0 SD_length_hrs SD_length_hrs],[75 77 77 75],[1 1 1]);
		%p1.FaceColor = [1 1 1];  % white face color
		hatchfill2(p1,'single','HatchAngle',45,'HatchDensity',65,'HatchLineWidth',1.5);
		r2 = rectangle('Position',[SD_length_hrs 75 12-SD_length_hrs 2]);
		r2.LineWidth = line_thick;
		hold off 
		
		%axB=gca;
		axB = HB.Children;
		axB.XLim = [0 12];
		axB.YLim = [75 axB.YLim(2)];
		axB.Box= 'off';
		axB.LineWidth=line_thick;
		axB.Title.String = [Sexes{1} ' and ' Sexes{2} ' (2-hr bins)'];
		axB.Title.FontWeight = 'normal';
		axB.XTick = [ 0 6 12];
		axB.YTick = [75:25:ymaxBOTH];
		axB.XLabel.String = 'Hour';
		axB.YLabel.String = 'Normalized NREM Delta';
		%axF.YLabel.FontWeight = 'bold';
		axB.FontSize = 16;
		l=legend([eb1m eb2m eb1f eb2f],[Sexes{1} ' WT'],[Sexes{1} ' Mut'],[Sexes{2} ' WT'],[Sexes{2} ' Mut']);
		l.Location = 'northwest';
		l.Box = 'off';
		l.FontSize = 14;


		% -------  Make a plot with 2 panels, each with 2 curves: WT M vs F and Mut M vs F -------
		HB=figure('Renderer', 'painters');
		HB.Name = 'NREM Power after SD (separated by sex) (2-hr bins)';
		
		% - WT M vs F -
		subplot(1,2,1)
		ymaxBOTH = max([max(Means_norm_power_2hrbins.WT.Female + SEMs_norm_power_2hrbins.WT.Female) max(Means_norm_power_2hrbins.WT.Male + SEMs_norm_power_2hrbins.WT.Male)]) + 5;
		

		% may not need these 5 lines
		% pm=plot(6:12,Means_norm_power.SD.WT.Female(6:12),'ko',6:12,Means_norm_power.SD.WT.Male(6:12),'k');
		% pm(1).LineWidth=line_thick;
		% pm(2).LineWidth=line_thick;
		% pm(1).MarkerSize = open_marker_size;
		% pm(2).MarkerSize = open_marker_size;

		hold on 
		eb1=errorbar(xvals,Means_norm_power_2hrbins.WT.Female,SEMs_norm_power_2hrbins.WT.Female,'o--','MarkerFaceColor','w','MarkerEdgeColor','k','Color','k');
		eb2=errorbar(xvals,Means_norm_power_2hrbins.WT.Male,  SEMs_norm_power_2hrbins.WT.Male,'o-','MarkerFaceColor','k','MarkerEdgeColor','k','Color','k');
		eb1.LineWidth = line_thick;
		eb2.LineWidth = line_thick;
		eb1.CapSize = 0;
		eb2.CapSize = 0;
		yline(100,'--','LineWidth',3);
		% r1 = rectangle('Position',[0 75 5 2]);
		% r1.FaceColor = 'red';
		% r1.LineWidth = 1;
		p1 = patch([0 0 SD_length_hrs SD_length_hrs],[75 77 77 75],[1 1 1]);
		hatchfill2(p1,'single','HatchAngle',45,'HatchDensity',65,'HatchLineWidth',1.5);
		r2 = rectangle('Position',[SD_length_hrs 75 12-SD_length_hrs 2]);
		r2.LineWidth = line_thick;
		
		% if WT_MvsF_logical is 1, add an asterisk above the data
		% for i=6:12
		% 	if WT_MvsF_logical(i-5)
		% 		text(i,max([Means_norm_power.SD.WT.Female(i)+SEMs_norm_power.SD.WT.Female(i),Means_norm_power.SD.Mut.Female(i)+SEMs_norm_power.SD.Mut.Female(i)])+3,'*','FontSize',30,'HorizontalAlignment','center');
		% 	end 
		% end 

		%axWTMvF=gca;
		axWTMvF = HB.Children;
		axWTMvF.XLim = [0 12];
		axWTMvF.YLim = [75 ymaxBOTH];
		axWTMvF.Box= 'off';
		axWTMvF.LineWidth=line_thick;
		axWTMvF.Title.String = ['WT ' Sexes{1} ' V ' Sexes{2} ' (2-hr bins)'];
		axWTMvF.Title.FontWeight = 'normal';
		axWTMvF.XTick = [ 0 6 12];
		axWTMvF.YTick = [75:25:ymaxBOTH];
		axWTMvF.XLabel.String = 'Hour';
		axWTMvF.YLabel.String = 'Normalized NREM Delta';
		%axF.YLabel.FontWeight = 'bold';
		axWTMvF.FontSize = 16;

		% axM.YLim = [75 ymaxBOTH];  % Set both M and F vertical axis to same value
		% axM.YTick = [75:25:ymaxBOTH];

		hold off
		l=legend([eb1 eb2],Sexes{2},Sexes{1});
		l.Location = 'northwest';
		l.Box = 'off';
		l.FontSize = 14;
		
		% - Mut M vs F -
		subplot(1,2,2)
		ymaxBOTH = max([max(Means_norm_power_2hrbins.Mut.Female + SEMs_norm_power_2hrbins.Mut.Female) max(Means_norm_power_2hrbins.Mut.Male + SEMs_norm_power_2hrbins.Mut.Male)]) + 5;
		

		% may not need these 5 lines
		% pm=plot(6:12,Means_norm_power.SD.WT.Female(6:12),'ko',6:12,Means_norm_power.SD.WT.Male(6:12),'k');
		% pm(1).LineWidth=line_thick;
		% pm(2).LineWidth=line_thick;
		% pm(1).MarkerSize = open_marker_size;
		% pm(2).MarkerSize = open_marker_size;

		hold on 
		eb3=errorbar(xvals,Means_norm_power_2hrbins.Mut.Female,SEMs_norm_power_2hrbins.Mut.Female,'o--','MarkerFaceColor','w','MarkerEdgeColor','r','Color','r');
		eb4=errorbar(xvals,Means_norm_power_2hrbins.Mut.Male,  SEMs_norm_power_2hrbins.Mut.Male,'o-','MarkerFaceColor','r','MarkerEdgeColor','r','Color','r');
		eb3.LineWidth = line_thick;
		eb4.LineWidth = line_thick;
		eb3.CapSize = 0;
		eb4.CapSize = 0;
		yline(100,'--','LineWidth',3);
		% r1 = rectangle('Position',[0 75 5 2]);
		% r1.FaceColor = 'red';
		% r1.LineWidth = 1;
		p1 = patch([0 0 SD_length_hrs SD_length_hrs],[75 77 77 75],[1 1 1]);
		hatchfill2(p1,'single','HatchAngle',45,'HatchDensity',65,'HatchLineWidth',1.5);
		r2 = rectangle('Position',[SD_length_hrs 75 12-SD_length_hrs 2]);
		r2.LineWidth = line_thick;
		
		% if WT_MvsF_logical is 1, add an asterisk above the data
		% for i=6:12
		% 	if WT_MvsF_logical(i-5)
		% 		text(i,max([Means_norm_power.SD.WT.Female(i)+SEMs_norm_power.SD.WT.Female(i),Means_norm_power.SD.Mut.Female(i)+SEMs_norm_power.SD.Mut.Female(i)])+3,'*','FontSize',30,'HorizontalAlignment','center');
		% 	end 
		% end 

		%axMutMvF=gca;
		axMutMvF = HB.Children(1);
		axMutMvF.XLim = [0 12];
		axMutMvF.YLim = [75 ymaxBOTH];
		axMutMvF.Box= 'off';
		axMutMvF.LineWidth=line_thick;
		axMutMvF.Title.String = ['Mut ' Sexes{1} ' V ' Sexes{2} ' (2-hr bins)'];
		axMutMvF.Title.FontWeight = 'normal';
		axMutMvF.XTick = [ 0 6 12];
		axMutMvF.YTick = [75:25:ymaxBOTH];
		axMutMvF.XLabel.String = 'Hour';
		axMutMvF.YLabel.String = 'Normalized NREM Delta';
		%axF.YLabel.FontWeight = 'bold';
		axMutMvF.FontSize = 16;
		hold off

		l=legend([eb3 eb4],Sexes{2},Sexes{1});
		l.Location = 'northwest';
		l.Box = 'off';
		l.FontSize = 14;



	end % end if female data present 


else  % BinSize choice
	error('You entered a BinSizeHrs that was not 1 or 2.  We currently only have code for 1 or 2.')
end 

% ---- End of making the SD figure ---------------------------------------------------
% ---------------------------------------------------------------------------------


if Plot_BL_too

	% --- Make the BL figure -------------------------------------------------------------------------------------------------
	% ------------------------------------------------------------------------------------------------------------------------
	line_thick = 2;
	filled_marker_size = 18;
	open_marker_size   = 6;
	ymaxM = max([max(Means_norm_power.BL.WT.Male(1:12) + SEMs_norm_power.BL.WT.Male(1:12)) max(Means_norm_power.BL.Mut.Male(1:12) + SEMs_norm_power.BL.Mut.Male(1:12))]) + 5;

	if FemaleDataPresent
		HBL=figure('Renderer', 'painters', 'Position', [600 400 650 320]);
		HBL.Name = 'NREM Power during Baseline LP';
		HBL.Renderer = 'painters';
	else
		HBL=figure('Renderer', 'painters', 'Position', [600 400 380 320]);
		HBL.Name = 'NREM Power during Baseline LP';
		HBL.Renderer = 'painters';
	end 

	
	if FemaleDataPresent
		subplot(1,2,1)    				% Male panel
	end  
	pm=plot(1:12,Means_norm_power.BL.WT.Male(1:12),'k.',1:12,Means_norm_power.BL.Mut.Male(1:12),'r.');
	pm(1).LineWidth=line_thick;
	pm(2).LineWidth=line_thick;
	pm(1).MarkerSize = filled_marker_size;
	pm(2).MarkerSize = filled_marker_size;
	hold on 
	eb1=errorbar(1:12,Means_norm_power.BL.WT.Male(1:12),SEMs_norm_power.BL.WT.Male(1:12),'k');
	eb2=errorbar(1:12,Means_norm_power.BL.Mut.Male(1:12),SEMs_norm_power.BL.Mut.Male(1:12),'r');
	eb1.LineWidth = line_thick;
	eb2.LineWidth = line_thick;
	eb1.CapSize = 0;
	eb2.CapSize = 0;
	yline(100,'--','LineWidth',3);
	% r1 = rectangle('Position',[0 75 5 2]);
	% r1.FaceColor = 'red';
	% r1.LineWidth = 1;
	r2 = rectangle('Position',[0 75 12 2]);
	r2.LineWidth = line_thick;
	% if Male_WTvsMut_logical is 1, add an asterisk above the data
	% Unfinished for BL data
	% for i=6:12
	% 	if Male_WTvsMut_logical(i-5)
	% 		text(i,max([Means_norm_power.BL.WT.Male(i)+SEMs_norm_power.BL.WT.Male(i),Means_norm_power.BL.Mut.Male(i)+SEMs_norm_power.BL.Mut.Male(i)])+3,'*','FontSize',30,'HorizontalAlignment','center');
	% 	end 
	% end 
	hold off
	%axM=gca;
	axM = HBL.Children;
	axM.XLim = [0 12];
	axM.YLim = [75 ymaxM];
	axM.Box= 'off';
	axM.LineWidth=line_thick;
	axM.Title.String = Sexes{1};
	axM.Title.FontWeight = 'normal';
	axM.XTick = [ 0 6 12];
	axM.YTick = [75:25:ymaxM];
	axM.XLabel.String = 'Hour';
	axM.YLabel.String = 'Normalized NREM Delta';
	axM.YLabel.FontWeight = 'bold';
	axM.YLabel.FontSize   = 14;
	axM.FontSize = 16;
	text(0.05, 0.95,{["NOTE: no stats have been done for these data,"], ["so no * or # appear automatically."]}, Units="normalized",VerticalAlignment="top", HorizontalAlignment="left")


	 
	if FemaleDataPresent
		subplot(1,2,2) 				% Female panel
		
		ymaxF = max([max(Means_norm_power.BL.WT.Female(1:12) + SEMs_norm_power.BL.WT.Female(1:12)) max(Means_norm_power.BL.Mut.Female(1:12) + SEMs_norm_power.BL.Mut.Female(1:12))]) + 5;
		ymaxBOTH = max([ymaxM ymaxF]);

		
		pm=plot(1:12,Means_norm_power.BL.WT.Female(1:12),'ko',1:12,Means_norm_power.BL.Mut.Female(1:12),'ro');
		pm(1).LineWidth=line_thick;
		pm(2).LineWidth=line_thick;
		pm(1).MarkerSize = open_marker_size;
		pm(2).MarkerSize = open_marker_size;

		hold on 
		eb1=errorbar(1:12,Means_norm_power.BL.WT.Female(1:12),SEMs_norm_power.BL.WT.Female(1:12),'o-','MarkerFaceColor','w','MarkerEdgeColor','k','Color','k');
		eb2=errorbar(1:12,Means_norm_power.BL.Mut.Female(1:12),SEMs_norm_power.BL.Mut.Female(1:12),'o-','MarkerFaceColor','w','MarkerEdgeColor','r','Color','r');
		eb1.LineWidth = line_thick;
		eb2.LineWidth = line_thick;
		eb1.CapSize = 0;
		eb2.CapSize = 0;
		yline(100,'--','LineWidth',3);
		% r1 = rectangle('Position',[0 75 5 2]);
		% r1.FaceColor = 'red';
		% r1.LineWidth = 1;
		p1 = patch([0 0 SD_length_hrs SD_length_hrs],[75 77 77 75],[1 1 1]);
		hatchfill2(p1,'single','HatchAngle',45,'HatchDensity',65,'HatchLineWidth',1.5);
		r2 = rectangle('Position',[SD_length_hrs 75 12-SD_length_hrs 2]);
		r2.LineWidth = line_thick;
		
		% if Female_WTvsMut_logical is 1, add an asterisk above the data
		% for i=6:12
		% 	if Female_WTvsMut_logical(i-5)
		% 		text(i,max([Means_norm_power.BL.WT.Female(i)+SEMs_norm_power.BL.WT.Female(i),Means_norm_power.BL.Mut.Female(i)+SEMs_norm_power.BL.Mut.Female(i)])+3,'*','FontSize',30,'HorizontalAlignment','center');
		% 	end 
		% end 

		%axF=gca;
		axF = HBL.Children(1);
		axF.XLim = [0 12];
		axF.YLim = [75 ymaxBOTH];
		axF.Box= 'off';
		axF.LineWidth=line_thick;
		axF.Title.String = Sexes{2};
		axF.Title.FontWeight = 'normal';
		axF.XTick = [ 0 6 12];
		axF.YTick = [75:25:ymaxBOTH];
		axF.XLabel.String = 'Hour';
		axF.YLabel.String = 'Normalized NREM Delta';
	%axF.YLabel.FontWeight = 'bold';
		axF.FontSize = 16;

		% axM.YLim = [75 ymaxBOTH];  % Set both M and F vertical axis to same value
		% axM.YTick = [75:25:ymaxBOTH];
		text(0.05, 0.95,'NOTE: no stats have been done for these data, \\ so no * or # appear automatically.', Units="normalized")

		hold off

		subplot(1,2,1)  % go back to the first subplot to set the axis limits
		axM.YLim = [75 ymaxBOTH];  % Set both M and F vertical axis to same value
		axM.YTick = [75:25:ymaxBOTH];

		%ax.YLabel.String = 'Normalized NREM Delta';
		% ax.YLabel.FontWeight = 'bold';
		% ax.YLabel.FontSize   = 14;
	end % end if female data present 


	% ---- End of making the BL figure ---------------------------------------------------
	% ---------------------------------------------------------------------------------
end  % end of if Plot_BL_too


% return the ANOVA tables
if BinSizeHrs == 1
	RM_ANOVA_tables.Male.WTvsMut   = RManovatbl_SleepPress_Male_WTVsMut;
	RM_ANOVA_tables.Female.WTvsMut = RManovatbl_SleepPress_Female_WTVsMut;
	RM_ANOVA_tables.MvsF.WTvsMut   = RManovatbl_SleepPress_MaleVsFemale_WTVsMut; 
	RM_ANOVA_tables.WT.MvsF        = RManovatbl_SleepPress_WT_MVsF;
	RM_ANOVA_tables.Mut.MvsF       = RManovatbl_SleepPress_Mut_MVsF;

	if exist('postHoc_SleepPress_Male_WTVsMut')   PostHocs_Recovery.Male.WTvsMut   = postHoc_SleepPress_Male_WTVsMut;   else PostHocs_Recovery.Male.WTvsMut = [];   end
	if exist('postHoc_SleepPress_Female_WTVsMut') PostHocs_Recovery.Female.WTvsMut = postHoc_SleepPress_Female_WTVsMut; else PostHocs_Recovery.Female.WTvsMut = []; end



elseif BinSizeHrs == 2
	RM_ANOVA_tables.TwoHrBins.Male.WTvsMut   = RManovatbl_SleepPress_2hrBins_Male_WTVsMut;
	RM_ANOVA_tables.TwoHrBins.Female.WTvsMut = RManovatbl_SleepPress_2hrBins_Female_WTVsMut;
	RM_ANOVA_tables.TwoHrBins.MvsF.WTvsMut   = RManovatbl_SleepPress_2hrBins_MaleVsFemale_WTVsMut; 
	RM_ANOVA_tables.TwoHrBins.WT.MvsF        = RManovatbl_SleepPress_2hrBins_WT_MVsF;
	RM_ANOVA_tables.TwoHrBins.Mut.MvsF       = RManovatbl_SleepPress_2hrBins_Mut_MVsF;

	if exist('postHoc_SleepPress_2hrBins_Male_WTVsMut')   PostHocs_Recovery.TwoHrBins.Male.WTvsMut   = postHoc_SleepPress_2hrBins_Male_WTVsMut;  else PostHocs_Recovery.TwoHrBins.Male.WTvsMut = []; end
	if exist('postHoc_SleepPress_2hrBins_Female_WTVsMut') PostHocs_Recovery.TwoHrBins.Female.WTvsMut = postHoc_SleepPress_2hrBins_Female_WTVsMut;else PostHocs_Recovery.TwoHrBins.Female.WTvsMut = []; end


else
	error('You entered a BinSizeHrs other than 1 or 2.  We only have code for 1 or 2.')
end 

% return the post-hoc tables
%RM_ANOVA_PostHoc_tables.Male.WTvsMut = 























