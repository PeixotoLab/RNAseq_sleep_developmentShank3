function [pvals] = spectral_analysis(options)
%
% USAGE: [pvals] = spectral_analysis(EEG_struct,SD_length_hrs,LegendLabels,EEG_bin_edges,EEGLowerLimit_Hz,Normalization,SeparateSpectralIntoLPDP,Sexes)
%
% where EEG_struct has the following field structure:  EEG_struct.WT.Male.BL.Wake
% where WT could be Mut, Male could be Female, BL could be SD, and Wake could be NREM or REM 
%
% This function produces figure panels plotting normalized EEG spectral power plotted with 
% frequency on the horizontal axis and % of total power on the vertical axis.  
%


arguments
	options.EEG_struct				 = [];
	options.SD_length_hrs			 = [];
	options.LegendLabels			 = []
	options.EEG_bin_edges			 = [];
	options.EEGLowerLimit_Hz		 = [];
	options.Normalization 			 = [];
	options.SeparateSpectralIntoLPDP = [];
	options.Sexes					 = [];
end 

EEG_struct 				 = options.EEG_struct;
SD_length_hrs 			 = options.SD_length_hrs;
LegendLabels 			 = options.LegendLabels;
EEG_bin_edges 			 = options.EEG_bin_edges;
EEGLowerLimit_Hz 		 = options.EEGLowerLimit_Hz;
Normalization  			 = options.Normalization;
SeparateSpectralIntoLPDP = options.SeparateSpectralIntoLPDP;
Sexes					 = options.Sexes;

% Do you have female data?
if isempty(EEG_struct.WT.(Sexes{2}).BL.Wake)
	FemaleDataPresent = 0;
else 
	FemaleDataPresent = 1;
end 


dx = unique(round(EEG_bin_edges(2:end)-EEG_bin_edges(1:end-1),4));  % the EEG bin size in Hz
if length(dx)>1 error('You have uneven EEG frequency bins.  Try again.');  end  


% First unpack the data structure into matrices (it's too confusing leaving it in structures)
% Male WT
EEG_Male_WT_BL_Wake = EEG_struct.WT.(Sexes{1}).BL.Wake; 
EEG_Male_WT_BL_NREM = EEG_struct.WT.(Sexes{1}).BL.NREM; 
EEG_Male_WT_BL_REM  = EEG_struct.WT.(Sexes{1}).BL.REM; 
EEG_Male_WT_SD_Wake = EEG_struct.WT.(Sexes{1}).SD.Wake(SD_length_hrs+1:24,:,:); 		% Remove the actual SD 
EEG_Male_WT_SD_NREM = EEG_struct.WT.(Sexes{1}).SD.NREM(SD_length_hrs+1:24,:,:); 
EEG_Male_WT_SD_REM  = EEG_struct.WT.(Sexes{1}).SD.REM(SD_length_hrs+1:24,:,:); 

% Male Mutants
EEG_Male_Mut_BL_Wake = EEG_struct.Mut.(Sexes{1}).BL.Wake; 
EEG_Male_Mut_BL_NREM = EEG_struct.Mut.(Sexes{1}).BL.NREM; 
EEG_Male_Mut_BL_REM  = EEG_struct.Mut.(Sexes{1}).BL.REM; 
EEG_Male_Mut_SD_Wake = EEG_struct.Mut.(Sexes{1}).SD.Wake(SD_length_hrs+1:24,:,:); 	% Remove the actual SD
EEG_Male_Mut_SD_NREM = EEG_struct.Mut.(Sexes{1}).SD.NREM(SD_length_hrs+1:24,:,:); 
EEG_Male_Mut_SD_REM  = EEG_struct.Mut.(Sexes{1}).SD.REM(SD_length_hrs+1:24,:,:); 

if FemaleDataPresent
	% Female WT
	EEG_Female_WT_BL_Wake = EEG_struct.WT.(Sexes{2}).BL.Wake; 
	EEG_Female_WT_BL_NREM = EEG_struct.WT.(Sexes{2}).BL.NREM; 
	EEG_Female_WT_BL_REM  = EEG_struct.WT.(Sexes{2}).BL.REM; 
	EEG_Female_WT_SD_Wake = EEG_struct.WT.(Sexes{2}).SD.Wake(SD_length_hrs+1:24,:,:);	% Remove the actual SD 
	EEG_Female_WT_SD_NREM = EEG_struct.WT.(Sexes{2}).SD.NREM(SD_length_hrs+1:24,:,:); 
	EEG_Female_WT_SD_REM  = EEG_struct.WT.(Sexes{2}).SD.REM(SD_length_hrs+1:24,:,:); 

	% Female Mutants
	EEG_Female_Mut_BL_Wake = EEG_struct.Mut.(Sexes{2}).BL.Wake; 
	EEG_Female_Mut_BL_NREM = EEG_struct.Mut.(Sexes{2}).BL.NREM; 
	EEG_Female_Mut_BL_REM  = EEG_struct.Mut.(Sexes{2}).BL.REM; 
	EEG_Female_Mut_SD_Wake = EEG_struct.Mut.(Sexes{2}).SD.Wake(SD_length_hrs+1:24,:,:); % Remove the actual SD
	EEG_Female_Mut_SD_NREM = EEG_struct.Mut.(Sexes{2}).SD.NREM(SD_length_hrs+1:24,:,:); 
	EEG_Female_Mut_SD_REM  = EEG_struct.Mut.(Sexes{2}).SD.REM(SD_length_hrs+1:24,:,:); 
end 


% First normalize by the "total state specific power across 24 hours of baseline" or sleepdep per animal (quotes from Lizzy's poster)
% perhaps put this in a function.  It's messy.  EEG_Normalized = normalize_EEG(EEG_struct) or something.  
% First compute the average per animal
% UPDATE 2.24:  After talking with Peter Achermann, normalize by area under the curve so the 
% area under the normalized curve is equal to 1.  
% These all used to be squeeze(mean(mean(EEG_Male_WT_BL_Wake,'omitnan'))); 
% UPDATE 3.28:  From Frank lab meeting: make the default to normalize to the power across all states (in BL?)
% 				and the option to normalize so area under the curve is 1.  

if strcmp(Normalization,'MeanPowerAllStates')
	All_States_ByAnimal.WT.BL.M  = [EEG_Male_WT_BL_Wake; EEG_Male_WT_BL_NREM; EEG_Male_WT_BL_REM];
	All_States_ByAnimal.WT.SD.M  = [EEG_Male_WT_SD_Wake; EEG_Male_WT_SD_NREM; EEG_Male_WT_SD_REM];
	All_States_ByAnimal.Mut.BL.M = [EEG_Male_Mut_BL_Wake; EEG_Male_Mut_BL_NREM; EEG_Male_Mut_BL_REM];
	All_States_ByAnimal.Mut.SD.M = [EEG_Male_Mut_SD_Wake; EEG_Male_Mut_SD_NREM; EEG_Male_Mut_SD_REM];

	% Male WT
	norm_by_animal_M_WT_BL_Wake = squeeze(mean(mean(All_States_ByAnimal.WT.BL.M,1,'omitnan'),2,'omitnan'));
	norm_by_animal_M_WT_BL_NREM = squeeze(mean(mean(All_States_ByAnimal.WT.BL.M,1,'omitnan'),2,'omitnan'));
	norm_by_animal_M_WT_BL_REM  = squeeze(mean(mean(All_States_ByAnimal.WT.BL.M,1,'omitnan'),2,'omitnan'));
	norm_by_animal_M_WT_SD_Wake = squeeze(mean(mean(All_States_ByAnimal.WT.SD.M,1,'omitnan'),2,'omitnan'));
	norm_by_animal_M_WT_SD_NREM = squeeze(mean(mean(All_States_ByAnimal.WT.SD.M,1,'omitnan'),2,'omitnan'));
	norm_by_animal_M_WT_SD_REM  = squeeze(mean(mean(All_States_ByAnimal.WT.SD.M,1,'omitnan'),2,'omitnan'));
	
	% Male Mut
	norm_by_animal_M_Mut_BL_Wake = squeeze(mean(mean(All_States_ByAnimal.Mut.BL.M,1,'omitnan'),2,'omitnan'));
	norm_by_animal_M_Mut_BL_NREM = squeeze(mean(mean(All_States_ByAnimal.Mut.BL.M,1,'omitnan'),2,'omitnan'));
	norm_by_animal_M_Mut_BL_REM  = squeeze(mean(mean(All_States_ByAnimal.Mut.BL.M,1,'omitnan'),2,'omitnan'));
	norm_by_animal_M_Mut_SD_Wake = squeeze(mean(mean(All_States_ByAnimal.Mut.SD.M,1,'omitnan'),2,'omitnan'));
	norm_by_animal_M_Mut_SD_NREM = squeeze(mean(mean(All_States_ByAnimal.Mut.SD.M,1,'omitnan'),2,'omitnan'));
	norm_by_animal_M_Mut_SD_REM  = squeeze(mean(mean(All_States_ByAnimal.Mut.SD.M,1,'omitnan'),2,'omitnan'));

	if FemaleDataPresent
		All_States_ByAnimal.WT.BL.F  = [EEG_Female_WT_BL_Wake;  EEG_Female_WT_BL_NREM;  EEG_Female_WT_BL_REM];
		All_States_ByAnimal.WT.SD.F  = [EEG_Female_WT_SD_Wake;  EEG_Female_WT_SD_NREM;  EEG_Female_WT_SD_REM];
		All_States_ByAnimal.Mut.BL.F = [EEG_Female_Mut_BL_Wake; EEG_Female_Mut_BL_NREM; EEG_Female_Mut_BL_REM];
		All_States_ByAnimal.Mut.SD.F = [EEG_Female_Mut_SD_Wake; EEG_Female_Mut_SD_NREM; EEG_Female_Mut_SD_REM];

		% Female WT
		norm_by_animal_F_WT_BL_Wake = squeeze(mean(mean(All_States_ByAnimal.WT.BL.F,1,'omitnan'),2,'omitnan'));
		norm_by_animal_F_WT_BL_NREM = squeeze(mean(mean(All_States_ByAnimal.WT.BL.F,1,'omitnan'),2,'omitnan'));
		norm_by_animal_F_WT_BL_REM  = squeeze(mean(mean(All_States_ByAnimal.WT.BL.F,1,'omitnan'),2,'omitnan'));
		norm_by_animal_F_WT_SD_Wake = squeeze(mean(mean(All_States_ByAnimal.WT.SD.F,1,'omitnan'),2,'omitnan'));
		norm_by_animal_F_WT_SD_NREM = squeeze(mean(mean(All_States_ByAnimal.WT.SD.F,1,'omitnan'),2,'omitnan'));
		norm_by_animal_F_WT_SD_REM  = squeeze(mean(mean(All_States_ByAnimal.WT.SD.F,1,'omitnan'),2,'omitnan'));
		
	% Female Mut
		norm_by_animal_F_Mut_BL_Wake = squeeze(mean(mean(All_States_ByAnimal.Mut.BL.F,1,'omitnan'),2,'omitnan'));
		norm_by_animal_F_Mut_BL_NREM = squeeze(mean(mean(All_States_ByAnimal.Mut.BL.F,1,'omitnan'),2,'omitnan'));
		norm_by_animal_F_Mut_BL_REM  = squeeze(mean(mean(All_States_ByAnimal.Mut.BL.F,1,'omitnan'),2,'omitnan'));
		norm_by_animal_F_Mut_SD_Wake = squeeze(mean(mean(All_States_ByAnimal.Mut.SD.F,1,'omitnan'),2,'omitnan'));
		norm_by_animal_F_Mut_SD_NREM = squeeze(mean(mean(All_States_ByAnimal.Mut.SD.F,1,'omitnan'),2,'omitnan'));
		norm_by_animal_F_Mut_SD_REM  = squeeze(mean(mean(All_States_ByAnimal.Mut.SD.F,1,'omitnan'),2,'omitnan'));
	end

elseif strcmp(Normalization,'AreaUnderCurve')
	% Male WT
	norm_by_animal_M_WT_BL_Wake = calc_normalization_factor(EEG_Male_WT_BL_Wake,dx);
	norm_by_animal_M_WT_BL_NREM = calc_normalization_factor(EEG_Male_WT_BL_NREM,dx);
	norm_by_animal_M_WT_BL_REM  = calc_normalization_factor(EEG_Male_WT_BL_REM,dx);
	norm_by_animal_M_WT_SD_Wake = calc_normalization_factor(EEG_Male_WT_SD_Wake,dx);
	norm_by_animal_M_WT_SD_NREM = calc_normalization_factor(EEG_Male_WT_SD_NREM,dx);
	norm_by_animal_M_WT_SD_REM  = calc_normalization_factor(EEG_Male_WT_SD_REM,dx);

	% Male Mut
	norm_by_animal_M_Mut_BL_Wake = calc_normalization_factor(EEG_Male_Mut_BL_Wake,dx);
	norm_by_animal_M_Mut_BL_NREM = calc_normalization_factor(EEG_Male_Mut_BL_NREM,dx);
	norm_by_animal_M_Mut_BL_REM  = calc_normalization_factor(EEG_Male_Mut_BL_REM,dx);
	norm_by_animal_M_Mut_SD_Wake = calc_normalization_factor(EEG_Male_Mut_SD_Wake,dx);
	norm_by_animal_M_Mut_SD_NREM = calc_normalization_factor(EEG_Male_Mut_SD_NREM,dx);
	norm_by_animal_M_Mut_SD_REM  = calc_normalization_factor(EEG_Male_Mut_SD_REM,dx);

	if FemaleDataPresent
		% Female WT
		norm_by_animal_F_WT_BL_Wake  = calc_normalization_factor(EEG_Female_WT_BL_Wake,dx);
		norm_by_animal_F_WT_BL_NREM  = calc_normalization_factor(EEG_Female_WT_BL_NREM,dx);
		norm_by_animal_F_WT_BL_REM   = calc_normalization_factor(EEG_Female_WT_BL_REM,dx);
		norm_by_animal_F_WT_SD_Wake  = calc_normalization_factor(EEG_Female_WT_SD_Wake,dx);
		norm_by_animal_F_WT_SD_NREM  = calc_normalization_factor(EEG_Female_WT_SD_NREM,dx);
		norm_by_animal_F_WT_SD_REM   = calc_normalization_factor(EEG_Female_WT_SD_REM,dx);

		% Female Mut
		norm_by_animal_F_Mut_BL_Wake = calc_normalization_factor(EEG_Female_Mut_BL_Wake,dx);
		norm_by_animal_F_Mut_BL_NREM = calc_normalization_factor(EEG_Female_Mut_BL_NREM,dx);
		norm_by_animal_F_Mut_BL_REM  = calc_normalization_factor(EEG_Female_Mut_BL_REM,dx);
		norm_by_animal_F_Mut_SD_Wake = calc_normalization_factor(EEG_Female_Mut_SD_Wake,dx);
		norm_by_animal_F_Mut_SD_NREM = calc_normalization_factor(EEG_Female_Mut_SD_NREM,dx);
		norm_by_animal_F_Mut_SD_REM  = calc_normalization_factor(EEG_Female_Mut_SD_REM,dx);
	end 

else
	error(['Invalid choice for Normalization.  Options are MeanPowerAllStates or AreaUnderCurve. You chose ',Normalization])
end  % end of choosing normalization


% Now reshape so we can use in bsxfun(@rdivide  ...)
% Male WT 
norm_matrix_Male_WT_BL_Wake = reshape(norm_by_animal_M_WT_BL_Wake,1,1,[]);   
norm_matrix_Male_WT_BL_NREM = reshape(norm_by_animal_M_WT_BL_NREM,1,1,[]);   
norm_matrix_Male_WT_BL_REM  = reshape(norm_by_animal_M_WT_BL_REM,1,1,[]);   
norm_matrix_Male_WT_SD_Wake = reshape(norm_by_animal_M_WT_SD_Wake,1,1,[]);   
norm_matrix_Male_WT_SD_NREM = reshape(norm_by_animal_M_WT_SD_NREM,1,1,[]);   
norm_matrix_Male_WT_SD_REM  = reshape(norm_by_animal_M_WT_SD_REM,1,1,[]);   

% Male Mut 
norm_matrix_Male_Mut_BL_Wake = reshape(norm_by_animal_M_Mut_BL_Wake,1,1,[]);   
norm_matrix_Male_Mut_BL_NREM = reshape(norm_by_animal_M_Mut_BL_NREM,1,1,[]);   
norm_matrix_Male_Mut_BL_REM  = reshape(norm_by_animal_M_Mut_BL_REM,1,1,[]);   
norm_matrix_Male_Mut_SD_Wake = reshape(norm_by_animal_M_Mut_SD_Wake,1,1,[]);   
norm_matrix_Male_Mut_SD_NREM = reshape(norm_by_animal_M_Mut_SD_NREM,1,1,[]);   
norm_matrix_Male_Mut_SD_REM  = reshape(norm_by_animal_M_Mut_SD_REM,1,1,[]);  

if FemaleDataPresent
	% Female WT 
	norm_matrix_Female_WT_BL_Wake = reshape(norm_by_animal_F_WT_BL_Wake,1,1,[]);   
	norm_matrix_Female_WT_BL_NREM = reshape(norm_by_animal_F_WT_BL_NREM,1,1,[]);   
	norm_matrix_Female_WT_BL_REM  = reshape(norm_by_animal_F_WT_BL_REM,1,1,[]);   
	norm_matrix_Female_WT_SD_Wake = reshape(norm_by_animal_F_WT_SD_Wake,1,1,[]);   
	norm_matrix_Female_WT_SD_NREM = reshape(norm_by_animal_F_WT_SD_NREM,1,1,[]);   
	norm_matrix_Female_WT_SD_REM  = reshape(norm_by_animal_F_WT_SD_REM,1,1,[]);   

	% Female Mut 
	norm_matrix_Female_Mut_BL_Wake = reshape(norm_by_animal_F_Mut_BL_Wake,1,1,[]);   
	norm_matrix_Female_Mut_BL_NREM = reshape(norm_by_animal_F_Mut_BL_NREM,1,1,[]);   
	norm_matrix_Female_Mut_BL_REM  = reshape(norm_by_animal_F_Mut_BL_REM,1,1,[]);   
	norm_matrix_Female_Mut_SD_Wake = reshape(norm_by_animal_F_Mut_SD_Wake,1,1,[]);   
	norm_matrix_Female_Mut_SD_NREM = reshape(norm_by_animal_F_Mut_SD_NREM,1,1,[]);   
	norm_matrix_Female_Mut_SD_REM  = reshape(norm_by_animal_F_Mut_SD_REM,1,1,[]); 
end 

% Next normalize the EEG 
% Male WT
EEG_Normalized_Male_WT_BL_Wake = bsxfun(@rdivide,EEG_Male_WT_BL_Wake,norm_matrix_Male_WT_BL_Wake);  % this is the normalized EEG
EEG_Normalized_Male_WT_BL_NREM = bsxfun(@rdivide,EEG_Male_WT_BL_NREM,norm_matrix_Male_WT_BL_NREM);  % this is the normalized EEG
EEG_Normalized_Male_WT_BL_REM  = bsxfun(@rdivide,EEG_Male_WT_BL_REM, norm_matrix_Male_WT_BL_REM);   % this is the normalized EEG
EEG_Normalized_Male_WT_SD_Wake = bsxfun(@rdivide,EEG_Male_WT_SD_Wake,norm_matrix_Male_WT_SD_Wake);  % this is the normalized EEG
EEG_Normalized_Male_WT_SD_NREM = bsxfun(@rdivide,EEG_Male_WT_SD_NREM,norm_matrix_Male_WT_SD_NREM);  % this is the normalized EEG
EEG_Normalized_Male_WT_SD_REM  = bsxfun(@rdivide,EEG_Male_WT_SD_REM, norm_matrix_Male_WT_SD_REM);   % this is the normalized EEG

% Male Mut
EEG_Normalized_Male_Mut_BL_Wake = bsxfun(@rdivide,EEG_Male_Mut_BL_Wake,norm_matrix_Male_Mut_BL_Wake);  % this is the normalized EEG
EEG_Normalized_Male_Mut_BL_NREM = bsxfun(@rdivide,EEG_Male_Mut_BL_NREM,norm_matrix_Male_Mut_BL_NREM);  % this is the normalized EEG
EEG_Normalized_Male_Mut_BL_REM  = bsxfun(@rdivide,EEG_Male_Mut_BL_REM, norm_matrix_Male_Mut_BL_REM);   % this is the normalized EEG
EEG_Normalized_Male_Mut_SD_Wake = bsxfun(@rdivide,EEG_Male_Mut_SD_Wake,norm_matrix_Male_Mut_SD_Wake);  % this is the normalized EEG
EEG_Normalized_Male_Mut_SD_NREM = bsxfun(@rdivide,EEG_Male_Mut_SD_NREM,norm_matrix_Male_Mut_SD_NREM);  % this is the normalized EEG
EEG_Normalized_Male_Mut_SD_REM  = bsxfun(@rdivide,EEG_Male_Mut_SD_REM, norm_matrix_Male_Mut_SD_REM);   % this is the normalized EEG

if FemaleDataPresent
	% Female WT
	EEG_Normalized_Female_WT_BL_Wake = bsxfun(@rdivide,EEG_Female_WT_BL_Wake,norm_matrix_Female_WT_BL_Wake);  % this is the normalized EEG
	EEG_Normalized_Female_WT_BL_NREM = bsxfun(@rdivide,EEG_Female_WT_BL_NREM,norm_matrix_Female_WT_BL_NREM);  % this is the normalized EEG
	EEG_Normalized_Female_WT_BL_REM  = bsxfun(@rdivide,EEG_Female_WT_BL_REM, norm_matrix_Female_WT_BL_REM);   % this is the normalized EEG
	EEG_Normalized_Female_WT_SD_Wake = bsxfun(@rdivide,EEG_Female_WT_SD_Wake,norm_matrix_Female_WT_SD_Wake);  % this is the normalized EEG
	EEG_Normalized_Female_WT_SD_NREM = bsxfun(@rdivide,EEG_Female_WT_SD_NREM,norm_matrix_Female_WT_SD_NREM);  % this is the normalized EEG
	EEG_Normalized_Female_WT_SD_REM  = bsxfun(@rdivide,EEG_Female_WT_SD_REM, norm_matrix_Female_WT_SD_REM);   % this is the normalized EEG

	% Female Mut
	EEG_Normalized_Female_Mut_BL_Wake = bsxfun(@rdivide,EEG_Female_Mut_BL_Wake,norm_matrix_Female_Mut_BL_Wake);  % this is the normalized EEG
	EEG_Normalized_Female_Mut_BL_NREM = bsxfun(@rdivide,EEG_Female_Mut_BL_NREM,norm_matrix_Female_Mut_BL_NREM);  % this is the normalized EEG
	EEG_Normalized_Female_Mut_BL_REM  = bsxfun(@rdivide,EEG_Female_Mut_BL_REM, norm_matrix_Female_Mut_BL_REM);   % this is the normalized EEG
	EEG_Normalized_Female_Mut_SD_Wake = bsxfun(@rdivide,EEG_Female_Mut_SD_Wake,norm_matrix_Female_Mut_SD_Wake);  % this is the normalized EEG
	EEG_Normalized_Female_Mut_SD_NREM = bsxfun(@rdivide,EEG_Female_Mut_SD_NREM,norm_matrix_Female_Mut_SD_NREM);  % this is the normalized EEG
	EEG_Normalized_Female_Mut_SD_REM  = bsxfun(@rdivide,EEG_Female_Mut_SD_REM, norm_matrix_Female_Mut_SD_REM);   % this is the normalized EEG
end 

% Finally, compute the means (and sems) so we can make the plots
% --- Means ---
if strcmp(Normalization,'AreaUnderCurve')
	NormFactor = 1;
elseif strcmp(Normalization,'MeanPowerAllStates')
	NormFactor = 100;
else
	error('You chose an invalid normalization factor.  Please choose AreadUnderCurve or MeanPowerAllStates.')
end 

% Male WT  			was 100*mean(mean .... 
Means.(Sexes{1}).WT.BL.Wake = NormFactor*mean(mean(EEG_Normalized_Male_WT_BL_Wake,1,'omitnan'),3,'omitnan');
Means.(Sexes{1}).WT.BL.NREM = NormFactor*mean(mean(EEG_Normalized_Male_WT_BL_NREM,1,'omitnan'),3,'omitnan');
Means.(Sexes{1}).WT.BL.REM  = NormFactor*mean(mean(EEG_Normalized_Male_WT_BL_REM,1,'omitnan'),3,'omitnan');
Means.(Sexes{1}).WT.SD.Wake = NormFactor*mean(mean(EEG_Normalized_Male_WT_SD_Wake,1,'omitnan'),3,'omitnan');
Means.(Sexes{1}).WT.SD.NREM = NormFactor*mean(mean(EEG_Normalized_Male_WT_SD_NREM,1,'omitnan'),3,'omitnan');
Means.(Sexes{1}).WT.SD.REM  = NormFactor*mean(mean(EEG_Normalized_Male_WT_SD_REM,1,'omitnan'),3,'omitnan');

% Male Mut
Means.(Sexes{1}).Mut.BL.Wake = NormFactor*mean(mean(EEG_Normalized_Male_Mut_BL_Wake,1,'omitnan'),3,'omitnan');
Means.(Sexes{1}).Mut.BL.NREM = NormFactor*mean(mean(EEG_Normalized_Male_Mut_BL_NREM,1,'omitnan'),3,'omitnan');
Means.(Sexes{1}).Mut.BL.REM  = NormFactor*mean(mean(EEG_Normalized_Male_Mut_BL_REM,1,'omitnan'),3,'omitnan');
Means.(Sexes{1}).Mut.SD.Wake = NormFactor*mean(mean(EEG_Normalized_Male_Mut_SD_Wake,1,'omitnan'),3,'omitnan');
Means.(Sexes{1}).Mut.SD.NREM = NormFactor*mean(mean(EEG_Normalized_Male_Mut_SD_NREM,1,'omitnan'),3,'omitnan');
Means.(Sexes{1}).Mut.SD.REM  = NormFactor*mean(mean(EEG_Normalized_Male_Mut_SD_REM,1,'omitnan'),3,'omitnan');

if FemaleDataPresent
	% Female WT
	Means.(Sexes{2}).WT.BL.Wake = NormFactor*mean(mean(EEG_Normalized_Female_WT_BL_Wake,1,'omitnan'),3,'omitnan');
	Means.(Sexes{2}).WT.BL.NREM = NormFactor*mean(mean(EEG_Normalized_Female_WT_BL_NREM,1,'omitnan'),3,'omitnan');
	Means.(Sexes{2}).WT.BL.REM  = NormFactor*mean(mean(EEG_Normalized_Female_WT_BL_REM,1,'omitnan'),3,'omitnan');
	Means.(Sexes{2}).WT.SD.Wake = NormFactor*mean(mean(EEG_Normalized_Female_WT_SD_Wake,1,'omitnan'),3,'omitnan');
	Means.(Sexes{2}).WT.SD.NREM = NormFactor*mean(mean(EEG_Normalized_Female_WT_SD_NREM,1,'omitnan'),3,'omitnan');
	Means.(Sexes{2}).WT.SD.REM  = NormFactor*mean(mean(EEG_Normalized_Female_WT_SD_REM,1,'omitnan'),3,'omitnan');

	% Female Mut
	Means.(Sexes{2}).Mut.BL.Wake = NormFactor*mean(mean(EEG_Normalized_Female_Mut_BL_Wake,1,'omitnan'),3,'omitnan');
	Means.(Sexes{2}).Mut.BL.NREM = NormFactor*mean(mean(EEG_Normalized_Female_Mut_BL_NREM,1,'omitnan'),3,'omitnan');
	Means.(Sexes{2}).Mut.BL.REM  = NormFactor*mean(mean(EEG_Normalized_Female_Mut_BL_REM,1,'omitnan'),3,'omitnan');
	Means.(Sexes{2}).Mut.SD.Wake = NormFactor*mean(mean(EEG_Normalized_Female_Mut_SD_Wake,1,'omitnan'),3,'omitnan');
	Means.(Sexes{2}).Mut.SD.NREM = NormFactor*mean(mean(EEG_Normalized_Female_Mut_SD_NREM,1,'omitnan'),3,'omitnan');
	Means.(Sexes{2}).Mut.SD.REM  = NormFactor*mean(mean(EEG_Normalized_Female_Mut_SD_REM,1,'omitnan'),3,'omitnan');
end 

% -- split into LP+DP? --
if SeparateSpectralIntoLPDP
	Means.(Sexes{1}).WT.BLFirst12hr.Wake = NormFactor*mean(mean(EEG_Normalized_Male_WT_BL_Wake(1:12,:,:),1,'omitnan'),3,'omitnan');
	Means.(Sexes{1}).WT.BLLast12hr.Wake  = NormFactor*mean(mean(EEG_Normalized_Male_WT_BL_Wake(13:24,:,:),1,'omitnan'),3,'omitnan');
	Means.(Sexes{1}).WT.BLFirst12hr.NREM = NormFactor*mean(mean(EEG_Normalized_Male_WT_BL_NREM(1:12,:,:),1,'omitnan'),3,'omitnan');
	Means.(Sexes{1}).WT.BLLast12hr.NREM  = NormFactor*mean(mean(EEG_Normalized_Male_WT_BL_NREM(13:24,:,:),1,'omitnan'),3,'omitnan');
	Means.(Sexes{1}).WT.BLFirst12hr.REM  = NormFactor*mean(mean(EEG_Normalized_Male_WT_BL_REM(1:12,:,:),1,'omitnan'),3,'omitnan');
	Means.(Sexes{1}).WT.BLLast12hr.REM   = NormFactor*mean(mean(EEG_Normalized_Male_WT_BL_REM(13:24,:,:),1,'omitnan'),3,'omitnan');
	Means.(Sexes{1}).WT.SDFirst12hr.Wake = NormFactor*mean(mean(EEG_Normalized_Male_WT_SD_Wake(1:12-SD_length_hrs,:,:),1,'omitnan'),3,'omitnan');
	Means.(Sexes{1}).WT.SDLast12hr.Wake  = NormFactor*mean(mean(EEG_Normalized_Male_WT_SD_Wake(12-SD_length_hrs+1:24-SD_length_hrs,:,:),1,'omitnan'),3,'omitnan');
	Means.(Sexes{1}).WT.SDFirst12hr.NREM = NormFactor*mean(mean(EEG_Normalized_Male_WT_SD_NREM(1:12-SD_length_hrs,:,:),1,'omitnan'),3,'omitnan');
	Means.(Sexes{1}).WT.SDLast12hr.NREM  = NormFactor*mean(mean(EEG_Normalized_Male_WT_SD_NREM(12-SD_length_hrs+1:24-SD_length_hrs,:,:),1,'omitnan'),3,'omitnan');
	Means.(Sexes{1}).WT.SDFirst12hr.REM  = NormFactor*mean(mean(EEG_Normalized_Male_WT_SD_REM(1:12-SD_length_hrs,:,:),1,'omitnan'),3,'omitnan');
	Means.(Sexes{1}).WT.SDLast12hr.REM   = NormFactor*mean(mean(EEG_Normalized_Male_WT_SD_REM(12-SD_length_hrs+1:24-SD_length_hrs,:,:),1,'omitnan'),3,'omitnan');

	% Male Mut
	Means.(Sexes{1}).Mut.BLFirst12hr.Wake = NormFactor*mean(mean(EEG_Normalized_Male_Mut_BL_Wake(1:12,:,:),1,'omitnan'),3,'omitnan');
	Means.(Sexes{1}).Mut.BLLast12hr.Wake  = NormFactor*mean(mean(EEG_Normalized_Male_Mut_BL_Wake(13:24,:,:),1,'omitnan'),3,'omitnan');
	Means.(Sexes{1}).Mut.BLFirst12hr.NREM = NormFactor*mean(mean(EEG_Normalized_Male_Mut_BL_NREM(1:12,:,:),1,'omitnan'),3,'omitnan');
	Means.(Sexes{1}).Mut.BLLast12hr.NREM  = NormFactor*mean(mean(EEG_Normalized_Male_Mut_BL_NREM(13:24,:,:),1,'omitnan'),3,'omitnan');
	Means.(Sexes{1}).Mut.BLFirst12hr.REM  = NormFactor*mean(mean(EEG_Normalized_Male_Mut_BL_REM(1:12,:,:),1,'omitnan'),3,'omitnan');
	Means.(Sexes{1}).Mut.BLLast12hr.REM   = NormFactor*mean(mean(EEG_Normalized_Male_Mut_BL_REM(13:24,:,:),1,'omitnan'),3,'omitnan');
	Means.(Sexes{1}).Mut.SDFirst12hr.Wake = NormFactor*mean(mean(EEG_Normalized_Male_Mut_SD_Wake(1:12-SD_length_hrs,:,:),1,'omitnan'),3,'omitnan');
	Means.(Sexes{1}).Mut.SDLast12hr.Wake  = NormFactor*mean(mean(EEG_Normalized_Male_Mut_SD_Wake(12-SD_length_hrs+1:24-SD_length_hrs,:,:),1,'omitnan'),3,'omitnan');
	Means.(Sexes{1}).Mut.SDFirst12hr.NREM = NormFactor*mean(mean(EEG_Normalized_Male_Mut_SD_NREM(1:12-SD_length_hrs,:,:),1,'omitnan'),3,'omitnan');
	Means.(Sexes{1}).Mut.SDLast12hr.NREM  = NormFactor*mean(mean(EEG_Normalized_Male_Mut_SD_NREM(12-SD_length_hrs+1:24-SD_length_hrs,:,:),1,'omitnan'),3,'omitnan');
	Means.(Sexes{1}).Mut.SDFirst12hr.REM  = NormFactor*mean(mean(EEG_Normalized_Male_Mut_SD_REM(1:12-SD_length_hrs,:,:),1,'omitnan'),3,'omitnan');
	Means.(Sexes{1}).Mut.SDLast12hr.REM   = NormFactor*mean(mean(EEG_Normalized_Male_Mut_SD_REM(12-SD_length_hrs+1:24-SD_length_hrs,:,:),1,'omitnan'),3,'omitnan');

	if FemaleDataPresent
		% Female WT
		Means.(Sexes{2}).WT.BLFirst12hr.Wake = NormFactor*mean(mean(EEG_Normalized_Female_WT_BL_Wake(1:12,:,:),1,'omitnan'),3,'omitnan');
		Means.(Sexes{2}).WT.BLLast12hr.Wake  = NormFactor*mean(mean(EEG_Normalized_Female_WT_BL_Wake(13:24,:,:),1,'omitnan'),3,'omitnan');
		Means.(Sexes{2}).WT.BLFirst12hr.NREM = NormFactor*mean(mean(EEG_Normalized_Female_WT_BL_NREM(1:12,:,:),1,'omitnan'),3,'omitnan');
		Means.(Sexes{2}).WT.BLLast12hr.NREM  = NormFactor*mean(mean(EEG_Normalized_Female_WT_BL_NREM(13:24,:,:),1,'omitnan'),3,'omitnan');
		Means.(Sexes{2}).WT.BLFirst12hr.REM  = NormFactor*mean(mean(EEG_Normalized_Female_WT_BL_REM(1:12,:,:),1,'omitnan'),3,'omitnan');
		Means.(Sexes{2}).WT.BLLast12hr.REM   = NormFactor*mean(mean(EEG_Normalized_Female_WT_BL_REM(13:24,:,:),1,'omitnan'),3,'omitnan');
		Means.(Sexes{2}).WT.SDFirst12hr.Wake = NormFactor*mean(mean(EEG_Normalized_Female_WT_SD_Wake(1:12-SD_length_hrs,:,:),1,'omitnan'),3,'omitnan');
		Means.(Sexes{2}).WT.SDLast12hr.Wake  = NormFactor*mean(mean(EEG_Normalized_Female_WT_SD_Wake(12-SD_length_hrs+1:24-SD_length_hrs,:,:),1,'omitnan'),3,'omitnan');
		Means.(Sexes{2}).WT.SDFirst12hr.NREM = NormFactor*mean(mean(EEG_Normalized_Female_WT_SD_NREM(1:12-SD_length_hrs,:,:),1,'omitnan'),3,'omitnan');
		Means.(Sexes{2}).WT.SDLast12hr.NREM  = NormFactor*mean(mean(EEG_Normalized_Female_WT_SD_NREM(12-SD_length_hrs+1:24-SD_length_hrs,:,:),1,'omitnan'),3,'omitnan');
		Means.(Sexes{2}).WT.SDFirst12hr.REM  = NormFactor*mean(mean(EEG_Normalized_Female_WT_SD_REM(1:12-SD_length_hrs,:,:),1,'omitnan'),3,'omitnan');
		Means.(Sexes{2}).WT.SDLast12hr.REM   = NormFactor*mean(mean(EEG_Normalized_Female_WT_SD_REM(12-SD_length_hrs+1:24-SD_length_hrs,:,:),1,'omitnan'),3,'omitnan');

		% Female Mut
		Means.(Sexes{2}).Mut.BLFirst12hr.Wake = NormFactor*mean(mean(EEG_Normalized_Female_Mut_BL_Wake(1:12,:,:),1,'omitnan'),3,'omitnan');
		Means.(Sexes{2}).Mut.BLLast12hr.Wake  = NormFactor*mean(mean(EEG_Normalized_Female_Mut_BL_Wake(13:24,:,:),1,'omitnan'),3,'omitnan');
		Means.(Sexes{2}).Mut.BLFirst12hr.NREM = NormFactor*mean(mean(EEG_Normalized_Female_Mut_BL_NREM(1:12,:,:),1,'omitnan'),3,'omitnan');
		Means.(Sexes{2}).Mut.BLLast12hr.NREM  = NormFactor*mean(mean(EEG_Normalized_Female_Mut_BL_NREM(13:24,:,:),1,'omitnan'),3,'omitnan');
		Means.(Sexes{2}).Mut.BLFirst12hr.REM  = NormFactor*mean(mean(EEG_Normalized_Female_Mut_BL_REM(1:12,:,:),1,'omitnan'),3,'omitnan');
		Means.(Sexes{2}).Mut.BLLast12hr.REM   = NormFactor*mean(mean(EEG_Normalized_Female_Mut_BL_REM(13:24,:,:),1,'omitnan'),3,'omitnan');
		Means.(Sexes{2}).Mut.SDFirst12hr.Wake = NormFactor*mean(mean(EEG_Normalized_Female_Mut_SD_Wake(1:12-SD_length_hrs,:,:),1,'omitnan'),3,'omitnan');
		Means.(Sexes{2}).Mut.SDLast12hr.Wake  = NormFactor*mean(mean(EEG_Normalized_Female_Mut_SD_Wake(12-SD_length_hrs+1:24-SD_length_hrs,:,:),1,'omitnan'),3,'omitnan');
		Means.(Sexes{2}).Mut.SDFirst12hr.NREM = NormFactor*mean(mean(EEG_Normalized_Female_Mut_SD_NREM(1:12-SD_length_hrs,:,:),1,'omitnan'),3,'omitnan');
		Means.(Sexes{2}).Mut.SDLast12hr.NREM  = NormFactor*mean(mean(EEG_Normalized_Female_Mut_SD_NREM(12-SD_length_hrs+1:24-SD_length_hrs,:,:),1,'omitnan'),3,'omitnan');
		Means.(Sexes{2}).Mut.SDFirst12hr.REM  = NormFactor*mean(mean(EEG_Normalized_Female_Mut_SD_REM(1:12-SD_length_hrs,:,:),1,'omitnan'),3,'omitnan');
		Means.(Sexes{2}).Mut.SDLast12hr.REM   = NormFactor*mean(mean(EEG_Normalized_Female_Mut_SD_REM(12-SD_length_hrs+1:24-SD_length_hrs,:,:),1,'omitnan'),3,'omitnan');
	end 

end % end of if SeparateSpectralIntoLPDP


% -- End of Means --

% -- SEMS --
% Male WT          was 100*std(mean ....
SEMS.(Sexes{1}).WT.BL.Wake = NormFactor*std(mean(EEG_Normalized_Male_WT_BL_Wake,1,'omitnan'),0,3,'omitnan')./sqrt(size(EEG_Normalized_Male_WT_BL_Wake,3));
SEMS.(Sexes{1}).WT.BL.NREM = NormFactor*std(mean(EEG_Normalized_Male_WT_BL_NREM,1,'omitnan'),0,3,'omitnan')./sqrt(size(EEG_Normalized_Male_WT_BL_NREM,3));
SEMS.(Sexes{1}).WT.BL.REM  = NormFactor*std(mean(EEG_Normalized_Male_WT_BL_REM,1,'omitnan'),0,3,'omitnan')./sqrt(size(EEG_Normalized_Male_WT_BL_REM,3));
SEMS.(Sexes{1}).WT.SD.Wake = NormFactor*std(mean(EEG_Normalized_Male_WT_SD_Wake,1,'omitnan'),0,3,'omitnan')./sqrt(size(EEG_Normalized_Male_WT_SD_Wake,3));
SEMS.(Sexes{1}).WT.SD.NREM = NormFactor*std(mean(EEG_Normalized_Male_WT_SD_NREM,1,'omitnan'),0,3,'omitnan')./sqrt(size(EEG_Normalized_Male_WT_SD_NREM,3));
SEMS.(Sexes{1}).WT.SD.REM  = NormFactor*std(mean(EEG_Normalized_Male_WT_SD_REM,1,'omitnan'),0,3,'omitnan')./sqrt(size(EEG_Normalized_Male_WT_SD_REM,3));

% Male Mut
SEMS.(Sexes{1}).Mut.BL.Wake = NormFactor*std(mean(EEG_Normalized_Male_Mut_BL_Wake,1,'omitnan'),0,3,'omitnan')./sqrt(size(EEG_Normalized_Male_Mut_BL_Wake,3));
SEMS.(Sexes{1}).Mut.BL.NREM = NormFactor*std(mean(EEG_Normalized_Male_Mut_BL_NREM,1,'omitnan'),0,3,'omitnan')./sqrt(size(EEG_Normalized_Male_Mut_BL_NREM,3));
SEMS.(Sexes{1}).Mut.BL.REM  = NormFactor*std(mean(EEG_Normalized_Male_Mut_BL_REM,1,'omitnan'),0,3,'omitnan')./sqrt(size(EEG_Normalized_Male_Mut_BL_REM,3));
SEMS.(Sexes{1}).Mut.SD.Wake = NormFactor*std(mean(EEG_Normalized_Male_Mut_SD_Wake,1,'omitnan'),0,3,'omitnan')./sqrt(size(EEG_Normalized_Male_Mut_SD_Wake,3));
SEMS.(Sexes{1}).Mut.SD.NREM = NormFactor*std(mean(EEG_Normalized_Male_Mut_SD_NREM,1,'omitnan'),0,3,'omitnan')./sqrt(size(EEG_Normalized_Male_Mut_SD_NREM,3));
SEMS.(Sexes{1}).Mut.SD.REM  = NormFactor*std(mean(EEG_Normalized_Male_Mut_SD_REM,1,'omitnan'),0,3,'omitnan')./sqrt(size(EEG_Normalized_Male_Mut_SD_REM,3));

if FemaleDataPresent
	% Female WT
	SEMS.(Sexes{2}).WT.BL.Wake = NormFactor*std(mean(EEG_Normalized_Female_WT_BL_Wake,1,'omitnan'),0,3,'omitnan')./sqrt(size(EEG_Normalized_Female_WT_BL_Wake,3));
	SEMS.(Sexes{2}).WT.BL.NREM = NormFactor*std(mean(EEG_Normalized_Female_WT_BL_NREM,1,'omitnan'),0,3,'omitnan')./sqrt(size(EEG_Normalized_Female_WT_BL_NREM,3));
	SEMS.(Sexes{2}).WT.BL.REM  = NormFactor*std(mean(EEG_Normalized_Female_WT_BL_REM,1,'omitnan'),0,3,'omitnan')./sqrt(size(EEG_Normalized_Female_WT_BL_REM,3));
	SEMS.(Sexes{2}).WT.SD.Wake = NormFactor*std(mean(EEG_Normalized_Female_WT_SD_Wake,1,'omitnan'),0,3,'omitnan')./sqrt(size(EEG_Normalized_Female_WT_SD_Wake,3));
	SEMS.(Sexes{2}).WT.SD.NREM = NormFactor*std(mean(EEG_Normalized_Female_WT_SD_NREM,1,'omitnan'),0,3,'omitnan')./sqrt(size(EEG_Normalized_Female_WT_SD_NREM,3));
	SEMS.(Sexes{2}).WT.SD.REM  = NormFactor*std(mean(EEG_Normalized_Female_WT_SD_REM,1,'omitnan'),0,3,'omitnan')./sqrt(size(EEG_Normalized_Female_WT_SD_REM,3));

	% Female Mut
	SEMS.(Sexes{2}).Mut.BL.Wake = NormFactor*std(mean(EEG_Normalized_Female_Mut_BL_Wake,1,'omitnan'),0,3,'omitnan')./sqrt(size(EEG_Normalized_Female_Mut_BL_Wake,3));
	SEMS.(Sexes{2}).Mut.BL.NREM = NormFactor*std(mean(EEG_Normalized_Female_Mut_BL_NREM,1,'omitnan'),0,3,'omitnan')./sqrt(size(EEG_Normalized_Female_Mut_BL_NREM,3));
	SEMS.(Sexes{2}).Mut.BL.REM  = NormFactor*std(mean(EEG_Normalized_Female_Mut_BL_REM,1,'omitnan'),0,3,'omitnan')./sqrt(size(EEG_Normalized_Female_Mut_BL_REM,3));
	SEMS.(Sexes{2}).Mut.SD.Wake = NormFactor*std(mean(EEG_Normalized_Female_Mut_SD_Wake,1,'omitnan'),0,3,'omitnan')./sqrt(size(EEG_Normalized_Female_Mut_SD_Wake,3));
	SEMS.(Sexes{2}).Mut.SD.NREM = NormFactor*std(mean(EEG_Normalized_Female_Mut_SD_NREM,1,'omitnan'),0,3,'omitnan')./sqrt(size(EEG_Normalized_Female_Mut_SD_NREM,3));
	SEMS.(Sexes{2}).Mut.SD.REM  = NormFactor*std(mean(EEG_Normalized_Female_Mut_SD_REM,1,'omitnan'),0,3,'omitnan')./sqrt(size(EEG_Normalized_Female_Mut_SD_REM,3));
end 

if SeparateSpectralIntoLPDP
	% Male WT          was 100*std(mean ....
	SEMS.(Sexes{1}).WT.BLFirst12hr.Wake = NormFactor*std(mean(EEG_Normalized_Male_WT_BL_Wake(1:12,:,:),1,'omitnan'),0,3,'omitnan')./sqrt(size(EEG_Normalized_Male_WT_BL_Wake,3));
	SEMS.(Sexes{1}).WT.BLLast12hr.Wake  = NormFactor*std(mean(EEG_Normalized_Male_WT_BL_Wake(13:24,:,:),1,'omitnan'),0,3,'omitnan')./sqrt(size(EEG_Normalized_Male_WT_BL_Wake,3));
	SEMS.(Sexes{1}).WT.BLFirst12hr.NREM = NormFactor*std(mean(EEG_Normalized_Male_WT_BL_NREM(1:12,:,:),1,'omitnan'),0,3,'omitnan')./sqrt(size(EEG_Normalized_Male_WT_BL_NREM,3));
	SEMS.(Sexes{1}).WT.BLLast12hr.NREM  = NormFactor*std(mean(EEG_Normalized_Male_WT_BL_NREM(13:24,:,:),1,'omitnan'),0,3,'omitnan')./sqrt(size(EEG_Normalized_Male_WT_BL_NREM,3));
	SEMS.(Sexes{1}).WT.BLFirst12hr.REM  = NormFactor*std(mean(EEG_Normalized_Male_WT_BL_REM(1:12,:,:),1,'omitnan'),0,3,'omitnan')./sqrt(size(EEG_Normalized_Male_WT_BL_REM,3));
	SEMS.(Sexes{1}).WT.BLLast12hr.REM   = NormFactor*std(mean(EEG_Normalized_Male_WT_BL_REM(13:24,:,:),1,'omitnan'),0,3,'omitnan')./sqrt(size(EEG_Normalized_Male_WT_BL_REM,3));
	SEMS.(Sexes{1}).WT.SDFirst12hr.Wake = NormFactor*std(mean(EEG_Normalized_Male_WT_SD_Wake(1:12-SD_length_hrs,:,:),1,'omitnan'),0,3,'omitnan')./sqrt(size(EEG_Normalized_Male_WT_SD_Wake,3));
	SEMS.(Sexes{1}).WT.SDLast12hr.Wake  = NormFactor*std(mean(EEG_Normalized_Male_WT_SD_Wake(12-SD_length_hrs+1:24-SD_length_hrs,:,:),1,'omitnan'),0,3,'omitnan')./sqrt(size(EEG_Normalized_Male_WT_SD_Wake,3));
	SEMS.(Sexes{1}).WT.SDFirst12hr.NREM = NormFactor*std(mean(EEG_Normalized_Male_WT_SD_NREM(1:12-SD_length_hrs,:,:),1,'omitnan'),0,3,'omitnan')./sqrt(size(EEG_Normalized_Male_WT_SD_NREM,3));
	SEMS.(Sexes{1}).WT.SDLast12hr.NREM  = NormFactor*std(mean(EEG_Normalized_Male_WT_SD_NREM(12-SD_length_hrs+1:24-SD_length_hrs,:,:),1,'omitnan'),0,3,'omitnan')./sqrt(size(EEG_Normalized_Male_WT_SD_NREM,3));
	SEMS.(Sexes{1}).WT.SDFirst12hr.REM  = NormFactor*std(mean(EEG_Normalized_Male_WT_SD_REM(1:12-SD_length_hrs,:,:),1,'omitnan'),0,3,'omitnan')./sqrt(size(EEG_Normalized_Male_WT_SD_REM,3));
	SEMS.(Sexes{1}).WT.SDLast12hr.REM   = NormFactor*std(mean(EEG_Normalized_Male_WT_SD_REM(12-SD_length_hrs+1:24-SD_length_hrs,:,:),1,'omitnan'),0,3,'omitnan')./sqrt(size(EEG_Normalized_Male_WT_SD_REM,3));

	% Male Mut
	SEMS.(Sexes{1}).Mut.BLFirst12hr.Wake = NormFactor*std(mean(EEG_Normalized_Male_Mut_BL_Wake(1:12,:,:),1,'omitnan'),0,3,'omitnan')./sqrt(size(EEG_Normalized_Male_Mut_BL_Wake,3));
	SEMS.(Sexes{1}).Mut.BLLast12hr.Wake  = NormFactor*std(mean(EEG_Normalized_Male_Mut_BL_Wake(13:24,:,:),1,'omitnan'),0,3,'omitnan')./sqrt(size(EEG_Normalized_Male_Mut_BL_Wake,3));
	SEMS.(Sexes{1}).Mut.BLFirst12hr.NREM = NormFactor*std(mean(EEG_Normalized_Male_Mut_BL_NREM(1:12,:,:),1,'omitnan'),0,3,'omitnan')./sqrt(size(EEG_Normalized_Male_Mut_BL_NREM,3));
	SEMS.(Sexes{1}).Mut.BLLast12hr.NREM  = NormFactor*std(mean(EEG_Normalized_Male_Mut_BL_NREM(13:24,:,:),1,'omitnan'),0,3,'omitnan')./sqrt(size(EEG_Normalized_Male_Mut_BL_NREM,3));
	SEMS.(Sexes{1}).Mut.BLFirst12hr.REM  = NormFactor*std(mean(EEG_Normalized_Male_Mut_BL_REM(1:12,:,:),1,'omitnan'),0,3,'omitnan')./sqrt(size(EEG_Normalized_Male_Mut_BL_REM,3));
	SEMS.(Sexes{1}).Mut.BLLast12hr.REM   = NormFactor*std(mean(EEG_Normalized_Male_Mut_BL_REM(13:24,:,:),1,'omitnan'),0,3,'omitnan')./sqrt(size(EEG_Normalized_Male_Mut_BL_REM,3));
	SEMS.(Sexes{1}).Mut.SDFirst12hr.Wake = NormFactor*std(mean(EEG_Normalized_Male_Mut_SD_Wake(1:12-SD_length_hrs,:,:),1,'omitnan'),0,3,'omitnan')./sqrt(size(EEG_Normalized_Male_Mut_SD_Wake,3));
	SEMS.(Sexes{1}).Mut.SDLast12hr.Wake  = NormFactor*std(mean(EEG_Normalized_Male_Mut_SD_Wake(12-SD_length_hrs+1:24-SD_length_hrs,:,:),1,'omitnan'),0,3,'omitnan')./sqrt(size(EEG_Normalized_Male_Mut_SD_Wake,3));
	SEMS.(Sexes{1}).Mut.SDFirst12hr.NREM = NormFactor*std(mean(EEG_Normalized_Male_Mut_SD_NREM(1:12-SD_length_hrs,:,:),1,'omitnan'),0,3,'omitnan')./sqrt(size(EEG_Normalized_Male_Mut_SD_NREM,3));
	SEMS.(Sexes{1}).Mut.SDLast12hr.NREM  = NormFactor*std(mean(EEG_Normalized_Male_Mut_SD_NREM(12-SD_length_hrs+1:24-SD_length_hrs,:,:),1,'omitnan'),0,3,'omitnan')./sqrt(size(EEG_Normalized_Male_Mut_SD_NREM,3));
	SEMS.(Sexes{1}).Mut.SDFirst12hr.REM  = NormFactor*std(mean(EEG_Normalized_Male_Mut_SD_REM(1:12-SD_length_hrs,:,:),1,'omitnan'),0,3,'omitnan')./sqrt(size(EEG_Normalized_Male_Mut_SD_REM,3));
	SEMS.(Sexes{1}).Mut.SDLast12hr.REM   = NormFactor*std(mean(EEG_Normalized_Male_Mut_SD_REM(12-SD_length_hrs+1:24-SD_length_hrs,:,:),1,'omitnan'),0,3,'omitnan')./sqrt(size(EEG_Normalized_Male_Mut_SD_REM,3));

	if FemaleDataPresent
		% Female WT
		SEMS.(Sexes{2}).WT.BLFirst12hr.Wake = NormFactor*std(mean(EEG_Normalized_Female_WT_BL_Wake(1:12,:,:),1,'omitnan'),0,3,'omitnan')./sqrt(size(EEG_Normalized_Female_WT_BL_Wake,3));
		SEMS.(Sexes{2}).WT.BLLast12hr.Wake  = NormFactor*std(mean(EEG_Normalized_Female_WT_BL_Wake(13:24,:,:),1,'omitnan'),0,3,'omitnan')./sqrt(size(EEG_Normalized_Female_WT_BL_Wake,3));
		SEMS.(Sexes{2}).WT.BLFirst12hr.NREM = NormFactor*std(mean(EEG_Normalized_Female_WT_BL_NREM(1:12,:,:),1,'omitnan'),0,3,'omitnan')./sqrt(size(EEG_Normalized_Female_WT_BL_NREM,3));
		SEMS.(Sexes{2}).WT.BLLast12hr.NREM  = NormFactor*std(mean(EEG_Normalized_Female_WT_BL_NREM(13:24,:,:),1,'omitnan'),0,3,'omitnan')./sqrt(size(EEG_Normalized_Female_WT_BL_NREM,3));
		SEMS.(Sexes{2}).WT.BLFirst12hr.REM  = NormFactor*std(mean(EEG_Normalized_Female_WT_BL_REM(1:12,:,:),1,'omitnan'),0,3,'omitnan')./sqrt(size(EEG_Normalized_Female_WT_BL_REM,3));
		SEMS.(Sexes{2}).WT.BLLast12hr.REM   = NormFactor*std(mean(EEG_Normalized_Female_WT_BL_REM(13:24,:,:),1,'omitnan'),0,3,'omitnan')./sqrt(size(EEG_Normalized_Female_WT_BL_REM,3));
		SEMS.(Sexes{2}).WT.SDFirst12hr.Wake = NormFactor*std(mean(EEG_Normalized_Female_WT_SD_Wake(1:12-SD_length_hrs,:,:),1,'omitnan'),0,3,'omitnan')./sqrt(size(EEG_Normalized_Female_WT_SD_Wake,3));
		SEMS.(Sexes{2}).WT.SDLast12hr.Wake  = NormFactor*std(mean(EEG_Normalized_Female_WT_SD_Wake(12-SD_length_hrs+1:24-SD_length_hrs,:,:),1,'omitnan'),0,3,'omitnan')./sqrt(size(EEG_Normalized_Female_WT_SD_Wake,3));
		SEMS.(Sexes{2}).WT.SDFirst12hr.NREM = NormFactor*std(mean(EEG_Normalized_Female_WT_SD_NREM(1:12-SD_length_hrs,:,:),1,'omitnan'),0,3,'omitnan')./sqrt(size(EEG_Normalized_Female_WT_SD_NREM,3));
		SEMS.(Sexes{2}).WT.SDLast12hr.NREM  = NormFactor*std(mean(EEG_Normalized_Female_WT_SD_NREM(12-SD_length_hrs+1:24-SD_length_hrs,:,:),1,'omitnan'),0,3,'omitnan')./sqrt(size(EEG_Normalized_Female_WT_SD_NREM,3));
		SEMS.(Sexes{2}).WT.SDFirst12hr.REM  = NormFactor*std(mean(EEG_Normalized_Female_WT_SD_REM(1:12-SD_length_hrs,:,:),1,'omitnan'),0,3,'omitnan')./sqrt(size(EEG_Normalized_Female_WT_SD_REM,3));
		SEMS.(Sexes{2}).WT.SDLast12hr.REM   = NormFactor*std(mean(EEG_Normalized_Female_WT_SD_REM(12-SD_length_hrs+1:24-SD_length_hrs,:,:),1,'omitnan'),0,3,'omitnan')./sqrt(size(EEG_Normalized_Female_WT_SD_REM,3));

		% Female Mut
		SEMS.(Sexes{2}).Mut.BLFirst12hr.Wake = NormFactor*std(mean(EEG_Normalized_Female_Mut_BL_Wake(1:12,:,:),1,'omitnan'),0,3,'omitnan')./sqrt(size(EEG_Normalized_Female_Mut_BL_Wake,3));
		SEMS.(Sexes{2}).Mut.BLLast12hr.Wake  = NormFactor*std(mean(EEG_Normalized_Female_Mut_BL_Wake(13:24,:,:),1,'omitnan'),0,3,'omitnan')./sqrt(size(EEG_Normalized_Female_Mut_BL_Wake,3));
		SEMS.(Sexes{2}).Mut.BLFirst12hr.NREM = NormFactor*std(mean(EEG_Normalized_Female_Mut_BL_NREM(1:12,:,:),1,'omitnan'),0,3,'omitnan')./sqrt(size(EEG_Normalized_Female_Mut_BL_NREM,3));
		SEMS.(Sexes{2}).Mut.BLLast12hr.NREM  = NormFactor*std(mean(EEG_Normalized_Female_Mut_BL_NREM(13:24,:,:),1,'omitnan'),0,3,'omitnan')./sqrt(size(EEG_Normalized_Female_Mut_BL_NREM,3));
		SEMS.(Sexes{2}).Mut.BLFirst12hr.REM  = NormFactor*std(mean(EEG_Normalized_Female_Mut_BL_REM(1:12,:,:),1,'omitnan'),0,3,'omitnan')./sqrt(size(EEG_Normalized_Female_Mut_BL_REM,3));
		SEMS.(Sexes{2}).Mut.BLLast12hr.REM   = NormFactor*std(mean(EEG_Normalized_Female_Mut_BL_REM(13:24,:,:),1,'omitnan'),0,3,'omitnan')./sqrt(size(EEG_Normalized_Female_Mut_BL_REM,3));
		SEMS.(Sexes{2}).Mut.SDFirst12hr.Wake = NormFactor*std(mean(EEG_Normalized_Female_Mut_SD_Wake(1:12-SD_length_hrs,:,:),1,'omitnan'),0,3,'omitnan')./sqrt(size(EEG_Normalized_Female_Mut_SD_Wake,3));
		SEMS.(Sexes{2}).Mut.SDLast12hr.Wake  = NormFactor*std(mean(EEG_Normalized_Female_Mut_SD_Wake(12-SD_length_hrs+1:24-SD_length_hrs,:,:),1,'omitnan'),0,3,'omitnan')./sqrt(size(EEG_Normalized_Female_Mut_SD_Wake,3));
		SEMS.(Sexes{2}).Mut.SDFirst12hr.NREM = NormFactor*std(mean(EEG_Normalized_Female_Mut_SD_NREM(1:12-SD_length_hrs,:,:),1,'omitnan'),0,3,'omitnan')./sqrt(size(EEG_Normalized_Female_Mut_SD_NREM,3));
		SEMS.(Sexes{2}).Mut.SDLast12hr.NREM  = NormFactor*std(mean(EEG_Normalized_Female_Mut_SD_NREM(12-SD_length_hrs+1:24-SD_length_hrs,:,:),1,'omitnan'),0,3,'omitnan')./sqrt(size(EEG_Normalized_Female_Mut_SD_NREM,3));
		SEMS.(Sexes{2}).Mut.SDFirst12hr.REM  = NormFactor*std(mean(EEG_Normalized_Female_Mut_SD_REM(1:12-SD_length_hrs,:,:),1,'omitnan'),0,3,'omitnan')./sqrt(size(EEG_Normalized_Female_Mut_SD_REM,3));
		SEMS.(Sexes{2}).Mut.SDLast12hr.REM   = NormFactor*std(mean(EEG_Normalized_Female_Mut_SD_REM(12-SD_length_hrs+1:24-SD_length_hrs,:,:),1,'omitnan'),0,3,'omitnan')./sqrt(size(EEG_Normalized_Female_Mut_SD_REM,3));
	end 
end 	% end of if SeparateSpectralIntoLPDP

%ts_W_Mut = tinv([0.025  0.975],size(WakeRel1hMut,3)-1); % 95% confidence interval


% Make the individual spectral power vs frequency figures
% BL
male_w_bl_fig_handle  = plot_spectral_power_vs_freq(Means_struct=Means,SEMS_struct=SEMS,sex=Sexes{1},BLorSD='BL',state='Wake',LegendLabels=LegendLabels,EEG_bin_edges=EEG_bin_edges,EEGLowerLimitHz=EEGLowerLimit_Hz,Normalization=Normalization);
male_nr_bl_fig_handle = plot_spectral_power_vs_freq(Means_struct=Means,SEMS_struct=SEMS,sex=Sexes{1},BLorSD='BL',state='NREM',LegendLabels=LegendLabels,EEG_bin_edges=EEG_bin_edges,EEGLowerLimitHz=EEGLowerLimit_Hz,Normalization=Normalization);
male_r_bl_fig_handle  = plot_spectral_power_vs_freq(Means_struct=Means,SEMS_struct=SEMS,sex=Sexes{1},BLorSD='BL',state='REM',LegendLabels=LegendLabels,EEG_bin_edges=EEG_bin_edges,EEGLowerLimitHz=EEGLowerLimit_Hz,Normalization=Normalization);

if FemaleDataPresent
	female_w_bl_fig_handle  = plot_spectral_power_vs_freq(Means_struct=Means,SEMS_struct=SEMS,sex=Sexes{2},BLorSD='BL',state='Wake',LegendLabels=LegendLabels,EEG_bin_edges=EEG_bin_edges,EEGLowerLimitHz=EEGLowerLimit_Hz,Normalization=Normalization);
	female_nr_bl_fig_handle = plot_spectral_power_vs_freq(Means_struct=Means,SEMS_struct=SEMS,sex=Sexes{2},BLorSD='BL',state='NREM',LegendLabels=LegendLabels,EEG_bin_edges=EEG_bin_edges,EEGLowerLimitHz=EEGLowerLimit_Hz,Normalization=Normalization);
	female_r_bl_fig_handle  = plot_spectral_power_vs_freq(Means_struct=Means,SEMS_struct=SEMS,sex=Sexes{2},BLorSD='BL',state='REM',LegendLabels=LegendLabels,EEG_bin_edges=EEG_bin_edges,EEGLowerLimitHz=EEGLowerLimit_Hz,Normalization=Normalization);
end 

% SD
male_w_sd_fig_handle  = plot_spectral_power_vs_freq(Means_struct=Means,SEMS_struct=SEMS,sex=Sexes{1},BLorSD='SD',state='Wake',LegendLabels=LegendLabels,EEG_bin_edges=EEG_bin_edges,EEGLowerLimitHz=EEGLowerLimit_Hz,Normalization=Normalization);
male_nr_sd_fig_handle = plot_spectral_power_vs_freq(Means_struct=Means,SEMS_struct=SEMS,sex=Sexes{1},BLorSD='SD',state='NREM',LegendLabels=LegendLabels,EEG_bin_edges=EEG_bin_edges,EEGLowerLimitHz=EEGLowerLimit_Hz,Normalization=Normalization);
male_r_sd_fig_handle  = plot_spectral_power_vs_freq(Means_struct=Means,SEMS_struct=SEMS,sex=Sexes{1},BLorSD='SD',state='REM',LegendLabels=LegendLabels,EEG_bin_edges=EEG_bin_edges,EEGLowerLimitHz=EEGLowerLimit_Hz,Normalization=Normalization);

if FemaleDataPresent
	female_w_sd_fig_handle  = plot_spectral_power_vs_freq(Means_struct=Means,SEMS_struct=SEMS,sex=Sexes{2},BLorSD='SD',state='Wake',LegendLabels=LegendLabels,EEG_bin_edges=EEG_bin_edges,EEGLowerLimitHz=EEGLowerLimit_Hz,Normalization=Normalization);
	female_nr_sd_fig_handle = plot_spectral_power_vs_freq(Means_struct=Means,SEMS_struct=SEMS,sex=Sexes{2},BLorSD='SD',state='NREM',LegendLabels=LegendLabels,EEG_bin_edges=EEG_bin_edges,EEGLowerLimitHz=EEGLowerLimit_Hz,Normalization=Normalization);
	female_r_sd_fig_handle  = plot_spectral_power_vs_freq(Means_struct=Means,SEMS_struct=SEMS,sex=Sexes{2},BLorSD='SD',state='REM',LegendLabels=LegendLabels,EEG_bin_edges=EEG_bin_edges,EEGLowerLimitHz=EEGLowerLimit_Hz,Normalization=Normalization);
end 

% Combine panels into one figure each for BL and sleep dep?  
combine_spectral_power_by_state_into_one_big_figure(FemaleDataPresent,Normalization,{'Spectral Baseline Figure','Spectral SD Figure'})


% ---- If separating by LP and DP, make 4 figures: BL LP, BL DP, SD LP, SD DP ---- 
if SeparateSpectralIntoLPDP
	% -- BL LP --
	male_w_bllp_fig_handle  = plot_spectral_power_vs_freq(Means_struct=Means,SEMS_struct=SEMS,sex=Sexes{1},BLorSD='BLFirst12hr',state='Wake',LegendLabels=LegendLabels,EEG_bin_edges=EEG_bin_edges,EEGLowerLimitHz=EEGLowerLimit_Hz,Normalization=Normalization);
	male_nr_bllp_fig_handle = plot_spectral_power_vs_freq(Means_struct=Means,SEMS_struct=SEMS,sex=Sexes{1},BLorSD='BLFirst12hr',state='NREM',LegendLabels=LegendLabels,EEG_bin_edges=EEG_bin_edges,EEGLowerLimitHz=EEGLowerLimit_Hz,Normalization=Normalization);
	male_r_bllp_fig_handle  = plot_spectral_power_vs_freq(Means_struct=Means,SEMS_struct=SEMS,sex=Sexes{1},BLorSD='BLFirst12hr',state='REM',LegendLabels=LegendLabels,EEG_bin_edges=EEG_bin_edges,EEGLowerLimitHz=EEGLowerLimit_Hz,Normalization=Normalization);

	if FemaleDataPresent
		female_w_bllp_fig_handle  = plot_spectral_power_vs_freq(Means_struct=Means,SEMS_struct=SEMS,sex=Sexes{2},BLorSD='BLFirst12hr',state='Wake',LegendLabels=LegendLabels,EEG_bin_edges=EEG_bin_edges,EEGLowerLimitHz=EEGLowerLimit_Hz,Normalization=Normalization);
		female_nr_bllp_fig_handle = plot_spectral_power_vs_freq(Means_struct=Means,SEMS_struct=SEMS,sex=Sexes{2},BLorSD='BLFirst12hr',state='NREM',LegendLabels=LegendLabels,EEG_bin_edges=EEG_bin_edges,EEGLowerLimitHz=EEGLowerLimit_Hz,Normalization=Normalization);
		female_r_bllp_fig_handle  = plot_spectral_power_vs_freq(Means_struct=Means,SEMS_struct=SEMS,sex=Sexes{2},BLorSD='BLFirst12hr',state='REM',LegendLabels=LegendLabels,EEG_bin_edges=EEG_bin_edges,EEGLowerLimitHz=EEGLowerLimit_Hz,Normalization=Normalization);
	end 

	% -- SD LP --
	male_w_sdlp_fig_handle  = plot_spectral_power_vs_freq(Means_struct=Means,SEMS_struct=SEMS,sex=Sexes{1},BLorSD='SDFirst12hr',state='Wake',LegendLabels=LegendLabels,EEG_bin_edges=EEG_bin_edges,EEGLowerLimitHz=EEGLowerLimit_Hz,Normalization=Normalization);
	male_nr_sdlp_fig_handle = plot_spectral_power_vs_freq(Means_struct=Means,SEMS_struct=SEMS,sex=Sexes{1},BLorSD='SDFirst12hr',state='NREM',LegendLabels=LegendLabels,EEG_bin_edges=EEG_bin_edges,EEGLowerLimitHz=EEGLowerLimit_Hz,Normalization=Normalization);
	male_r_sdlp_fig_handle  = plot_spectral_power_vs_freq(Means_struct=Means,SEMS_struct=SEMS,sex=Sexes{1},BLorSD='SDFirst12hr',state='REM',LegendLabels=LegendLabels,EEG_bin_edges=EEG_bin_edges,EEGLowerLimitHz=EEGLowerLimit_Hz,Normalization=Normalization);

	if FemaleDataPresent
		female_w_sdlp_fig_handle  = plot_spectral_power_vs_freq(Means_struct=Means,SEMS_struct=SEMS,sex=Sexes{2},BLorSD='SDFirst12hr',state='Wake',LegendLabels=LegendLabels,EEG_bin_edges=EEG_bin_edges,EEGLowerLimitHz=EEGLowerLimit_Hz,Normalization=Normalization);
		female_nr_sdlp_fig_handle = plot_spectral_power_vs_freq(Means_struct=Means,SEMS_struct=SEMS,sex=Sexes{2},BLorSD='SDFirst12hr',state='NREM',LegendLabels=LegendLabels,EEG_bin_edges=EEG_bin_edges,EEGLowerLimitHz=EEGLowerLimit_Hz,Normalization=Normalization);
		female_r_sdlp_fig_handle  = plot_spectral_power_vs_freq(Means_struct=Means,SEMS_struct=SEMS,sex=Sexes{2},BLorSD='SDFirst12hr',state='REM',LegendLabels=LegendLabels,EEG_bin_edges=EEG_bin_edges,EEGLowerLimitHz=EEGLowerLimit_Hz,Normalization=Normalization);
	end 

	combine_spectral_power_by_state_into_one_big_figure(FemaleDataPresent,Normalization,{'Spectral Baseline LP Figure','Spectral SD LP Figure'})

	% -- BL DP --
	male_w_bldp_fig_handle  = plot_spectral_power_vs_freq(Means_struct=Means,SEMS_struct=SEMS,sex=Sexes{1},BLorSD='BLLast12hr',state='Wake',LegendLabels=LegendLabels,EEG_bin_edges=EEG_bin_edges,EEGLowerLimitHz=EEGLowerLimit_Hz,Normalization=Normalization);
	male_nr_bldp_fig_handle = plot_spectral_power_vs_freq(Means_struct=Means,SEMS_struct=SEMS,sex=Sexes{1},BLorSD='BLLast12hr',state='NREM',LegendLabels=LegendLabels,EEG_bin_edges=EEG_bin_edges,EEGLowerLimitHz=EEGLowerLimit_Hz,Normalization=Normalization);
	male_r_bldp_fig_handle  = plot_spectral_power_vs_freq(Means_struct=Means,SEMS_struct=SEMS,sex=Sexes{1},BLorSD='BLLast12hr',state='REM',LegendLabels=LegendLabels,EEG_bin_edges=EEG_bin_edges,EEGLowerLimitHz=EEGLowerLimit_Hz,Normalization=Normalization);

	if FemaleDataPresent
		female_w_bldp_fig_handle  = plot_spectral_power_vs_freq(Means_struct=Means,SEMS_struct=SEMS,sex=Sexes{2},BLorSD='BLLast12hr',state='Wake',LegendLabels=LegendLabels,EEG_bin_edges=EEG_bin_edges,EEGLowerLimitHz=EEGLowerLimit_Hz,Normalization=Normalization);
		female_nr_bldp_fig_handle = plot_spectral_power_vs_freq(Means_struct=Means,SEMS_struct=SEMS,sex=Sexes{2},BLorSD='BLLast12hr',state='NREM',LegendLabels=LegendLabels,EEG_bin_edges=EEG_bin_edges,EEGLowerLimitHz=EEGLowerLimit_Hz,Normalization=Normalization);
		female_r_bldp_fig_handle  = plot_spectral_power_vs_freq(Means_struct=Means,SEMS_struct=SEMS,sex=Sexes{2},BLorSD='BLLast12hr',state='REM',LegendLabels=LegendLabels,EEG_bin_edges=EEG_bin_edges,EEGLowerLimitHz=EEGLowerLimit_Hz,Normalization=Normalization);
	end 

	% -- SD DP --
	male_w_sddp_fig_handle  = plot_spectral_power_vs_freq(Means_struct=Means,SEMS_struct=SEMS,sex=Sexes{1},BLorSD='SDLast12hr',state='Wake',LegendLabels=LegendLabels,EEG_bin_edges=EEG_bin_edges,EEGLowerLimitHz=EEGLowerLimit_Hz,Normalization=Normalization);
	male_nr_sddp_fig_handle = plot_spectral_power_vs_freq(Means_struct=Means,SEMS_struct=SEMS,sex=Sexes{1},BLorSD='SDLast12hr',state='NREM',LegendLabels=LegendLabels,EEG_bin_edges=EEG_bin_edges,EEGLowerLimitHz=EEGLowerLimit_Hz,Normalization=Normalization);
	male_r_sddp_fig_handle  = plot_spectral_power_vs_freq(Means_struct=Means,SEMS_struct=SEMS,sex=Sexes{1},BLorSD='SDLast12hr',state='REM',LegendLabels=LegendLabels,EEG_bin_edges=EEG_bin_edges,EEGLowerLimitHz=EEGLowerLimit_Hz,Normalization=Normalization);

	if FemaleDataPresent
		female_w_sddp_fig_handle  = plot_spectral_power_vs_freq(Means_struct=Means,SEMS_struct=SEMS,sex=Sexes{2},BLorSD='SDLast12hr',state='Wake',LegendLabels=LegendLabels,EEG_bin_edges=EEG_bin_edges,EEGLowerLimitHz=EEGLowerLimit_Hz,Normalization=Normalization);
		female_nr_sddp_fig_handle = plot_spectral_power_vs_freq(Means_struct=Means,SEMS_struct=SEMS,sex=Sexes{2},BLorSD='SDLast12hr',state='NREM',LegendLabels=LegendLabels,EEG_bin_edges=EEG_bin_edges,EEGLowerLimitHz=EEGLowerLimit_Hz,Normalization=Normalization);
		female_r_sddp_fig_handle  = plot_spectral_power_vs_freq(Means_struct=Means,SEMS_struct=SEMS,sex=Sexes{2},BLorSD='SDLast12hr',state='REM',LegendLabels=LegendLabels,EEG_bin_edges=EEG_bin_edges,EEGLowerLimitHz=EEGLowerLimit_Hz,Normalization=Normalization);
	end 

	combine_spectral_power_by_state_into_one_big_figure(FemaleDataPresent,Normalization,{'Spectral Baseline DP Figure','Spectral SD DP Figure'})





end % end of if SeparateSpectralIntoLPDP


% Now do the stats to compare the spectral curves? 
pvals = [];