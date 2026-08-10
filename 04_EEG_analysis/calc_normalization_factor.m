function norm_by_animal = calc_normalization_factor(EEG_Matrix,dx)
%
% USAGE: norm_by_animal = calc_normalization_factor(EEG_Matrix,dx)
%
%
%  This function calculates the area under the spectral curve (power vs freq)
%  For each animal individually.  It uses the trapeziodal rule. 
%
% INPUTS: 	EEG_Matrix 	A 24xMxN matrix of spectral power where M is the number of frequency bins and N=# of animals
%
%			dx:  The spacing between frequencies in the EEG (0.25 Hz or 1 Hz for instance)

% Handle the case where no data are present.  
if isempty(EEG_Matrix)
	norm_by_animal = [];
	return
end


trap_rule = inline('dx*sum(v)-dx/2*v(1)-dx/2*v(end)','v','dx');


for i=1:size(EEG_Matrix,3)
	m = mean(EEG_Matrix(:,:,i),1,'omitnan');
	norm_by_animal(i) = trap_rule(m,dx);
end

norm_by_animal = norm_by_animal';
