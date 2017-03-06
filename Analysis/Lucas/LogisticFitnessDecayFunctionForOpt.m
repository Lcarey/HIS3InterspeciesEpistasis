function [ fitness_diff_real_pred  , predicted_fitness ] = LogisticFitnessDecayFunctionForOpt( fitness_function_parameters__npos__deltaG_vect , xdata_for_all_variants , fitness_for_all_variants )
%% fitness_diff_real_pred = LogisticFitnessDecayFunctionForOpt( fitness_function_parameters__npos__deltaG_vect , xdata_for_all_variants , fitness_for_all_variants )
% 
% fitness_function_parameters    : a vector of paramters for the fitness function (eg: sigmoid)
%
% deltaG_matrix      : the MxN matrix of ddG values (fitness penalties) for each AA at each position
%                      N positions that vary, for M (20) possible amino acids at each pos. 0 for no fitness defect at a give position / AA
%
% xdata_for_all_variants   :  an NxV matrix w/values (1-20) for each AA at each of the N positions
%
% fitness_for_all_variants   : a V length vector with the fitness values [0 1] for each variant
%
%
% LBC 2017

%fitness_for_all_variants = NaN( 1e3 , 1);
%deltaG_matrix = randi(10,20,5) ; % 15 positions for each of 20 AAs
%deltaG_matrix = (deltaG_matrix./max(max(deltaG_matrix))).^10 ;
n_possible_aas = 20 ;
k  = fitness_function_parameters__npos__deltaG_vect(1) ;
L  = fitness_function_parameters__npos__deltaG_vect(2) ;
x0 = fitness_function_parameters__npos__deltaG_vect(3) ;
n_positions = fitness_function_parameters__npos__deltaG_vect(4) ;
deltaG_vect = reshape( fitness_function_parameters__npos__deltaG_vect(5:end) , [] , 1) ;
n_variants = numel(fitness_for_all_variants);

% % % % 
% % % % n_variants = numel(fitness_for_all_variants);
% % % % n_positions = size(deltaG_matrix,2);
% % % % n_possible_aas = size(deltaG_matrix,1);
% % % % 
% % % % deltaG_vect = deltaG_matrix(:);

%   xdata_for_all_variants = randi( 20 , n_variants , n_positions );


sum_delta_Gs  = arrayfun( @(I)sum( deltaG_vect( ((0:(n_positions-1)) .* n_possible_aas)+xdata_for_all_variants(I,:))) , 1:n_variants) ;

predicted_fitness = L ./ ( 1 + exp( -k .* (sum_delta_Gs-x0) ) ) ;

fitness_diff_real_pred = sum( abs(predicted_fitness - fitness_for_all_variants) );

end

% y = L / ( 1 + exp( -k*(x-x0) ) )