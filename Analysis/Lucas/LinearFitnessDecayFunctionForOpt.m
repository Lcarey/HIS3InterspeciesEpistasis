function [ predicted_fitness ] = LinearFitnessDecayFunctionForOpt( init_m_b_ddGvect , aa_for_all_variants_vect  )
%% [ predicted_fitness ] = LinearFitnessDecayFunctionForOpt( init_a_b_ddGvect , aa_for_all_variants_vect  )  
%
% init_x0_k_L_ddGvect    : a vector of paramters for the fitness function (eg: sigmoid)
%             followed by delta_G_vect_init
%
% deltaG_matrix      : the MxN matrix of ddG values (fitness penalties) for each AA at each position
%                      N positions that vary, for M (20) possible amino acids at each pos. 0 for no fitness defect at a give position / AA
%
% aa_for_all_variants_vect   :  vectorized NxV matrix w/values (1-20) for each AA at each of the N positions
%
%
% LBC 2017

% global aa_for_all_variants_vect ; 
% global fitness_for_all_variants_vect ; 

m = init_m_b_ddGvect(1) ;
b  = init_m_b_ddGvect(2) ;
global n_possible_aas  ;


deltaG_matrix = reshape( init_m_b_ddGvect(4:end) , n_possible_aas , [] ) ;
deltaG_vect = deltaG_matrix(:);  % this is what we want to learn
n_positions = size(deltaG_matrix,2); 

aa_for_all_variants = reshape( aa_for_all_variants_vect , [] ,n_positions) ;
n_variants = size(aa_for_all_variants,1);

sum_delta_Gs  = arrayfun( @(I)sum( deltaG_vect( ((0:(n_positions-1)) .* n_possible_aas)+aa_for_all_variants(I,:))) , 1:n_variants) ;


predicted_fitness = b + (m .* sum_delta_Gs)  ; 

end