%% load data
S12 = readtable('S12_exponent_s_gridfit_bybest10000_aa_avg_with_dist.csv');
S12.fitness = S12.s ;
T = S12(~S12.stop,:);
clear 'S12' ;
%% convert aa_seq_short to numeric vector
%%
aa_single_letter_codes = 'ACDEFGHIKLMNPQRSTVWY' ;
n_positions = length(T.aa_seq_short{1}) ; 
n_possible_aa = length(aa_single_letter_codes) ; 

n_variants = height(T) ;
varmat = zeros(n_positions , n_possible_aa);

real_x = zeros( n_variants , 1);

variants_fitness = T.fitness ; 
variants_seqs    = cell( n_variants , 1);

for VarI = 1:n_variants
    variants_seqs{VarI} = cell2mat(arrayfun( @(I)T.aa_seq_short{VarI}(I) == aa_single_letter_codes , 1:n_positions,'UniformOutput',false)');
end

%
Y = variants_fitness ; 
X = cell2mat(cellfun( @(X)X(:) , variants_seqs,'UniformOutput',false)')' ;

mdl = fitglm( X , Y);
fprintf('%0.0e\t%0.02f\n' , n_variants , mdl.Rsquared.Ordinary);
%% 
fitness_coeffs = reshape( mdl.Coefficients.Estimate(2:end) , [] , n_possible_aa)
fitness_pvals = reshape( mdl.Coefficients.pValue(2:end) , [] , n_possible_aa) 


%%
ltr = 'I';
ltr_pos = find(aa_single_letter_codes==ltr);
%figure; hold on;  grid on ;
%boxplot( T.fitness , cellfun( @(X)X(12)=='I' , T.aa_seq_short));
%set(gca,'xticklabel',{'I at 12','not'})
for I = 1:n_positions
idx = cellfun( @(X)X(I)==ltr , T.aa_seq_short);
%idx = cellfun( @(X)X(2)=='C' , T.aa_seq_short);
%idx = T.dist_Scer<3 ; 
figure; hold on;  grid on ;
%ecdf( T.fitness(idx)); ecdf( T.fitness(~idx));
boxplot(T.fitness,~idx);
set(gca,'xticklabel',{['I at ' num2str(I)] 'not'});
%legend({['I at ' num2str(I)] 'not'},'location','nw')
txt = sprintf( 'Coef=%0.02f p=%0.02f' , fitness_coeffs(I,ltr_pos) , fitness_pvals(I,ltr_pos)  );
title(txt)
end