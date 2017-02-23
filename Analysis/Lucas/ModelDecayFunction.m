% k = 10 ; % steepness
% L = 1 ; % max Y
% x0 = 2 ; % midpoint of curve
% logistic = @(X) L / ( 1 + exp( -k*(X-x0) ) ) ; 
% % y = L / ( 1 + exp( -k*(x-x0) ) )

%%
fitness_function_parameters = [ 10 1 2] ;
n_positions = 10 ; 
n_possible_aas = 20 ; 

n_variants = 1e4 ;

% fitmat contains fitness defects for each AA at each pos
%  0 == no fitness defect
deltaG_matrix = randi(10,n_possible_aas,n_positions) ; % 15 positions for each of 20 AAs
deltaG_matrix = (deltaG_matrix./max(max(deltaG_matrix))).^10 ;
deltaG_vect = deltaG_matrix(:);
%%
xdata_for_all_variants = randi( 20 , n_variants , n_positions ) ;

sum_delta_Gs  = arrayfun( @(I)sum( deltaG_vect( ((0:(n_positions-1)) .* n_possible_aas)+xdata_for_all_variants(I,:))) , 1:n_variants) ;

k = fitness_function_parameters(1) ; % steepness
L = fitness_function_parameters(2) ; % max Y
x0 = fitness_function_parameters(3) ; % midpoint of curve
fitness_for_all_variants = L ./ ( 1 + exp( -k .* (sum_delta_Gs-x0) ) ) ;

fitness_function_parameters__npos__deltaG_vect = [fitness_function_parameters  , n_positions , deltaG_vect' ];

[ fitness_diff_real_pred  , predicted_fitness ] = LogisticFitnessDecayFunctionForOpt( ...
    fitness_function_parameters__npos__deltaG_vect , xdata_for_all_variants , fitness_for_all_variants )

lb = zeros( numel(fitness_function_parameters__npos__deltaG_vect) , 1);
lb(4) = n_positions ;
ub = ones( numel(fitness_function_parameters__npos__deltaG_vect) , 1);
ub(4) = n_positions ;
ub(1:3) = [100 10 10] ;

x = lsqcurvefit(LogisticFitnessDecayFunctionForOpt , fitness_function_parameters__npos__deltaG_vect , xdata_for_all_variants , fitness_for_all_variants , ...
    lb , ub )

%%
real_x = zeros( n_variants , 1);

amino_acids = 1:20 ;
fitness_penalties = zeros( numel(amino_acids) , 1);
fitness_penalties(1) = 1   ;

variants_fitness = NaN( n_variants , 1);
variants_seqs    = cell( n_variants , 1);

for I = 1:n_variants
    variants_seqs{I} = varmat ;
    nmut = randi(n_positions) ; 
    idx = randsample( n_positions, nmut) ;
    variants_seqs{I}(idx , 1 ) = 1 ;
    for J = find(variants_seqs{I}(:,1)==0)'
        variants_seqs{I}(J,randi(n_possible_aa))=1;
    end
 %   if  sum(variants_seqs{I}(:,3))   ; % for testing 'epistasis'
        real_x(I) = sum( fitmat( logical(variants_seqs{I})) ) ;
  %  end
    variants_fitness(I) = 1-logistic( real_x(I) );
end

%
Y = variants_fitness ; 
X = cell2mat(cellfun( @(X)X(:) , variants_seqs,'UniformOutput',false)')' ;

figure ; 
plot( real_x , variants_fitness,'ok');
xlabel('real x')
ylabel('variants fitness')


% %%
% variants_seqs{ variants_fitness == min(variants_fitness) }
% fitmat
% variants_seqs{ find(variants_fitness == max(variants_fitness),1) }
% 

%% 
[strongest_fitness_variant_r,strongest_fitness_variant_c]=find((fitmat==max(fitmat(:)))) ;
idx = arrayfun( @(I) variants_seqs{I}(strongest_fitness_variant_r,strongest_fitness_variant_c)==1 , 1:n_variants);
figure; hold on ; grid on; 
boxplot( variants_fitness , ~idx);
set(gca,'xticklabel',{'has worst AA' 'no'});
ylabel('fitness')
%histogram( variants_fitness(~idx),'Normalization','Probability')
%histogram( variants_fitness(idx),'Normalization','Probability')

%%
%figure; 
%plot( fitmat(:) ,'ok');
toplot = mdl.Coefficients.Estimate(2:end) ; % (1) is Intercept
fitvals = round(fitmat(:)*10) ; 
pVals = mdl.Coefficients.pValue(2:end);
figure;  hold on; 
boxplot( toplot(pVals<10.05) , fitvals(pVals<10.05) ) ;
%boxplot( toplot ,  fitmat(:)  ) ;

%set(gca,'xticklabel',{'not important' , 'Important' })
ylabel('Model coefficient')

%% find best fit to sigmoid
rbest = 0 ; 
parbest = [0 0 0];
for k = 0:0.1:2
    for L = 0:0.1:2
        for x0 = 0:0.1:2
            ypred = 1 - ( L ./ ( 1 + exp( -k.*(real_x-x0) ) ) ) ; 
            r = corr(ypred,variants_fitness);
            if(r>rbest)
                parbest = [k L x0];
                rbest = r;
            end
        end
    end
end
parbest  , rbest


%% Cross-validation
Y = variants_fitness ; 
X = cell2mat(cellfun( @(X)X(:) , variants_seqs,'UniformOutput',false)')' ;

Indices = crossvalind('Kfold', n_variants, 10);
fit_coeffs_lin = NaN( n_positions , n_possible_aa , 10);
fit_coeffs_log = NaN( n_positions , n_possible_aa , 10);

for I = 1:10
    test_idx = Indices==I ; 
    mdl_logistic = fitglm( X(~test_idx,:) , Y(~test_idx) ,'Link','logit');
    mdl_nologit  = fitglm( X(~test_idx,:) , Y(~test_idx) );
    
    y_logistic_pred = predict( mdl_logistic , X(test_idx,:)); 
    y_nologit_pred  = predict( mdl_nologit , X(test_idx,:)); 
    
    rl = corr(y_logistic_pred,Y(test_idx));
    rnol = corr(y_nologit_pred,Y(test_idx) );
    
    fit_coeffs_lin( : ,  : , I) = reshape( mdl_nologit.Coefficients.Estimate(2:end) , [] , n_possible_aa) ;
    fit_coeffs_log( : ,  : , I) = reshape( mdl_logistic.Coefficients.Estimate(2:end) , [] , n_possible_aa) ; 
    fprintf('%0.02f\t%0.02f\n' , corr( fitmat(:) , reshape(fit_coeffs_log(:,:,I),[],1) ) ,  corr( fitmat(:) , reshape(fit_coeffs_lin(:,:,I),[],1) ) ); 

    fprintf('%0.0e\t%0.02f,%0.02f\t%0.02f,%0.02f\n' , n_variants , mdl_logistic.Rsquared.Ordinary , rl^2 , mdl_nologit.Rsquared.Ordinary , rnol^2);
end
fit_coeffs_lin = median(fit_coeffs_lin,3) ;
fit_coeffs_log = median(fit_coeffs_log,3) ;

corr( fitmat(:) , fit_coeffs_log(:))
corr( fitmat(:) , fit_coeffs_lin(:))