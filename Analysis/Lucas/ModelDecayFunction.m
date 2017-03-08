% k = 10 ; % steepness
% L = 1 ; % max Y
% x0 = 2 ; % midpoint of curve
% logistic = @(X) L / ( 1 + exp( -k*(X-x0) ) ) ; 
% % y = L / ( 1 + exp( -k*(x-x0) ) )

%%   
n_positions = 10 ; 
n_possible_aas = 20 ; 

n_variants = 3e3 ;

% fitmat contains fitness defects for each AA at each pos
%  0 == no fitness defect
deltaG_matrix = random('uniform',0,1,n_possible_aas,n_positions) ; % fitness defect caused by each AA at each position
%deltaG_matrix = (deltaG_matrix./max(max(deltaG_matrix))).^10 ;LSQNONLIN

deltaG_vect = deltaG_matrix(:);  % this is what we want to learn

aa_for_all_variants = randi( 20 , n_variants , n_positions ) ;

sum_delta_Gs  = arrayfun( @(I)sum( deltaG_vect( ((0:(n_positions-1)) .* n_possible_aas)+aa_for_all_variants(I,:))) , 1:n_variants) ;

% generated a matrix w/the fitness cost of each variant (rows -> AAs) at each position (column)
%   the sum of the fitness costs at each position is the sum_delta_Gs
%   from this sum we can calculate the actual fitness using a sigmoid function

%
k=3 ;  % steepness
x0 = 2.5 ;  % midpoint
L = 1 ; % max Y 
logistic_function = @(x0,k,L,xdata) L ./ (1+exp( (-1.*k).*(xdata-x0) ) ) ; % x0 , k , L
logistic_function_onevar = @(x,xdata) x(3) ./ (1+exp( (-1.*x(2)).*(xdata-x(1)) ) ) ; % x0 , k , L

fitness_for_all_variants_logistic = 1-logistic_function(x0,k,L,sum_delta_Gs) ;
fitness_for_all_variants_linear = 1-( sum_delta_Gs./max(sum_delta_Gs)) ; 

fitness_for_all_variants_logistic_noise = random('normal',0,0.00000005,size(fitness_for_all_variants_logistic)) + fitness_for_all_variants_logistic ;
max_theo_r2 = ( corr( fitness_for_all_variants_logistic_noise(:) , fitness_for_all_variants_logistic(:)) )^2 

%fitness_for_all_variants_logistic = fitness_for_all_variants_logistic_noise  ; 
% figure; hold on; 
% plot( sum_delta_Gs , fitness_for_all_variants_logistic ,'ok' , 'DisplayName' ,'Logistic')
% plot( sum_delta_Gs , fitness_for_all_variants_linear ,'ob', 'DisplayName' ,'Linear')
% 
% xlabel('Sum delta Gs')
% ylabel('Fitness')
% 
% log_parm_pred = lsqcurvefit( logistic_function_onevar , [x0 k L] , sum_delta_Gs , fitness_for_all_variants_linear ,[ -100 -100 -100], [100 100 100]);
% fit_pred = logistic_function_onevar( log_parm_pred , sum_delta_Gs ); 
% r2 = rsquare( fitness_for_all_variants_linear , fit_pred)
% 

%% Cross-validation

fit_opts = optimset( 'Diagnostics' , 'off' , 'Display' , 'off', 'PlotFcn' , [] ,'UseParallel' , false ) ;

init_x0_k_L_ddGvect = [0 0 0 deltaG_vect'] ; 
lb = zeros( size(init_x0_k_L_ddGvect));
ub = ones( size(init_x0_k_L_ddGvect));
ub(1:3) = 10 ;

Indices = crossvalind('Kfold', n_variants, 10);
R = dataset();
R.fit_ddGvect = cell(10,1);
R.Train = cell(10,1);
R.Test = cell(10,1);
R.train_rmse = NaN(10,1);
R.train_r2 = NaN(10,1) ;
R.test_rmse = NaN(10,1);
R.test_r2 = NaN(10,1) ;
R.runtime = R.test_r2 ; 
R.x0 = R.test_r2 ; 
R.k = R.test_r2 ; 
R.L = R.test_r2 ; 
R.FitOutputS = cell(10,1);

for I = 1:10
    R.Train{I} = find( Indices ~= I );
    R.Test{I} = find( Indices == I );
end

for I = 1:10
    tic; 
    aa = aa_for_all_variants( R.Train{I} , :);
    fit = fitness_for_all_variants_logistic_noise(  R.Train{I} );
    [fit_x0_k_L_ddGvect , resnorm , residual , exitflag , output] = lsqcurvefit( @LogisticFitnessDecayFunctionForOpt , init_x0_k_L_ddGvect , aa(:)  , fit , lb , ub , fit_opts ) ;
    R.FitOutputS{I} = output ;
    R.runtime(I) = toc; 
    R.fit_ddGvect{I} = fit_x0_k_L_ddGvect(4:end) ;
    R.x0(I) = fit_x0_k_L_ddGvect(1);
    R.k(I) = fit_x0_k_L_ddGvect(2);
    R.L(I) = fit_x0_k_L_ddGvect(3);
    R.x0(I) = fit_x0_k_L_ddGvect(1);

    pred_fit_train = LogisticFitnessDecayFunctionForOpt( fit_x0_k_L_ddGvect , aa(:)  ) ; 
    [ R.train_r2(I) , R.train_rmse(I)] = rsquare( fit , pred_fit_train ) ; 

    aa_test = aa_for_all_variants( R.Test{I} , :);
    fit_test = fitness_for_all_variants_logistic_noise(  R.Test{I} );
    pred_fit_test = LogisticFitnessDecayFunctionForOpt( fit_x0_k_L_ddGvect , aa_test(:)  ) ; 
    [ R.test_r2(I) , R.test_rmse(I)] = rsquare( fit_test , pred_fit_test )

end



% % % % %% we need to fit: 
% % % % %  x0,k,L & deltaG_vect
% % % % %  
% % % % %  given: fitness & aa_seq for each variant
% % % % %  
% % % % %  lsqcurvefit( function_to_sum_aa_dGs_and_do_logistic_on_the_sum , init_x0_k_L_ddGvect , aa_sequence_vector , measured_fitness , lb , ub )
% % % % 
% % % % % function [ fitness_diff_real_pred] = LogisticFitnessDecayFunctionForOpt( ...
% % % % %    init_x0_k_L_ddGvect , aa_for_all_variants_vect , fitness_for_all_variants_vect )
% % % % 
% % % % init_x0_k_L_ddGvect = [0 0 0 deltaG_vect'] ; 
% % % % lb = zeros( size(init_x0_k_L_ddGvect));
% % % % ub = ones( size(init_x0_k_L_ddGvect));
% % % % ub(1:3) = 100 ;
% % % % 
% % % % % lsqcurvefit
% % % % fit_opts = optimset( 'Diagnostics' , 'on' , 'Display' , 'iter-detailed', 'PlotFcn' , @optimplotx ) ;
% % % % fit_opts = optimset( 'Diagnostics' , 'off' , 'Display' , 'off', 'PlotFcn' , [] ) ;
% % % % 
% % % % fitness_for_all_variants_logistic_noise = random('normal',0,0.2,size(fitness_for_all_variants_logistic)) + fitness_for_all_variants_logistic ;
% % % % 
% % % % tic; 
% % % % [fit_x0_k_L_ddGvect , resnorm , residual , exitflag , output] = lsqcurvefit( @LogisticFitnessDecayFunctionForOpt ...
% % % %     , init_x0_k_L_ddGvect , aa_for_all_variants(:)  , fitness_for_all_variants_logistic_noise , lb , ub , fit_opts) ;
% % % % toc
% % % % %% compare model fit params to input params
% % % % txt = sprintf('corr = %0.02f' , corr(deltaG_vect , fit_x0_k_L_ddGvect(4:end)'  ) );
% % % % figure; 
% % % % subplot(2,1,1);
% % % % hold on; plot( deltaG_vect , fit_x0_k_L_ddGvect(4:end) ,'ok','DisplayName',txt); 
% % % % xlabel('Fit');ylabel('Input');title('ddG');grid on ;xlim([0 1]);ylim(xlim);line(xlim,xlim);
% % % % legend('location','best');
% % % % subplot(2,1,2);
% % % % bar( [x0 k L ; fit_x0_k_L_ddGvect(1:3)]' ); 
% % % % legend({'Orig' 'Fit'})
% % % % set(gca,'xtick',1:3);set(gca,'xticklabel',{'x0' 'k' 'L'})
% % % % 
% % % % 
