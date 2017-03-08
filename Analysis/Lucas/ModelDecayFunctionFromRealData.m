%% load data
DataDir = '~/Develop/HIS3InterspeciesEpistasis/Data/';
T = readtable([ DataDir 'S12_scaled_info.csv']);
T = T( T.middle & T.nogap & ~T.stop , :) ; % no crap
T = T( logical(T.nat_lib) , :);
T = T( 1:10000 , :);
%% setup variables from data
global n_possible_aas ; 
n_positions = T.len(1) ;
all_aas_all_positions = cell2mat(T.aa_seq(:)) ; all_aas_all_positions = unique(all_aas_all_positions(:) )  ;
n_possible_aas = numel( all_aas_all_positions );
n_variants = height(T) ; 

MapAA2I = containers.Map( arrayfun(@(X){X},all_aas_all_positions) ,1:n_possible_aas) ;
aa_for_all_variants = NaN( n_variants , n_positions );
for I = 1:n_variants
    aa_for_all_variants(I,:) = arrayfun( @(X)MapAA2I(X) ,  T.aa_seq{I});
end

fitness_for_all_variants = ( 1-T.s)'  ; % XDATA & function output have to be same size & shape vector

%% fit model w/ cross-validation

fit_opts = optimset( 'Diagnostics' , 'off' , 'Display' , 'off', 'PlotFcn' , [] ,'UseParallel' , false ) ;

deltaG_matrix_init = random('uniform',0,1,n_possible_aas,n_positions) ; % fitness defect caused by each AA at each position
deltaG_vect_init = deltaG_matrix_init(:);  % this is what we want to learn

init_x0_k_L_ddGvect = [0 0 0 deltaG_vect_init'] ; 
lb = zeros( size(init_x0_k_L_ddGvect));
ub = ones( size(init_x0_k_L_ddGvect));
ub(1:3) = 100 ;

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
    fit = fitness_for_all_variants(  R.Train{I} ) ; 
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
    fit_test = fitness_for_all_variants(  R.Test{I} );
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
