%% load data
start_runtime = char(datetime) ;
DataDir = '~/Develop/HIS3InterspeciesEpistasis/Data/';
T = readtable([ DataDir 'S12_scaled_info.csv']);
T = T( T.middle & T.nogap & ~T.stop , :) ; % no crap
T = T( logical(T.nat_lib) , :);
T = T( 1:2000 , :);
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

%% fit logistic model w/ cross-validation

fit_opts = optimset( 'Diagnostics' , 'off' , 'Display' , 'off', 'PlotFcn' , [] ,'UseParallel' , true ) ;

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

R.pred_fit_train = cell(10,1);
R.pred_fit_test = cell(10,1);
R.fit_train = cell(10,1);
R.fit_test = cell(10,1);

for I = 1:10
    R.Train{I} = find( Indices ~= I );
    R.Test{I} = find( Indices == I );
end

for I = 1:10
    fprintf('.');
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

    pred_fit_train = LogisticFitnessDecayFunctionForOpt( fit_x0_k_L_ddGvect , aa(:)  ) ; 
    [ R.train_r2(I) , R.train_rmse(I)] = rsquare( fit , pred_fit_train ) ; 

    aa_test = aa_for_all_variants( R.Test{I} , :);
    fit_test = fitness_for_all_variants(  R.Test{I} );
    pred_fit_test = LogisticFitnessDecayFunctionForOpt( fit_x0_k_L_ddGvect , aa_test(:)  ) ; 
    [ R.test_r2(I) , R.test_rmse(I)] = rsquare( fit_test , pred_fit_test )

    
    R.pred_fit_train{I} = pred_fit_train ; 
    R.pred_fit_test{I} = pred_fit_test ;
    R.fit_train{I} = fit ;
    R.fit_test{I} = fit_test ; 
end

R
save( ['R_logistic' start_runtime '.mat'] ,'R')


%% fit linear model w/ cross-validation

fit_opts = optimset( 'Diagnostics' , 'off' , 'Display' , 'off', 'PlotFcn' , [] ,'UseParallel' , true ) ;

deltaG_matrix_init = random('uniform',0,1,n_possible_aas,n_positions) ; % fitness defect caused by each AA at each position
deltaG_vect_init = deltaG_matrix_init(:);  % this is what we want to learn

init_x0_k_L_ddGvect = [0 0 0 deltaG_vect_init'] ; 
lb = zeros( size(init_x0_k_L_ddGvect));
ub = ones( size(init_x0_k_L_ddGvect));
ub(1:2) = 1000 ;
ub(1:2) = -1000 ;

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

R.pred_fit_train = cell(10,1);
R.pred_fit_test = cell(10,1);
R.fit_train = cell(10,1);
R.fit_test = cell(10,1);

for I = 1:10
    R.Train{I} = find( Indices ~= I );
    R.Test{I} = find( Indices == I );
end

for I = 1:10
    fprintf('.');
    tic; 
    aa = aa_for_all_variants( R.Train{I} , :);
    fit = fitness_for_all_variants(  R.Train{I} ) ; 
    [fit_x0_k_L_ddGvect , resnorm , residual , exitflag , output] = lsqcurvefit( @LinearFitnessDecayFunctionForOpt , init_x0_k_L_ddGvect , aa(:)  , fit , lb , ub , fit_opts ) ;
    R.FitOutputS{I} = output ;
    R.runtime(I) = toc; 
    R.fit_ddGvect{I} = fit_x0_k_L_ddGvect(4:end) ;
    R.x0(I) = fit_x0_k_L_ddGvect(1);
    R.k(I) = fit_x0_k_L_ddGvect(2);
    R.L(I) = fit_x0_k_L_ddGvect(3);

    pred_fit_train = LinearFitnessDecayFunctionForOpt( fit_x0_k_L_ddGvect , aa(:)  ) ; 
    [ R.train_r2(I) , R.train_rmse(I)] = rsquare( fit , pred_fit_train ) ; 

    aa_test = aa_for_all_variants( R.Test{I} , :);
    fit_test = fitness_for_all_variants(  R.Test{I} );
    pred_fit_test = LinearFitnessDecayFunctionForOpt( fit_x0_k_L_ddGvect , aa_test(:)  ) ; 
    [ R.test_r2(I) , R.test_rmse(I)] = rsquare( fit_test , pred_fit_test )

    
    R.pred_fit_train{I} = pred_fit_train ; 
    R.pred_fit_test{I} = pred_fit_test ;
    R.fit_train{I} = fit ;
    R.fit_test{I} = fit_test ; 
end

R
save( ['R_linear' start_runtime '.mat'] ,'R')