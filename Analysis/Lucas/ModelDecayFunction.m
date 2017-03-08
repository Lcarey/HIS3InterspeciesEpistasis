%%  setup simulated data
global n_possible_aas ; 
n_positions = 5 ; 
n_possible_aas = 10 ; 

n_variants = 500 ;

% fitmat contains fitness defects for each AA at each pos
%  0 == no fitness defect
deltaG_matrix = random('uniform',0,1,n_possible_aas,n_positions) ; % fitness defect caused by each AA at each position
%deltaG_matrix = (deltaG_matrix./max(max(deltaG_matrix))).^10 ;LSQNONLIN

deltaG_vect = deltaG_matrix(:);  % this is what we want to learn

aa_for_all_variants = randi( n_possible_aas , n_variants , n_positions ) ;

sum_delta_Gs  = arrayfun( @(I)sum( deltaG_vect( ((0:(n_positions-1)) .* n_possible_aas)+aa_for_all_variants(I,:))) , 1:n_variants) ;

% generated a matrix w/the fitness cost of each variant (rows -> AAs) at each position (column)
%   the sum of the fitness costs at each position is the sum_delta_Gs
%   from this sum we can calculate the actual fitness using a sigmoid function

%
k=10 ;  % steepness
x0 = 2.5 ;  % midpoint
L = 1 ; % max Y 
logistic_function = @(x0,k,L,xdata) L ./ (1+exp( (-1.*k).*(xdata-x0) ) ) ; % x0 , k , L
logistic_function_onevar = @(x,xdata) x(3) ./ (1+exp( (-1.*x(2)).*(xdata-x(1)) ) ) ; % x0 , k , L

fitness_for_all_variants_logistic = 1-logistic_function(x0,k,L,sum_delta_Gs) ;
fitness_for_all_variants_linear = 1-( sum_delta_Gs./max(sum_delta_Gs)) ; 

fitness_for_all_variants_logistic_noise = random('normal',0,0.1,size(fitness_for_all_variants_logistic)) + fitness_for_all_variants_logistic ;
fitness_for_all_variants_linear_noise = random('normal',0,0.05,size(fitness_for_all_variants_linear)) + fitness_for_all_variants_linear ;

fitness_for_all_variants_noise = fitness_for_all_variants_logistic_noise ;

max_theo_r2 = ( corr( fitness_for_all_variants_logistic_noise(:) , fitness_for_all_variants_logistic(:)) )^2 
max_theo_r2 = ( corr( fitness_for_all_variants_linear_noise(:) , fitness_for_all_variants_linear(:)) )^2 

%%
figure; hold on; 
plot( sum_delta_Gs , fitness_for_all_variants_logistic ,'ok' , 'DisplayName' ,'Logistic')
plot( sum_delta_Gs , fitness_for_all_variants_linear ,'ob', 'DisplayName' ,'Linear')

xlabel('Sum delta Gs')
ylabel('Fitness')


%% Cross-validation w/logistic

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
    fit = fitness_for_all_variants_noise(  R.Train{I} );
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
    fit_test = fitness_for_all_variants_noise(  R.Test{I} );
    pred_fit_test = LogisticFitnessDecayFunctionForOpt( fit_x0_k_L_ddGvect , aa_test(:)  ) ; 
    [ R.test_r2(I) , R.test_rmse(I)] = rsquare( fit_test , pred_fit_test )

    R.pred_fit_train{I} = pred_fit_train ; 
    R.pred_fit_test{I} = pred_fit_test ;
    R.fit_train{I} = fit ;
    R.fit_test{I} = fit_test ; 
end

%% Linear

fit_opts = optimset( 'Diagnostics' , 'off' , 'Display' , 'off', 'PlotFcn' , [] ,'UseParallel' , false ) ;

init_m_b_ddGvect = [0 0 0 deltaG_vect'] ; 
lb = zeros( size(init_m_b_ddGvect));
ub = ones( size(init_m_b_ddGvect));
ub(1:2) = 1000 ;
lb(1:2) = -1000 ;

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
R.FitOutputS = cell(10,1);

for I = 1:10
    R.Train{I} = find( Indices ~= I );
    R.Test{I} = find( Indices == I );
end

for I = 1:10
    tic; 
    aa = aa_for_all_variants( R.Train{I} , :);
    fit = fitness_for_all_variants_noise(  R.Train{I} );
    [fit_x0_k_L_ddGvect , resnorm , residual , exitflag , output] = lsqcurvefit( @LinearFitnessDecayFunctionForOpt , init_m_b_ddGvect , aa(:)  , fit , lb , ub , fit_opts ) ;
    R.FitOutputS{I} = output ;
    R.runtime(I) = toc; 
    R.fit_ddGvect{I} = fit_x0_k_L_ddGvect(4:end) ;
    R.m(I) = fit_x0_k_L_ddGvect(1);
    R.b(I) = fit_x0_k_L_ddGvect(2);


    pred_fit_train = LinearFitnessDecayFunctionForOpt( fit_x0_k_L_ddGvect , aa(:)  ) ; 
    [ R.train_r2(I) , R.train_rmse(I)] = rsquare( fit , pred_fit_train ) ; 

    aa_test = aa_for_all_variants( R.Test{I} , :);
    fit_test = fitness_for_all_variants_noise(  R.Test{I} );
    pred_fit_test = LinearFitnessDecayFunctionForOpt( fit_x0_k_L_ddGvect , aa_test(:)  ) ; 
    [ R.test_r2(I) , R.test_rmse(I)] = rsquare( fit_test , pred_fit_test )

    R.pred_fit_train{I} = pred_fit_train ; 
    R.pred_fit_test{I} = pred_fit_test ;
    R.fit_train{I} = fit ;
    R.fit_test{I} = fit_test ; 
end


%%
figure;
subplot(2,2,1)
plot( R.fit_train{1} , R.pred_fit_train{1} ,'ok');
subplot(2,2,2)
plot( R.fit_test{1} , R.pred_fit_test{1} ,'ok');
subplot(2,2,3)
plot( R.fit_train{2} , R.pred_fit_train{2} ,'ok');
subplot(2,2,4)
plot( R.fit_test{2} , R.pred_fit_test{2} ,'ok');

