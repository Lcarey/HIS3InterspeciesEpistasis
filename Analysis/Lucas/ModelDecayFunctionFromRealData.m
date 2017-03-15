function ModelDecayFunctionFromRealData( data_file_name , N_variants_to_fit , ParamID )

% flags to control behavior
RUN_LINEAR_FLAG = false   ; 
ONLY_NATLIB_FLAG = true  ; 
ONLY_LIB_FLAG  = true  ; 
NO_STOP_FLAG = true ; 

%% setup paths
addpath(genpath('~/Develop/matlab/'));
addpath(genpath('~/Develop/Matlab/'));
addpath(genpath('~/Develop/HIS3InterspeciesEpistasis/'));


% if numel(gcp('nocreate'))==0 , 	parpool local , end
%% load data
DataDir = '~/Develop/HIS3InterspeciesEpistasis/Data/';
T = readtable([ DataDir data_file_name ] );
T = T( logical(T.nogap)  , :) ; % can't handle seq's w/varying length
T = T( T.len == mode(T.len) , :); % can't handle seq's w/varying length

if NO_STOP_FLAG
	T = T( T.middle & T.nogap & ~T.stop & ~T.nonsense , :) ; 
end
if ONLY_NATLIB_FLAG
	T = T( logical(T.nat_lib) , :);
end
if ONLY_LIB_FLAG
	T = T( logical(T.lib) , :);
end

N_variants_to_fit = min(N_variants_to_fit,height(T));
T = T( 1:N_variants_to_fit , :);

output_file_basename = [ regexprep( char(datetime) , ' ' ,'') '_' num2str(N_variants_to_fit) '_' data_file_name '_' ] ;;
if exist('ParamID','var') , % optional argument for tagging various runs
	if isnumeric(ParamID)
		ParamID = num2str(ParamID);
	end
	output_file_basename = [ ParamID '_' output_file_basename ];
end

%% setup variables from data
global n_possible_aas ; 
n_positions = T.len(1) ;
all_aas_all_positions = cell2mat(T.aa_seq(:)) ; all_aas_all_positions = unique(all_aas_all_positions(:) )  ;
n_possible_aas = numel( all_aas_all_positions );
n_variants = height(T) ; 

%  map AA letters to numbers
MapAA2I = containers.Map(  arrayfun(@(X){X},['A':'Z' '_' ])  , 1:27 ) ; % all AAs + stop
n_possible_aas = 27 ;
% MapAA2I = containers.Map( arrayfun(@(X){X},all_aas_all_positions) ,1:n_possible_aas) ; % only AAs in this lib. harder to map back

aa_for_all_variants = NaN( n_variants , n_positions );
for I = 1:n_variants
    aa_for_all_variants(I,:) = arrayfun( @(X)MapAA2I(X) ,  T.aa_seq{I});
end

fitness_for_all_variants = ( T.s)'  ; % XDATA & function output have to be same size & shape vector

%% fit logistic model w/ cross-validation

fit_opts = optimset( 'Diagnostics' , 'off' , 'Display' , 'off', 'PlotFcn' , [] ,'UseParallel' , false ) ;
kfold = 10 ; 
Indices = crossvalind('Kfold', n_variants, kfold);
R = dataset();
R.fit_ddGvect = cell(kfold,1);
R.Train = cell(kfold,1);
R.Test = cell(kfold,1);
R.train_rmse = NaN(kfold,1);
R.train_r2 = NaN(kfold,1) ;
R.test_rmse = NaN(kfold,1);
R.test_r2 = NaN(kfold,1) ;
R.runtime = R.test_r2 ; 
R.x0 = R.test_r2 ; 
R.k = R.test_r2 ; 
R.L = R.test_r2 ; 
R.FitOutputS = cell(kfold,1);

R.pred_fit_train = cell(kfold,1);
R.pred_fit_test = cell(kfold,1);
R.fit_train = cell(kfold,1);
R.fit_test = cell(kfold,1);
R.train_sumddG = cell(kfold,1);
R.test_sumddG = cell(kfold,1);
R.aa_for_all_variants = cell(kfold,1);
R.all_aas_all_positions = cell(kfold,1);

for I = 1:kfold
    R.Train{I} = find( Indices ~= I );
    R.Test{I} = find( Indices == I );
	R.aa_for_all_variants{I} = uint8(aa_for_all_variants);
	R.all_aas_all_positions{I} = all_aas_all_positions ;
end

for I = 1:kfold
	
	deltaG_vect_init = repmat( 0.5 , n_possible_aas,n_positions);
	%deltaG_matrix_init = random('uniform',0.01,0.99,n_possible_aas,n_positions) ; % fitness defect caused by each AA at each position
	deltaG_vect_init = deltaG_matrix_init(:);  % this is what we want to learn

	init_x0_k_L_ddGvect = [0 0 0 deltaG_vect_init'] ; 
	lb = ones( size(init_x0_k_L_ddGvect)) * 0 ; % (or -1 or -10)
	ub = ones( size(init_x0_k_L_ddGvect)) * 1  ;% ( or 10)
	ub(1:3) = [10 100 10] ;	lb(1:3) = [0 0 0] ;

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

    [pred_fit_train , R.train_sumddG{I}] = LogisticFitnessDecayFunctionForOpt( fit_x0_k_L_ddGvect , aa(:)  ) ; 
    [ R.train_r2(I) , R.train_rmse(I)] = rsquare( fit , pred_fit_train ) ; 

    aa_test = aa_for_all_variants( R.Test{I} , :);
    fit_test = fitness_for_all_variants(  R.Test{I} );
    [pred_fit_test , R.test_sumddG{I}] = LogisticFitnessDecayFunctionForOpt( fit_x0_k_L_ddGvect , aa_test(:)  ) ; 
    [ R.test_r2(I) , R.test_rmse(I)] = rsquare( fit_test , pred_fit_test ) ;

    
    R.pred_fit_train{I} = pred_fit_train ; 
    R.pred_fit_test{I} = pred_fit_test ;
    R.fit_train{I} = fit ;
    R.fit_test{I} = fit_test ; 

	fprintf('%d\t%fsec\n' , I , R.runtime(I) );
	save( ['R_logistic_' output_file_basename '_' num2str(I) '.mat'] ,'R') ;
	delete(  ['R_logistic_' output_file_basename '_' num2str(I-1) '.mat'] );
end

R
delete(  ['R_logistic_' output_file_basename '_' num2str(kfold) '.mat'] );
save( ['R_logistic_' output_file_basename '.mat'] ,'R')


%% fit linear model w/ cross-validation
if RUN_LINEAR_FLAG

fit_opts = optimset( 'Diagnostics' , 'off' , 'Display' , 'off', 'PlotFcn' , [] ,'UseParallel' , false ) ;

deltaG_matrix_init = random('uniform',0,1,n_possible_aas,n_positions) ; % fitness defect caused by each AA at each position
deltaG_vect_init = deltaG_matrix_init(:);  % this is what we want to learn

%  init_x0_k_L_ddGvect = [0 0 0 deltaG_vect_init'] ; 
%  lb = zeros( size(init_x0_k_L_ddGvect));
%  ub = ones( size(init_x0_k_L_ddGvect));
%  ub(1:2) = 100 ;
%  ub(1:2) = -100 ;

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
    fprintf(2,'.');
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
    [ R.test_r2(I) , R.test_rmse(I)] = rsquare( fit_test , pred_fit_test ) ;

    
    R.pred_fit_train{I} = pred_fit_train ; 
    R.pred_fit_test{I} = pred_fit_test ;
    R.fit_train{I} = fit ;
    R.fit_test{I} = fit_test ; 

	fprintf(2,'%d\t%fsec\tR^2=%0.02f %0.02f\n' , I , R.runtime(I) ,R.test_r2(I),R.train_r2(I));
	save( ['R_linear_' output_file_basename  '_' num2str(I)  '.mat'] ,'R')
end

R
save( ['R_linear_' output_file_basename '.mat'] ,'R')

end
