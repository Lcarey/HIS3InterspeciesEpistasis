function [R , T ] = ModelDecayFunctionFromRealData( data_file_name  ,  varargin )
%% [R , T ] = ModelDecayFunctionFromRealData( 7 ,'NO_SINGLES_FROM_ANY_OTHER',true,'ONLY_NATLIB_FLAG',false,'ONLY_LIB_FLAG',true);
%
% LBC March 2017


p = inputParser;
addRequired(p,'data_file_name');
addOptional(p,'N_variants_to_fit',999999999,@isnumeric)
addOptional(p,'RUN_LINEAR_FLAG',0,@islogical)
addOptional(p,'ONLY_NATLIB_FLAG',1,@islogical)
addOptional(p,'ONLY_LIB_FLAG',1,@islogical)
addOptional(p,'NO_STOP_FLAG',1,@islogical)
addOptional(p,'ParamID','XXXX',@ischar)
addOptional(p,'KFOLD','10',@isnumeric)

addOptional(p,'NO_SINGLES_FROM_ANY_OTHER',0,@islogical)


parse(p,data_file_name,varargin{:});

MapAA2I = containers.Map(  arrayfun(@(X){X},['A':'Z' '_' ])  , 1:27 ) ; % all AAs + stop


%% setup paths
addpath(genpath('~/Develop/matlab/'));
addpath(genpath('~/Develop/Matlab/'));
addpath(genpath('~/Develop/HIS3InterspeciesEpistasis/'));


% if numel(gcp('nocreate'))==0 , 	parpool local , end
%% load data
% data_file_name can be a segment # or file name
DataDir = '~/Develop/HIS3InterspeciesEpistasis/Data/';
if isnumeric(data_file_name)
    data_file_name = [ 'S' num2str(data_file_name) '_scaled_info.csv'];
end
T = readtable([ DataDir data_file_name ] );

T = T( logical(T.nogap)  , :) ; % can't handle seq's w/varying length
T = T( T.len == mode(T.len) , :); % can't handle seq's w/varying length

if p.Results.NO_STOP_FLAG
	T = T( T.middle & T.nogap & ~T.stop & ~T.nonsense , :) ; 
end
if p.Results.ONLY_NATLIB_FLAG
	T = T( logical(T.nat_lib) , :);
end
if p.Results.ONLY_LIB_FLAG
	T = T( logical(T.lib) , :);
end



%%
N_variants_to_fit = min(p.Results.N_variants_to_fit,height(T));
T = T( 1:N_variants_to_fit , :);

output_file_basename = [ regexprep( char(datetime) , ' ' ,'') '_' num2str(N_variants_to_fit) '_' data_file_name '_' ] ;;
if exist('ParamID','var')  % optional argument for tagging various runs
	if isnumeric(p.ParamID)
		p.Results.ParamID = num2str(p.Results.ParamID);
	end
	output_file_basename = [ p.Results.ParamID '_' output_file_basename ];
end

%% shrink sequences down to the variable part
[ T  , columns_that_vary , uniq_aas ] = EpistasisLBC__ShortedAAseqsOnlyPositionsThatVary( T );
T.columns_that_vary = repmat( uint8(find(columns_that_vary)) , height(T) , 1) ; 


%%  there are VERY few seqs that are 1 away from any other
if p.Results.NO_SINGLES_FROM_ANY_OTHER
    T.aa_seq_variable_numeric = NaN( height(T) , numel(T.aa_seq_variable{1}));
    for I = 1:height(T)
        T.aa_seq_variable_numeric(I,:) = arrayfun( @(X)MapAA2I(X) ,  T.aa_seq_variable{I});
    end
    stalling = 0 ;
    oht = height(T);
    while stalling<1e3
        for I = randi( height(T) , 1)
            neighbors = pdist2( T.aa_seq_variable_numeric , T.aa_seq_variable_numeric(I,:) ,'Hamming')*numel(T.aa_seq_variable{1}) ;
            T = T( neighbors>1 | neighbors==0,:);
        end
        %   fprintf('%d\n' , oht - height(T));
        if  (oht - height(T)) == 0 % no change
            stalling = stalling + 1;
        end
        oht = height(T);
    end
    
end

%% setup variables from data
global n_possible_aas ; 
n_positions = sum(columns_that_vary) ;
all_aas_all_positions = cell2mat(T.aa_seq_variable(:)) ; all_aas_all_positions = unique(all_aas_all_positions(:) )  ;
n_possible_aas = numel( all_aas_all_positions );
n_variants = height(T) ; 

%  map AA letters to numbers
n_possible_aas = 27 ;
% MapAA2I = containers.Map( arrayfun(@(X){X},all_aas_all_positions) ,1:n_possible_aas) ; % only AAs in this lib. harder to map back

aa_for_all_variants = NaN( n_variants , n_positions );
for I = 1:n_variants
    aa_for_all_variants(I,:) = arrayfun( @(X)MapAA2I(X) ,  T.aa_seq_variable{I});
end

fitness_for_all_variants = ( T.s)'  ; % XDATA & function output have to be same size & shape vector

%% fit logistic model w/ cross-validation

fit_opts = optimset( 'Diagnostics' , 'off' , 'Display' , 'off', 'PlotFcn' , [] ,'UseParallel' , false ) ;
kfold = p.Results.KFOLD ; 
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
	
	deltaG_matrix_init = repmat( 0.5 , n_possible_aas,n_positions);
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

	fprintf('%d\t%d\t%fsec\t%0.02f\n' , I , height(T) , R.runtime(I) , R.test_r2(I) );
	save( ['R_logistic_' output_file_basename '_' num2str(I) '.mat'] ,'R') ;
	if(I>1) , delete(  ['R_logistic_' output_file_basename '_' num2str(I-1) '.mat'] ); end
end

R
delete(  ['R_logistic_' output_file_basename '_' num2str(kfold) '.mat'] );
save( ['R_logistic_' output_file_basename '.mat'] ,'R')



end
