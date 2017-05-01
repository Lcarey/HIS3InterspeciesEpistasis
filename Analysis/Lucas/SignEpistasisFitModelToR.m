function  [ G , R]  = SignEpistasisFitModelToR( R , fitnessthreshold , min_N_minor_on_which_to_predict )
%% G = SignEpistasisFitModelToR( R , fitnessthreshold , min_N_minor_on_which_to_predict)
% given an R dataset (for measuring sign epistais)
%
% fit a model to it and measure sign epistais
%
% LBC April 2017

% calculate the frequency w/major & w/minor fitness effects
if ~exist('fitnessthreshold','var')
    fitnessthreshold = 0.2 ;
end
if ~exist('min_N_minor_on_which_to_predict','var')
    min_N_minor_on_which_to_predict = 50  ;
end

R.HasMajorSignFitEffect = false(length(R),1);
R.HasMinorSignFitEffect = false(length(R),1);
R.MajorSignFitEffectSize = NaN(length(R),1);
R.MinorSignFitEffectSize = NaN(length(R),1);
IDs = strcat(R.Perm,num2str(R.VarPos)) ;
uids = unique(IDs)' ;
for I = uids
    idx = strcmp(IDs,I) ;
    
    if ( sum(R.FitImpact(idx)>fitnessthreshold) > sum(R.FitImpact(idx) < (-1*fitnessthreshold) )  )% more positive
        R.HasMajorSignFitEffect(idx) = R.FitImpact(idx) > fitnessthreshold ;
        R.HasMinorSignFitEffect(idx) = R.FitImpact(idx) < (-1*fitnessthreshold) ;
    else
        R.HasMajorSignFitEffect(idx) = R.FitImpact(idx) < (-1*fitnessthreshold) ;
        R.HasMinorSignFitEffect(idx) = R.FitImpact(idx) > fitnessthreshold ;
    end
end
R.MajorSignFitEffectSize(R.HasMajorSignFitEffect) = R.FitImpact(R.HasMajorSignFitEffect) ;
R.MinorSignFitEffectSize(R.HasMinorSignFitEffect) = R.FitImpact(R.HasMinorSignFitEffect) ;

G = grpstats( R , {'VarPos' 'Perm'} , {'mean' 'sum'} ,'DataVars' , { 'MinorSignFitEffect' 'MajorSignFitEffect' 'MajorSignFitEffectSize' 'MinorSignFitEffectSize' }) ;


%% predict

%%
G.ValidationAcc_weighted = NaN(length(G),1);
G.ValidationAcc_noweights = NaN(length(G),1);
G.AUC_weighted = NaN(length(G),1);
G.AUC_noweights = NaN(length(G),1);
G.Classifier_weighted = cell(length(G),1);
G.Classifier_noweights = cell(length(G),1);

% run predictor for all substuttions w/enought minor effect backgrounds
for GI = find( G.sum_MinorSignFitEffect >  min_N_minor_on_which_to_predict )'
    Q = R( R.VarPos==G.VarPos(GI) & strcmp(R.Perm,G.Perm{GI}),:);
    % generate sparsevect for prediction
    SparseVect = cell( length(Q.Seq1) , 1);
    seq_length = length(Q.Seq1{1});
    a = cellfun( @(X) arrayfun( @(I) [ num2str(I)  X(I)] , 1:seq_length ,'UniformOutput',false) ...
        , Q.Seq1 ,'UniformOutput',false);
    a = vertcat(a{:}) ;
    
    unique_variants = unique(a(:)) ;
    for I = 1:length(Q.Seq1)
        SparseVect{I} = ismember(unique_variants , a(I,:))'  ;
    end
    idx_to_classify = Q.MajorSignFitEffect | Q.MinorSignFitEffect ;
    W = ones( length(Q) , 1) ;
    X = double(cell2mat(SparseVect)) ;
    X = X(idx_to_classify,:);
    Y = Q.MinorSignFitEffect(idx_to_classify) ;
    W(Q.MinorSignFitEffect) = sum(Q.MajorSignFitEffect)/sum(Q.MinorSignFitEffect) ;
    W = W(idx_to_classify);
    classificationKNN = fitcknn( X , Y , 'Weights' , W ,   ...
        'Distance', 'Hamming', ...
        'CategoricalPredictors','all' , ...
        'Exponent', [], ...
        'NumNeighbors', 10, ...
        'DistanceWeight', 'SquaredInverse', ...
        'Standardize', false, ...
        'ClassNames', [0; 1]);
    
    % Perform cross-validation
    partitionedModel = crossval(classificationKNN, 'KFold', 10);
    validationAccuracy = 1 - kfoldLoss(partitionedModel, 'LossFun', 'ClassifError') ;
    
    [~,~,~,AUC,~] = perfcurve( Y, partitionedModel.kfoldPredict , 1 ) ;
    
    G.ValidationAcc_weighted(GI) = validationAccuracy  ;
    G.AUC_weighted(GI) = AUC ;
    G.Classifier_weighted{GI} = partitionedModel  ;
    
    
    % w/out weights
    classificationKNN = fitcknn( X , Y  ,   ...
        'Distance', 'Hamming', ...
        'CategoricalPredictors','all' , ...
        'Exponent', [], ...
        'NumNeighbors', 10, ...
        'DistanceWeight', 'SquaredInverse', ...
        'Standardize', false, ...
        'ClassNames', [0; 1]);
    
    % Perform cross-validation
    partitionedModel = crossval(classificationKNN, 'KFold', 10);
    validationAccuracy = 1 - kfoldLoss(partitionedModel, 'LossFun', 'ClassifError') ;
    
    G.ValidationAcc_noweights(GI) = validationAccuracy  ;
    G.AUC_noweights(GI) = AUC ;
    G.Classifier_noweights{GI} = partitionedModel  ;
    
end

end
