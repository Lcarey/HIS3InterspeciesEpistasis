%% Fig__SigmoidCantPredictHighDimEpiInteractions.m
%  scatter plot of how poorly the model does (residuals) vs dimensionality of sign epistasis
%  figure 4f in the intial submission. 
% LBC

% load data
DATADIR = '~/Develop/HIS3InterspeciesEpistasis/' ; 
load( [ DATADIR 'Analysis/Lucas/Data/SignEpi_Ronly.mat' ] )

BigG = readtable([ DATADIR 'Analysis/Lucas/Data/SignEpiPairs.xlsx' ] );


RESID = readtable( [ DATADIR 'Analysis/Katya/NN/residuals/S2.csv' ] );
RESID.PredictedFitness = RESID.predictedMinusObserved + RESID.observed  ;

%% which segment to work with
SegN = 2 ; 

R = s(SegN).R ;
% clear 's';

%% calculate the predicted fitness impact for each pair of seqs in R
%   each pair of seqs is a 1AA substitution (396,439 pairs)
%   
%  G will have the avg residual between actual and predict fitness impact
%   this is very very very slow
tic ; 
R.PredFitSeq1 = NaN( length(R)  , 1);
R.PredFitSeq2 = NaN( length(R)  , 1);
for I = 1:length(R)
    R.PredFitSeq1(I) = RESID.PredictedFitness( strcmp( RESID.aa_seq , R.Seq1{I})) ;
    R.PredFitSeq2(I) = RESID.PredictedFitness( strcmp( RESID.aa_seq , R.Seq2{I})) ;
end
R.PredFitImpact = R.PredFitSeq1 - R.PredFitSeq2 ;

R.AbsDiffPredReal_FitImpact = abs( R.PredFitImpact - R.FitImpact );
R.SqrDiffPredReal_FitImpact = ( R.PredFitImpact - R.FitImpact ).^2;

G = grpstats( R ,{'VarPos' 'Perm'} ,{'mean' 'median' 'std'} ,'DataVars' , {'SqrDiffPredReal_FitImpact' 'AbsDiffPredReal_FitImpact'});
toc
%% Count how often  each substitution exhibits sign epistasis. 
BigG = BigG( BigG.SegN == SegN,:);
BigG.HasSignEpi = BigG.pBon < 0.05 ;
BigGG = grpstats( BigG , {'VarPos' 'Perm' 'SubPos' 'SubPerm'} ,'sum' ,'DataVars','HasSignEpi');
BigGG = grpstats( BigGG , {'VarPos' 'Perm' } ,'sum' ,'DataVars','sum_HasSignEpi');

%% calculate correlation

Q = innerjoin( BigGG , dataset2table(G) ,'Key',{'VarPos' 'Perm'});

X = Q.sum_sum_HasSignEpi ;
Y = Q.mean_SqrDiffPredReal_FitImpact ; 
fh = figure('units','centimeters','position',[5 5 3.5 4]);
plot( X , Y  ,'ok','MarkerFaceColor',[.7 .7 .7])
[c,p] = corr(X,Y);
[cS,pS] = corr(X,Y,'Type','Spearman')
xlim([0 15])
ylim([0.04 .11])
%xlabel('# of substitutions showing sign epistasis per site')
%ylabel('Squared residual of predicted fitness impact vs measured fitness impact')