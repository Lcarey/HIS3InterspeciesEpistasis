%% Fig__SigmoidCantPredictHighDimEpiInteractions.m
load SignEpi_Ronly.mat
RESID = readtable('/Users/lcarey/Develop/HIS3InterspeciesEpistasis/Analysis/Katya/NN/residuals/S2.csv');
RESID.PredictedFitness = RESID.predictedMinusObserved + RESID.observed  ;

R = s(2).R ;
clear 's';

%%
R.PredFitSeq1 = NaN( length(R)  , 1);
R.PredFitSeq2 = NaN( length(R)  , 1);
for I = 1:length(R)
    R.PredFitSeq1(I) = RESID.PredictedFitness( strcmp( RESID.aa_seq , R.Seq1{I})) ;
    R.PredFitSeq2(I) = RESID.PredictedFitness( strcmp( RESID.aa_seq , R.Seq2{I})) ;
end
R.PredFitImpact = R.PredFitSeq1 - R.PredFitSeq2 ;
%%
R.AbsDiffPredReal_FitImpact = abs( R.PredFitImpact - R.FitImpact );
R.SqrDiffPredReal_FitImpact = ( R.PredFitImpact - R.FitImpact ).^2;

G = grpstats( R ,{'VarPos' 'Perm'} ,{'mean' 'median' 'std'} ,'DataVars' , {'SqrDiffPredReal_FitImpact' 'AbsDiffPredReal_FitImpact'});

%%
BigG = readtable('SignEpiPairs.xlsx');
BigG = BigG( BigG.SegN == 2,:);
BigG.HasSignEpi = BigG.pBon < 0.05 ;
BigGG = grpstats( BigG , {'VarPos' 'Perm' 'SubPos' 'SubPerm'} ,'sum' ,'DataVars','HasSignEpi');
BigGG = grpstats( BigGG , {'VarPos' 'Perm' } ,'sum' ,'DataVars','sum_HasSignEpi');

%%
Q = join( BigGG , dataset2table(G) ,'Key',{'VarPos' 'Perm'});

%%
X = Q.sum_sum_HasSignEpi ;
fh = figure('units','centimeters','position',[5 5 3.5 4]);
plot( X , Q.mean_SqrDiffPredReal_FitImpact  ,'ok','MarkerFaceColor',[.7 .7 .7])
[c,p] = corr(  Q.mean_SqrDiffPredReal_FitImpact , X );
[cS,pS] = corr(  Q.mean_SqrDiffPredReal_FitImpact , X ,'Type','Spearman')
xlim([0 15])
ylim([0.04 .11])
%xlabel('# of substitutions showing sign epistasis per site')
%ylabel('Squared residual of predicted fitness impact vs measured fitness impact')