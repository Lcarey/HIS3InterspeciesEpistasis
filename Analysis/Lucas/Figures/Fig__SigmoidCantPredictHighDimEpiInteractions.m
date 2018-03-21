function Q = Fig__SigmoidCantPredictHighDimEpiInteractions( SegN , Npairs , pTHRESHOLD )
%% Q = Fig__SigmoidCantPredictHighDimEpiInteractions( SegN , Npairs , pTHRESHOLD)
%  scatter plot of how poorly the model does (residuals) vs dimensionality of sign epistasis
%  figure 4f in the intial submission. 
% 
%
% LBC

% load data
%  SignEpi_Ronly is a struct array of tables. Each table has the FitnessImpact for each substitution. eg: the first five subsitutions in segment 2:
%  s(2).R(1:5,:)
%
%    Seq1                                  Seq2                                  Fit1       Fit2       VarPos    Perm        FitImpact
%    'SDRAFAFIDLGLQREKVGDLSCEMIPHF'        'SGRAFAFIDLGLQREKVGDLSCEMIPHF'         1.1222      1.012    2         'DG'          0.11023
%    'SDRAFAFIDLGLQREKVGDLSCEMIPHI'        'SGRAFAFIDLGLQREKVGDLSCEMIPHI'        0.82869    0.97377    2         'DG'         -0.14508
%    'SDRAFAFIDLGLQREKVGDLSCEMIPHL'        'SGRAFAFIDLGLQREKVGDLSCEMIPHL'        0.90375     0.9998    2         'DG'        -0.096051
%    'SDRAFAFIDLGLQREKVGDLSCEMIPHV'        'SGRAFAFIDLGLQREKVGDLSCEMIPHV'        0.46428     1.0583    2         'DG'         -0.59407
%    'SDRAFAFIDLGLQREKVGDLSCEMITHV'        'SGRAFAFIDLGLQREKVGDLSCEMITHV'        0.74974     1.0626    2         'DG'         -0.31288

% BigG (SignEpiPairs.xlsx) has the p values for each pair of substitutions saying if VarPos-Perm is affected by SubPos-SubPerm
%  SegN	VarPos	Perm	SubPos	SubPerm	p	logodds
%  2	6	CS	8	IT	2.0444E-246	0.009041752
%  2	6	AC	8	IT	2.1484E-126	0.031832934
%  2	8	TV	28	FV	1.4223E-115	0.003212437

% RESID has predicted & observed fitness for each genotype
 
DATADIR = '~/Develop/HIS3InterspeciesEpistasis/' ; 
load( [ DATADIR 'Analysis/Lucas/Data/SignEpi_Ronly.mat' ] )

BigG = readtable([ DATADIR 'Analysis/Lucas/Data/SignEpiPairs.xlsx' ] );

RESID = readtable( [ DATADIR 'Analysis/Katya/NN/residuals/S' num2str(SegN) '.csv' ] );
RESID.PredictedFitness = RESID.predictedMinusObserved + RESID.observed  ;

% which segment to work with

R = s(SegN).R ;
% if passed Npairs to speed things up.
if exist('Npairs','var')
if Npairs < size(R,1)
    R = R( randsample( size(R,1) , Npairs) , : );
end
end


% clear 's';

%% calculate the predicted fitness impact for each pair of seqs in R
%   each pair of seqs is a 1AA substitution (396,439 pairs)
%   
%  G will have the avg residual between actual and predict fitness impact
%   this is very very very slow
R.PredFitSeq1 = NaN( length(R)  , 1);
R.PredFitSeq2 = NaN( length(R)  , 1);
tic ; 
for I = 1:size(R,1)
    R.PredFitSeq1(I) = RESID.PredictedFitness( strcmp( RESID.aa_seq , R.Seq1{I})) ;
    R.PredFitSeq2(I) = RESID.PredictedFitness( strcmp( RESID.aa_seq , R.Seq2{I})) ;
end
toc
R.PredFitImpact = R.PredFitSeq1 - R.PredFitSeq2 ;
R.AbsDiffPredReal_FitImpact = abs( R.PredFitImpact - R.FitImpact );
R.SqrDiffPredReal_FitImpact = ( R.PredFitImpact - R.FitImpact ).^2;

G = grpstats( R ,{'VarPos' 'Perm'} ,{'mean' 'median' 'std'} ,'DataVars' , {'SqrDiffPredReal_FitImpact' 'AbsDiffPredReal_FitImpact'});
save( [ '~/Downloads/Fig__SigmoidCantPredictHighDimEpiInteractions__G__' num2str(SegN) '.mat'] , 'G' , 'BigG' ,'SegN'  );
%% Count how often  each substitution exhibits sign epistasis. 
BigG = BigG( BigG.SegN == SegN,:);
if ~exist( 'pTHRESHOLD','var')
    BigG.HasSignEpi = BigG.pBon < 0.05 ;
else
    BigG.HasSignEpi = BigG.p < pTHRESHOLD ;
end

BigGG = grpstats( BigG , {'VarPos' 'Perm' 'SubPos' 'SubPerm'} ,'sum' ,'DataVars','HasSignEpi');
BigGG = grpstats( BigGG , {'VarPos' 'Perm' } ,'sum' ,'DataVars','sum_HasSignEpi');

%% calculate correlation

Q = innerjoin( BigGG , dataset2table(G) ,'Key',{'VarPos' 'Perm'});
Q.SegN = repmat( SegN , height(Q) , 1) ; 
Q.DateTimeForThisAnalysis = repmat( datetime() , height(Q) , 1) ; 

X = Q.sum_sum_HasSignEpi ;
Y = Q.mean_SqrDiffPredReal_FitImpact ; 


fh = figure('units','centimeters','position',[5 5 5 5]);
plot( X , Y  ,'ok','MarkerFaceColor',[.7 .7 .7])
[c,p] = corr(X,Y);
[cS,pS] = corr(X,Y,'Type','Spearman')
%xlim([0 15])
%ylim([0.04 .11])
txt = sprintf('%d %0.02f %0.02e' , SegN , c , p) ;
title( txt )
%xlabel('# of substitutions showing sign epistasis per site')
%ylabel('Squared residual of predicted fitness impact vs measured fitness impact')
