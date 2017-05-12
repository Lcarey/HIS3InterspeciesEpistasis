%% Given two fit genotypes, what is the fitness & fitness potential distribution of measured intermediate genotypes? 
%    if the two anchor genotypes are near the edge of the sigmoid, then
%    intermediate states will be unfit
%
% LBC May 11, 2017
%

%% load data
SegN = 7 ;
cd('/Users/lcarey/Desktop/HIS3scratch/CliffCausesUnfitIntermediateStates/');

T = readtable( [ '~/Develop/HIS3InterspeciesEpistasis/Data/S' num2str(SegN) '_scaled_info_v2.csv'] , 'FileType','text','Delimiter','\t');
T = T( logical(T.nat_lib) & logical(T.nogap) & logical(T.middle) & ~logical(T.stop) & ~logical(T.nonsense) ,:);
[ T.aa_seq_variable ,   columns_that_vary ,  uniq_aas ] = EpistasisLBC__ShortedAAseqsOnlyPositionsThatVary( T ) ;
%T.SparseVect = SparseVect ; 
%[ SparseVect , unique_variants]  = EpistasisLBC__FindVariableRegion_GenSparseMatrix(  aa_seq_variable ) ;
NNFP = readtable( [ '~/Develop/HIS3InterspeciesEpistasis/Analysis/Katya/NN/residuals/S' num2str(SegN) '.csv'] , 'FileType','text','Delimiter',',');
T = innerjoin( T , NNFP ,'Key','aa_seq');
T.predicted_Fitness = T.predictedMinusObserved + T.observed ; 
%% calculate hamming distance matrix

MapAA2I = containers.Map(  arrayfun(@(X){X},['A':'Z' '_' ])  , uint8(1:27) ) ; % all AAs + stop
%MapI2AA = containers.Map( uint8(1:27) ,   arrayfun(@(X){X},['A':'Z' '_' ])  ) ; % all AAs + stop
T.aa_num = cell( height(T) , 1);
for I = 1:height(T)
    T.aa_num{I} = arrayfun( @(X)MapAA2I(X) , T.aa_seq_variable{I});
end

distance_matrix = squareform(pdist( cell2mat(T.aa_num) ,'hamming')) .* length(T.aa_num{1}) ;

%% pick random pairs
% used for picking anchor genotypes
FitThreshold = 0.95 ; 
UnFitThreshold = 0.4 ;


%% pick fit pairs
idx_fit = find(T.s >= FitThreshold)  ; 
N_Random_Fit_Pairs = 3e4;
chosen_fit_pairs = NaN( N_Random_Fit_Pairs , 2);
for I = 1:N_Random_Fit_Pairs
    chosen_fit_pairs(I,:) = sort(randsample( idx_fit , 2)') ; 
end
chosen_fit_pairs = unique( chosen_fit_pairs , 'rows') ; 

% Find Intermediate States
DS = dataset();
DS.IntermediateStates = cell( nrows(chosen_fit_pairs) , 1);
DS.IntermediateStates_idx = cell( nrows(chosen_fit_pairs) , 1);
DS.Seq1_idx = NaN( nrows(chosen_fit_pairs) , 1);
DS.Seq2_idx = NaN( nrows(chosen_fit_pairs) , 1);
DS.IntermediateStates_fitness = cell( nrows(chosen_fit_pairs) , 1);
DS.IntermediateStates_FitnessPotential = cell( nrows(chosen_fit_pairs) , 1);
DS.Seq1_FP = NaN( nrows(chosen_fit_pairs) , 1);
DS.Seq2_FP = NaN( nrows(chosen_fit_pairs) , 1);
DS.HD =  NaN( nrows(chosen_fit_pairs) , 1);

for I = 1:length(DS.IntermediateStates)
    s1 = T.aa_seq_variable{ chosen_fit_pairs(I,1) } ;
    s2 = T.aa_seq_variable{ chosen_fit_pairs(I,2) } ;
    DS.Seq1_idx(I) = find( strcmp( T.aa_seq_variable , s1) ) ;
    DS.Seq2_idx(I) = find( strcmp( T.aa_seq_variable , s2) ) ;
    IntermediateStates  = ExpandSeqAlign( s1 , s2 );
    DS.IntermediateStates{I} = IntermediateStates( ~ismember(IntermediateStates,s1) & ~ismember(IntermediateStates,s2) ) ;
    DS.IntermediateStates_idx{I} = find( ismember( T.aa_seq_variable , DS.IntermediateStates{I} ) );
    DS.IntermediateStates_fitness{I} = T.s(  DS.IntermediateStates_idx{I} );
    
    DS.IntermediateStates_FitnessPotential{I} = T.fitnessPotential(  DS.IntermediateStates_idx{I} );
    DS.Seq1_FP(I) = T.fitnessPotential( DS.Seq1_idx(I)) ;
    DS.Seq2_FP(I) = T.fitnessPotential( DS.Seq2_idx(I)) ;
    
    DS.HD(I) = distance_matrix(  DS.Seq1_idx(I) ,  DS.Seq2_idx(I));
end


DS.ParentStates_MeanFitPot = mean([ DS.Seq1_FP DS.Seq2_FP ] , 2) ;
DS.ParentStates_MaxFitPot = max([ DS.Seq1_FP DS.Seq2_FP ] , [], 2) ;
DS.ParentStates_MaxFitPot = max([ DS.Seq1_FP DS.Seq2_FP ] , [], 2) ;

DS.IntermediateStates_PctUnfit = cellfun( @(X) mean( X <= UnFitThreshold )* 100 , DS.IntermediateStates_fitness) ; 
DS.IntermediateStates_PctFit = cellfun( @(X) mean( X >= FitThreshold )* 100 , DS.IntermediateStates_fitness) ; 
DS.IntermediateStates_PctIntfit = 100 - DS.IntermediateStates_PctFit - DS.IntermediateStates_PctUnfit ;

DS = sortrows(DS,'IntermediateStates_PctUnfit','descend');
DS = DS( ~isnan(DS.IntermediateStates_PctUnfit) ,:);
save( [ 'CliffCausesUnfitIntermediateStates__Seg' num2str(SegN) ] ,'DS'  )
%%
%% Pretty Pictures

Q = DS( cellfun(@length,DS.IntermediateStates_fitness) > 100 ,:);


figname = 'IntermediateStatesFitnessCliff.eps';
delete(figname);
idxs_to_plot = unique([ find(Q.IntermediateStates_PctUnfit>20)' 1:10:length(Q)]) ;
idxs_to_plot = randsample( idxs_to_plot , 3);
X = T.fitnessPotential ; 
Y = T.s ;
% a few examples
for I = idxs_to_plot
    fh = figure('units','centimeters','position',[5 5 5 5]);
    hold on ; 
    dscatter( X , Y );
    colormap(flipud(gray(100)));
    axis tight;
    plot( X( Q.IntermediateStates_idx{I} ) , Y( Q.IntermediateStates_idx{I} ) , '.','Color',[.5 0 .5]);
    plot( X( Q.Seq1_idx(I)) , Y( Q.Seq1_idx(I) ) ,'o','MarkerFaceColor','b','Color','b');
    plot( X( Q.Seq2_idx(I)) , Y( Q.Seq2_idx(I) ) ,'o','MarkerFaceColor','r','Color','r');
    set(gca,'xtick',[]);
    set(gca,'ytick',[]);
    xlabel('Fitness potential');
    ylabel('Fitness');
    
    
    line( xlim , [FitThreshold FitThreshold],'LineStyle','--','Color',[.7 .7 .7])
    line( xlim , [UnFitThreshold UnFitThreshold],'LineStyle','--','Color',[.7 .7 .7])
    
    title( sprintf('%0.0f%% %0.0f%%  %0.0f%% ' , Q.IntermediateStates_PctUnfit(I)  , Q.IntermediateStates_PctIntfit(I) , Q.IntermediateStates_PctFit(I)  ) );
    print('-dpsc2',figname,'-append');
    close
end

%%
for I = 0:10
    idx = DS.HD>I;
    if sum(idx)>10
         c1=corr( DS.ParentStates_MeanFitPot(idx) , DS.IntermediateStates_PctUnfit(idx) );
         c2=corr( DS.ParentStates_MaxFitPot(idx) , DS.IntermediateStates_PctUnfit(idx) );
         c3=corr( mean( [ DS.ParentStates_MaxFitPot(idx)   DS.ParentStates_MeanFitPot(idx)] , 2) , DS.IntermediateStates_PctUnfit(idx) );
         fprintf('%d\t%0.03f\t%0.03f\t%0.03f\n' , I , c1,c2 , c3);
    end
end
%% summary figures;
fh = figure('units','centimeters','position',[5 5 5 5 ]);
idx = DS.HD == 7 ; 
plot( DS.ParentStates_MeanFitPot(idx) , DS.IntermediateStates_PctUnfit(idx) ,'ok','MarkerFaceColor',[.7 .7 .7]);
title( sprintf('Pearson = %0.02f' , corr( DS.ParentStates_MeanFitPot(idx) , DS.IntermediateStates_PctUnfit(idx) )))
axis tight;
xlabel('Mean fitness potential')
ylabel('% unfit intermediates')
print('-dpsc2', [ num2str(SegN) '.eps']);
close;
%
fh = figure('units','centimeters','position',[5 5 5 5 ]);
G = grpstats( DS , 'HD' ,'mean' , 'DataVars' , 'IntermediateStates_PctUnfit');
bar(G.HD , G.mean_IntermediateStates_PctUnfit ,'FaceColor',[.7 .7 .7])
xlim([1.5 7.5])
ylabel('% unfit intermediate')
xlabel('Disance between pairs')
print('-dpsc2', [ num2str(SegN) '.eps'],'-append');
close;
