%% Fig_NNModelWeightsCorrelateEvolDistance
% do the weights asigned to each AA state (eg: segment 7 , pos 0 , F)
% correlate w/the evoultionary distance to when that state was first observed? 
%
%
ALN = readtable('/Users/lcarey/Desktop/HIS3scratch/EvolDistance/HIS3_alignment_500species.fas.tab','FileType','text','Delimiter','\t','ReadVariableNames',false);
DISTS = readtable('/Users/lcarey/Desktop/HIS3scratch/EvolDistance/HIS3_alignment_500species.fas.tab.pairwise.proportion','FileType','text','Delimiter',',','ReadVariableNames',false);
DISTS.Properties.VariableNames = {'Sp1' 'Sp2' 'Dist'} ;
scer_aligned_seq = ALN.Var2{ regexpcmp( ALN.Var1 , 'cerevisiae')};
PosMapVect = MapPositionsFromAlignmentWithInsertionsBackIntoRefSeqCoords( scer_aligned_seq ) ;

%%
idx1 = find(regexpcmp(DISTS.Sp1,'cerevi'));
idx2 = find(regexpcmp(DISTS.Sp2,'cerevi'));

D = dataset();
D.Species = vertcat( DISTS.Sp2(idx1) , DISTS.Sp1(idx2));
D.Distances = vertcat( DISTS.Dist(idx1) , DISTS.Dist(idx2));
