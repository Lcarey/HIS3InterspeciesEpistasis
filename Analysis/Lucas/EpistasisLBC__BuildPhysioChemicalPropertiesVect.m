function T = EpistasisLBC__BuildPhysioChemicalPropertiesVect( T )
%% T = EpistasisLBC__BuildPhysioChemicalPropertiesVect( T )
% 
% given output T from  EpistasisLBC__LoadData_FindVariableRegion_GenSparseMatrix()
%   build a vector for each variant that has all physiochemical property scores 
%
%
% LBC March 207

%% load data
T = EpistasisLBC__LoadData_FindVariableRegion_GenSparseMatrix(7);
T = T(: ,  ~regexpcmp( T.Properties.VariableNames ,'^dist_') );
T = T(: ,  ~regexpcmp( T.Properties.VariableNames ,'^t.*_fr') );

%%
AAI = readtable('~/Google Drive/CareyLab/ExternalData/Atchley05/Atchley05.xlsx','Sheet','Table2')
AAI_AAs = AAI.AminoAcid ;
AAI_mat = table2array( AAI(:,2:end));
%% How to score aa_seq_variable for each of the five Factors in AAI
len_seq = numel( T.aa_seq_variable{1}) ;
n_factors = size(AAI_mat,2) ;

T.physio_chem_vector = cellfun( @(X) ...
    cell2mat( arrayfun( @(SeqI) arrayfun( @(FactI)AAI_mat(ismember(AAI_AAs , X(SeqI)) ,FactI)   , 1:n_factors) , 1:len_seq,'UniformOutput',false))...
    , T.aa_seq_variable ,'UniformOutput',false) ;


%%

distmat = squareform(pdist( cell2mat(T.physio_chem_vector) , 'seuclidean'))  ;
fitmat = squareform(pdist(T.s , 'euclidean'))  ;
corr(fitmat(:),distmat(:))
%%
WT = 'IHALAKHSGWSLIVECIGDLHIDDHHTTED'
WTvar = WT(T.columns_that_vary(1,:))
WT_physio_chem_vector = cell2mat( arrayfun( @(SeqI) arrayfun( @(FactI)AAI_mat(ismember(AAI_AAs , WTvar(SeqI)) ,FactI)   , 1:n_factors) , 1:len_seq,'UniformOutput',false)) ;

distmat = pdist2( WT_physio_chem_vector , cell2mat(T.physio_chem_vector) )  ;
fitmat = pdist2(T.s(strcmp(T.aa_seq_variable,WTvar)) , T.s  )  ;
corr(fitmat(:),distmat(:))
% % 
% % %% PCA
% % [coeff,score,latent,tsquared,explained,mu] = pca(cell2mat(T.SparseVect) ,'NumComponents',10);
% % %[coeff,score,latent,tsquared,explained,mu] = pca(cell2mat(T.physio_chem_vector) ,'NumComponents',10);
% % 
% % 
% % figure; plot(cumsum(explained),'-ok');xlim([0 25]); grid on; set(gca,'xtick',1:2:100)
% % ylabel('% variance explained')
% % xlabel('# of principal coponents')
% % title('binary sparse vector')
% % 
% % figure; scatter3(score(:,1),score(:,2),score(:,3),30,T.s,'filled');
% % colorbar;
% % xlabel('PC1');ylabel('PC2');zlabel('PC3')
% % title('binary sparse vector')
