function EpistasisLBC__tSNE_Visualize( Segment , n_dims , initial_dims , perplexity , theta , max_iter)
%% try different tSNE paramters
tic;
s = struct('Segment' , Segment , 'n_dims',n_dims,'datetime',datetime,'initial_dims',initial_dims,'perplexity',perplexity,'theta',theta,'max_iter',max_iter);


%% load table
T = EpistasisLBC__LoadData_FindVariableRegion_GenSparseMatrix( Segment );

%% add physical params
%%
try
	AAI = readtable('~/Google Drive/CareyLab/ExternalData/Atchley05/Atchley05.xlsx','Sheet','Table2');
catch
	AAI = readtable('~/single_cell_behavior/Data/Atchley05/Atchley05.xlsx','Sheet','Table2');
end
AAI_AAs = AAI.AminoAcid ;
AAI_mat = table2array( AAI(:,2:end));
%% How to score aa_seq_variable for each of the five Factors in AAI
len_seq = numel( T.aa_seq_variable{1}) ;
n_factors = size(AAI_mat,2) ;

T.physio_chem_vector = cellfun( @(X) ...
    cell2mat( arrayfun( @(SeqI) arrayfun( @(FactI)AAI_mat(ismember(AAI_AAs , X(SeqI)) ,FactI)   , 1:n_factors) , 1:len_seq,'UniformOutput',false))...
	    , T.aa_seq_variable ,'UniformOutput',false) ;

%% run tSNE
distance_metrics = {'euclidean' 'squaredeuclidean' 'seuclidean' 'cityblock' 'minkowski' 'chebychev' 'mahalanobis' 'cosine' 'correlation' 'spearman'};

s.M = cell2mat( T.physio_chem_vector  );
s.Y = fast_tsne(s.M, n_dims , initial_dims, perplexity , theta , 'svd', max_iter ); 
for I = 1:numel(distance_metrics)
	s.( ['Dist_M_' distance_metrics{I}] ) = squareform(pdist(s.M , distance_metrics{I})  ) ;
	s.( ['Y_' distance_metrics{I}] ) = fast_tsne( s.( ['Dist_M_' distance_metrics{I}] ) , n_dims , initial_dims, perplexity , theta , 'svd', max_iter ); 
end

output_fname = sprintf('tSNEVisualize_%d_%d_%d_%d_%0.02f_%d_%s.mat' , ...
	Segment , n_dims , initial_dims , perplexity , theta , max_iter , char(s.datetime) ) 

s.runtime = toc ;
save( output_fname ,'s');



% Y = fast_tsne(M, n_dims , initial_dims, 30, 0.5, 'svd', 1000); % pass it data



