function [ T , unique_variants]  = EpistasisLBC__FindVariableRegion_GenSparseMatrix( T )
%% [ T , unique_variants]  = EpistasisLBC__FindVariableRegion_GenSparseMatrix( T )
% find variable regions 
% generate SparseVect

% find the variable part
[ T  , ~ , ~ ] = EpistasisLBC__ShortedAAseqsOnlyPositionsThatVary( T );

% generate sparse matrix


SparseVect = cell( height(T) , 1);

seq_length = length(T.aa_seq_variable{1});
a = cellfun( @(X) arrayfun( @(I) [ num2str(I)  X(I)] , 1:seq_length ,'UniformOutput',false) ...
    , T.aa_seq_variable ,'UniformOutput',false);
a = vertcat(a{:}) ;
unique_variants = unique(a(:)) ; 
for I = 1:height(T)
    SparseVect{I} = ismember(unique_variants , a(I,:))'  ;
end

T.SparseVect = SparseVect ; 