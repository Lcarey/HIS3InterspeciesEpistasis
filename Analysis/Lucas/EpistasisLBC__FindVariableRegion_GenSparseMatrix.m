function [ SparseVect , unique_variants]  = EpistasisLBC__FindVariableRegion_GenSparseMatrix(  aa_seq_variable )
%% [ SparseVect , unique_variants]  = EpistasisLBC__FindVariableRegion_GenSparseMatrix(  aa_seq_variable )
% find variable regions 
% generate SparseVect

% find the variable part
%  should run this independently!
%   [ aavar  , ~ , ~ ] = EpistasisLBC__ShortedAAseqsOnlyPositionsThatVary( T );
%   T.aavar = T.aavar 
%   [ SparseVect , unique_variants]  = EpistasisLBC__FindVariableRegion_GenSparseMatrix( aavar )
% generate sparse matrix


SparseVect = cell( length(aa_seq_variable) , 1);

seq_length = length(aa_seq_variable{1});
a = cellfun( @(X) arrayfun( @(I) [ num2str(I)  X(I)] , 1:seq_length ,'UniformOutput',false) ...
    , aa_seq_variable ,'UniformOutput',false);
a = vertcat(a{:}) ;
unique_variants = unique(a(:)) ; 
for I = 1:length(aa_seq_variable)
    SparseVect{I} = ismember(unique_variants , a(I,:))'  ;
end

end