function [ aa_seq_variable ,   columns_that_vary ,  uniq_aas ] = EpistasisLBC__ShortedAAseqsOnlyPositionsThatVary( T )
%% [ T  columns_that_vary uniq_aas ] = EpistasisLBC__ShortedAAseqsOnlyPositionsThatVary( T )
% 
% given  a table already read in
%   add in two additional columns showing which positions in the AA seq
%   vary. This used to be aa_seq_short
%
% the results of this will differ by how T is filtered, so if you want to
% filter (eg: nat_lib , no-stop, etc). Filter first!
%
% LBC

if ~istable(T)
    error( '[ T  columns_that_vary uniq_aas ] = EpistasisLBC__ShortedAAseqsOnlyPositionsThatVary( T )');
end

seqmat = cell2mat(T.aa_seq) ;
columns_that_vary = false( 1 , size(seqmat,2) ) ;
uniq_aas = cell( 1 , size(seqmat,2) ) ;

for I = 1:numel(columns_that_vary)
    uniq_aas{I} = unique( seqmat(:,I));
    if numel(uniq_aas{I})>1
        columns_that_vary(I)=true;
    end
end

aa_seq_variable = cellfun( @(X)X(columns_that_vary) , T.aa_seq  , 'UniformOutput' , false ) ; 

end