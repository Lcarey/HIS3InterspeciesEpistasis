function all_transition_states = ExpandSeqAlign( seq1 , seq2)
%% cell_array_all_transition_states = ExpandSeqAlign( seq1 , seq2)
% given two sequences of the same length return a cell array of all intermediate sequences
% LBC January 2017

mat = ff2n(length(seq1))+1 ;
seqs = [seq1 ; seq2] ;

all_transition_states = cell(nrows(mat),1);
l = length(seq1);
for I = 1:nrows(mat)
    all_transition_states{I} = buildseq( seqs , mat(I,:) , l) ;
end

all_transition_states = unique(all_transition_states);

end

function this_seq = buildseq( seqs , r , l)
    this_seq = '';
    for posI = 1:l
        this_seq = [ this_seq , seqs( r(posI) , posI)];
    end
end