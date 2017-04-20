function s = EpistasisLBC__CompareUnfitFreqBetweenFit_VS_EntireLibrary(   SegRange , N_Pairs_Fast_to_measure , save_struct_file_name )
%%  s = EpistasisLBC__CompareUnfitFreqBetweenFit_VS_EntireLibrary( SegRange  , N_Pairs_Fast_to_measure , save_struct_file_name )
%
%
% LBC 2017


if ~exist('fast_fit_cutoff','var')
    fast_fit_cutoff = 0.4 ;
end
if ~exist('slow_fit_cutoff','var')
    slow_fit_cutoff = 0.2  ;
end

if ~exist('SegRange','var')
    SegRange = 1:12 ; 
end

if ~exist('N_Pairs_Fast_to_measure','var')
    N_Pairs_Fast_to_measure = 1000 ; 
end


MapAA2I = containers.Map(  arrayfun(@(X){X},['A':'Z' '_' ])  , uint8(1:27) ) ; % all AAs + stop

s = struct();

DataDir = '~/Develop/HIS3InterspeciesEpistasis/Data/';


% to compute dist of each seq from wild-type
WTseqs = readtable( [ DataDir '../Data_Small_Tables/wt_seq.csv'] ,'Format','%s%s','ReadVariableNames',true,'HeaderLines',0);
WTseqs.Var1 = str2double( regexprep(WTseqs.Var1 , '^S',''));
WTseqs.Properties.VariableNames = {'SegN' 'wt_aa_seq'};

for SegN = SegRange
    thisone = length(s)+1; if(length(s)==1) , thisone=1; end
    s(thisone).SegN = SegN ;
    T = readtable([ DataDir filesep num2str(SegN) '.tab' ],'FileType','text') ;
    
    T = T( ~logical(T.stop) & ~logical(T.nonsense)  ,:);
    T = T( logical(T.middle) ,:);
    T = T( logical(T.lib) ,:);
    T = T( T.len == mode(T.len)  , :);
    T = T( logical(T.nat_lib) , :);
    
    wtseq = WTseqs.wt_aa_seq{ WTseqs.SegN == SegN } ;
    T.DistFromScer = cellfun( @(X) HammingDistance(X,wtseq) , T.aa_seq  );
    
    
    fast_seqs =  find( T.s >=  fast_fit_cutoff ) ;
    
    %% find all columns that vary
    aa_num = cellfun( @(A) arrayfun( @(X)MapAA2I(X),A) , T.aa_seq ,'UniformOutput' , false );
    aa_num_mat = cell2mat(aa_num) ;
    n_aa_per_col = NaN( 27 , length(T.aa_seq{1}) ) ;
    for I = 1:size(n_aa_per_col,2)
        [a,b]=count_unique(aa_num_mat(:,I));
        n_aa_per_col(a,I) = b ;
    end
    cols_with_variation = arrayfun( @(I)sum(n_aa_per_col(:,I)>1) , 1:size(n_aa_per_col,2)) > 1  ;
    T.aa_seq_varies = cellfun( @(X) X(cols_with_variation) , T.aa_seq,'UniformOutput',false);
    
    
    % choose N random pairs of fast sequences
    pairs = [ randsample(fast_seqs,N_Pairs_Fast_to_measure*10,true) randsample(fast_seqs,N_Pairs_Fast_to_measure*10,true)];
    pairs = pairs( pairs(:,1)~=pairs(:,2) ,:);
    pairs = pairs( 1:N_Pairs_Fast_to_measure , :);
    N_Pairs_Fast_to_measure = size(pairs,1);
    
    R = dataset();
    R.FitnessDistributions = cell( N_Pairs_Fast_to_measure , 1);
    R.HammingDistances = NaN( N_Pairs_Fast_to_measure , 1);
    R.NUnfitMeasured = NaN( N_Pairs_Fast_to_measure , 1);
    R.NFitMeasured = NaN( N_Pairs_Fast_to_measure , 1);
    R.NIntMeasured = NaN( N_Pairs_Fast_to_measure , 1);
    R.StatesMeasured = NaN( N_Pairs_Fast_to_measure , 1);
    R.Seq1idx = NaN( N_Pairs_Fast_to_measure , 1);
    R.Seq2idx = NaN( N_Pairs_Fast_to_measure , 1);
    
    R.Seq1 = cell( N_Pairs_Fast_to_measure , 1);
    R.Seq2 = cell( N_Pairs_Fast_to_measure , 1);
    R.ShortSeq1 = cell( N_Pairs_Fast_to_measure , 1);
    R.ShortSeq2 = cell( N_Pairs_Fast_to_measure , 1);
    R.IntermediateStatesIDX = cell( N_Pairs_Fast_to_measure , 1);
    
    tic;
    for I = 1:N_Pairs_Fast_to_measure
        seq1 = T.aa_seq_varies{pairs(I,1)}; R.Seq1idx(I) = pairs(I,1) ;
        seq2 = T.aa_seq_varies{pairs(I,2)}; R.Seq2idx(I) = pairs(I,2) ;
        R.Seq1{I} = T.aa_seq{pairs(I,1)} ;
        R.Seq2{I} = T.aa_seq{pairs(I,2)} ;
        R.ShortSeq1{I} = seq1 ;
        R.ShortSeq2{I} = seq2 ;
        all_transition_states = ExpandSeqAlign( seq1 , seq2);
        idx = find( ismember( T.aa_seq_varies , all_transition_states) );
        R.FitnessDistributions{I} = T.s(idx) ;
        R.IntermediateStatesIDX{I} = idx ;
        
        R.NUnfitMeasured(I) = sum( R.FitnessDistributions{I} <= slow_fit_cutoff);
        R.NFitMeasured(I) = sum( R.FitnessDistributions{I} >= fast_fit_cutoff);
        R.NIntMeasured(I) = sum( R.FitnessDistributions{I} < fast_fit_cutoff &  R.FitnessDistributions{I} > slow_fit_cutoff);
        
        R.HammingDistances(I) = HammingDistance( seq1 , seq2 ) ;
        R.StatesMeasured(I) = numel(idx)  ;
        if (mod(I,20)==0) , fprintf('.') , end
        if (mod(I,100)==0) , fprintf(' ') , end
        if (mod(I,1000)==0) , fprintf(' %0.0f%%\n',I/N_Pairs_Fast_to_measure*100) , end
    end
    toc
    R.pct_unfit  = ( R.NUnfitMeasured ) ./ R.StatesMeasured * 100 ;
    R.pct_states_measured = R.StatesMeasured ./ (2.^R.HammingDistances) .* 100 ;
    
    s(thisone).T = T ;
    s(thisone).R = R ;
end

% optionally save the structure
if exist('save_struct_file_name','var')
    save(save_struct_file_name , 's');
end

end
