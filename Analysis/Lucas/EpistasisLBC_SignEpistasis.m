%% build R
s = struct();
for SegN = 1:12
    tic;
    % natlib
    T  = EpistasisLBC__LoadData( SegN  ,'ONLY_NATLIB_FLAG',true, 'ONLY_MIDDLE_FLAG',true,'NO_STOP_FLAG',true);
    
    %everything. basically the same result
    % T  = EpistasisLBC__LoadData( I  ,'ONLY_NATLIB_FLAG',false, 'ONLY_MIDDLE_FLAG',false,'NO_STOP_FLAG',true);
    
  %  T  = EpistasisLBC__FindVariableRegion_GenSparseMatrix( T ) ;
    
    fitness = T.s ;
    aaseqs = T.aa_seq  ;
    
    
    % given this seq: YTSYVTV
    % what is the effect of introducing the Y into position 1
    %  so, given all seqs that have a Y in position 1
    %   find all seqs that are 0 edits away i
    
    % given V at pos 7   &  A at pos 7
    % eg: YTSYVTV & YTSYVTA
    %  this is one pair. calc fit diff
    %   for all other pairs that go V->A at pos 7
    %    calc fit diff
    
    % mat = double(cell2mat(T.SparseVect)) ;
    % distance_mat = squareform( pdist(mat,'hamming') ) * size(mat,2)  ;
    
    %
    Seq1 = cell(1e6,1);
    Seq2 = cell(1e6,1);
    Fit1 = NaN(1e6,1);
    Fit2 = NaN(1e6,1);
    VarPos = NaN(1e6,1);
    Perm = cell(1e6,1);
    c  = 0 ;
    for this_variant_position = 1:length(aaseqs{1})
        unique_aas_at_this_position = unique(cellfun( @(X)X(this_variant_position) , aaseqs)) ;
        permutations = nchoosek(unique_aas_at_this_position,2) ;
        for pI = 1:size(permutations,1)
            idx_with_first  = find( cellfun( @(X)X(this_variant_position) , aaseqs) == permutations(pI,1)) ;
            for I = 1:numel(idx_with_first)
                this_seq = aaseqs{idx_with_first(I)} ;
                next_seq = this_seq ;
                next_seq(this_variant_position) = permutations(pI,2) ;
                idx_of_pair = strcmp( aaseqs , next_seq)  ;
                if sum(idx_of_pair)==1
                    c = c + 1 ;
                    Seq1{c} = aaseqs{ idx_with_first(I) };
                    Seq2{c} = aaseqs{ idx_of_pair } ;
                    Fit1(c) = fitness1( idx_with_first(I) ) ;
                    Fit2(c) = fitness( idx_of_pair );
                    VarPos(c) = this_variant_position ;
                    Perm{c} = permutations(pI,:) ;
                end
            end
        end
    end
    
    %
    R = dataset();
    R.Seq1 = Seq1 ;
    R.Seq2 = Seq2 ;
    R.Fit1 = Fit1 ;
    R.Fit2 = Fit2 ;
    R.VarPos = VarPos ;
    R.Perm = Perm ;    
    R.FitImpact = R.Fit1 - R.Fit2 ;
    R = R( ~isnan(R.VarPos) ,:); % in case we pre-allocated too much space
    fprintf('Seg %d took %0.01f min\n' , SegN , toc/60) ;
	s(SegN).R = R ; 
	save( [ '~/Desktop/HIS3scratch/SignEpi/SignEpi_R' num2str(SegN) '.mat' ] , 'R' );
    save('~/Desktop/HIS3scratch/SignEpi/SignEpi_Ronly.mat' , 's');
end

% %% remememer to remove excess pre-allocated space
% for I = 1:12
%     s(I).R =     s(I).R( ~isnan(s(I).R.VarPos) ,:);  
%     length(s(I).R)
% end
%% Fit to model ;
% load('~/Desktop/HIS3scratch/SignEpi/SignEpi_Ronly.mat');

% find sign epi
for SegN = 1:12
    SegN
    [ G , ~]  = SignEpistasisFitModelToR( s(SegN).R , 0.4  ) ;
    G.Properties.ObsNames = [] ; 
    G.SegN = repmat( SegN , length(G) , 1);
    if SegN==1
        BigG = G;
    else
        BigG = vertcat(BigG,G);
    end
end

save('~/Desktop/HIS3scratch/SignEpi/SignEpi_BigG_04.mat' , 'BigG' );
