%% build R
s = struct();
for SegN = 1:12
    % natlib
    T  = EpistasisLBC__LoadData( SegN  ,'ONLY_NATLIB_FLAG',true, 'ONLY_MIDDLE_FLAG',true,'NO_STOP_FLAG',true);
    
    %everything. basically the same result
    % T  = EpistasisLBC__LoadData( I  ,'ONLY_NATLIB_FLAG',false, 'ONLY_MIDDLE_FLAG',false,'NO_STOP_FLAG',true);
    
    T  = EpistasisLBC__FindVariableRegion_GenSparseMatrix( T ) ;
    
    fitness = T.s ;
    aaseqs = T.aa_seq  ;
    
    
    %% given this seq: YTSYVTV
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
    Seq1 = cell(1e5,1);
    Seq2 = cell(1e5,1);
    Fit1 = NaN(1e5,1);
    Fit2 = NaN(1e5,1);
    VarPos = NaN(1e5,1);
    Perm = cell(1e5,1);
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
                    Fit1(c) = fitness( idx_with_first(I) ) ;
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
    R = R( 1:(c-1) , :);
    
    R.FitImpact = R.Fit1 - R.Fit2 ;
    R.FitImpactThreshGood2 = (R.Fit1 - R.Fit2) >= 0.2   ;
    R.FitImpactThreshBad2  = (R.Fit1 - R.Fit2) <= -0.2   ;
    R.FitImpactThreshGood3 = (R.Fit1 - R.Fit2) >= 0.3   ;
    R.FitImpactThreshBad3  = (R.Fit1 - R.Fit2) <= -0.3   ;
    
    
end


%% Fit to model ;
load('~/Desktop/HIS3scratch/SignEpi/SignEpiWithKNNModel.mat');
%% find sign epi
for SegN = 1:12
    SegN
    [ G , ~]  = SignEpistasisFitModelToR( s(SegN).R , 0.2 , 50 ) ;
    G.Properties.ObsNames = [] ; 
    if SegN==1
        BigG = G;
    else
        BigG = vertcat(BigG,G);
    end
end

save('~/Desktop/HIS3scratch/SignEpi/SignEpiWithKNNModel_BigG.mat' , 'BigG' );

%% how common is sign epi? 
N=50;
BigG.HasEnoughBigEffects = BigG.sum_MinorSignFitEffect > 2*N | BigG.sum_MajorSignFitEffect > 2*N ;
BigG.HasSignEpiByCount   = BigG.sum_MinorSignFitEffect > N & BigG.sum_MajorSignFitEffect > N ;
BigG.HasSignEpiByAUC     = BigG.AUC_noweights > 0.6   ;
BigG.HasSignEpiByCountAndAUC = BigG.HasSignEpiByCount & BigG.HasSignEpiByAUC ; 

fprintf('All subs = %d\n' , length(BigG))
fprintf('enough w/big effectsd = %d (%0.02f%%)\n' , sum(BigG.HasEnoughBigEffects) ,  100*sum(BigG.HasEnoughBigEffects)/length(BigG))
fprintf('enough w/both signs = %d (%0.02f%%)\n' , sum(BigG.HasSignEpiByCount), 100*sum(BigG.HasSignEpiByCount)/length(BigG))
fprintf('pass AUC = %d (%0.02f%%)\n' , sum(BigG.HasSignEpiByAUC), 100*sum(BigG.HasSignEpiByAUC)/length(BigG))
fprintf('pass AUC & enough w/both signs = %d (%0.02f%%)\n' , sum(BigG.HasSignEpiByCountAndAUC), 100*sum(BigG.HasSignEpiByCountAndAUC)/length(BigG))

%%
fh = figure('units','centimeters','position',[5 5 7 7]);
hold on; 
histogram( 100*BigG.mean_MinorSignFitEffect( BigG.HasSignEpiByCountAndAUC  ) , 10 ,'Normalization','Count','EdgeColor','k','FaceColor','k')
xlabel('% of genontypes with minor sign')
ylabel('# of substitutions')
axis tight ;
xlim([0 004.4])
set(gca,'xtick',0:1:100)
print('-dpsc2','SignEpistasis.eps','-append')
close ;

%%
