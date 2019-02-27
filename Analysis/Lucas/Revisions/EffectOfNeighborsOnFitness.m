%% load synonymous variant data
T = readtable( '~/Develop/HIS3InterspeciesEpistasis/Data/synonymous_variants_rescaled_data.tab' , 'FileType','text','Delimiter','\t');

%% get data for the segment of interest
results_cell_vect = cell( 12 , 2) ; 
for SegN = 1:12
    
    NT = T( T.SegN == SegN , :);
    NT.fitness_nt = NT.s ; NT.s = [] ; 
    
    % rescale fitness
    AA = readtable( [ '~/Develop/HIS3InterspeciesEpistasis/Data/S' num2str(SegN) '_scaled_info_v2.csv' ]  , 'FileType','text','Delimiter','\t');
    Q = innerjoin(AA(  : ,{'aa_seq' 'size' 's'}) , NT(:,{'aa_seq','fitness_nt' 'seq'}) , 'Key','aa_seq');
    near_wt_aa_3 = AA.aa_seq(AA.dist_Scer<=3 & AA.nat_lib) ; % fitness for near-WT strains
    mf_near_wt_fitness = modefit(Q.fitness_nt( ismember(Q.aa_seq , near_wt_aa_3) )) ; 
    ns_95 = prctile(NT.fitness_nt(regexpcmp(NT.aa_seq,'_')) , 95 ) ; 
    
    NT.fitness_nt_rescaled = NT.fitness_nt - ns_95 ; 
    NT.fitness_nt_rescaled( NT.fitness_nt_rescaled < 0 ) = 0 ;
    NT.fitness_nt_rescaled = NT.fitness_nt_rescaled ./ (mf_near_wt_fitness - ns_95) ; 
   
    
    HIGH_FITNESS_THRESHOLD = 0.6 ;
    basename =  [ '~/Downloads/nonsense_NT_genotypes_neighbor_high_fitness_effect_on_fitness__Segment_' num2str(SegN)  ]

    NT.nt_seq_lengths = cellfun(@length , NT.seq) ;
    keep_idx =  NT.nt_seq_lengths == mode( NT.nt_seq_lengths )  ;
    nonsense_idx = regexpcmp(NT.aa_seq,'_') ; 
    NS = NT( nonsense_idx & keep_idx , :);
    HF = NT( NT.fitness_nt_rescaled > HIGH_FITNESS_THRESHOLD & ~nonsense_idx & keep_idx , :);
    
    % calculate the # of neighbors that are the closest distance
    tic ;
    nearest_neighbor_count_vect_1 = NaN( height(NS) , 1) ;
    nearest_neighbor_distance_vect_1 = NaN( height(NS) , 1) ;
    nearest_neighbor_count_vect_10 = NaN( height(NS) , 1) ;
    nearest_neighbor_distance_vect_10 = NaN( height(NS) , 1) ;
    ns_nt_seqs = NS.seq ;
    hf_nt_seqs = HF.seq ;
    tic;
    parfor I = 1:numel(ns_nt_seqs)
        ntseq = ns_nt_seqs{I} ;
        dv = cellfun( @(X)HammingDistance( X , ntseq  ) , hf_nt_seqs );
        
        % nearest_neighbor count == 1? or count == 10 ? 
        [all_distances,all_distances_counts]=count_unique(dv) ;
        nearest_neighbor_distance_vect_10(I) =  min(all_distances(all_distances_counts >= 10))  ; % min count == 10 
        nearest_neighbor_distance_vect_1(I) = min(dv) ; % min count == 1
 
        nearest_neighbor_count_vect_10(I) = sum( dv == nearest_neighbor_distance_vect_10(I) ) ;
        nearest_neighbor_count_vect_1(I) = sum( dv == nearest_neighbor_distance_vect_1(I) ) ;

    end
    toc
    
    % calculate summary statistics
    R = table();
    R.fitness = NS.fitness_nt_rescaled ;
    R.nearest_neighbor_distance_1 = nearest_neighbor_distance_vect_1 ;
    R.nearest_neighbor_count_1 = nearest_neighbor_count_vect_1 ;
    R.nearest_neighbor_distance_10 = nearest_neighbor_distance_vect_10 ;
    R.nearest_neighbor_count_10 = nearest_neighbor_count_vect_10 ;

    results_cell_vect{SegN,1} = R ; % save a copy of R

    R.nearest_neighbor_distance_1( R.nearest_neighbor_distance_1>5) = 7 ;
    R.nearest_neighbor_distance_10( R.nearest_neighbor_distance_10>5) = 7 ;


    % what fraction taken up by each
    R.nearest_neighbor_distance = R.nearest_neighbor_distance_1 ; 
    G = grpstats(R,'nearest_neighbor_distance' , 'mean' , 'DataVars' , 'fitness' );
    G = G( G.GroupCount >= 3 , :) ;
    G.pct  = G.GroupCount ./ sum(G.GroupCount) * 100 ;
    for I = 1:height(G)
        m = bootstrp( 1000 , @mean , R.fitness( R.nearest_neighbor_distance == G.nearest_neighbor_distance(I) ));
        G.mean_fitness(I) =  mean(m);
        G.std_fitness(I) =  std(m) ;
    end
    
    results_cell_vect{SegN,2} = G ; % save a copy of G
    delete( '~/Downloads/results_cell_vect.mat'  ) ;
    save('~/Downloads/results_cell_vect.mat' , 'results_cell_vect') ; 
    
    
    save([ basename '.mat'] , 'R' , 'G' , 'SegN' , 'basename' , 'AA' , 'NT' , 'results_cell_vect' ) ; 

end

%% plot all segments into a single figure
Q = table();
for SegN = 1:12
    basename =  [ '~/Downloads/nonsense_NT_genotypes_neighbor_high_fitness_effect_on_fitness__Segment_' num2str(SegN)  ] ;
    load([ basename '.mat'],'R') ;  
    Q = vertcat(Q,R);
end
R = Q ; 
    G = grpstats(R,'nearest_neighbor_distance' , 'mean' , 'DataVars' , 'fitness' );
    G = G( G.GroupCount >= 3 , :) ;
    G.pct  = G.GroupCount ./ sum(G.GroupCount) * 100 ;
    for I = 1:height(G)
        m = bootstrp( 1000 , @mean , R.fitness( R.nearest_neighbor_distance == G.nearest_neighbor_distance(I) ));
        G.mean_fitness(I) =  mean(m);
        G.std_fitness(I) =  std(m) ;
    end

  
%% figure
    figname =   '~/Downloads/nonsense_NT_genotypes_neighbor_high_fitness_effect_on_fitness_3.png' ;
    G.nearest_neighbor_distance(end) = 6.5  ; 
    xl = arrayfun(@(I){num2str(I)},G.nearest_neighbor_distance) ;
    xl{end} = '>5' ;
 
    figure('Position',[0 0 250 300]) ;
    
    subplot(2,1,1)
    hold on ;
    bar( G.nearest_neighbor_distance , G.pct , 'FaceColor',[.8 .8 .8]);
    set(gca,'xtick',G.nearest_neighbor_distance )  ;
    set(gca,'xticklabel',xl)
    ylabel('% of genotypes')
    xlabel('NT to closest high fitness genotype')
    ylim([0 50])   
    set(gca,'ytick',0:25:100)
    
    subplot(2,1,2)
    hold on ;
    bar( G.nearest_neighbor_distance , G.mean_fitness );
    errorbar( G.nearest_neighbor_distance , G.mean_fitness , G.std_fitness , 'ok' , 'LineWidth' , 2) ; 
    set(gca,'xtick',G.nearest_neighbor_distance )  ;
    set(gca,'xticklabel',xl)
    ylabel('Mean fitness')
    xlabel('NT to closest high fitness genotype')
   print('-dpng',figname , '-r600');
    close ;


    
%% replot figures
for SegN = 1:12
    basename =  [ '~/Downloads/nonsense_NT_genotypes_neighbor_high_fitness_effect_on_fitness__Segment_' num2str(SegN)  ] ;
    load([ basename '.mat'],'G') ;  
    
    % figure
    xl = arrayfun(@(I){num2str(I)},G.nearest_neighbor_distance) ;
    xl{end} = '>5' ;
 
    figure('Position',[0 0 250 300]) ;
    
    subplot(2,1,1)
    hold on ;
    bar( G.nearest_neighbor_distance , G.pct , 'FaceColor',[.8 .8 .8]);
    set(gca,'xtick',G.nearest_neighbor_distance )  ;
    set(gca,'xticklabel',xl)
    ylabel('% of genotypes')
    xlabel('Edit distance to closest high fitness genotype')
    title(['Segment ' num2str(SegN) ]);
    ylim([0 100])   
    set(gca,'ytick',0:25:100)
    
    subplot(2,1,2)
    hold on ;
    bar( G.nearest_neighbor_distance , G.mean_fitness );
    errorbar( G.nearest_neighbor_distance , G.mean_fitness , G.std_fitness , 'ok' , 'LineWidth' , 2) ; 
    set(gca,'xtick',G.nearest_neighbor_distance )  ;
    set(gca,'xticklabel',xl)
    ylabel('Mean fitness')
    xlabel('Edit distance to closest high fitness genotype')
  %  ylim([0 0.2]) ;     set(gca,'ytick',0:0.1:1);

    print('-dpng',[ basename '.png'] , '-r600');
    close ;
end




%% 
% % % %% figure;
% % % %     NS = NS(1:I,:);
% % % w = cell(4,1) ;
% % % for wI = 1:numel(w)
% % %     idx =  NS.(genvarname(['Nd' num2str(wI)]))  > 1 & NS.(genvarname(['Nd' num2str(wI-1)])) == 0  ;
% % %     w{wI} = NS.s( idx ) ;
% % % end
% % %
% % % figure; hold on ;
% % % X =  [ 1:numel(w) numel(w)  + 1.25 ] ;
% % % bar(  X , vertcat( cellfun(@nanmean,w) , nanmean(NS.s( NS.Nlt5 == 0 ))  ) )
% % % for wI = 1:numel(w)
% % %     m = bootstrp(1000, @mean, w{wI});
% % %     errorbar( X(wI) , mean(m) , std(m),'o-k','LineWidth',3);
% % % end
% % % m = bootstrp(1000, @mean,  NS.s( NS.Nlt5 == 0 ) ) ;
% % % errorbar( X(end) , mean(m) , std(m),'o-k','LineWidth',3);
% % % set(gca,'xtick', X )
% % % set(gca,'xticklabel' , {'1' '2' '3' '4' '>5'})
% % % ylabel('Mean fitness of nonsense genotypes')
% % % xlabel('Edit distance to closest high fitness genotype')
% % % xlim([0.5 wI+1.75])
% % % grid on ;
% % %
% % %
% % % %%
% % % THRESH = 0.2 ;
% % % figure; hold on ;
% % % X =  [ 1:numel(w) numel(w)  + 1.25 ] ;
% % % bar(  X , vertcat( cellfun(@(X)mean(X>THRESH),w) , mean(NS.s( NS.Nlt5 == 0 ) > THRESH) ) * 100 ) ;
% % %
% % % for wI = 1:numel(w)
% % %     m = bootstrp(1000, @(X)mean(X>THRESH), w{wI}) * 100 ;
% % %     errorbar( X(wI) , mean(m) , std(m),'o-k','LineWidth',3);
% % % end
% % % m = bootstrp(1000, @(X)mean(X>THRESH)   , NS.s( NS.Nlt5 == 0 )) * 100 ;
% % % errorbar( X(end) , mean(m) , std(m),'o-k','LineWidth',3);
% % %
% % % set(gca,'xtick', X )
% % % set(gca,'xticklabel' , {'1' '2' '3' '4' '>5'})
% % % ylabel('% of genotypes with fitness > 0.2')
% % % xlabel('Edit distance to closest high fitness genotype')
% % % xlim([0.5 wI+1.75])
% % % grid on ;
