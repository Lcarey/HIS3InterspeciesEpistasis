%%
% nat - 1 if all mutations in genotype appear in some species. Based on uniprot alignment with ~500 species. Is a subject to slightly change because Fedya doesn't like alignment
% lib - 1 if all mutations in genotype were in library design
% nat_lib nat&lib
% nogap - 1 if no gaps or stop codons
% shift - 1 if there is frame shift 
% stop - 1 if there is stop codon
% middle - 1 if there are no mutations in middle (non-mutated) part
% size - how many nucleotide genotypes have this amino acid genotype


%% get all data files
WD = '~/Develop/HIS3InterspeciesEpistasis/';
% S12 = readtable([ WD '/Data/' 'S12_exponent_s_gridfit_bybest10000_aa_avg_with_dist.csv' ]);
% S12.fitness = S12.s ;
library_data_files  = dir( [  WD '/Data/' '*.csv']);

%% For each library, find AA seqs for each species
% plot fitness distribution for all libraries
figure; hold on ;

for I = 1:numel(library_data_files)
    T = readtable(  [  library_data_files(I).folder filesep library_data_files(I).name ]);
    %T = T( T.middle & T.nogap & T.nat_lib , :);
    
    % find WT sequences and the indexes into the table for each seq's WT
    wt_species_names = T.Properties.VariableNames ; 
    wt_species_names = wt_species_names( regexpcmp(wt_species_names , '^dist_[A-Z][a-z]'));
    wt_sequences_idx = NaN(1,numel(wt_species_names));
    for J = 1:numel(wt_species_names)
        idx = find( T.(wt_species_names{J}) == 0  & T.nat_lib & T.middle & T.nogap );
        if ~isempty(idx)
            wt_sequences_idx(J) = idx ;
        end
    end
    % remove species for which we don't have a sequence
    wt_species_names = wt_species_names(~isnan(wt_sequences_idx));
    wt_sequences_idx = wt_sequences_idx(~isnan(wt_sequences_idx));

    wt_specices_names = regexprep( wt_species_names , 'dist_' ,'');
    
    subplot(3,4,I) ; hold on; 
    [f,x] = ecdf( T.s);  plot(x,f,'-k','LineWidth',3);  % ALL
    [f,x] = ecdf( T.s(logical(T.stop)));  plot(x,f,'-r','LineWidth',3);  % stop
    [f,x] = ecdf( T.s( T.middle & T.nogap & ~logical(T.stop) ));   plot(x,f,'-b','LineWidth',3);  % natural
    [f,x] = ecdf( T.s( T.middle & T.nogap & ~logical(T.stop) & T.nat_lib ));  plot(x,f,'-c','LineWidth',3);  % natural

    for J = 1:numel(wt_sequences_idx)
        line( [ T.s( wt_sequences_idx(J) )  T.s( wt_sequences_idx(J) ) ],ylim);
        text( T.s( wt_sequences_idx(J) ) , random('uniform',.1,.9,1) , wt_specices_names{J});
    end
    xlim([-2 1.5])
    xlabel('Fitness (s)')
    grid on ;
    t = regexp( library_data_files(I).name , '_' ,'split');
    title( t{1} );
    

    library_data_files(I).wt_sequences_idx = wt_sequences_idx ;
    library_data_files(I).wt_sequences_distances = pdist( string(T.aa_seq(wt_sequences_idx)) , @HammingDistance) ;
    library_data_files(I).wt_sequences_short_distances = pdist( string(T.aa_seq_short(wt_sequences_idx)) , @HammingDistance) ;

    wt_sequences_idx = wt_sequences_idx(T.s(wt_sequences_idx)>0) ;
    library_data_files(I).wt_sequences_idx_fit = wt_sequences_idx ;
    library_data_files(I).wt_sequences_distances_fit = pdist( string(T.aa_seq(wt_sequences_idx)) , @HammingDistance) ;
    library_data_files(I).wt_sequences_short_distances_fit = pdist( string(T.aa_seq_short(wt_sequences_idx)) , @HammingDistance) ;

end