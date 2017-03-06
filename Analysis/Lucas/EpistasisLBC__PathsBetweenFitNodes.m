%%
% nat - 1 if all mutations in genotype appear in some species. Based on uniprot alignment with ~500 species. Is a subject to slightly change because Fedya doesn't like alignment
% lib - 1 if all mutations in genotype were in library design
% nat_lib nat&lib
% nogap - 1 if no gaps or stop codons
% shift - 1 if there is frame shift 
% stop - 1 if there is stop codon
% middle - 1 if there are no mutations in middle (non-mutated) part
% size - how many nucleotide genotypes have this amino acid genotype


%% load data
WD = '~/Develop/HIS3InterspeciesEpistasis/';
% S12 = readtable([ WD '/Data/' 'S12_exponent_s_gridfit_bybest10000_aa_avg_with_dist.csv' ]);
% S12.fitness = S12.s ;

T = readtable([ WD '/Data/' 'S3_exponent_s_gridfit_bybest10000_aa_avg_with_dist.csv' ]);
T.fitness = T.s ;

%% Determine fitness threshold using nonsense mutants
Y = T.fitness ;
Y(Y<-2) = -2 ;
Y(Y > 1) = 1 ;

xl = -2.1:.05:1.1 ;

fh = figure('units','centimeters','position',[5 5 10 10]);

subplot(2,1,2)
hold on; 
ecdf( Y( logical(T.stop)) );
ecdf( Y( ~logical(T.stop)));
legend({'+stop' 'no stop'},'location','nw')
xlim([-2 1])
grid on ;
set(gca,'xtick',-2:.5:1)
set(gca,'ytick',0:.1:1)


subplot(2,1,1)
hold on; 
histogram( Y( logical(T.stop)) , xl ,'Normalization','Probability');
histogram( Y( ~logical(T.stop)) , xl,'Normalization','Probability');
legend({'+stop' 'no stop'},'location','nw')
axis tight;

for I = [95 97.5 99 99.5 99.9 99.99]
    v = prctile( Y( logical(T.stop)) , I);
    fprintf('%0.02f%% of nonsense left w/a fitness threshold of %0.02f\n' , 100-I , v);
end
    
% Declare fitness threhsold 
fit_cutoff = -0.25 ; % less than this is unfit

%% remove all nonsense
T = T( ~T.stop , :);
T = T( T.len == mode(T.len) , :);

%% find WT sequences and the indexes into the table for each seq's WT
wt_species_names = T.Properties.VariableNames ; 
wt_species_names = wt_species_names( regexpcmp(wt_species_names , '^dist_[A-Z][a-z]'));
wt_sequences_idx = NaN(1,numel(wt_species_names));
for I = 1:numel(wt_species_names)
    idx = find( T.(wt_species_names{I}) == 0  & T.nat_lib );
    if ~isempty(idx)
        wt_sequences_idx(I) = idx ;
    end
end
% remove species for which we don't have a sequence
wt_species_names = wt_species_names(~isnan(wt_sequences_idx));
wt_sequences_idx = wt_sequences_idx(~isnan(wt_sequences_idx));

wt_specices_names = regexprep( wt_species_names , 'dist_' ,'')


%% Plot fitness of WT genotypes
figure ; hold on ;
[f,x] = ecdf( T.fitness);
plot(x,f,'-k','LineWidth',3);
for I = 1:numel(wt_sequences_idx)
    line( [ T.fitness( wt_sequences_idx(I) )  T.fitness( wt_sequences_idx(I) ) ],ylim);
    text( T.fitness( wt_sequences_idx(I) ) , random('uniform',0,1,1) , wt_specices_names{I});
end


%% remove any unfit WT sequences

keep = T.fitness( wt_sequences_idx) > fit_cutoff ;
wt_sequences_idx = wt_sequences_idx(keep);
wt_specices_names = wt_specices_names(keep) ;

WTseqsMap = containers.Map( wt_specices_names , wt_sequences_idx);



%%
fast_seqs =  T.aa_seq_short( wt_sequences_idx ) ;
fast_seqs =  T.aa_seq_short( T.fitness >  fit_cutoff ) ;
Scer_seq = T.aa_seq_short( wt_sequences_idx(strcmp( wt_specices_names , 'Scer'))) ; 
J = find( ismember( fast_seqs , Scer_seq) , 1 );

IntermediateIDX = cell( numel(fast_seqs) , 1);
ParentIDX = cell( numel(fast_seqs) , 1);

FitnessDistributions = cell( numel(fast_seqs) , 1);
HammingDistances = NaN( numel(fast_seqs) , 1);
FractionIntermediateUnfit = NaN( numel(fast_seqs) , 1);
NumIntermediateStatesMeasured = NaN( numel(fast_seqs) , 1);
for I = 1:numel(fast_seqs)
    seq1 = fast_seqs{J};
    seq2 = fast_seqs{I};
    all_transition_states = ExpandSeqAlign( seq1 , seq2);
    idx = find( ismember( T.aa_seq_short , all_transition_states) );
    FractionIntermediateUnfit(I) = mean( T.fitness(idx)<fit_cutoff);
    HammingDistances(I) = HammingDistance( seq1 , seq2 ) ;
    FitnessDistributions{I} = T.fitness(idx) ;
    NumIntermediateStatesMeasured(I) = numel(idx) - 2  ;
    IntermediateIDX{I} = idx ; 
%     if numel(idx)>5
%         fh = figure('units','centimeters','position',[5 5 5 5]);
%         hold on ;
%         ecdf(T.fitness)
%         ecdf(T.fitness(idx))
%         legend({'all' sprintf('%d %0.01f%%',numel(idx),fraction_intermediate_unfit*100)},'location','nw');
%         title_string = sprintf('%d x %d %s x %s' , I , J , seq1 , seq2);
%         title(title_string);
%         xlabel('Fitness')
%         ylabel('Fraction of seqs')
%         set(gcf,'PaperPosition',[0 0 5 5]);
%         print('-dpsc2',figname,'-append');
%         close;
%     end
    if (mod(I,100)==0) , fprintf('.') , end;
    if (mod(I,1000)==0) , fprintf('\n') , end;
end


%%
R = dataset();
R.FractionIntermediateUnfit = FractionIntermediateUnfit;
R.HammingDistances = HammingDistances;
R.FitnessDistributions = FitnessDistributions ;
R.NumIntermediateStatesMeasured = NumIntermediateStatesMeasured ;
R.IntermediateIDX = IntermediateIDX;

R = R( ~isnan(R.FractionIntermediateUnfit),:);

%%

figure; 
boxplot(R.FractionIntermediateUnfit*100,R.HammingDistances,'notch','on')
%xlim([3.5 10.5])
xlabel('Distance between pairs')
ylabel('% unfit (s < 0.45)')
set(gca,'ytick',0:5:50);
grid on; 
ylim([-1 25])
% 
% figure; 
% boxplot(cellfun( @(X)mean(X<0) , R.FitnessDistributions)*100,R.HammingDistances,'notch','on')
% xlim([3.5 10.5])
% xlabel('Distance between pairs')
% ylabel('% unfit (s < 0)')
% set(gca,'ytick',0:5:50);
% grid on; 
% ylim([-1 15])

%%
idx = find( ismember( T.aa_seq , all_transition_states) );
histogram(T.fitness(idx))
