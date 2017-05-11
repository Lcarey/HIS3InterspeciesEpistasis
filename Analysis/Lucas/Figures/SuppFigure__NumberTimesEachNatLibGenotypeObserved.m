%% Supp figure
% # of times each genotype is observed in a fit genetic background
% First, each of extant amino acid states was found in at least 50 different genotypes with wild-type fitness in our library (Supplementary Figure 4). 
% LBC May 1, 2017

%% load data
SparseVects = cell( 12 , 1);
Fitnesses = cell(12,1);
for SegN = 1:12
    T  = EpistasisLBC__LoadData( SegN  ,'ONLY_NATLIB_FLAG',true, 'ONLY_MIDDLE_FLAG',true,'NO_STOP_FLAG',false);
    T  = EpistasisLBC__FindVariableRegion_GenSparseMatrix( T ) ;
    SparseVects{SegN} = cell2mat( T.SparseVect);
    Fitnesses{SegN} = T.s ; 
end

%% 
is_fit_threshold = 0.4 ; 
is_unfit_threshold = 0.2 ; 
num_times_variant_observed_genotype_is_fit = NaN(0);
num_times_variant_observed_genotype_is_unfit = NaN(0);

for SegN = 1:12
    num_times_variant_observed_genotype_is_fit = [num_times_variant_observed_genotype_is_fit sum( SparseVects{SegN}( Fitnesses{SegN} > is_fit_threshold , :) ) ] ;
    num_times_variant_observed_genotype_is_unfit = [num_times_variant_observed_genotype_is_unfit sum( SparseVects{SegN}( Fitnesses{SegN} < is_unfit_threshold , :) ) ] ;
end

%%
fh = figure('units','centimeters','position',[5 5 7.5 7]);
histogram( log10(num_times_variant_observed_genotype_is_fit)  , 25 ,'FaceColor',[.7 .7 .7] ) 
set(gca,'xtick',1:5)
set(gca,'xticklabel',[10 1e2 1e3 1e4 1e5])
xlabel('# of times observed in fit genotypes')
ylabel('# of extant amino acid states')
xlim([1.5 4.5])
print('-dpng','# of times each genotype is observed in a fit genetic background.png','-r600');
print('-dpsc2','# of times each genotype is observed in a fit genetic background.eps');

%%
fh = figure('units','centimeters','position',[5 5 7.5 7]);
hold on ;
histogram( log10(num_times_variant_observed_genotype_is_unfit+5)  , 25 ,'FaceColor',[.7 .7 .7] ) 
set(gca,'xtick',[ log10(5) 2:5] )
set(gca,'xticklabel',[0 1e2 1e3 1e4 1e5])
xlabel('# of times observed in unfit genotypes')
ylabel('# of extant amino acid states')
xlim([0 4.7])
print('-dpng','# of times each genotype is observed in a UNfit genetic background.png','-r600');
print('-dpsc2','# of times each genotype is observed in a UNfit genetic background.eps');

%%
histogram2(  log10(num_times_variant_observed_genotype_is_fit) , log10(num_times_variant_observed_genotype_is_unfit+5)  )