%% load data
DataDir = '~/Develop/HIS3InterspeciesEpistasis/Data/';
T = readtable([ DataDir 'S7_scaled_info.csv']) ;
T = T( T.len == mode(T.len)  , :);


%%
figure; hold on; 
ecdf(T.s);
ecdf(T.s(logical(T.stop)));
ecdf(T.s(~logical(T.stop)));
ecdf(T.s(~logical(T.stop) & logical(T.nat)))
ecdf(T.s(~logical(T.stop) & logical(T.lib)))
ecdf(T.s(~logical(T.stop) & logical(T.nat_lib)))
ecdf(T.s(~logical(T.stop) & logical(T.nat_lib) & logical(T.middle)) );
xlabel('Fitness')
title('S7')
set(gca,'xtick',0:.05:1)
grid on ;
set(gca,'ytick',0:.1:1)
legend({'All' '+stop' '-stop' '-stop nat'  '-stop lib'  '-stop nat-lib' '-stop nat-lib & middle' } ,'location','best')
%%
figure;hold on; 

ecdf(T.s(logical(T.stop)));
T = T( ~logical(T.stop)  , :);

T = T( logical(T.nat_lib)  , :);

ecdf(T.s);

legend({'stop' 'nat_lib'})
grid on; set(gca,'xtick',0:0.1:1)


%%
%%
% % % % MapAA2I = containers.Map(  arrayfun(@(X){X},['A':'Z' '_' ])  , uint8(1:27) ) ; % all AAs + stop
% % % % MapI2AA = containers.Map( uint8(1:27) ,   arrayfun(@(X){X},['A':'Z' '_' ])  ) ; % all AAs + stop
% % % % 
% % % % aa_num = cellfun( @(A) arrayfun( @(X)MapAA2I(X),A) , T.aa_seq ,'UniformOutput' , false ); 

variation_in_seq_mat = mean(cell2mat(aa_num) == aa_num{1} ) ~= 1  ; % these are the columns for which the AA seqs are not all exactly the same

T.aa_seq_varies = cellfun( @(X) X(variation_in_seq_mat) , T.aa_seq,'UniformOutput',false); 
%%

T = T( : , {'aa_seq' 'aa_seq_varies'  's'});
fast_fit_cutoff = 0.425 ; 
slow_fit_cutoff = 0.25  ; 

fast_seqs =  T.aa_seq_varies( T.s >=  fast_fit_cutoff ) ;
slow_seqs =  T.aa_seq_varies( T.s <=  slow_fit_cutoff ) ;

% choose N random pairs of fast sequences
N_Pairs_Fast_to_measure = 1e4 ; 
pairs = randi( numel(fast_seqs) , N_Pairs_Fast_to_measure*2 , 2) ;
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


for I = 1:N_Pairs_Fast_to_measure
    seq1 = fast_seqs{pairs(I,1)};
    seq2 = fast_seqs{pairs(I,2)};
    all_transition_states = ExpandSeqAlign( seq1 , seq2);
    idx = find( ismember( T.aa_seq_varies , all_transition_states) );
    R.FitnessDistributions{I} = T.s(idx) ;    
    R.NUnfitMeasured(I) = sum( R.FitnessDistributions{I} <= slow_fit_cutoff);
    R.NFitMeasured(I) = sum( R.FitnessDistributions{I} >= fast_fit_cutoff);
    R.NIntMeasured(I) = sum( R.FitnessDistributions{I} < fast_fit_cutoff &  R.FitnessDistributions{I} > slow_fit_cutoff);

    R.HammingDistances(I) = HammingDistance( seq1 , seq2 ) ;
    R.StatesMeasured(I) = numel(idx)  ;
    if (mod(I,100)==0) , fprintf('.') , end;
    if (mod(I,1000)==0) , fprintf(' %0.0f%%\n',I/N_Pairs_Fast_to_measure*100) , end;
end

R.pct_unfit  = ( R.NUnfitMeasured ) ./ R.StatesMeasured * 100 ; 
R.pct_states_measured = R.StatesMeasured ./ (2.^R.HammingDistances) .* 100 ;
%%

figure;

subplot(2,2,1)
boxplot( R.pct_unfit , R.HammingDistances ,'notch','on')
xlabel('Distance between pairs')
ylabel(sprintf('%% unfit (s < %0.2f)' , slow_fit_cutoff ) )
set(gca,'ytick',0:5:100);
grid on; 
title('median')

subplot(2,2,2); hold on;
bar(unique(R.HammingDistances), grpstats(R.pct_unfit, R.HammingDistances))
xlabel('Distance between pairs')
ylabel(sprintf('%% unfit (s < %0.2f)' , slow_fit_cutoff ) )
set(gca,'ytick',0:1:100);
set(gca,'xtick',1:20)
grid on; 
title('mean')


subplot(2,2,3); hold on;
bar(unique(R.HammingDistances), grpstats(  ( R.NIntMeasured ) ./ R.StatesMeasured * 100 , R.HammingDistances))
xlabel('Distance between pairs')
ylabel(sprintf('%% (s>%0.2f & s<%0.02f)' , slow_fit_cutoff , fast_fit_cutoff) )
set(gca,'ytick',0:5:100);
set(gca,'xtick',1:20)
grid on; 
title('mean')


subplot(2,2,4); hold on;
bar(unique(R.HammingDistances), grpstats(  ( R.NFitMeasured ) ./ R.StatesMeasured * 100 , R.HammingDistances))
xlabel('Distance between pairs')
ylabel(sprintf('%% s>=%0.02f)' , fast_fit_cutoff) )
set(gca,'ytick',0:5:100);
set(gca,'xtick',1:20)
grid on; 
title('mean')

subtitle('S7 all same length')