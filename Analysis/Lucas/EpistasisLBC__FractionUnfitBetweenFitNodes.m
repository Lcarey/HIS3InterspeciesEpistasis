function [ R  , T ] = EpistasisLBC__FractionUnfitBetweenFitNodes( fitness_file_csv  , N_Pairs_Fast_to_measure , fast_fit_cutoff ,  slow_fit_cutoff )
%%  [ R , T ] = EpistasisLBC__FractionUnfitBetweenFitNodes( fitness_file_csv  , N_Pairs_Fast_to_measure , fast_fit_cutoff ,  slow_fit_cutoff )
%
%
% LBC 2017
%% load data
DataDir = '~/Develop/HIS3InterspeciesEpistasis/Data/';

if isnumeric(fitness_file_csv)
    fitness_file_csv = sprintf('S%d_scaled_info.csv',fitness_file_csv) ;
    T = readtable([ DataDir filesep  fitness_file_csv ]) ;
elseif exist( fitness_file_csv , 'file')
    T = readtable(fitness_file_csv) ;
elseif exist( [ DataDir filesep fitness_file_csv ] , 'file')
    T = readtable([ DataDir filesep fitness_file_csv ]) ;
else
   error([ 'Cant find ' fitness_file_csv]);
end

if ~exist('fast_fit_cutoff','var')
    fast_fit_cutoff = 0.4 ; 
end
if ~exist('slow_fit_cutoff','var')
    slow_fit_cutoff = 0.25  ; 
end

base_file_name = regexp( fitness_file_csv , filesep , 'split')
base_file_name = regexprep( base_file_name{end} , '\....' , '') ; 

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
set(gca,'xtick',0:.05:1)
grid on ;
set(gca,'ytick',0:.1:1)
legend({'All' '+stop' '-stop' '-stop nat'  '-stop lib'  '-stop nat-lib' '-stop nat-lib & middle' } ,'location','best');
xlim([0 0.55])

line([fast_fit_cutoff fast_fit_cutoff],ylim,'LineStyle','--','Color','k')
line([slow_fit_cutoff slow_fit_cutoff],ylim,'LineStyle','--','Color','k')

title(base_file_name)

set(gcf,'PaperPosition',[0 0 5 5])
print('-dpng', [ base_file_name '_dist.png'] , '-r300');
close; 

%% Must remove stop codons for results to work!
T = T( ~logical(T.stop) & ~logical(T.nonsense)  ,:);
T = T( logical(T.middle) ,:);
T = T( logical(T.lib) ,:);


T = T( T.len == mode(T.len)  , :);

T = T( logical(T.nat_lib) , :);


MapAA2I = containers.Map(  arrayfun(@(X){X},['A':'Z' '_' ])  , uint8(1:27) ) ; % all AAs + stop
MapI2AA = containers.Map( uint8(1:27) ,   arrayfun(@(X){X},['A':'Z' '_' ])  ) ; % all AAs + stop

%% find all columns that vary
%   shrink the sequence down to only those columns
aa_num = cellfun( @(A) arrayfun( @(X)MapAA2I(X),A) , T.aa_seq ,'UniformOutput' , false ); 
aa_num_mat = cell2mat(aa_num) ; 
n_aa_per_col = NaN( 27 , length(T.aa_seq{1}) ) ; 
for I = 1:size(n_aa_per_col,2)
    [a,b]=count_unique(aa_num_mat(:,I));
    n_aa_per_col(a,I) = b ;
end
cols_with_variation = arrayfun( @(I)sum(n_aa_per_col(:,I)>1) , 1:size(n_aa_per_col,2)) > 1  ; 

%variation_in_seq_mat = mean(aa_num_mat == aa_num{1} ) ~= 1    ; % these are the columns for which the AA seqs are not all exactly the same

T.aa_seq_varies = cellfun( @(X) X(cols_with_variation) , T.aa_seq,'UniformOutput',false); 
%%

T = T( : , {'aa_seq' 'aa_seq_varies'  's'});

fast_seqs =  find( T.s >=  fast_fit_cutoff ) ;
slow_seqs =  find( T.s <=  slow_fit_cutoff ) ;

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
    if (mod(I,20)==0) , fprintf('.') , end;
    if (mod(I,100)==0) , fprintf(' ') , end;
    if (mod(I,1000)==0) , fprintf(' %0.0f%%\n',I/N_Pairs_Fast_to_measure*100) , end;
end
toc
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

subtitle( sprintf('%s nat-lib & middle & nostop (%d)' , base_file_name, height(T) ) )

set(gcf,'PaperPosition',[0 0 10 7])
print('-dpng', [ base_file_name '.png'] , '-r300');
close; 


save( [ base_file_name '.mat']  , 'T' , 'R')

end
