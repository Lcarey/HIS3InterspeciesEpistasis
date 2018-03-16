%% generate a table with fitness & FitPot data for synonymous variants 
%     for Laura
% LBC June 2017
SegN = 7 ;

%% Load data
SYN = readtable('~/Develop/HIS3InterspeciesEpistasis/Data/synonymous_variants_rescaled_data.tab' ,'FileType','text');
SYN = SYN( SYN.SegN == SegN  , :);
FITPOT = readtable( ['~/Develop/HIS3InterspeciesEpistasis/Analysis/Katya/NN/residuals/S' num2str(SegN) '.csv']);


%% create a joined table and save output
T = innerjoin( FITPOT , SYN , 'Key' , 'aa_seq');
writetable( T , ['~/Downloads/SynonymousVariantsData__' num2str(SegN) '.tab'] ,'FileType','text','Delimiter','\t');

%% a few example plots

% plot the sigmoid
figure;
plot( T.fitnessPotential , T.observed ,'.k')
xlabel('Fitness Potential')
ylabel('Fitness (measured)')
title( sprintf( 'Segment %d' , SegN) )
axis tight ;
set(gca,'xtick',-20:1:20)
set(gca,'ytick',0:.1:2)
grid on; 

%% find an AA seq with lots of synonymous variants, and look at the fitness distribution of the syn variants
G = grpstats( T , 'aa_seq' , 'median' ,  'DataVars' , 's');
G = sortrows( G  , 'GroupCount' ,'descend');

most_common_group_nt_seq_index = strcmp( T.aa_seq , G.aa_seq{1} ) ;

figure;
histogram( T.s( most_common_group_nt_seq_index) , 0:.01:0.55 )
xlabel('Fitness')
ylabel('# of synonymous variants')
title( T.aa_seq( find(most_common_group_nt_seq_index,1)))