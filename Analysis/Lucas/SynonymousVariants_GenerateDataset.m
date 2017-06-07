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