%% supplementary figure showing
%  that WT seqs from other species have WT growth rates when put on a
%  plasmid in S cer
%
%
%Furthermore, the epistatic interactions must be intragenic because when the entire 
% sequence of His3 from other species is transfected into baker´s yeast there was no effect 
% on fitness (Figure [SFullWTs]).
%
% SupFig_FullHIS3otherSpeciesHaveWTgrowthrates
%  LBC April 2017

%% load data
T = readtable('~/Develop/HIS3InterspeciesEpistasis/Data_Small_Tables/control.csv','FileType','text','Delimiter','\t');
T = T( cellfun(@length,T.aa_seq) < 10 & cellfun(@length,T.aa_seq) > 1 ,:);
T = T(: , [3:end]);
T = sortrows(T , 'aa_seq' );

T.aa_seq = regexprep( T.aa_seq , '(\d)$' , '_$1');


%% make figure
figure; hold on; 
bar( T.mean ,'FaceColor',[.7 .7 .7]);

idx = find(regexpcmp(lower(T.aa_seq),'wt'));
bar( idx ,  T.mean(idx) ,'FaceColor','r');


idx = find(regexpcmp(lower(T.aa_seq),'scer'));
bar( idx ,  T.mean(idx) ,'FaceColor','b');


errorbar( 1:height(T) , T.mean , T.sem , 'ok');
ylabel('Mean growth rate (hr^{-1})');
set(gca,'xtick',1:height(T))
set(gca,'xticklabel',T.aa_seq)
xtickangle(45);
axis tight;
set(gca,'ytick',0:0.1:1);
ylim([0 0.55])

line(xlim,[0.2 0.2],'LineStyle','--','Color','k')
line(xlim,[0.4 0.4],'LineStyle','--','Color','k')
line(xlim,[0.325 0.325],'LineStyle','--','Color','c')
