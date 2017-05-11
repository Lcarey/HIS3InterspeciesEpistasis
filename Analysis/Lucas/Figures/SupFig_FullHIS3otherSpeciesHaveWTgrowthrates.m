%% supplementary figure showing
%  that WT seqs from other species have WT growth rates when put on a
%  plasmid in S cer
%
%
%Furthermore, the epistatic interactions must be intragenic because when the entire 
% sequence of His3 from other species is transfected into baker´s yeast there was no effect 
% on fitness (Figure [SFullWTs]).
%
% Second, the entire His3 coding sequence from the extant species fully complements a His3 deletion in  S. cerevisiae (Supplementary Figure 5). 
% SupFig_FullHIS3otherSpeciesHaveWTgrowthrates
%
%  LBC April 2017

%% load data
T = readtable('~/Develop/HIS3InterspeciesEpistasis/Data_Small_Tables/control.csv','FileType','text','Delimiter','\t');
T = T( cellfun(@length,T.aa_seq) < 10 & cellfun(@length,T.aa_seq) > 1 ,:);
T = T(: , [3:end]);
T = sortrows(T , 'aa_seq' );

T.aa_seq = regexprep( T.aa_seq , '(\d)$' , '_$1');


%% make figure
fh = figure('units','centimeters','position',[5 5 17 10]);
hold on; 
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
set(gca,'ytick',0:0.05:1);
ylim([0 0.55])
set(gca,'ygrid','on')
print('-dpsc2','Full His3 other species have WT growth rates v1.eps')
print('-dpng','Full His3 other species have WT growth rates v1.png','-r600')
close ; 
% line(xlim,[0.2 0.2],'LineStyle','--','Color','k')
% line(xlim,[0.4 0.4],'LineStyle','--','Color','k')
% line(xlim,[0.325 0.325],'LineStyle','--','Color','c')

%% version 2
T = readtable('~/Develop/HIS3InterspeciesEpistasis/Data_Small_Tables/control.csv','FileType','text','Delimiter','\t');
T = T( ~regexpcmp( T.Comments , 'Has multiple peaks') ,:);
T = sortrows(T,'mean','descend');

fh = figure('units','centimeters','position',[5 5 10 10]);
hold on; 

wtidx = strcmp(T.x_library ,'wt');
stopidx = regexpcmp(T.aa_seq,'_');
histogram( T.mean , 100 , 'FaceColor','k')
histogram( T.mean(wtidx) , 10  , 'FaceColor','r')
histogram( T.mean(stopidx) , 20  , 'FaceColor','b')
ylabel('# of genotypes')
xlabel('Mean growth rate (hr^{-1})');
legend({'all' 'WT' 'nonsense'},'location','nw')
xlim( [ 0 max(T.mean)*1.025])
print('-dpsc2','Full His3 other species have WT growth rates v2.eps')
print('-dpng','Full His3 other species have WT growth rates v2.png','-r600')
close ; 