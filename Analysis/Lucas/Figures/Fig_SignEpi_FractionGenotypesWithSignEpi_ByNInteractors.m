%% Figure 6
% Fig_SignEpi_FractionGenotypesWithSignEpi_ByNInteractors
% We find that 40% of all considered sites exhibit sign epistasis (Figure 6) 
%
%
% ecdf()
% compare # of interacting sights that determine sign epi
%
% LBC May 11, 2017
WORKDIR =  '/Users/lcarey/Desktop/HIS3scratch/SignEpi/';
load( [ WORKDIR 'SignEpi_Pairs_04.mat' ] )
%%

R.HasSign = R.pBon < 0.05 ;
G = grpstats( R , {'SegN' 'VarPos' },'sum','DataVars','HasSign');

figure;
[x,c]=count_unique(G.sum_HasSign) ;

% 0 , 1-5 6-10 11-15 16-20
    cnts = [ c(x==0) sum(c(x>0 & x<6)) sum(c(x>=6 & x<=10)) ...
    sum(c(x>10 )) ]
bar( cnts , 'FaceColor',[.75 .75 .75])
set(gca,'xtick',1:5)
set(gca,'xticklabel',{'0' '1-5' '6-10' '>10'})
xlabel('# of interacting pairs per site')
ylabel('# of sites with sign epistasis')
xlim([0.5 4.5])

fprintf('%d / %d (%0.02f%%)\n '  , sum(G.sum_HasSign>0) , length(G) , mean(G.sum_HasSign>0) *100)

%% how many have recirpocal sign epi