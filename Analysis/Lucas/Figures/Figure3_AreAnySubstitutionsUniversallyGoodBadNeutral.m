% Fig_AreAnySubstitutionsUniversallyGoodBadNeutral
% Remarkably, not a single amino acid state was found to be universally deleterious or universally neutral.
%
% May 15, 2017
% load the pair-wise FitImpacts for each sub

load('~/Desktop/HIS3scratch/SignEpi/SignEpi_Ronly.mat');

%% write out data to .tab files
for SegN = 1:12
    writetable( dataset2table(s(SegN).R) , [num2str(SegN) '.tab'] , 'FileType','text','Delimiter','\t');
    system( [ 'gzip --fast ' num2str(SegN) '.tab'] )
end



%%
figname = 'How many substitutions are neutral, only good or only bad Main Fig ';
threshold_change = 0.4 ;
min_FRACTION_with_change =  0.02 ;

for SegN = 1:12
    Q = s(SegN).R ;
    Q.BadFI = Q.FitImpact <= (-1 * threshold_change) ;
    Q.GoodFI = Q.FitImpact >= threshold_change ;
    Q.SegN = repmat(SegN , length(Q) , 1);
    if SegN == 1
        G = grpstats( Q , {'VarPos' 'Perm' 'SegN'} ,{'sum' 'mean'},'DataVars',{'GoodFI' 'BadFI'} );
    else
        G = vertcat( G ,  grpstats( Q , {'VarPos' 'Perm'  'SegN'} ,{'sum' 'mean'},'DataVars',{'GoodFI' 'BadFI'} ));
    end
end
G.PctNeutral = (1-(( G.sum_GoodFI + G.sum_BadFI) ./ G.GroupCount) ) * 100 ; 
G.PctGood    = ((( G.sum_GoodFI ) ./ G.GroupCount) ) * 100 ; 
G.PctBad     = ((( G.sum_BadFI) ./ G.GroupCount) ) * 100 ; 

vn = {'VarPos' 'Perm' 'SegN' 'GroupCount' 'PctNeutral' 'PctGood' 'PctBad'};
writetable( dataset2table( G( : , vn)) , 'PctGoodBadNeutral.tab','FileType','text','Delimiter','\t');
%%
X = double( G( : , {'PctNeutral' 'PctGood' 'PctBad'}));
%X = X( G.SegN ~= 9 ,:);
X0 = log10(X+1);
%[~,o] = sortrows( geomean( X(:,2:3) , 2) );
X0 = X0(o,:);
X0(:,1) = round(X0(:,1)*1)/1 ;
%X0 = round(X0*100)/100 ; 
X0 = sortrows( X0 , 1:3 ,{'ascend' 'descend' 'descend'});

figure;
imagesc(X0);
set(gca,'xtick',1:3);
set(gca,'xticklabel',{'% neutral' '% good' '% bad'})

set(gca,'ytick',[])
ylabel('Genotypes')
%title('log10( % )')
colormap( flipud(winter(100)) )
colorbar('Ticks' , [0:.5:3])
%colorbar('Ticks' , [0:.5:3] ,'TickLabels' , [0 25 50 75 100])

%%
X = sortrows( X , 1:3 ,{'ascend' 'descend' 'descend'});
X = sortrows( X , 2:3 ,{'descend' 'descend'});

figure; hold on; 
for I = 1:nrows(X)
    plot([0 X(I,3)],[I I],'-r');
    plot([X(I,3) X(I,3)+X(I,1) ],[I I],'-','Color',[.7 .7 .7]);
    plot([100-X(I,2) 100],[I I],'-b');
end
axis tight;
ylabel('Genotypes')
xlabel('% observed')

%%
figure;
gscatter(G.PctGood,G.PctBad,G.SegN)
xlabel('% backgrounds increase > 0.4')
ylabel('% backgrounds decrease > -0.4')
%T = clusterdata( X , 'distance' , 'correlation'  , 'linkage' ,'centroid' , 'maxclust',10)

% % 
% % 
% % %%
% % G.PercentNeutral =  G.sum_BadFI <= min_N_with_change & G.sum_GoodFI <= min_N_with_change
% % G.OnlyNeutral    =  G.sum_BadFI <= min_N_with_change & G.sum_GoodFI <= min_N_with_change ;
% % G.NeutralAndGood =  G.sum_BadFI <= min_N_with_change & G.sum_GoodFI <= min_N_with_change ;
% % 
% % 
% % fh = figure('units','centimeters','position',[5 5 4 6]);
% % data = 100 * [ mean( G.sum_BadFI <= min_N_with_change & G.sum_GoodFI <= min_N_with_change) ...
% %     mean( G.sum_BadFI <= min_N_with_change & G.sum_GoodFI > min_N_with_change) ...
% %     mean( G.sum_BadFI > min_N_with_change & G.sum_GoodFI <= min_N_with_change) ...
% %     mean( G.sum_BadFI > min_N_with_change & G.sum_GoodFI > min_N_with_change)  ] ;
% % bar(data ,'FaceColor','r','FaceAlpha',0.5)
% % set(gca,'xtick',[])
% % %set(gca,'xticklabel',{'neutral', '+' '-' '+ & -'})
% % axis tight ;
% % xlim([.5 4.5])
% % ylabel('% of extant substitutions')
% % %xlabel('Effect of substitution')
% % 
% % print('-dpsc2',sprintf( '%s %0.02f %0.02f.eps' , figname ,threshold_change , min_FRACTION_with_change) ,'-append');
% % close;

