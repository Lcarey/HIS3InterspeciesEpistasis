%% no correlation between # syn sites and distance from S cer


figname = '/Users/lcarey/Dropbox/Pokusaeva17/Figures/SupplementaryFigures/Correlation_#SynVariants_DistFromScer_';
for SegN = 1:12
    T = readtable( sprintf('~/Develop/HIS3InterspeciesEpistasis/Data/S%d_scaled_info_v2.csv',SegN),'FileType','text','Delimiter','\t');

    idx = logical(T.nat_lib) & logical(T.nogap) & logical(T.middle); 
    ymax = 2+max(max( prctile(T.size,90) , max(T.size(T.dist_Scer==0))));    

    figure;
    subplot(2,1,1)
    boxplot( T.size(idx) , T.dist_Scer(idx),'notch','on')
    ylim([ 0 ymax ] )
    xlabel('Distance from S cer')
    ylabel('# syn variants')
    title([ num2str(SegN) ' nat lib' ])
    xl=xlim();
    grid on ;
 
    ymax = 2+max(max( prctile(T.size,90) , max(T.size(T.dist_Scer==0))));    
    subplot(2,1,2)
    boxplot( T.size  , T.dist_Scer ,'notch','on')
    ylim([ 0 ymax ] )
    xlabel('Distance from S cer')
    ylabel('# syn variants')  
    title([ num2str(SegN) ' all' ])
    xlim(xl)
    grid on; 
    print('-dpng',[ figname num2str(SegN) '.png'] );
    close;
end