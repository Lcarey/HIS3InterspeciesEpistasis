%% Fitness_LogRatio_vs_OurCalculation_natlib
% Our calculation assigns genotypes that drop out of the competition to the maximum likely fitness based on the initial population size (see methods). That said, the two methods are highly correlated (Extended Data Figure X (Fitness_LogRatio_vs_OurCalculation_natlib)). 
% 
% April 2018 LBC

%% load data
DATADIR = '~/Develop/HIS3InterspeciesEpistasis/Data/' ;
C = cell(12,1);
for SegN = 1:12
    T = readtable( [DATADIR 'S' num2str(SegN) '_scaled_info_v2.csv'] , 'FileType','text','Delimiter','\t');
    idx = T.nat_lib & ~T.stop & T.nogap & ~T.nonsense & T.t2_fr>0 & T.t1_fr>0 ;
    C{SegN} = T(idx,:);
end

%% plot

fh = figure('units','centimeters','position',[5 5 25 20]);
for SegN = 1:12
    subplot(3,4,SegN)
    S = C{SegN}.s ; 
    LR = log2( C{SegN}.t2_fr  ./ C{SegN}.t0_fr) ;
    idx = ~isinf(LR) & ~isinf(S) ; 
    S = S(idx);
    LR = LR(idx); 
    [c(SegN),p] = corr(S,LR,'Rows','Complete');
    sh = scatter( S , LR , 20 , 'k' ,'MarkerFaceColor',[.7 .7 .7],'MarkerEdgeAlpha',0.1,'MarkerFaceAlpha',0.15);
    text(0.1,4, sprintf('corr = %0.02f' , c(SegN)))
    xlim([0 1.2])
    set(gca,'xtick',0:0.2:2)
    ylim([-5 5])
    set(gca,'ytick',-10:2:10)
    title( sprintf( 'segment %d' , SegN))
end
FIGDIR = '~/Google Drive/Private_CareyLab/Manuscripts/Pokuseeva__HIS3/His3_Figures_March2018/LBC/' ; 
print('-dpng',[ FIGDIR 'Fitness_LogRatio_vs_OurCalculation_natlib.png'] ,'-r300');
close ; 