%%  The statistical inference procedure used here is quite odd in the context of the literature.
%       It is based on a least squares fit to an exponential growth model. 
%       More typical in the literature is to consider log enrichment ratios or 
%          more generally log ratios of counts (e.g. Hietpas et al. 2011).
% 
% LBC January 2017

%% load data
clear all ; 
DATADIR = '~/Develop/HIS3InterspeciesEpistasis/Data/' ; 
FIGDIR = '~/Google Drive/CareyLab/Projects/Finished or Retired/HIS3/OLD/HIS3_Vika/Figures_For_Revisions/' ; 
for I = 1:12
    s(I).T = readtable([DATADIR 'S' num2str(I) '_scaled_info_v2.csv'],'Delimiter','\t','FileType','text'); 
end
%%
figname = 'Fitness_LogRatio_vs_OurCalculation_3syn' ; 
c = NaN(12,2);
figure; 
set(gca,'Visible', 'off') ;
ha = tight_subplot( 4 , 3 );
fontsize = 15 ; 
for I = 1:12
    T = s(I).T ; 
    eps = 0.01 ; 
    lr = log2( (T.t2_fr + eps) ./ T.t0_fr  );

    % nat_lib only or all?
    idx = T.nat_lib & ~T.nonsense & (T.len==mode(T.len)) & ~isnan(T.s) & ~isnan(lr) & ~isinf(T.s) & ~isinf(lr) ;
 %   idx = ~isnan(T.s) & ~isnan(lr) & ~isinf(T.s) & ~isinf(lr) ;
    idx = ~isnan(T.s) & ~isnan(lr) & ~isinf(T.s) & ~isinf(lr) & T.size>2 ; 
    
    axes( ha(I));
    sh = scatter(lr(idx) , T.s(idx) , 5 , 'k' ,'MarkerFaceColor',[.7 .7 .7]...
        ,'MarkerEdgeAlpha',0.1,'MarkerFaceAlpha',0.1);
    c(I,1) = corr(lr(idx),T.s(idx),'rows','complete') ;
    c(I,2) = corr(lr,T.s,'rows','complete') ;
    xlim([-6 2]);
    ylim([0 1]);
    text(-5.5 ,  0.93 , sprintf('seg %d',I) ,'FontSize',fontsize);
    text(-5.5 ,  0.78 , sprintf('\x3c1 = %0.02f',c(I,1)) ,'FontSize',fontsize);
end
set(ha(1:9),'XTickLabel','');
set(ha([2 3 5 6 8 9 11 12]),'YTickLabel','');
print('-dpsc2',[ FIGDIR figname ]);
print('-dpng2',[ FIGDIR figname ] ,'-r600');
