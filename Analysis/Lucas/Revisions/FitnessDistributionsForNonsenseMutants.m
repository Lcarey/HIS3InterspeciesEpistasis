%% show that nonsense mutatants have low fitness
%
% LBC January 2017


DATADIR = '~/Develop/HIS3InterspeciesEpistasis/Data/' ; 
FIGDIR = '~/Google Drive/CareyLab/Projects/Finished or Retired/HIS3/OLD/HIS3_Vika/Figures_For_Revisions/' ;

%% load data
for I = 1:12
    s(I).T = readtable([DATADIR 'S' num2str(I) '_scaled_info_v2.csv'],'Delimiter','\t','FileType','text'); 
end
%% plot data
figname = 'FitnessDistributionsForNonsenseMutants_NoFilter';

figure;
set(gca,'Visible', 'off') ;

clrs = get(gca,'ColorOrder');
ha = tight_subplot( 4 , 3 );
totals_nonsense  = NaN(0);
totals_nearWT  = NaN(0);

for I = 1:12
        T = s(I).T ; 
        T.s(T.s>1) = 1 ; 
      %   T = T( T.t0_fr > prctile(T.t0_fr,50) , :); %filter by abundance
      %   T = T( T.size > 2 ,:); % filter by #syn
        idx_nonsense = logical(T.nonsense) ; 
        idx_natlib = T.nat_lib & ~T.nonsense & T.nogap & (T.len==mode(T.len)) ; 
        idx_wt1 = idx_natlib & T.dist_Scer < 2 ; 
        idx_wt0 = idx_natlib & T.dist_Scer < 1 ; 
       
        totals_nonsense = vertcat( totals_nonsense , T.s(idx_nonsense));
        totals_nearWT = vertcat( totals_nearWT , T.s(idx_wt1));
        
        axes( ha(I));
        xlim([0 1]) ; 
        hold on; 
     	line([0.8 0.8],ylim,'LineStyle','--','Color',[.7 .7 .7])
        line([0.4 0.4],ylim,'LineStyle','--','Color',[.7 .7 .7])
 
        [f,x]=ecdf(T.s);
        plot(x,f,'-','LineWidth',3,'Color',clrs(1,:))
        
        if sum(idx_wt1)>0
            [f,x]=ecdf(T.s(idx_wt1));
            plot(x,f,'-','LineWidth',3,'Color',[.5 .5 .5])
        end
        
        if sum(idx_wt0)>0
            [f,x]=ecdf(T.s(idx_wt0));
            plot(x,f,'-','LineWidth',3,'Color','k')
        end
        
        if sum(idx_nonsense)>0
            [f,x]=ecdf(T.s(idx_nonsense));
        	plot(x,f,'-','LineWidth',3,'Color',clrs(2,:))
            text( 0.1,  0.8 , sprintf('%0.01f%% of nonsense < 0.4',  mean(T.s(idx_nonsense)<0.4)*100) );
        end
        
       
        text( 0.1,  0.9 , sprintf('seg %d',I) );
        
end
set(ha(1:9),'XTickLabel','');
set(ha([2 3 5 6 8 9 11 12]),'YTickLabel','');
set(ha(10:12),'Xtick',0:.2:1)
set(ha(10:12),'XtickLabel',0:.2:1)
set(ha([1 4 7 10]),'Ytick',0:.2:1)
set(ha([1 4 7 10]),'YtickLabel',0:.2:1)
print('-dpng',[FIGDIR figname],'-r300');

fprintf('%0.04f%% nonsense < 0.4\n' , mean( totals_nonsense < 0.4) * 100 )
fprintf('%0.04f%% near WT > 0.4\n' , mean( totals_nearWT > 0.4) * 100 )
fprintf('%0.04f%% near WT > 0.8\n' , mean( totals_nearWT > 0.8) * 100 )
