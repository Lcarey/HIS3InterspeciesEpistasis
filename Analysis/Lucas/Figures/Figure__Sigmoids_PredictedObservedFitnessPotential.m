%% Figure : sigmoids showing fitness potential vs fitness, and R^2
figname = '/Users/lcarey/Dropbox/Pokusaeva17/Figures/Figure__Sigmoids_Predicted_vs_Observed.png';
delete(figname);
RDIR = '~/Develop/HIS3InterspeciesEpistasis/Analysis/Katya/NN/residuals/';
fh = figure('units','centimeters','position',[5 5 40 10]);
hold on;
for SegN = 1:12
    subplot(2,12,SegN); hold on; 
    R = readtable( [ RDIR 'S' num2str(SegN) '.csv'] ,'FileType','text','Delimiter',',');
    T = readtable( sprintf('~/Develop/HIS3InterspeciesEpistasis/Data/S%d_scaled_info_v2.csv',SegN),'FileType','text','Delimiter','\t');
    T = T( logical(T.nat_lib) & logical(T.nogap) & logical(T.middle) ,:);
    T = innerjoin( R , T(:,{'aa_seq','s' 't0_fr' 'size' }) ,'Key','aa_seq');
    
    dscatter( T.fitnessPotential   , T.observed , 'msize' , 5 );
    plot( T.fitnessPotential   , T.predictedMinusObserved+T.observed,'.k','MarkerSize',1)
    axis tight;
    ylim([0  1.2])
    set(gca,'xtick',[]); set(gca,'ytick',[]); 
    %  set(gca,'xtick',-100:2:100);  set(gca,'xticklabel',100:-2:-100);
    % xlabel('Predicted fitness potential');  ylabel('Measured fitness')
    title( sprintf( 'Seg %d' , SegN ) )
    
    subplot(2,12,SegN+12) ; hold on; 
    dscatter(  (T.predictedMinusObserved+T.observed)   , T.observed , 'msize' , 5 )
    line([0 2],[0 2],'Color','k','LineStyle','-','LineWidth',1);
    xlim([ 0 1]); % set(gca,'xtick',0:0.2:10);
    ylim([0  1.2]); % set(gca,'ytick',-0:0.2:10);
   set(gca,'xtick',[]); set(gca,'ytick',[])

     %   xlabel('Predicted fitness')
     %   ylabel('Measured fitness')
    title( sprintf( 'R^2=%0.02f' , rsquare( T.observed , T.predictedMinusObserved+T.observed) ) )
    
    
end
print('-dpng',figname,'-r600');


%%   Misc figures for data exploration here (not for paper) %% 

%%
clrs = parula(3);
figname = '/Users/lcarey/Dropbox/Pokusaeva17/Figures/SupplementaryFigures/Figure__Predicted_vs_Observed_cutoff';
RDIR = '~/Develop/HIS3InterspeciesEpistasis/Analysis/Katya/NN/residuals/';
for SegN = 1:12
    R = readtable( [ RDIR 'S' num2str(SegN) '.csv'] ,'FileType','text','Delimiter',',');
    T = readtable( sprintf('~/Develop/HIS3InterspeciesEpistasis/Data/S%d_scaled_info_v2.csv',SegN),'FileType','text','Delimiter','\t');
    T = T( logical(T.nat_lib) & logical(T.nogap) & logical(T.middle) ,:);
    T = innerjoin( R , T(:,{'aa_seq','s' 't0_fr' 'size' }) ,'Key','aa_seq');
    
    switch_sign = 1 ;
    if corr(T.fitnessPotential , T.s) > 0
        switch_sign = -1;
    end
    
    
    
    rsq1 = NaN(100,2);
    for I = 1:100
        idx = T.size>=I;
        rsq1(I,1) = rsquare( T.observed(idx) , T.predictedMinusObserved(idx)+T.observed(idx));
        rsq1(I,2) =  sum(idx) / length(idx) * 100 ;
    end
    % fh = figure('units','centimeters','position',[5 5 15 10]);
    % [ax,h1,h2]=plotyy( 1:length(rsq1) , rsq1(:,1) , 1:length(rsq1) , rsq1(:,2));
    % ax(1).YLabel.String = 'R^2'  ;
    % ax(2).YLabel.String =  '% kept' ;
    % xlabel('min # syn variants');
    % title([ 'Segment ' num2str(SegN) ])
    % grid on ;
    
    
    rsq2 = NaN(1000,2);
    ft = linspace(0,100,1000);
    for I = 1:numel(ft)
        idx = T.t0_fr>=ft(I);
        rsq2(I,1) = rsquare( T.observed(idx) , T.predictedMinusObserved(idx)+T.observed(idx));
        rsq2(I,2) =  sum(idx) / length(idx) * 100 ;
    end
    
    rsq3 = NaN(100,2);
    ft = linspace(0,100,100);
    ntries = 100 ;
    for I = 1:numel(ft)
        N = sum(T.t0_fr>=ft(I)) ;
        jxn = NaN(ntries,1);
        for J = 1:ntries
            idx = randsample( height(T) , N );
            jxn(J) = rsquare( T.observed(idx) , T.predictedMinusObserved(idx)+T.observed(idx));
        end
        rsq3(I,1) = median(jxn);
        rsq3(I,2) =  N / height(T) * 100 ;
    end
    % fh = figure('units','centimeters','position',[5 5 15 10]);
    % [ax,h1,h2]=plotyy( ft, rsq2(:,1) , ft , rsq2(:,2));
    % ax(1).YLabel.String = 'R^2'  ;
    % ax(2).YLabel.String =  '% kept' ;
    % xlabel('min t0_{fr}');
    % title([ 'Segment ' num2str(SegN) ])
    % grid on ;
    fh = figure('units','centimeters','position',[5 5 4 4]);
    hold on;
    plot( rsq3(:,2) , rsq3(:,1),'-','LineWidth',3,'Color',[.7 .7 .7])
    plot( rsq2(:,2) , rsq2(:,1),'-','LineWidth',3,'Color',clrs(2,:))
    plot( rsq1(:,2) , rsq1(:,1),'-','LineWidth',3,'Color',clrs(1,:))
    
    
    %legend({'t0_{fr}' '# syn'})
    title([ 'segment ' num2str(SegN) ])
    grid on;
    ylim([0.4 1])
    xlim([0 100])
    xlabel('% kept')
    ylabel('R^2')
    set(gca,'xtick',0:25:100)
    set(gca,'ytick',0:0.1:1)
    print('-dpng2',[figname '_' num2str(SegN)]);
    close;
    
end

%%
%%
figname = '/Users/lcarey/Dropbox/Pokusaeva17/Figures/SupplementaryFigures/Figure__Predicted_vs_Observed_cutoff';
RDIR = '~/Develop/HIS3InterspeciesEpistasis/Analysis/Katya/NN/residuals/';
SegN = 5
R = readtable( [ RDIR 'S' num2str(SegN) '.csv'] ,'FileType','text','Delimiter',',');
T = readtable( sprintf('~/Develop/HIS3InterspeciesEpistasis/Data/S%d_scaled_info_v2.csv',SegN),'FileType','text','Delimiter','\t');
T = T( logical(T.nat_lib) & logical(T.nogap) & logical(T.middle) ,:);
T = innerjoin( R , T(:,{'aa_seq','s' 't0_fr' 'size' }) ,'Key','aa_seq');

switch_sign = 1 ;
if corr(T.fitnessPotential , T.s) > 0
    switch_sign = -1;
end
for I = 1:100
    idx = T.size>=I;
    fh = figure('units','centimeters','position',[5 5 4 4]);
    X =  (T.predictedMinusObserved(idx)+T.observed(idx))   ;
    Y=T.observed(idx)  ;
    dscatter( X , Y, 'msize' , 5 )
    xlabel('predicted')
    ylabel('measured')
    title( sprintf('%d %0.0f%% %0.02f' , I , mean(idx)*100 ,  corr(X,Y) ) )
    axis tight;
    print('-dpsc2',[ figname '_S5.eps'],'-append');
    close;
end
