% %% EpistasisLBC__CompareUnfitFreqBetweenFit_VS_EntireLibrary Figures
% % generate all data on laptop (takes ~24 hrs)
% N_Pairs_Fast_to_measure = 5e4 ;
% for I = 1:12
%     q = EpistasisLBC__CompareUnfitFreqBetweenFit_VS_EntireLibrary(   I , N_Pairs_Fast_to_measure );
%     save( sprintf('Segment%d.mat' , I ) ,  'q');
% end

%%  make inidividual segment figures
figbasename = 'CompareUnfitFreqBetweenFit_VS_EntireLibrary_';
segment_files = dir('Segment*.mat');
for SegFileI = 1:numel(segment_files)
    load( [ segment_files(SegFileI).folder filesep segment_files(SegFileI).name]);
    T = q.T ;
    R = q.R ;
    SegN = q.SegN ;
    figname = [ figbasename num2str(SegN) '.eps'];
    delete(figname);
    
    % Fitness for variants as a fucntion of distance from S cer
    clrs = hot(max(T.DistFromScer)+5) ;
    ph = NaN(0);
    fh = figure('units','centimeters','position',[5 5 7 10 ]);
    hold on ;
    for I = 1:max(T.DistFromScer)
        [f,x] = ecdf(T.s( T.DistFromScer==I));
        ph(I) = plot(x,f,'-','Color',clrs(I,:),'DisplayName',sprintf('Dist_{Scer}=%d' , I),'LineWidth',2);
    end
 %   legend( ph ,  'location','nw')
    ylabel('Fraction of variants')
    xlabel('Fitness')
    xlim([0 0.47])
    print('-dpsc2',figname,'-append');
    close;
    
    % Fitness of intermediate states as a function of distance between the pairs
    clrs = hot(max(R.HammingDistances)+5) ;
    ph = NaN(0);
    fh = figure('units','centimeters','position',[5 5 7 10 ]);
    hold on ;
    for I = 1:max(R.HammingDistances)
        [f,x] = ecdf(  cell2mat(R.FitnessDistributions(  R.HammingDistances == I ) ) );
        ph(I) = plot(x,f,'-','Color',clrs(I,:),'DisplayName',sprintf('Dist_{pairs}=%d' , I),'LineWidth',2);
    end
%    legend( ph ,  'location','nw')
    ylabel('Fraction of variants')
    xlabel('Fitness')
    xlim([0 0.47])
    print('-dpsc2',figname,'-append');
    close;
    
    %  Direct comparison of the two
    N = max(R.HammingDistances);
    between_pairs = NaN(N,1);
    from_Scer  = NaN(N,1);
    for I = 1:N
        between_pairs(I) = 100 * mean( cell2mat(R.FitnessDistributions(  R.HammingDistances == I )) < 0.2);
        from_Scer(I) = 100 * mean(T.s( T.DistFromScer==I)  < 0.2);
    end
    
    fh = figure('units','centimeters','position',[5 5 7 5 ]);
    bh = bar( [between_pairs ,   from_Scer  ] );
    bh(1).FaceColor=[0 0 1] ; bh(1).EdgeColor=[0 0 1] ;
    bh(2).FaceColor=[1 0 0] ;  bh(2).EdgeColor=[1 0 0] ;
    bh(1).BarWidth = 1 ;
    bh(2).BarWidth = 1 ;
    ylabel('% unfit variants')
    xlabel('Distance from Scer or between fit pairs')
  %  legend({'between fit pairs' 'from S. cer' },'location','sw')
    set(gca,'Ygrid','on')
    axis tight
    set(gca,'ytick',0:5:100);
    if( max(max(between_pairs),max(from_Scer)) < 0.10)
        set(gca,'ytick',0:1:100);
    end
    xlim([0.9 N+.4])
    print('-dpsc2',figname,'-append');
    close;    
    %
    
    fh = figure('units','centimeters','position',[5 5 7 5 ]);
    Y = ((from_Scer-between_pairs) ./ from_Scer)  ;  %  as a percent
    Y =  from_Scer - between_pairs  ;  % just the difference
    plot( 1:N ,  Y ,'ok','MarkerFaceColor',[.7 .7 .7])
    %xlim([0.9 7.5])
    set(gca,'xtick',0:100)
    set(gca,'ytick',-100:10:100)
    if (range(Y)<30)
        set(gca,'ytick',-100:5:100)
    end
    set(gca,'Ygrid','on')
    axis tight; 
    ylabel('% intersegnental epistasis')
    xlabel('Mutational distance')

    % for the exponential fit idea
%     Y = 100* ((from_Scer-between_pairs) ./ from_Scer)' ;
%     X = [ 2:numel(Y) 220] ; % His3 is 220aa long
%     Y = [ Y(2:end) 0] ;
%     
    print('-dpsc2',figname,'-append');
    close;    
end

%% make one figure combining data for all segments

figname = 'CompareUnfitFreqBetweenFit_VS_EntireLibrary.eps';
delete(figname);
segment_files = dir('Segment*.mat');

dist_at_which_to_calc = 7 ; 
fast_thresh  = 0.4  ;
slow_thresh  = 0.2 ; 

DS = dataset();
DS.SegN = NaN( 12*15 , 1);
DS.pct_unfit_between_pairs = NaN( 12*15 , 1);
DS.pct_unfit_from_Scer   = NaN( 12*15 , 1);
DS.dist_at_which_to_calc = NaN(12*15,1);
c = 0;
for I = 1:numel(segment_files)
    load( [ segment_files(I).folder filesep segment_files(I).name]);
    for dist_at_which_to_calc = 1:15 ; 
        c = c+1 ; 
        DS.SegN(c) = q.SegN ;   
        DS.pct_unfit_between_pairs(c) = 100*mean( cell2mat(q.R.FitnessDistributions(  q.R.HammingDistances == dist_at_which_to_calc )) < slow_thresh );
        DS.pct_unfit_from_Scer(c) = 100*mean(q.T.s( q.T.DistFromScer==dist_at_which_to_calc)  < slow_thresh );    
        DS.dist_at_which_to_calc(c) = dist_at_which_to_calc;
    end
end
DS.difference  =  DS.pct_unfit_from_Scer - DS.pct_unfit_between_pairs ;
%DS.difference  = 100* ((DS.pct_unfit_from_Scer-DS.pct_unfit_between_pairs) ./ DS.pct_unfit_from_Scer)  ;

DS = DS( ~isnan(DS.difference) ,:);
%% figure
% fh = figure('units','centimeters','position',[5 5 9 5 ]);
% bar( double(DS(:,2:end))' )
% set(gca,'xticklabel' , {'between pairs' 'from S. cer' 'difference'})
% ylabel('difference')
% set(gca,'Ygrid','on')
% title('% unfit at dist=7')
values_to_plot = {'difference' 'pct_unfit_from_Scer' 'pct_unfit_between_pairs'} ;


G = grpstats( DS , 'dist_at_which_to_calc' ,{'mean' 'median' 'std'} )
DS = DS( DS.dist_at_which_to_calc < 12 ,:);
nu = numel(unique(DS.dist_at_which_to_calc));

d=0.15 ; 

fh = figure('units','centimeters','position',[5 5 7 7 ]); 
hold on ;


bh = boxplot(  DS.difference , DS.dist_at_which_to_calc ,'Positions', (1:nu)  ,'PlotStyle','compact','Symbol','','Color',[.75 .75 .75],'LabelOrientation','horizontal' )


bh = boxplot(  DS.pct_unfit_from_Scer , DS.dist_at_which_to_calc ,'Positions', (1:nu)-d  ,'PlotStyle','compact','Symbol','','Color',[1 0.75 0.75 ],'LabelOrientation','horizontal' )

bh = boxplot(  DS.pct_unfit_between_pairs , DS.dist_at_which_to_calc ,'Positions', (1:nu)+d  ,'PlotStyle','compact','Symbol','','Color',[0.75 0.75 1 ],'LabelOrientation','horizontal' )

plot( G.dist_at_which_to_calc , G.median_difference ,'.k','MarkerSize',15)
[fitresult1, gof1] = fit( G.dist_at_which_to_calc, G.median_difference , fittype( 'power1' ) );
plot( 0:.1:100  , fitresult1(  0:.1:100 ) ,'-k')

plot( (1:nu)-d  , G.median_pct_unfit_from_Scer ,'.','Color',[1 0 0 ],'MarkerSize',15)
[fitresult2, gof2] = fit( G.dist_at_which_to_calc, G.median_pct_unfit_from_Scer , fittype( 'power1' ) );
plot( 0:.1:100  , fitresult2(  0:.1:100 ) ,'-','Color',[1 0 0 ])


plot( (1:nu)+d  , G.median_pct_unfit_between_pairs ,'.','Color',[0 0 1 ],'MarkerSize',15)
[fitresult3, gof3] = fit( G.dist_at_which_to_calc, G.median_pct_unfit_between_pairs , fittype( 'power1' ) );
plot( 0:.1:100  , fitresult3(  0:.1:100 ) ,'-','Color',[0 0 1 ])


set(gca,'Ygrid','on')
ylabel('Amount of epsitasis')
set(gca,'ytick',[0   10:10:100 ] );
set(gca,'yticklabel',  arrayfun(@(X)sprintf('%d%%',X) , get(gca,'ytick') ,'UniformOutput',false) );
ylim([-1 40])
xlim([1.3 11.6])
%set(gca,'xtick',1:2:20) ; set(gca,'xticklabel',1:2:20)
xlabel('Mutational distance')
%print('-dpng','Inter-segmental epistasis increases faster than intrasegmental 1.png','-r600');
%close;

fprintf('Diff b=%0.03f R^2=%0.02f\n' , fitresult1.b , gof1.rsquare);
fprintf('From Scer b=%0.03f R^2=%0.02f\n' , fitresult2.b , gof2.rsquare); 
fprintf('From pairs b=%0.03f R^2=%0.02f\n' , fitresult3.b , gof3.rsquare); 

%% Fit all segment seperatly
usegs = unique(DS.SegN);
R = dataset();
R.SegN = NaN(12,1);
for I = 1:numel(usegs)
    idx = ( DS.SegN == usegs(I)  ) ;
    R.SegN(I) = usegs(I) ;
    [R.fitresultPAIR{I}, R.gofPAIR{I}] = fit( DS.dist_at_which_to_calc(idx), DS.pct_unfit_between_pairs(idx) , fittype( 'power1' ) );
    [R.fitresultSCER{I}, R.gofSCER{I}] = fit( DS.dist_at_which_to_calc(idx), DS.pct_unfit_from_Scer(idx) , fittype( 'power1' ) );
    [R.fitresultDIFF{I}, R.gofDIFF{I}] = fit( DS.dist_at_which_to_calc(idx), DS.difference(idx) , fittype( 'power1' ) );
    R.diff_b(I) = R.fitresultDIFF{I}.b ;
    R.pair_b(I) = R.fitresultPAIR{I}.b ;
    R.scer_b(I) = R.fitresultSCER{I}.b ;
    R.diff_r2(I) = R.gofDIFF{I}.rsquare ;
    R.pair_r2(I) = R.gofPAIR{I}.rsquare ;
    R.scer_r2(I) = R.gofSCER{I}.rsquare ;    
end

% only decent fits
R = R( R.scer_r2>0.5,:) ; 

figure; hold on; 
subplot(1,2,1);
%plot( R.diff_b , R.scer_b ,'ok','MarkerFaceColor',[.7 .7 .7]);
text( R.diff_b , R.scer_b , arrayfun(@num2str, R.SegN) )
xlabel('difference')
ylabel('S cer ')
xlim([0 6])
ylim(xlim)
line(xlim,xlim)

subplot(1,2,2);
%plot( R.diff_b , R.pair_b ,'ok','MarkerFaceColor',[.7 .7 .7]);
text( R.diff_b , R.pair_b , arrayfun(@num2str, R.SegN) )
xlabel('difference')
ylabel('Pairs')

suptitle('b exponent from fit')
xlim([0 6])
ylim(xlim)
line(xlim,xlim)
