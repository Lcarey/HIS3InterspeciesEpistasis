% %% EpistasisLBC__CompareUnfitFreqBetweenFit_VS_EntireLibrary Figures
% % generate all data on laptop (takes ~24 hrs)
% N_Pairs_Fast_to_measure = 5e4 ;
% for I = 1:12
%     q = EpistasisLBC__CompareUnfitFreqBetweenFit_VS_EntireLibrary(   I , N_Pairs_Fast_to_measure );
%     save( sprintf('Segment%d.mat' , I ) ,  'q');
% end

%%
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
        between_pairs(I) = mean( cell2mat(R.FitnessDistributions(  R.HammingDistances == I )) < 0.2);
        from_Scer(I) = mean(T.s( T.DistFromScer==I)  < 0.2);
    end
    
    fh = figure('units','centimeters','position',[5 5 7 5 ]);
    bh = bar( 100* [between_pairs ,   from_Scer  ] );
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
    plot( 1:N ,  100* ((from_Scer-between_pairs) ./ from_Scer)  ,'ok','MarkerFaceColor',[.7 .7 .7])
    %xlim([0.9 7.5])
    set(gca,'xtick',0:100)
    set(gca,'ytick',-100:10:100)
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

