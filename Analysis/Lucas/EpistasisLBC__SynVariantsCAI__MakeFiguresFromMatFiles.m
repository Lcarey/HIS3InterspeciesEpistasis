%% EpistasisLBC__SynVariantsCAI__MakeFiguresFromMatFiles

%% load mat files from cluster
%  EpistasisLBC__SynVariantsCAI()
MATDIR = 'outliersfilt_1std/';
matfiles = dir( [MATDIR 'SynVariantsCAI_*.mat']);
for I = 1:numel(matfiles)
    load( [MATDIR matfiles(I).name] );
    matfiles(I).LibN = s.segment ;
    matfiles(I).G = s.G ;
    matfiles(I).source_file_name = s.source_file_name ;
    matfiles(I).remove_outliers_method = s.remove_outliers_method ;
    matfiles(I).remove_outliers_value = s.remove_outliers_value ;
end
clear 's';
matfiles = sortStruct( matfiles , 'LibN',1)
%%
MinN = 10 ;
minGC  = 5;
min_CAI_range = 0.0 ;
ws = 20 ;

figure ; hold on;

clrs = parula(12);

ph = NaN(12,1);
for SegI = 1:12
    r_CAI = cellfun(@range,  matfiles(SegI).G.tAI_all);
    idx = matfiles(SegI).G.GroupCount >= minGC & r_CAI >= min_CAI_range ;
    fit_range = linspace( min(matfiles(SegI).G.median_s) , max(matfiles(SegI).G.median_s), 200);
    V = 1+(ws/2):(numel(fit_range)-(ws/2)) ;
    results = NaN( numel(V) , 3);
    for I = V
        idxF = matfiles(SegI).G.median_s >= fit_range(I-ws/2) & matfiles(SegI).G.median_s <= fit_range(I+ws/2)  ;
        C = ( matfiles(SegI).G.tAI_C(idx & idxF) );
        if numel(C > MinN)
            results(I,1) = mean(C);
            [~,results(I,2)] = ttest(C);
            results(I,3) = fit_range(I);
        end
        % if numel(C)>10
        % plot( fit_range(I) , mean(C) ,'o','Color',clrs(SegI,:));
        % if ttest(C)
        %     plot( fit_range(I) , mean(C) ,'p','Color',clrs(SegI,:));
        %  end
        % end
    end
    if sum(results(:,2)<0.05) > 20
       
    X = results(:,3);
    Y = smooth( X , results(:,1));
    ph(SegI) =  plot( X , Y  ,'-','DisplayName',num2str(SegI),'LineWidth',1);
    %ph(SegI) =  plot( X , Y  ,'-','Color',clrs(SegI,:),'DisplayName',num2str(SegI),'LineWidth',1);

    plot( results(results(:,2)<0.01,3) , results( results(:,2)<0.01 ,1)  ,'o','Color',get(ph(SegI),'Color'),'DisplayName',num2str(SegI) ,'MarkerFaceColor',get(ph(SegI),'Color'));
    
    end
end

ph = ph( ~isnan(ph));

xlim([-1 1.5])
xlabel('Fitness')
ylabel('Fitss vs tAI mean correlation')
ylim([-0.05 0.15])
grid on;
line(xlim,[0 0],'LineStyle','-','Color',[.5 .5 .5])
lh =legend(ph,'location','best')
title(MATDIR)
%% plot examples w/the most syn variants
SegI = 8;
G = matfiles(SegI).G ;
G = sortrows(G,{'GroupCount'},'descend');
G.r_CAI = cellfun(@range, G.tAI_all);
G = G( G.r_CAI>0.1  & G.GroupCount> 15 ,:);
G = G(G.median_s<0.75 & G.median_s>-0.75,:);


figname = 'Corrs.eps';
delete(figname)
for I = 1:height(G)
    Y = G.fit_all{I} ;
    X = G.tAI_all{I} ;
    
    idx = Y > ( median(Y)-2*std(Y)) & Y < ( median(Y)+2*std(Y)) ;
    Y = Y(idx);
    X = X(idx);
    [c,p] = corr(X,Y);
    [cS,pS] = corr(X,Y,'Type','Spearman');
    
    if p<0.2
    figure; 
    subplot(2,1,1)
    hold on; 
    plot(X,Y','ok','MarkerFaceColor',[.75 .75 .75]);

    [xData, yData] = prepareCurveData( X, Y );
    ft = fittype( 'poly1' );
    opts = fitoptions( 'Method', 'LinearLeastSquares' );
    opts.Robust = 'LAR';opts.Robust = 'off';
    [fitresult, gof] = fit( xData, yData, ft, opts );
    plot(fitresult)
    legend('off')
    title( sprintf('C=%0.02f , p=%0.02f , std=%0.02f' , c , p , std(X)) )
    ylabel('Fitness')
    xlabel('Codon Bias (tAI)')
    grid on;

    subplot(2,1,2);
    [grps,edgs]=discretize(X,4);
    boxplot( Y , grps );
    set(gca,'xtick',1:4);
    set(gca,'xticklabel',arrayfun(@(I)mean([ edgs(I) edgs(I+1) ]) , 1:4))
    ylabel('Fitness')
    xlabel('Codon Bias (tAI)')  
    
    set(gcf,'PaperPosition',[0 0 4 4]);
    print('-dpsc2',figname,'-append');
    close;
    end
end
    
%% testing
%plot( matfiles(10).G.nTE_C , log10(matfiles(10).G.nTE_P),'ok')
I=7;
v_min_fit = NaN(0);
v_max_fit = NaN(0);
v_min_CAI_range = NaN(0);
v_minGC = NaN(0);
v_mean_nTE_C = NaN(0);
v_mean_tAI_C = NaN(0);
v_N = NaN(0);

%r_CAI = cellfun(@range,  matfiles(I).G.nTE_all); %
r_CAI = cellfun(@range,  matfiles(I).G.tAI_all);
% tAI better
for minGC = 1:10
    for min_CAI_range = 0:0.005:0.05 ;
        for min_fit = -1:0.1:0
            for max_fit = 0:0.1:1
                if max_fit > min_fit
                    idx1 = r_CAI >= min_CAI_range  & matfiles(I).G.GroupCount >= minGC ;
                    idx2 = matfiles(I).G.median_s >= min_fit & matfiles(I).G.median_s <= max_fit ;
                    idx = idx1 & idx2;
                    if sum(idx)>100
                        v_minGC(end+1) = minGC ;
                        v_min_CAI_range(end+1) = min_CAI_range ;
                        v_mean_nTE_C(end+1) = mean( matfiles(I).G.nTE_C(idx));
                        v_mean_tAI_C(end+1) = mean( matfiles(I).G.tAI_C(idx));
                        v_N(end+1) = sum(idx);
                        v_min_fit(end+1) = min_fit ;
                        v_max_fit(end+1) = max_fit ;
                    end
                end
            end
        end
    end
end

R = dataset();
R.minGC = v_minGC';
R.min_CAI_range = v_min_CAI_range';
R.mean_nTE_C = v_mean_nTE_C';
R.mean_tAI_C = v_mean_tAI_C';
R.max_fit = v_max_fit' ;
R.min_fit = v_min_fit' ;
R.N = v_N';

