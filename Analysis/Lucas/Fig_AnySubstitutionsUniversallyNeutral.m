% The observed fitness landscape of His3 must involve
% the interaction of many sites because 95% of substitutions between
% extant amino acid states had a strong effect on fitness (>|0.2|) in at
% least 42 different genetic backgrounds (Figure X).
%
%
% LBC May 04, 2017
%% load subs effect data
load('~/Desktop/HIS3scratch/SignEpi/SignEpi_Ronly.mat');


%%  figure version 2
XtoCnt = 42 ;
figname = 'PctSubstitutionsWithAnyEfectEver_v2.eps';
delete(figname);

for MinLargeEffectSize = [0.05 0.1 0.2 0.3 0.4]
    
    N=100;
    
    for SegN = 1:12
        R = s(SegN).R  ;
        R.Pos = R.FitImpact > MinLargeEffectSize ;
        R.Neg = R.FitImpact < (-1 * MinLargeEffectSize) ;
        R.Neutral = R.FitImpact > (-1 * MinLargeEffectSize)  & R.FitImpact < MinLargeEffectSize ;
        R.SegN = repmat(SegN,length(R),1);
        if SegN==1
            G = grpstats( R, {'SegN' 'VarPos' 'Perm'} , 'sum' ,'DataVars' , {'Pos' 'Neg' 'Neutral'}) ;
        else
            G = vertcat( G , grpstats( R, {'SegN'  'VarPos' 'Perm'} , 'sum' ,'DataVars' , {'Pos' 'Neg' 'Neutral'}));
        end
    end
    
    
    cnt = NaN(12,N);
    for ThreshN = 1:N
        for SegN = 1:12
            cnt(SegN,ThreshN) = 100 *  sum( (G.sum_Pos >= ThreshN | G.sum_Neg >= ThreshN) & G.SegN==SegN ) / sum(G.SegN==SegN) ;
        end
        
    end
    
    fh = figure('units','centimeters','position',[5 5 10 10]);
    hold on;
    clrs = parula(12) ;
    for I = 1:12
        plot( 1:N , cnt(I,:),'Color',clrs(I,:) ,'LineWidth',3,'DisplayName',num2str(I))
    end
    ylabel('% of substitutions that effect fitness')
    xlabel('Minimum # of backgrounds with effect')
    grid on;
    set(gca,'ytick',0:5:100)
    set(gca,'xtick',0:10:N)
    title( sprintf('Min effect size = %0.02f' , MinLargeEffectSize))
    print('-dpsc2',figname,'-append');
    close;
    
    fh = figure('units','centimeters','position',[5 5 10 10]);
    bar( cnt(:,XtoCnt) ,'FaceColor',[.7 .7 .7])
    ylabel('% of substitutions that effect fitness')
    xlabel('Segment')
    grid on;
    set(gca,'ytick',0:5:100)
    set(gca,'xtick',1:12)
    title( sprintf('Min effect size = %0.02f ' , MinLargeEffectSize));
    axis tight;
    print('-dpsc2',figname,'-append');
    close;
    
end

%%  figure version 1
figname = 'PctSubstitutionsWithAnyEfectEver_v1.eps';
delete(figname)
efs =  [0.05 0.1 0.2 0.3 0.4] ;
cnts = NaN(12,100);


for I = 1:numel(efs)
    MinLargeEffectSize = efs(I);
    for SegN = 1:12
        R = s(SegN).R  ;
        R.Pos = R.FitImpact > MinLargeEffectSize ;
        R.Neg = R.FitImpact < (-1 * MinLargeEffectSize) ;
        R.Neutral = R.FitImpact > (-1 * MinLargeEffectSize)  & R.FitImpact < MinLargeEffectSize ;
        R.SegN = repmat(SegN,length(R),1);
        if SegN==1
            G = grpstats( R, {'SegN' 'VarPos' 'Perm'} , 'sum' ,'DataVars' , {'Pos' 'Neg' 'Neutral'}) ;
        else
            G = vertcat( G , grpstats( R, {'SegN'  'VarPos' 'Perm'} , 'sum' ,'DataVars' , {'Pos' 'Neg' 'Neutral'}));
        end
    end
    
    
    %
    cnts = NaN(3,100);
    
    for ThreshN = 1:100
        cnt(1,ThreshN) = sum( G.sum_Pos >= ThreshN ) ;
        cnt(2,ThreshN) = sum( G.sum_Neg >= ThreshN ) ;
        cnt(3,ThreshN) = sum( G.sum_Pos >= ThreshN | G.sum_Neg >= ThreshN ) ;
        
    end
    
    fh = figure('units','centimeters','position',[5 5 10 10]);
    plot( 1:100 ,  (100 * cnt) / length(G) ,'LineWidth',3)
    
    
    ylabel('% of substitutions that can be')
    xlabel('Minimum # of backgrounds with effect')
    legend({ 'strong positive' 'strong negative' 'strong + or -' },'location','sw')
    grid on;
    set(gca,'ytick',0:5:100)
    set(gca,'xtick',0:10:100)
    
    title( sprintf('Min effect size = %0.02f' , MinLargeEffectSize))
    print('-dpsc2',figname,'-append');
    close;
    
end