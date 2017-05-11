% The observed fitness landscape of His3 must involve
% the interaction of many sites because 95% of substitutions between
% extant amino acid states had a strong effect on fitness (>|0.2|) in at
% least 42 different genetic backgrounds (Figure X).
%
%
% LBC May 04, 2017
%% load subs effect data
load('~/Desktop/HIS3scratch/SignEpi/SignEpi_Ronly.mat');



%%  figure version 1
figname = 'PctSubstitutionsWithAnyEfectEver_v1.eps';
delete(figname)
efs =  [0.2 0.25 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0 1.1 1.25] ;
N = 100 ; 
threshes = linspace(0,0.1,N);

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
    pct = NaN(3,N);
    for tI = 1:N
        pct(1,tI) = 100 * sum( G.sum_Pos ./ G.GroupCount >= threshes(tI) ) ./ length(G)  ;
        pct(2,tI) = 100 * sum( G.sum_Neg ./ G.GroupCount >= threshes(tI) ) ./ length(G);
        pct(3,tI) = 100 * sum( G.sum_Pos ./ G.GroupCount >= threshes(tI) | G.sum_Neg  ./ G.GroupCount  >= threshes(tI) ) ./ length(G);
    end
    
    fh = figure('units','centimeters','position',[5 5 10 10]);
    plot( 100*threshes , pct'  ,'LineWidth',3)
    
    
    ylabel('% of substitutions that can be')
    xlabel('Minimum % of backgrounds with effect')
    legend({ 'strong positive' 'strong negative' 'strong + or -' },'location','sw')
    grid on;
    set(gca,'ytick',0:5:100)
    set(gca,'xtick',0:100)
    
    title( sprintf('Min effect size = %0.02f' , MinLargeEffectSize))
    print('-dpsc2',figname,'-append');
    close;
    
end