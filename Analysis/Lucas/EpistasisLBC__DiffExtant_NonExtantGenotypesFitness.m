function [ R  , T ] = EpistasisLBC__DiffExtant_NonExtantGenotypesFitness(    )
%%  [ R , T ] = EpistasisLBC__DiffExtant_NonExtantGenotypesFitness(   )
%
%  .  Approximately an equal fraction of genotypes, those composed of amino acids from fit
%    genotypes and those incorporating other amino acid states, show low fitness (Figure[F11] ),
%
% LBC 2017

fh = figure('units','centimeters','position',[5 5 30 15 ]);
hold on ;

R = dataset();
R.SegN = (1:12)';

R.f_NOSTOP_NOTEXTANT = NaN(12,1);
R.f_NOSTOP_EXTANT = NaN(12,1);


fast_fit_cutoff = 0.4  ;
slow_fit_cutoff = 0.2 ;

for SegN = 1:12
    subplot(2,6,SegN);
    hold on ;
    DataDir = '~/Develop/HIS3InterspeciesEpistasis/Data/';
    T = readtable([ DataDir filesep num2str(SegN) '.tab' ],'FileType','text') ;
    
    
    ph = NaN(0);
    stop_idx = logical(T.stop) | logical(T.nonsense) ;
    
    
    
    [f,x]=ecdf(T.s(stop_idx));
    ph(1)=plot(x,f,'-k','LineWidth',2,'DisplayName','nonsense') ;
    [f,x]=ecdf(T.s(~stop_idx));
    ph(2)=plot(x,f,'-','LineWidth',2,'DisplayName','extant') ;
    [f,x]=ecdf(T.s(~stop_idx & ~T.nat));
    ph(3)=plot(x,f,'-','LineWidth',2,'DisplayName','not extant') ;
    
 %   xlabel('Fitness')
 %   ylabel('Fraction of genotypes')
    
    grid on ;
    set(gca,'ytick',0:.2:1)
    set(gca,'yticklabel', ['0' arrayfun( @(X)regexprep( sprintf('%0.01f',X) , '^0' , '') , .2:.2:.8 ,'UniformOutput' , false ) '1' ] )
 
    set(gca,'xtick',0:.1:1)
    set(gca,'xticklabel', ['0' arrayfun( @(X)regexprep( sprintf('%0.01f',X) , '^0' , '') , .1:.1:.9 ,'UniformOutput' , false ) '1' ] )
    
    xlim([-0.01 0.49])
    
    line([fast_fit_cutoff fast_fit_cutoff],ylim,'LineStyle','--','Color','k')
    line([slow_fit_cutoff slow_fit_cutoff],ylim,'LineStyle','--','Color','k')
    
    %legend( ph  ,'location','se');
    
    title(SegN)
    
    
    
    R.f_NOSTOP_NOTEXTANT(SegN) = mean(T.s( ~stop_idx & ~T.nat)) ;
    R.f_NOSTOP_EXTANT(SegN) = mean(T.s( ~stop_idx & T.nat)) ;
end

set(gcf,'PaperPosition',[0 0 15 5])
print('-dpsc2', 'F11__seg_fit_dist.eps' );
%close;

%%
% G,CLR,SYM,SIZ,DOLEG,XNAM,YNAM)
clrs = parula(13);
clrs = clrs(1:12,:);
fh = figure('units','centimeters','position',[5 5 8 8  ]);
xlim([0 0.5])
ylim(xlim)
hold on ;
gh = gscatter( R.f_NOSTOP_EXTANT , R.f_NOSTOP_NOTEXTANT , R.SegN , clrs , 'o+' , 10);
for I = 1:numel(gh)
    set(gh(I),'LineWidth',3)
end
xlim([0 0.5])
ylim(xlim)
set(gca,'xtick',0:.1:1)
set(gca,'ytick',0:.1:1)
line(xlim,xlim,'Color',[.7 .7 .7])
lh = legend(gh,'location','nw');
set(lh,'box','on','color',[1 1 1])

xlabel({'Extant amino acids' 'mean library fitness' })
ylabel({'mean library fitness ' 'Non-extant aminio acids'})

set(gcf,'PaperPosition',[0 0 5 5])

