%% figure for sign epistasis
% after runing EpistasisLBC_SignEpistasis

%% an example of sign epistasis
%  SegN    Pos    Perm        SubPos    SubPerm     p              logodds    X                   pBon       
%    2       28     'LV'        8         'TV'        4.9191e-104    39.6       [2×2 double]        3.2751e-100
threshold = 0.0000001  ; 
ridx = 1 ; 
biggidx =  BigG.SegN==R.SegN(ridx) & BigG.VarPos==R.VarPos(ridx) & strcmp(BigG.Perm,R.Perm{ridx}) ;
FitImpact = BigG.FitImpact{biggidx} ;
Seqs = BigG.Seq1{biggidx} ;
aas_at_subpos = cellfun( @(X)X( R.SubPos(ridx) ) , Seqs) ; 
uaas = unique(aas_at_subpos);

data = NaN(0);
for I = 1:numel(uaas)
data(1,I) = sum( FitImpact > threshold & aas_at_subpos ==  uaas(I) ) / sum(aas_at_subpos ==  uaas(I)) ;
data(2,I) = sum( FitImpact < -1*threshold & aas_at_subpos ==  uaas(I) )  / sum(aas_at_subpos ==  uaas(I)) ;
data(1,I) = sum( FitImpact > threshold & aas_at_subpos ==  uaas(I) )  / sum(aas_at_subpos ==  uaas(I)) ;
data(2,I) = sum( FitImpact < -1*threshold & aas_at_subpos ==  uaas(I) )  / sum(aas_at_subpos ==  uaas(I)) ;
data(1,I) = sum( FitImpact > threshold & aas_at_subpos ==  uaas(I) )  / sum(aas_at_subpos ==  uaas(I)) ;
data(2,I) = sum( FitImpact < -1*threshold & aas_at_subpos ==  uaas(I) )  / sum(aas_at_subpos ==  uaas(I)) ;
end
data = 100 * data ; 

txt  =  sprintf( '%s -> %s fitness impact' , BigG.Perm{biggidx}(1) , BigG.Perm{biggidx}(2) ) ;

% V1
figure ;
bar( data(1,:) ,'FaceColor',[.75 .75 .75]) ;
ylabel('% of the time the effect is +')
axis tight; 
set(gca,'ygrid','on')
set(gca,'ytick',0:10:100)
set(gca,'xticklabel',uaas)
xlabel( txt ) 

%% V2
clrs = get(gca,'ColorOrder') ;
figure ;
bh = bar( data ) ;
bh(1).FaceColor = clrs(1,:);
bh(2).FaceColor = clrs(2,:);
bh(3).FaceColor = clrs(3,:);

ylabel('% of the time the effect is')
set(gca,'xtick',[.775 1 1.225 1.775 2 2.225])
set(gca,'xticklabel', [ uaas ])

axis tight; 
set(gca,'ygrid','on')
set(gca,'ytick',0:10:100)
title( txt ) 
% 
% ylim([ 0 5 ])
% set(gca,'ytick',0:5)
% set(gca,'yticklabel',[0 1 2 3 4 round(max(data(:)))] )
% 
% 
% figure; 
% bar( log2( data(1,:) ./ data(2,:)) , 'FaceColor' , [.7 .7 .7])
% set(gca,'xticklabel',uaas)
% xlabel('AA at position 8')
% ylabel({txt 'log2( + / - )'})
%% how common is sign epi? 
N=50;
BigG.HasEnoughBigEffects = BigG.sum_MinorSignFitEffect > 2*N | BigG.sum_MajorSignFitEffect > 2*N ;
BigG.HasSignEpiByCount   = BigG.sum_MinorSignFitEffect > N & BigG.sum_MajorSignFitEffect > N ;
BigG.HasSignEpiByAUC     = BigG.AUC_noweights > 0.6   ;
BigG.HasSignEpiByCountAndAUC = BigG.HasSignEpiByCount & BigG.HasSignEpiByAUC ; 

fprintf('All subs = %d\n' , length(BigG))
fprintf('enough w/big effectsd = %d (%0.02f%%)\n' , sum(BigG.HasEnoughBigEffects) ,  100*sum(BigG.HasEnoughBigEffects)/length(BigG))
fprintf('enough w/both signs = %d (%0.02f%%)\n' , sum(BigG.HasSignEpiByCount), 100*sum(BigG.HasSignEpiByCount)/length(BigG))
fprintf('pass AUC = %d (%0.02f%%)\n' , sum(BigG.HasSignEpiByAUC), 100*sum(BigG.HasSignEpiByAUC)/length(BigG))
fprintf('pass AUC & enough w/both signs = %d (%0.02f%%)\n' , sum(BigG.HasSignEpiByCountAndAUC), 100*sum(BigG.HasSignEpiByCountAndAUC)/length(BigG))

%%
fh = figure('units','centimeters','position',[5 5 7 7]);
hold on; 
histogram( 100*BigG.mean_MinorSignFitEffect( BigG.HasSignEpiByCountAndAUC  ) , 10 ,'Normalization','Count','EdgeColor','k','FaceColor','k')
xlabel('% of genontypes with minor sign')
ylabel('# of substitutions')
axis tight ;
xlim([0 004.4])
set(gca,'xtick',0:1:100)
print('-dpsc2','SignEpistasis.eps','-append')
close ;

%%
