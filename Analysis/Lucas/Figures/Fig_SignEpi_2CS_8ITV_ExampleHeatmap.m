
%%
% The substitution C6S in segment two has an opposite effect on fitness depending on amino acid at state 8 (I, V or T), and the substitution I8T exhibits sign epistasis depending on the amino acid at site 28 (F,I,V or L). These interactions can be represented by a graph, whereby nodes represent a pair of extant amino acid states at a specific site and nodes are connected by edges if strong sign epistasis has been detected between them (C6S - I8T - I28F) for the example in the previous sentence). The dimensionality of the fitness landscape can then be determined by the dimensionality of the graph of all pairs of sites under sign epistasis (Figure 7). 
% example of sign epistasis. 
% May 10, 2017
% LBC
load SignEpi_Ronly.mat
%%
DS = s(2).R( strcmp( s(2).R.Perm , 'CS')  & s(2).R.VarPos == 6 , :);
%%
DS.LargePosEffect = DS.FitImpact > 0.4 ;
DS.LargeNegEffect = DS.FitImpact < -0.4 ;
DS.Neutral = ~DS.LargePosEffect & ~DS.LargeNegEffect  ;
DS.next_pos = cellfun( @(X) X(8) , DS.Seq1) ; 
G = grpstats( DS , 'next_pos' ,{'sum' 'mean'} ,'DataVars',{'LargePosEffect' 'LargeNegEffect' 'Neutral'})


%%

clrs = gray(100);
clrs = flipud(clrs);
figure;
data = [ G.sum_LargeNegEffect  G.sum_LargePosEffect ]' ;
%data = [85	1305	165 ; 389	54	195] ; % from .X in SignEpi table
imagesc(log10(data ))
set(gca,'xtick',1:3)
set(gca,'xticklabel',{'8 I' '8 T' '8 V'})
xlabel('Background')
set(gca,'ytick',1:2)
set(gca,'yticklabel',{'- effect' '+ effect'})
colormap(clrs)
colorbar( 'Ticks',[log10(50) 2 2.5 3]  ,'TickLabels',[50 100 500 1000])

%%
next_pos = cellfun( @(X) X(8) , DS.Seq1) ; 
unext_pos = unique(next_pos);

figure; 
hold on ;
[f,x] = ecdf( DS.FitImpact );
plot(x,f,'-k','LineWidth',2,'DisplayName','all');
for l = unique(next_pos)'
	idx = next_pos == l ;
    [f,x] = ecdf( DS.FitImpact(idx) );
    plot(x,f,'-','LineWidth',2,'DisplayName', sprintf('%s (N=%d)' , l , sum(idx) ) );
end
xlabel('Fit w/C  -  Fit w/S  (+ = C better)')
ylabel('Fraction of genotypes')
xlim( [ -0.5 0.5])
set(gca,'xtick',-10:.1:10)
set(gca,'ytick',0:0.05:1)
grid on ;
line([ -.4 -.4] , ylim)
line([ .4 .4] , ylim)
legend('location','nw')



%%
clrs = gray(100);
clrs = flipud(clrs);
figure;
data = [85	1305	165 ; 389	54	195] ; % from .X in SignEpi table
imagesc(log10(data ))
set(gca,'xtick',1:3)
set(gca,'xticklabel',{'8 I' '8 T' '8 V'})
xlabel('Background')
set(gca,'ytick',1:2)
set(gca,'yticklabel',{'+ effect' '- effect'})
colormap(clrs)
colorbar( 'Ticks',[log10(50) 2 2.5 3]  ,'TickLabels',[50 100 500 1000])
%title('Effect of a C->S substitution at position 6')
