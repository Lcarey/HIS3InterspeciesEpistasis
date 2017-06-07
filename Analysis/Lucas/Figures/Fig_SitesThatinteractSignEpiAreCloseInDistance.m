%% Figure to show that residues w/sign epi 
%
%
cd('/Users/lcarey/Desktop/HIS3scratch/SignEpi/StructureStuff--PairwiseDistancesAndSignEpi')

%% ecdf
%fns = {'03.repSign.res' '03.SignEpiPairs.FT.res' '03.SignEpiPairs.PB.res' '03.SignEpiPairs.TF.res'};
fns = {'03.SignEpiPairs.FT.res' '03.SignEpiPairs.TF.res' '03.repSign.res'  };
names = {'no sign' 'sign' 'recip sign'};
fh = figure('units','centimeters','position',[5 5 6 6]);
hold on ;
clrs = get(gca,'ColorOrder') ;
for fi = 1:numel(fns)
     T = readtable( fns{fi} ,'FileType','text');
     [f,x] = ecdf( T.Var3);
     plot( x , f*100 , 'LineWidth',2 , 'Color' , clrs(fi,:) , ...
         'DisplayName' ,  names{fi});
end
legend('location','best');
xlabel('Pairwise distances (Å)')
ylabel('Cummulative % of pairs')
axis tight; 
set(gca,'ytick',0:10:100)
set(gca,'xtick',0:5:100)
grid on 
xlim([0 39])
%% ecdf log2(Anstroms)
ph = NaN(0);
fh = figure('units','centimeters','position',[5 5 6 6]);
hold on ;
clrs = get(gca,'ColorOrder') ;
for fi = 1:numel(fns)
     T = readtable( fns{fi} ,'FileType','text');
     [f,x] = ecdf( log2(T.Var3) );
     ph(fi) = plot( x , f*100 , 'LineWidth',2 , 'Color' , clrs(fi,:) , ...
         'DisplayName' ,  names{fi});
end
legend('location','nw');
xlabel('Pairwise distances log2(Å)')
ylabel('Cummulative % of pairs')
axis tight; 
set(gca,'ytick',0:10:100)
set(gca,'xtick',0:1:100)
grid on ;
xlim([1.9 5])
%%
data = NaN(1e6,numel(fns));
for fi = 1:numel(fns)
     T = readtable( fns{fi} ,'FileType','text');
     data( 1:height(T) , fi) = T.Var3 ; 
end
data = log2(data);
figure; 
bh = boxplot( data ,'notch','on' ,'PlotStyle','Compact','Color',[.85 .85 .85])
ylabel('Pairwise distances log2(Å)')
set(gca,'xticklabel', {'no sign' 'sign' 'repi sign'})
for fi = 1:3
    set(bh(1,fi),'Color', get(ph(fi),'Color'))
    set(bh(2,fi),'Color', get(ph(fi),'Color'))
end
[~,p12]=ttest2( (data(:,1)) , (data(:,2)))
[~,p13]=ttest2( (data(:,1)) , (data(:,3)))
[~,p23]=ttest2( (data(:,2)) , (data(:,3)))

%%
NoSign = data( ~isnan(data(:,1)),1);
Sign = data(~isnan(data(:,2)),2);
RepicSign = data(~isnan(data(:,3)),3);
NoSignL2 = log2(NoSign);
NoSignL10 = log10(NoSign);


figure; hold on; 
[parmhat,parmci] = lognfit(NoSign);
errorbar( 1 , parmhat(1) , parmci(1,1)-parmhat(1) , parmci(2,1)-parmhat(1) ,'ok')

[parmhat,parmci] = lognfit(Sign)
errorbar( 2 , parmhat(1) , parmci(1,1)-parmhat(1) , parmci(2,1)-parmhat(1) ,'ok')

[parmhat,parmci] = lognfit(RepicSign)
errorbar( 3 , parmhat(1) , parmci(1,1)-parmhat(1) , parmci(2,1)-parmhat(1) ,'ok')

set(gca,'xtick',1:3)
set(gca,'xticklabel', {'no sign' 'sign' 'repi sign'})
axis tight; 
grid on ;
xlim([0 4])
