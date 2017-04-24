%% calculate error using synonymous variants
SynDir = '~/Develop/HIS3InterspeciesEpistasis/Data/';
T = readtable( [ SynDir 'syndata.tab'] ,'ReadVariableNames',false , 'Format','%d%s%f' ,'FileType','text');
T.Properties.VariableNames = {'SegN' 'aa_seq' 'fitness'};
T = T( ~isnan(T.fitness) ,:);
%% rescale to [0 1] for each segment
for I = 1:12
    idx = (T.SegN == I);
    v = T.fitness(idx);
    v(v < prctile(v,2.5)) = prctile(v,2.5) ;
    v = v + abs(min(v));
    v(v > prctile(v,98.5)) = prctile(v,98.5) ;
    v = v ./ max(v) ;
    T.fitness(idx) = v ;
end

%% find threshold for unfit strains -- 99% of w/stop are here or lower
nonsense_threshold_FAST = NaN(12,1);
nonsense_threshold_SLOW = NaN(12,1);
gapidx = regexpcmp( T.aa_seq ,'_');  % same as STOP & NONSENSE

for I = 1:12
    idx = ( T.SegN==I) ;
    nonsense_threshold_SLOW(I) = prctile( T.fitness( idx & gapidx) , 90);
    nonsense_threshold_FAST(I) = prctile( T.fitness( idx & gapidx) , 99.9);
    
    figure('units','centimeters','position',[5 5 7 7]);
    hold on ;
    rectangle( 'Position', [nonsense_threshold_FAST(I) 0 , 1 , 1] , 'FaceColor',[0.9 1 .9],'EdgeColor',[.9 1 .9])
    rectangle( 'Position', [0 0 ,  nonsense_threshold_SLOW(I) , 1] , 'FaceColor',[1 .9 .9],'EdgeColor',[1 .9 .9])

    [f,x]=ecdf( T.fitness(idx & ~gapidx)); plot(x,f,'-k','LineWidth',3)
    [f,x]=ecdf( T.fitness(idx & gapidx));plot(x,f,'-r','LineWidth',3)
%   line( [ nonsense_threshold_SLOW(I) nonsense_threshold_SLOW(I)],ylim,'Color','g')
%    line( [ nonsense_threshold_FAST(I) nonsense_threshold_FAST(I)],ylim,'Color','b')
    xlabel('Fitness')
    ylabel('Fraction of variants')
    xlim([0 1])
    set(gca,'xtick',0:.2:1)
end
[   nonsense_threshold_SLOW nonsense_threshold_FAST]


%%
G = grpstats( T , {'SegN' 'aa_seq'} ,{ 'median' 'mean' 'std' } ,'DataVars','fitness');
A = grpstats( T , 'SegN' ,{ 'median' 'mean' 'max' 'min'} ,'DataVars','fitness');

G2 = G( G.GroupCount >= 5 ,:);
T2 = T( ismember(  T.aa_seq , G2.aa_seq ) ,:);

DS = outerjoin( T2 , G2 , 'Key','aa_seq' , 'MergeKeys' , true);

%% how often is the median fast (> nonsense_threshold_99) but
%   a syn var slow ( < nonsense_threshold_90) ?
% and vice - versa
DS.MedianIsFast = NaN(height(DS),1);
DS.MedianIsSlow = NaN(height(DS),1);
DS.SynVarIsFast = NaN(height(DS),1);
DS.SynVarIsSlow = NaN(height(DS),1);
data = NaN( 12 , 4);
for I = 1:12
    idx = ( DS.SegN_T2==I) ;
    DS.MedianIsFast(idx) = DS.median_fitness(idx) >= nonsense_threshold_FAST(I) ;
    DS.MedianIsSlow(idx) = DS.median_fitness(idx) <= nonsense_threshold_SLOW(I) ;
    DS.SynVarIsFast(idx) = DS.fitness(idx) >= nonsense_threshold_FAST(I) ;
    DS.SynVarIsSlow(idx) = DS.fitness(idx) <= nonsense_threshold_SLOW(I) ;

    data( I , 1) = mean(DS.MedianIsFast(idx) & ~DS.SynVarIsSlow(idx));
    data( I , 2) = mean(DS.MedianIsSlow(idx) & ~DS.SynVarIsFast(idx));
    data( I , 3) = mean(DS.MedianIsFast(idx) & DS.SynVarIsSlow(idx));
    data( I , 4) = mean(DS.MedianIsSlow(idx) & DS.SynVarIsFast(idx));

end
%%
figure('units','centimeters','position',[5 5 6 7]);
bh = bar( data(:,3)*100  )
xlabel('Segment')
title('AA fit & SynVar unfit' )
ylabel('% of NT seqs')
set(gca,'Ygrid','on')
set(gca,'yscale','log')
set(gca,'ytick',[1e-3 1e-2 1e-1 1 2 5 10]);
set(gca,'yticklabel',[0.001 0.01 0.1 1 2  5 10])
ylim([0.1 10])
figure('units','centimeters','position',[5 5 6 7]);
bh = bar( data(:,4)*100  )
xlabel('Segment')
title('AA unfit & SynVar fit' )
ylabel('% of NT seqs')
set(gca,'Ygrid','on')
set(gca,'yscale','log')
set(gca,'ytick',[1e-3 1e-2 1e-1 1 5 10]);
set(gca,'yticklabel',[0.001 0.01 0.1 1 5 10])
ylim([0.001 5])
%% bar plot of distance to contain 90%, 95% , 99% of variants
figure;
hold on;
X  = ([0.75 0.90 0.95] )' ;
Y = NaN( numel(X),12);

for SegN = 1:12
    idx = find(DS.SegN_T2 == SegN);
    if ~isempty(idx)
        data = abs(DS.median_fitness(idx)-DS.fitness(idx) ) ;
        [f,x]=ecdf(data);
        for I = 1:numel(X)
            Y(I,SegN) =  x(find( f>=X(I)  ,1, 'first') );
        end
     end
end
bh = bar(1:3,Y , 'LineWidth',0.01) ; 
set(gca,'xtick',1:3)
set(gca,'xticklabel',X)
set(gca,'xscale','lin')
set(gca,'yscale','lin')
axis tight;
set(gca,'Ygrid','on')
set(gca,'ytick',0:.1:1)
xlabel('% of synonymous variants')
ylabel('distance from median')
legend( arrayfun(@num2str,1:12,'UniformOutput',false)  , 'location','nw');


%% Plot a few examples
G3 = G( G.GroupCount >= 25 ,:);
G3 = sortrows(G3,'GroupCount','descend');

G3 = sortrows(G3,'std_fitness','ascend');
G3 = G3( G3.GroupCount >= 150 ,:);

G3 = G( G.GroupCount >= 50 & G.median_fitness < 0.3 ,:);
G3 = sortrows(G3,'GroupCount','descend');
G3 = sortrows(G3,'std_fitness','ascend');

%%
data =  T.fitness( strcmp(T.aa_seq , G3.aa_seq{1}))  ; 
fh = figure('units','centimeters','position',[5 5 8 7]);
hold on; 
title(G3.aa_seq{1})

h = histogram(data ,25 ,'Normalization','count','FaceColor','r','EdgeColor','r')
rectangle( 'Position' , [median(data)-0.1 , 0 , 0.2 , max(h.Values)] ,'FaceColor',[.97 .97 .97],'EdgeColor',[.9  .9  .9 ])

h = histogram(data ,25 ,'Normalization','count','FaceColor','r','EdgeColor','r')

[f,x]=ecdf(data ) ;
plot(x,f.*max(h.Values),'-k','LineWidth',2)

grid on ;
line([median(data) median(data)] , [ 0  max(h.Values) ])
axis tight;
xlabel('Fitness of synonymous variants')
ylabel('# of synonymous variants')
set(gca,'ytick',0:5:100)
set(gca,'xtick',0:.2:1)
xlim([0 1])