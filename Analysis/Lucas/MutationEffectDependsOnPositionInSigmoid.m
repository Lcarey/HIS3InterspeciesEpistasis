%% The effect of mutations is context dependent based on where the starting point is in the fitness potential sigmoid

DataDir = '~/Desktop/HIS3scratch/March20Fits/';
load([DataDir 'R_logistic_20-Mar-201716:11:07_4315_S7_scaled_info.csv_.mat'] );

fitness = R.fit_train{1}';
sumDDG  = 1-R.train_sumddG{1}';
figure;  hold on; 
dscatter(sumDDG , fitness);
grid on ;

%%

figure; hold on; 
boxplot( fitness , (round(sumDDG*20)./20) ,'notch','on');
set(gca,'xtick',1:100);set(gca,'xticklabel',bins);
xlabel('predicted fitness potential')
ylabel('Measured fitness')
%%


%%
bins = 0:0.05:1 ; 
[Y,E] = discretize( sumDDG , bins ) ;


aa_for_all_variants = double( R.aa_for_all_variants{1}( R.Train{1} , :) );
n_variants = size(aa_for_all_variants,1) ; 

all_fitnesses = NaN(0);
neighbor_ddG_bin = NaN(0);
this_ddG = NaN(0);
all_fitnesses_diffs = NaN(0);
all_this_fitness = NaN(0);
for I = 1:n_variants
    nearest_neighbors_idx = find( pdist2( aa_for_all_variants(I,:) , aa_for_all_variants,'Hamming')*numel(aa_for_all_variants(1,:))==1 );
    this_fitness = fitness(I);
    all_fitnesses_diffs = vertcat(all_fitnesses ,  fitness(nearest_neighbors_idx) - this_fitness );
    all_fitnesses = vertcat(all_fitnesses ,  fitness(nearest_neighbors_idx) );
    this_ddG = vertcat( this_ddG , repmat( sumDDG(I) , numel(nearest_neighbors_idx) , 1) );
    all_this_fitness = vertcat( all_this_fitness , repmat( this_fitness , numel(nearest_neighbors_idx) , 1) );

    neighbor_ddG_bin = vertcat( neighbor_ddG_bin , repmat( Y(I) , numel(nearest_neighbors_idx) , 1) );
end
DS = dataset();
DS.all_fitnesses = all_fitnesses;
DS.neighbor_ddG_bin = neighbor_ddG_bin ;
DS.all_fitnesses_diffs = all_fitnesses_diffs ; 
DS.this_ddG = this_ddG ; 
DS.all_this_fitness = all_this_fitness;
G = grpstats(DS , 'neighbor_ddG_bin' ,{'mean'  'std'}, 'DataVars' , 'all_fitnesses');

% figure; hold on; 
% bar( DS.neighbor_ddG_bin , DS.mean_all_fitnesses );
% errorbar( DS.neighbor_ddG_bin , DS.mean_all_fitnesses  , DS.std_all_fitnesses ,'ok')

figure; hold on; 
boxplot( all_fitnesses , neighbor_ddG_bin ,'notch','on');
set(gca,'xtick',1:100);set(gca,'xticklabel',bins);
xlabel('predicted fitness potential')
ylabel('change in fitness 1 edit away')
line(xlim,[0 0 ],'LineStyle','--','Color',[.7 .7 .7])



%%
idx = DS.this_ddG>0.35 & DS.this_ddG<0.45 ;
mdl = fitglm( [DS.this_ddG(idx) DS.all_this_fitness(idx) ] , DS.all_fitnesses_diffs(idx) ,'interactions','VarNames',{'ddG' 'fit' 'fittdiff'});
figure; hold on; 
dscatter( DS.this_ddG(idx) , mdl.Residuals.LinearPredictor(idx) );
plot(1-R.train_sumddG{1},R.pred_fit_train{1},'.k')
set(gca,'xtick',0:0.05:1)
grid on;
%%
figure; hold on ; 
ybins = unique(Y);
clrs = parula(numel(ybins));
for I = 1:numel(ybins)
    data = all_fitnesses( neighbor_ddG_bin==ybins(I));
    if ~isempty(data)
    [f,x]=ecdf(data);
    plot(x,f,'-','DisplayName',num2str(bins(I)),'LineWidth',3,'Color',clrs(I,:));
    end
end
legend('location','best')
xlabel('change in fitness 1 edit away')
ylabel('Fraction of variants')