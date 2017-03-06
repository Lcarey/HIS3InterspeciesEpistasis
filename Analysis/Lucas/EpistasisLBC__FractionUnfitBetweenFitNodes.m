%% load data
S12 = readtable('S12_exponent_s_gridfit_bybest10000_aa_avg_with_dist.csv');
S12.fitness = S12.s ;
%%


%%
Y = S12.s ;
Y(Y<-2) = -2 ;
Y(Y > 1) = 1 ;

xl = -2.1:.05:1.1 ;

fh = figure('units','centimeters','position',[5 5 10 10]);



subplot(2,1,2)
hold on; 
ecdf( Y( logical(S12.stop)) );
ecdf( Y( ~logical(S12.stop)));
legend({'+stop' 'no stop'},'location','nw')
xlim([-2 1])
grid on ;
set(gca,'xtick',-2:.5:1)
set(gca,'ytick',0:.1:1)


subplot(2,1,1)
hold on; 
histogram( Y( logical(S12.stop)) , xl ,'Normalization','Probability');
histogram( Y( ~logical(S12.stop)) , xl,'Normalization','Probability');
legend({'+stop' 'no stop'},'location','nw')
axis tight;

%%


%%
T = S12( : , {'aa_seq' 'aa_seq_short' 'fitness'});
fit_cutoff = 0.45 ; 
fast_seqs =  T.aa_seq_short( T.fitness >  0.6 ) ;
slow_seqs =  T.aa_seq_short( T.fitness <  0 ) ;
figname = 'fraction_unfit_intermediate_states.eps';
delete(figname);
J = 1; 

FitnessDistributions = cell( numel(fast_seqs) , 1);
HammingDistances = NaN( numel(fast_seqs) , 1);
FractionIntermediateUnfit = NaN( numel(fast_seqs) , 1);
NumIntermediateStatesMeasured = NaN( numel(fast_seqs) , 1);
for I = 2:numel(fast_seqs)
    seq1 = fast_seqs{J};
    seq2 = fast_seqs{I};
    all_transition_states = ExpandSeqAlign( seq1 , seq2);
    idx = find( ismember( T.aa_seq_short , all_transition_states) );
    
    FractionIntermediateUnfit(I) = mean( T.fitness(idx)<fit_cutoff);
    HammingDistances(I) = HammingDistance( seq1 , seq2 ) ;
    FitnessDistributions{I} = T.fitness(idx) ;
    NumIntermediateStatesMeasured(I) = numel(idx) - 2  ;
%     if numel(idx)>5
%         fh = figure('units','centimeters','position',[5 5 5 5]);
%         hold on ;
%         ecdf(T.fitness)
%         ecdf(T.fitness(idx))
%         legend({'all' sprintf('%d %0.01f%%',numel(idx),fraction_intermediate_unfit*100)},'location','nw');
%         title_string = sprintf('%d x %d %s x %s' , I , J , seq1 , seq2);
%         title(title_string);
%         xlabel('Fitness')
%         ylabel('Fraction of seqs')
%         set(gcf,'PaperPosition',[0 0 5 5]);
%         print('-dpsc2',figname,'-append');
%         close;
%     end
    if (mod(I,100)==0) , fprintf('.') , end;
    if (mod(I,1000)==0) , fprintf('\n') , end;
end


%%
R = dataset();
R.FractionIntermediateUnfit = FractionIntermediateUnfit;
R.HammingDistances = HammingDistances;
R.FitnessDistributions = FitnessDistributions ;
R.NumIntermediateStatesMeasured = NumIntermediateStatesMeasured ;
R = R( ~isnan(R.FractionIntermediateUnfit),:);

%%

figure; 
boxplot(R.FractionIntermediateUnfit*100,R.HammingDistances,'notch','on')
%xlim([3.5 10.5])
xlabel('Distance between pairs')
ylabel('% unfit (s < 0.45)')
set(gca,'ytick',0:5:50);
grid on; 
ylim([-1 25])
% 
% figure; 
% boxplot(cellfun( @(X)mean(X<0) , R.FitnessDistributions)*100,R.HammingDistances,'notch','on')
% xlim([3.5 10.5])
% xlabel('Distance between pairs')
% ylabel('% unfit (s < 0)')
% set(gca,'ytick',0:5:50);
% grid on; 
% ylim([-1 15])

%%
idx = find( ismember( T.aa_seq , all_transition_states) );
histogram(T.fitness(idx))
