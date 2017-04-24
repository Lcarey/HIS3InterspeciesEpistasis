%%   [F9]Figure with the distribution of genotype fitnesses
% We found that in many His3 segments a substantial fraction of combinations of natural amino acids
%   lead to genotypes with low fitness, many indistinguishable from lethal frameshift mutations (Figure[F9] ). 

all_nonsense_fitness = cell(12,1);
all_natlib_fitness = cell(12,1);
all_other_fitness = cell(12,1);

threshold_for_figures = 99 ; 
nonsense_fitness_thresholds = [99 99.5 99.75 99.9 100] ; 
DS = dataset();
DS.SegN = NaN(100,1);
DS.nonsense_fitness_threshold_pct = NaN(100,1);
DS.nonsense_fitness_threshold_val = NaN(100,1);
DS.PctNatLib_ThisBad = NaN(100,1) ;
DS.PctOther_ThisBad = NaN(100,1) ;

c = 0;
for I = 1:12
    T = readtable( [ '~/Develop/HIS3InterspeciesEpistasis/Data/' num2str(I) '_all.tab'] ,'FileType','text');
    nonsense_fitness = T.s( T.stop  & T.nonsense) ;
    natlib_fitness = T.s( ~T.stop  & ~T.nonsense & T.middle & T.nogap & T.nat_lib) ;
    other_fitness = T.s( ~T.stop  & ~T.nonsense & T.middle & T.nogap & ~T.nat_lib) ;

    for nfti = 1:numel(nonsense_fitness_thresholds)
        c = c+1;
        DS.SegN(c) = I ; 
        DS.nonsense_fitness_threshold_pct(c) = nonsense_fitness_thresholds(nfti);
        DS.nonsense_fitness_threshold_val(c) = prctile( nonsense_fitness ,  nonsense_fitness_thresholds(nfti)) ; 
        DS.PctNatLib_ThisBad(c) = 100 * mean( natlib_fitness < DS.nonsense_fitness_threshold_val(c));
        DS.PctOther_ThisBad(c) = 100 * mean( other_fitness < DS.nonsense_fitness_threshold_val(c));
        if nonsense_fitness_thresholds(nfti) == threshold_for_figures
            all_nonsense_fitness{I} = nonsense_fitness ;
            all_natlib_fitness{I} = natlib_fitness ;
            all_other_fitness{I}  = other_fitness ; 
        end

    end
end
DS = DS( ~isnan(DS.SegN) ,:);

%%
for I = 1:12
    thresh_val = DS.nonsense_fitness_threshold_val( DS.SegN==I & DS.nonsense_fitness_threshold_pct==threshold_for_figures) ; 

    figure('units','centimeters','position',[5 5 7 7]);
    hold on ;
    line([ thresh_val thresh_val] , [ 0 1] ,'Color','c','LineWidth',2,'LineStyle','--')
    [f,x]=ecdf( all_natlib_fitness{I});
    plot(x,f,'-r','LineWidth',3)
    [f,x]=ecdf( all_nonsense_fitness{I});
    plot(x,f,'-k','LineWidth',3)
    [f,x]=ecdf( all_other_fitness{I});
    plot(x,f,'-b','LineWidth',3)
    xlabel('Fitness')
    ylabel('Fraction of variants')
    xlim([0 1])
    set(gca,'xtick',0:.2:1)
    grid on; 
    axis tight;
    title( sprintf( 'segment %d' , I))
end

        
%%
figure('units','centimeters','position',[5 5 7 7]);
gscatter(  DS.nonsense_fitness_threshold_pct , DS.PctNatLib_ThisBad , DS.SegN);
set(gca,'xtick',nonsense_fitness_thresholds)
ylabel('% of variants w/fit < threshold')
xlabel('threshold for unfit (as % of nonsense mutants)')
%%
figure('units','centimeters','position',[5 5 7 7]);
bar( DS.PctNatLib_ThisBad( DS.nonsense_fitness_threshold_pct == 99) ,'FaceColor',[0.7 0 0])
axis tight;
xlabel('Segment')
ylabel('% of natural variants that are unfit')
set(gca,'ytick',0:10:100)
ylim([0 100])

figure('units','centimeters','position',[5 5 7 7]);
bar( DS.PctOther_ThisBad( DS.nonsense_fitness_threshold_pct == 99)  ,'FaceColor',[0 0 0.5])
axis tight;
xlabel('Segment')
ylabel('% of non-nat variants that are unfit')
set(gca,'ytick',0:10:100)
ylim([0 100])

%title('99% threshold')