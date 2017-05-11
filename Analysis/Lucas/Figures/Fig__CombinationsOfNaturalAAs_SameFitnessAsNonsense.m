%%   [F9]Figure with the distribution of genotype fitnesses
% We found that in many His3 segments a substantial fraction of combinations of natural amino acids
%   lead to genotypes with low fitness, many indistinguishable from lethal frameshift mutations (Figure[F9] ).
%
%
% In X out of 12 His3 segments a substantial fraction of combinations of extant amino acid states led to genotypes with low fitness, including that indistinguishable from lethal frameshift mutations (Figure 3).
%% In X out of 12 His3 segments a substantial fraction of combinations of extant
% amino acid states led to genotypes with low fitness, including that
% indistinguishable from lethal frameshift mutations (Figure 3).


all_nonsense_fitness = cell(12,1);
all_natlib_fitness = cell(12,1);
all_other_fitness = cell(12,1);
all_middlevar_fitness = cell(12,1);
all_other_pos_var_fitness = cell(12,1);


nonsense_fitness_thresholds = [95 99 99.5 99.75 99.9 100] ;
DS = dataset();
DS.SegN = NaN(100,1);
DS.nonsense_fitness_threshold_pct = NaN(100,1);
DS.nonsense_fitness_threshold_val = NaN(100,1);
DS.PctNatLib_ThisBad = NaN(100,1) ;
DS.PctOther_ThisBad = NaN(100,1) ;
DS.PctMiddleVar_ThisBad = NaN(100,1) ;

threshold_for_figures = 99 ;

c = 0;
for I = 1:12
    T = readtable( [ '~/Develop/HIS3InterspeciesEpistasis/Data/S' num2str(I) '_scaled_info_v2.csv'] ,'FileType','text','Delimiter','\t');
    nonsense_fitness = T.s( T.stop  & T.nonsense) ;
    natlib_fitness = T.s( ~T.stop  & ~T.nonsense & T.middle & T.nogap & T.nat_lib) ;
    other_fitness = T.s( ~T.stop  & ~T.nonsense & T.middle & T.nogap & ~T.nat_lib & T.lib ) ;
    other_pos_var_fitness = T.s( T.nogap & ~T.nonsense & ~T.lib & ~T.stop & ~T.nat ) ;
    middlevar_fitness = T.s( ~T.stop  & ~T.nonsense & ~T.middle & T.nogap & ~T.nat_lib & T.lib ) ;
    for nfti = 1:numel(nonsense_fitness_thresholds)
        c = c+1;
        DS.SegN(c) = I ;
        DS.nonsense_fitness_threshold_pct(c) = nonsense_fitness_thresholds(nfti);
        DS.nonsense_fitness_threshold_val(c) = prctile( nonsense_fitness ,  nonsense_fitness_thresholds(nfti)) ;
        DS.PctNatLib_ThisBad(c) = 100 * mean( natlib_fitness < DS.nonsense_fitness_threshold_val(c));
        DS.PctOther_ThisBad(c) = 100 * mean( other_fitness < DS.nonsense_fitness_threshold_val(c));
        DS.PctMiddleVar_ThisBad(c) = 100 * mean( middlevar_fitness < DS.nonsense_fitness_threshold_val(c));
        
        if nonsense_fitness_thresholds(nfti) == threshold_for_figures
            all_nonsense_fitness{I} = nonsense_fitness ; 
            all_natlib_fitness{I} = natlib_fitness ;
            all_other_fitness{I}  = other_fitness ;
            all_middlevar_fitness{I} = middlevar_fitness ;
            all_other_pos_var_fitness{I} = other_pos_var_fitness ; 
        end
        
    end
end
DS = DS( ~isnan(DS.SegN) ,:);


%%
basefigname = '~/Dropbox/Pokusaeva17/Figures/NatLibGivesUnfitCombinations__FitnessDistributionsPerSegment';


clrs = colormap();
nat_lib_color = 'r' ;
other_color   = clrs ;
%% ECDFS
for SegN = 1:12
    thresh_val = DS.nonsense_fitness_threshold_val( DS.SegN==SegN & DS.nonsense_fitness_threshold_pct==threshold_for_figures) ;
    figure('units','centimeters','position',[5 5 6 6]);
    hold on ;
    
    line([ thresh_val thresh_val] , [ 0 2] ,'Color',[.7 .7 .7],'LineWidth',1,'LineStyle','-.')
 
     
    if ~isempty(all_other_pos_var_fitness{SegN})
    [f,x]=ecdf( all_other_pos_var_fitness{SegN});
    ph_other = plot(x,f,'-','LineWidth',3,'Color','g') ;
    end
    
    [f,x]=ecdf( all_natlib_fitness{SegN});
    ph_natlib = plot(x,f,'-','LineWidth',3,'Color','r') ;


    if ~isempty(all_other_fitness{SegN})
    [f,x]=ecdf( all_other_fitness{SegN});
    ph_other = plot(x,f,'-','LineWidth',3,'Color','b') ;
    end    
    
    %     [f,x]=ecdf( all_middlevar_fitness{I});
    %     ph_midvar = plot(x,f,'-','LineWidth',3)    ;
    
    [f,x]=ecdf( all_nonsense_fitness{SegN});
    plot(x,f,'-k','LineWidth',3)
    
    ylabel('Fraction of variants')
    xlabel('Fitness')
    xlim([0 1.3])
    ylim([0 1])
    set(gca,'xtick',0:.2:10)
    set(gca,'ytick',0:.2:1)
    grid on;
    title( sprintf( 'segment %d' , SegN))
    print('-dpng',[basefigname '_ecdfs_' num2str(SegN) '.png']);
    close all  ;
end
%% histograms
delete([basefigname '_hists.eps'])

for I = 1:12
    thresh_val = DS.nonsense_fitness_threshold_val( DS.SegN==I & DS.nonsense_fitness_threshold_pct==threshold_for_figures) ;
    figure('units','centimeters','position',[5 5 7 7]);
    hold on ;
   
   h1 = histogram( all_natlib_fitness{I} ,0:.05:2 ,'FaceColor',nat_lib_color,'FaceAlpha',1,'EdgeColor','k');
   h0 = histogram( all_nonsense_fitness{I}  ,0:.05:2 ,'FaceColor','k','FaceAlpha',0.4);
   line([ thresh_val thresh_val] , ylim ,'Color',[.7 .7 .7],'LineWidth',1,'LineStyle','-.')
    
    
    xlabel('Fitness')
    ylabel('# of variants (\times1000)')
    xlim([0 1.3])
    set(gca,'xtick',0:.25:10)
    set(gca,'ytick',[0 2e3 5e3 1e4 2e4 3e4 5e4 6e4 7e4 8e4 9e4 1e5 ])
    set(gca,'yticklabel',[0 2 5 10:10:100])

    [f,x]=ecdf( all_nonsense_fitness{I});
    plot(x , f*max(ylim) ,'-k');

    
    [f,x]=ecdf( all_natlib_fitness{I});
    plot(x , f*max(ylim) ,'-','Color',nat_lib_color);
    
    title( sprintf( 'segment %d' , I))    
    
    
    print('-dpsc2',[basefigname '_hists.eps'],'-append');
    close all  ;
end


%%
 
figname = [ basefigname '_Summary_FractionUnfit_per_segment.eps'];
delete(figname);


figure('units','centimeters','position',[5 5 7 7]);
data = [  DS.PctNatLib_ThisBad( DS.nonsense_fitness_threshold_pct == 99) DS.PctOther_ThisBad( DS.nonsense_fitness_threshold_pct == 99) ] ;
bh = bar(data) ;
set(bh(1),'FaceColor',nat_lib_color,'EdgeColor',nat_lib_color);
set(bh(2),'FaceColor',other_color,'EdgeColor',other_color);
xlim([.5 12.5])
xlabel('Segment')
ylabel('% of variants that are unfit')
set(gca,'ytick',0:10:100)
ylim([0 100]) 
print('-dpsc2',figname,'-append');
close all  

figure('units','centimeters','position',[5 5 7 7]);
bar( DS.PctNatLib_ThisBad( DS.nonsense_fitness_threshold_pct == 99) ,'FaceColor', nat_lib_color )
axis tight;
xlabel('Segment')
ylabel('% of extant variants that are unfit')
set(gca,'ytick',0:10:100)
ylim([0 100])
print('-dpsc2',figname,'-append');
close all  ;


figure('units','centimeters','position',[5 5 7 7]);
bar( DS.PctOther_ThisBad( DS.nonsense_fitness_threshold_pct == 99)  ,'FaceColor',other_color )
axis tight;
xlabel('Segment')
ylabel('% of non-extant variants that are unfit')
set(gca,'ytick',0:10:100)
ylim([0 100])
print('-dpsc2',figname,'-append');
close all  

figure('units','centimeters','position',[5 5 7 7]);
hold on;
bh = bar(  DS.PctOther_ThisBad( DS.nonsense_fitness_threshold_pct == 99)  ,'FaceColor',other_color,'EdgeColor',other_color) ;
bh = bar(  DS.PctNatLib_ThisBad( DS.nonsense_fitness_threshold_pct == 99)  ,'FaceColor',nat_lib_color,'EdgeColor',nat_lib_color) ;
xlim([.5 12.5])
xlabel('Segment')
ylabel('% of variants that are unfit')
set(gca,'ytick',0:10:100)
set(gca,'xtick',1:12)
ylim([0 100]) 
print('-dpsc2',figname,'-append');
close all  

