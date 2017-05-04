%%   [F9]Figure with the distribution of genotype fitnesses
% We found that in many His3 segments a substantial fraction of combinations of natural amino acids
%   lead to genotypes with low fitness, many indistinguishable from lethal frameshift mutations (Figure[F9] ).
%
%
%% In X out of 12 His3 segments a substantial fraction of combinations of extant
% amino acid states led to genotypes with low fitness, including that
% indistinguishable from lethal frameshift mutations (Figure 3).
all_nonsense_fitness = cell(12,1);
all_natlib_fitness = cell(12,1);
all_other_fitness = cell(12,1);
all_middlevar_fitness = cell(12,1);

threshold_for_figures = 99 ; % as a % of nonsense?  or absolute thresh


nonsense_fitness_thresholds = [99 99.5 99.75 99.9 100] ;
DS = dataset();
DS.SegN = NaN(100,1);
DS.nonsense_fitness_threshold_pct = NaN(100,1);
DS.nonsense_fitness_threshold_val = NaN(100,1);
DS.PctNatLib_ThisBad = NaN(100,1) ;
DS.PctOther_ThisBad = NaN(100,1) ;
DS.PctMiddleVar_ThisBad = NaN(100,1) ;

c = 0;
for I = 1:12
    T = readtable( [ '~/Develop/HIS3InterspeciesEpistasis/Data/' num2str(I) '_all.tab'] ,'FileType','text');
    nonsense_fitness = T.s( T.stop  & T.nonsense) ;
    natlib_fitness = T.s( ~T.stop  & ~T.nonsense & T.middle & T.nogap & T.nat_lib) ;
    other_fitness = T.s( ~T.stop  & ~T.nonsense & T.middle & T.nogap & ~T.nat_lib & T.lib ) ;
    middlevar_fitness = T.s( ~T.stop  & ~T.nonsense & ~T.middle & T.nogap & ~T.nat_lib) ;

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
        end
        
    end
end
DS = DS( ~isnan(DS.SegN) ,:);

%%
figname = 'Distributions_per_segment_cdf.eps';
delete(figname);
for I = 1:12
    thresh_val = DS.nonsense_fitness_threshold_val( DS.SegN==I & DS.nonsense_fitness_threshold_pct==threshold_for_figures) ;
    figure('units','centimeters','position',[5 5 7 7]);
    hold on ;
    
    line([ thresh_val thresh_val] , [ 0 2] ,'Color',[.7 .7 .7],'LineWidth',1,'LineStyle','-.')
    
    [f,x]=ecdf( all_natlib_fitness{I});
    ph_natlib = plot(x,f,'-','LineWidth',3) ;
    
    if ~isempty(all_other_fitness{I})
    [f,x]=ecdf( all_other_fitness{I});
    ph_other = plot(x,f,'-','LineWidth',3) ;
    end
    
    %     [f,x]=ecdf( all_middlevar_fitness{I});
    %     ph_midvar = plot(x,f,'-','LineWidth',3)    ;
    
    [f,x]=ecdf( all_nonsense_fitness{I});
    plot(x,f,'-k','LineWidth',3)
    
    xlabel('Growth rate (hr^{-1})')
    ylabel('Fraction of variants')
    xlim([0 0.5])
    set(gca,'xtick',0:.1:1)
    ylim([0 1])
    set(gca,'ytick',0:.1:1)
    grid on;
    title( sprintf( 'segment %d' , I))
    print('-dpsc2',figname,'-append');
    close all  ;
end

figname = 'Distributions_per_segment_pdf.eps';
delete(figname);
for I = 1:12
    thresh_val = DS.nonsense_fitness_threshold_val( DS.SegN==I & DS.nonsense_fitness_threshold_pct==threshold_for_figures) ;
    figure('units','centimeters','position',[5 5 7 7]);
    hold on ;
    
    
    [f,x]=ksdensity( all_natlib_fitness{I});
    ph_natlib = plot(x,f,'-','LineWidth',3) ;
    
    if ~isempty(all_other_fitness{I})
    [f,x]=ksdensity( all_other_fitness{I});
    ph_other = plot(x,f,'-','LineWidth',3) ;
    end
    
    %     [f,x]=ksdensity( all_middlevar_fitness{I});
    %     ph_midvar = plot(x,f,'-','LineWidth',3)    ;
    %
    [f,x]=ksdensity( all_nonsense_fitness{I});
    plot(x,f,'-k','LineWidth',3)
    
    axis tight; 
    line([ thresh_val thresh_val] , ylim ,'Color',[.7 .7 .7],'LineWidth',1,'LineStyle','-.')
    
    
    xlabel('Growth rate (hr^{-1})')
    ylabel('% of variants')
    xlim([0 0.5])
    set(gca,'xtick',0:.1:1)
    grid on;
    title( sprintf( 'segment %d' , I))
    
    try
    nat_lib_color = ph_natlib.Color ;
    other_color   = ph_other.Color ; 
    catch
    end
    
    print('-dpsc2',figname,'-append');
    
 
    close all  ;
end


% %%
% figure('units','centimeters','position',[5 5 7 7]);
% gscatter(  DS.nonsense_fitness_threshold_pct , DS.PctNatLib_ThisBad , DS.SegN);
% set(gca,'xtick',nonsense_fitness_thresholds)
% ylabel('% of variants w/fit < threshold')
% xlabel('threshold for unfit (as % of nonsense mutants)')
%%

figname = 'FractionUnfit_per_segment.eps';
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
%%
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

