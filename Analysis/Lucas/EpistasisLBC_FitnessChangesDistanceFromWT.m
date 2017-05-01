%% load data
%WT = readtable('~/Develop/HIS3InterspeciesEpistasis/Data_Small_Tables/wt_seq.csv','Delimiter','\t')
s = struct();
for SegN = 1:12
    fn = ['~/Develop/HIS3InterspeciesEpistasis/Data/S' num2str(SegN) '_scaled_info.csv'] ;
    T = readtable(fn,'FileType','text');
    T = T( T.middle & T.nogap & T.nat  ,:);
    
    fit_cutoff = 0.28  ;
    T.isfit = T.s > fit_cutoff  ;
 
    s(SegN).T = T ;
    
end

%%
for SegN = 1:12
    T = s(SegN).T ; 
       vn = T.Properties.VariableNames( regexpcmp( T.Properties.VariableNames , '^dist_[A-Z]'));
    vn = vn( ~strcmp(vn,'dist_Scer'));
    figure; hold on;
    clrs = parula(numel(vn));
    for I = 1:numel(vn)
        G = grpstats( T , vn{I}, 'mean' ,'DataVars' , {'isfit'});
        mdl = fitglm( G.(vn{I}) ,  G.mean_isfit) ;
        txt = sprintf( '%s m=%0.02f'  ,  regexprep( vn{I} ,'dist_' ,'')  , 10*mdl.Coefficients.Estimate(2) );
        plot( G.(vn{I}) , 100*G.mean_isfit ,'-','Color',clrs(I,:) ,'DisplayName' ,txt,'LineWidth',2);
        
    end
    
    G = grpstats( T , 'dist_Scer', 'mean' ,'DataVars' , {'isfit'});
    mdl = fitglm( G.dist_Scer ,  G.mean_isfit) ;
    txt = sprintf( '%s m=%0.02f'  ,  'Scer'  , 10*mdl.Coefficients.Estimate(2) );
    plot( G.dist_Scer , 100*G.mean_isfit ,'-k' ,'DisplayName' ,txt,'LineWidth',2);
    
    legend('location','sw')
    ylabel('% fit genotypes')
    xlabel('Distance from this species')
    title(SegN)
    
    ylim([0 100]);
    
end

%%
SegN = 7 ; 
T = s(SegN).T ;
 vn = T.Properties.VariableNames( regexpcmp( T.Properties.VariableNames , '^dist_[A-Z]'));
T( any( table2array(T(:,vn))==0 , 2) ,:)

%% 
fitness = NaN(0);
SegN = NaN(0);
for I = 1:12
    idx = find( any( table2array(T(:,vn))==0 , 2));
    fitness = vertcat( fitness , T.s(idx) );
    SegN = vertcat( SegN , repmat( I , numel(idx) , 1));
end
R = dataset();
R.fitness = fitness ; 
R.SegN = SegN ; 
    

