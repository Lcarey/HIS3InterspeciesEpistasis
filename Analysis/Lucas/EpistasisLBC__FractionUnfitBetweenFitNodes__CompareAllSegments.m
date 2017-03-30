%% EpistasisLBC__FractionUnfitBetweenFitNodes__CompareAllSegments
% given a folder w/a single .mat file for each segment
% graph % unfit vs fraction unfit intermediate states
%
%

matfiles = dir( 'IntermediateStates__March22/S*scaled_info*.mat');

FAST=0.40 ;
SLOW=0.20 ;

DS = dataset;
DS.pct_unfit = NaN( numel(matfiles) , 1);
DS.max_unfit_int = NaN( numel(matfiles) , 1);
DS.pct_fit = NaN( numel(matfiles) , 1);
DS.pct_int = NaN( numel(matfiles) , 1);
DS.Segment = NaN( numel(matfiles) , 1);

for I = 1:numel(matfiles)
    load( [ matfiles(I).folder filesep matfiles(I).name ]);
    G = grpstats( R  , 'HammingDistances', 'mean', 'DataVars' , 'pct_unfit' );
    DS.pct_unfit(I) = mean( T.s < SLOW) * 100;
    DS.pct_fit(I) = mean( T.s > FAST) * 100;
    DS.pct_int(I) = mean( T.s > SLOW & T.s < FAST) * 100;
    DS.max_unfit_int(I) = G.mean_pct_unfit(end) ;
    
    metadata = regexp(  matfiles(I).name  ,'_' ,'split') ;
    DS.Segment(I) = str2double( metadata{1}(2:end));
end
DS = sortrows(DS,'Segment')
%%
figure; 

figure;
plot( DS.max_unfit_int  , DS.pct_unfit ,'.k' );
text( DS.max_unfit_int  , DS.pct_unfit , num2str(DS.Segment) ,'FontWeight','bold');

title( corr(DS.max_unfit_int  , DS.pct_unfit))
ylabel('% unfit in library')
xlabel('% unfit intermediate states @ max HD')
grid on; 
xlim([0 70])
ylim(xlim)
line(xlim,ylim)


%%
mdl = fitglm( DS ,'linear' , 'PredictorVars'  , {'pct_unfit' 'pct_int' 'pct_fit'} , 'ResponseVar' ,'max_unfit_int')