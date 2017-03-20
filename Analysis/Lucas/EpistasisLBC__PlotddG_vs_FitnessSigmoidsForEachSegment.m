%% EpistasisLBC__PlotddG_vs_FitnessSigmoidsForEachSegment.m
% given a set of .mat files generated by ModelDecayFunctionFromRealData()
% produce a dscatter() plot of sum(ddG) ( called fitness potential ) vs measured fitness
%
% visual representation of how well the model does
%
% LBC

global n_possible_aas  ;
n_possible_aas = 27 ;

ModelDecayFunctionFromRealData_output_filenames = dir('./27/R_logistic*-27*mat')

fit_threshold = 0.4 ;
unfit_threshold = 0.2 ;

for I = 1:numel(ModelDecayFunctionFromRealData_output_filenames)
    title_split_str = regexp( ModelDecayFunctionFromRealData_output_filenames(I).name ,'_' ,'split');
    ModelDecayFunctionFromRealData_output_filenames(I).segment_name =  title_split_str{ regexpcmp( title_split_str , '^S\d')} ;
    ModelDecayFunctionFromRealData_output_filenames(I).segment_N = str2double(regexprep(ModelDecayFunctionFromRealData_output_filenames(I).segment_name,'S',''));
end
ModelDecayFunctionFromRealData_output_filenames = sortStruct(ModelDecayFunctionFromRealData_output_filenames,'segment_N',1);


%% one figure for each file. one subplot for each x-validation
for I = 1:numel(ModelDecayFunctionFromRealData_output_filenames)
    load( [ ModelDecayFunctionFromRealData_output_filenames(I).folder filesep ModelDecayFunctionFromRealData_output_filenames(I).name ] );
    R = R( ~isnan(R.L) , :) ;
    figure;
    for J = 1:length(R.L)
        subplot(4,5,J)
        hold on;
        test_idx = R.Test{J};
        x0_k_L_ddGvect = [ R.x0(J) R.k(J) R.L(J) R.fit_ddGvect{J}] ;
        [ predicted_fitness , sum_delta_Gs ] = LogisticFitnessDecayFunctionForOpt( x0_k_L_ddGvect , double(R.aa_for_all_variants{1})  ) ;
        
        dscatter( 1-sum_delta_Gs(test_idx)' , R.fit_test{J}' )
        plot( 1-sum_delta_Gs , predicted_fitness , '.k');
        set(gca,'ytick',0:0.1:1);
        set(gca,'xtick',0:.1:0.5);
        xlim([0 0.5]);
        ylim([0 0.5])
        text(0,0.05,  sprintf('R^2=%0.02f' , R.test_r2(J)) ,'Color','r' )
        grid on ;
        
        subplot(4,5,J+10)
        dscatter( R.pred_fit_test{J}' , R.fit_test{J}' )
        xlim([0 0.5]);ylim(xlim);
        line(xlim,xlim,'Color','r');
        set(gca,'xtick',0:.1:0.5);
        set(gca,'ytick',0:0.1:1);
        grid on ;
        
        miss_classification  = sum( (R.fit_test{J} >= fit_threshold) &  (R.pred_fit_test{J} <= unfit_threshold) | (R.fit_test{J} <= unfit_threshold) &  (R.pred_fit_test{J} >= fit_threshold)) / numel(R.fit_test{J})   ;
        miss_classification_fit = sum( (R.fit_test{J} >= fit_threshold) &  (R.pred_fit_test{J} <= unfit_threshold)) /  sum(R.fit_test{J} >= fit_threshold) ;
        miss_classification_unfit = sum( (R.fit_test{J} <= unfit_threshold) &  (R.pred_fit_test{J} >= fit_threshold)) /  sum(R.fit_test{J} <= unfit_threshold) ;
        
        title( sprintf( '%0.0f%%,%0.0f%%,%0.0f%%' , miss_classification*100,miss_classification_fit*100 , miss_classification_unfit*100))
    end
    title_split_str = regexp( ModelDecayFunctionFromRealData_output_filenames(I).name ,'_' ,'split');
    subtitle( [ title_split_str{ regexpcmp( title_split_str , '^S\d')} ' N=' title_split_str{ find(regexpcmp( title_split_str , '^S\d'))-1} ]   );
    print('-dpng' , [ ModelDecayFunctionFromRealData_output_filenames(I).name '.png' ] , '-r600');
    close ;
end


%% one figure. one panel for each segment
figure;
for I = 1:numel(ModelDecayFunctionFromRealData_output_filenames)
    load( [ ModelDecayFunctionFromRealData_output_filenames(I).folder filesep ModelDecayFunctionFromRealData_output_filenames(I).name ] );
    R = R( ~isnan(R.L) , :) ;
    
    subplot(2,12,I)
    if(I==1),ylabel('Measured Fitness'),end
    hold on;
    miss_classification = NaN(0);
    miss_classification_fit = NaN(0);
    miss_classification_unfit = NaN(0);
    predicted_fitness_test = NaN(0);
    measured_fitness_test = NaN(0);
    pred_ddG_test = NaN(0);
    r2_test = NaN(0);
    for J = 1:length(R.L)
        
        test_idx = R.Test{J};
        x0_k_L_ddGvect = [ R.x0(J) R.k(J) R.L(J) R.fit_ddGvect{J}] ;
        [ predicted_fitness , sum_delta_Gs ] = LogisticFitnessDecayFunctionForOpt( x0_k_L_ddGvect , double(R.aa_for_all_variants{1})  ) ;
        
        measured_fitness_test = vertcat( measured_fitness_test ,  R.fit_test{J}' );
        predicted_fitness_test = vertcat( predicted_fitness_test ,  R.pred_fit_test{J}' );
        pred_ddG_test = vertcat( pred_ddG_test ,  R.test_sumddG{J}' );
        r2_test = vertcat( r2_test ,  R.test_r2(J)' );
        x0 = vertcat(x0 ,  R.x0(J)' );
        miss_classification(end+1)  = sum( (R.fit_test{J} >= fit_threshold) &  (R.pred_fit_test{J} <= unfit_threshold) | (R.fit_test{J} <= unfit_threshold) &  (R.pred_fit_test{J} >= fit_threshold)) / numel(R.fit_test{J})   ;
        miss_classification_fit(end+1) = sum( (R.fit_test{J} >= fit_threshold) &  (R.pred_fit_test{J} <= unfit_threshold)) /  sum(R.fit_test{J} >= fit_threshold) ;
        miss_classification_unfit(end+1) = sum( (R.fit_test{J} <= unfit_threshold) &  (R.pred_fit_test{J} >= fit_threshold)) /  sum(R.fit_test{J} <= unfit_threshold) ;
        
        
    end
    
    %find midpoint
    [ ~ , o] = sort(measured_fitness_test,'ascend');
    midpt_idx = find( measured_fitness_test(o)>0.25 ,1 , 'first');
    midpt = 1-pred_ddG_test(o(midpt_idx));
    
    
    
    miss_classification   = median(miss_classification)*100 ;
    miss_classification_fit  = median(miss_classification_fit)*100 ;
    miss_classification_unfit = median(miss_classification_unfit)*100 ;
    
    dscatter( 1-pred_ddG_test , measured_fitness_test );
    plot( 1-sum_delta_Gs , predicted_fitness , '.k');
    title( sprintf('%s R^2=%0.02f' , ModelDecayFunctionFromRealData_output_filenames(I).segment_name ,median(r2_test)));
  
    set(gca,'ytick',0:0.1:1);
    set(gca,'xtick',0:.25:1);
    xlim([0 0.5]);
    ylim([0 0.5])
    grid on ;
    line([midpt midpt],ylim,'LineStyle','-','Color','r')
    
    subplot(2,12,I+12); hold on;
    if(I==1),ylabel('Measured Fitness'),end
    dscatter( predicted_fitness_test , measured_fitness_test )
    title( sprintf( '%0.0f%%,%0.0f%%,%0.0f%%' , miss_classification,miss_classification_fit , miss_classification_unfit))
    
    xlim([0 0.5]);ylim(xlim);
    set(gca,'xtick',0:.25:1);
    set(gca,'ytick',0:0.1:1);
    grid on;
    
end



