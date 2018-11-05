ER = [0.0009	0.0006	0.0006	0.0009	0.0007	0.0004	0.0004	0.0004	0.0005	0.0005	0.0003	0.0005] ; 
R2biorep = [0.96	0.95	0.87	0.94	0.86	0.96	0.97	0.96	0.58	0.90	0.95	0.95] ; 
FalseFit = [0.02	0.11	0.4	0.11	0.35	0.13	0.03	0.18	0.8	0.01	0.22	0.03] ; 
FalseUnfit = [1	3.7	1.8	0.1	1.8	2.2	0.5	0.5	17.6	2	3.6	1.4]  ;
R2model = [.77 .70 .80 .92 .44 .81 .98 .91 .02 .01 .87 .63] ; 

mat = vertcat( ER*100 , R2biorep , FalseFit , FalseUnfit , R2model);
idx = [1:8 10:12] ; 

[c,p] = corr(mat(:,idx)') ; 
lbls = {'R^2 bio reps' 'FDR FalseFit' 'FDR FalseUnfit' 'R^2 model'} ;
%%
fh = figure('units','centimeters','position',[5 5 12 12]);
for I = 2:5
    subplot(2,2,I-1)
    hold on ;
   % plot( mat(1,:) , mat(I,:),'ok');
    plot( mat(1,idx) , mat(I,idx),'ok','MarkerFaceColor',[.7 .7 .7]);
    txt = sprintf( 'corr = %0.02f p=%0.02f' , c(I,1) , p(I,1));
    title(txt)
    xlabel('Sequencing error rate (%)')
    ylabel(lbls{I-1})
end

print('-dpng','ErrorRateIsNotPredictiveOfAnything.png','-r300');