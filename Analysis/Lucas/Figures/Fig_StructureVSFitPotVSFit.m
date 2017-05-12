% Indeed, across all segments, predicted protein stability correlates better with fitness potential than with fitness (Supplementary Figure 10). 
% compare the correlation of predicted ddG
%   vs 
%  fitness & fitness potential
% May 11, 2017

figure; 
for SegN = 1:12
    
ddG = readtable('~/Develop/HIS3InterspeciesEpistasis/Analysis/Sasha/rosetta_runs/run-170503-results.csv','FileType','text','Delimiter','\t');
ddG.ddG = median( table2array(ddG(:,3:5)) , 2) ; 
FP = readtable([ '~/Develop/HIS3InterspeciesEpistasis/Analysis/Katya/NN/residuals/S' num2str(SegN) '.csv'],'FileType','text','Delimiter',',');
T  = EpistasisLBC__LoadData( SegN  ,'ONLY_NATLIB_FLAG',true, 'ONLY_MIDDLE_FLAG',true,'NO_STOP_FLAG',true);
T = join( T(:,{'aa_seq' 's'}) , FP ,'Key','aa_seq');
T = innerjoin( T  ,ddG(:,[1 7]) ,'Key','aa_seq');
T.Fit = T.s;
T.fitPoten = T.fitnessPotential ;
vn = {'Fit' 'fitPoten' 'ddG' } ;

subplot(4,3,SegN)
%dscatter( T.ddG , T.fitPoten );
%title( sprintf( '%d corr=%0.02f' , SegN , corr( T.fitPoten , T.ddG )) )
%dscatter( T.ddG , T.fitPoten );
%title( sprintf( '%d corr=%0.02f' , SegN , corr( T.fitPoten , T.ddG )) )
dscatter( T.ddG , T.Fit );
title( sprintf( '%d corr=%0.02f' , SegN , corr( T.Fit , T.ddG )) )

axis tight;
%fh = figure('units','centimeters','position',[5 5 10 10]);
%cph = corrplot( table2array(T(:,vn)) ,'VarNames',vn)
%print('-dpsc2',figname,'-append');
%close;

end