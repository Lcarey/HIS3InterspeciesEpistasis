%  EpistasisLBC_Figure_HowManyDimensionsToDisplayMutationalLandscape_PCA.m 
%% how many dimensions are necessary to display mutational landscape of each segment? 
% calculated using linear PCA (no free params, and exact)
% hamming distance between mutations
%
% LBC March 2017
s=struct();
for I = 1:12
    tic;
    % natlib
    T  = EpistasisLBC__LoadData( I  ,'ONLY_NATLIB_FLAG',true, 'ONLY_MIDDLE_FLAG',true,'NO_STOP_FLAG',true);
    
    %everything. basically the same result
    % T  = EpistasisLBC__LoadData( I  ,'ONLY_NATLIB_FLAG',false, 'ONLY_MIDDLE_FLAG',false,'NO_STOP_FLAG',true);

    T  = EpistasisLBC__FindVariableRegion_GenSparseMatrix( T ) ; 
    s(I).T = T;
    s(I).Segment = I ; 
    [~,~,~,~, s(I).explained,~] = pca( double(cell2mat(s(I).T.SparseVect)),'Centered',false);
    fprintf('%d took %0.01f min\n' , I , toc/60); tic;
end


 
%%
% s = sortStruct(s,'Segment');


fewest_dim = min(cellfun(@length,{s.explained})) ; 
rmat = NaN(  fewest_dim+1, numel(s));

for I = 1:numel(s)
    s(I).thresh99 = find(cumsum(s(I).explained)>99,1,'first') ; 
end

for I = 1:numel(s)
    rmat(:,I) = vertcat( 0 , cumsum(s(I).explained(1:fewest_dim)));
end
fh = figure('units','centimeters','position',[5 5 8 6]);
hold on; 
imagesc(rmat);
colorbar
for I = 1:numel(s)
    line([I-0.5 I+0.5] , [ s(I).thresh99 s(I).thresh99] ,'Color','k')
end
ylabel('# of principal components')
xlabel('Segment')
axis tight;

%%
% % % % 
% % % % 
% % % % %%
% % % % for I = 1:12
% % % %     s(I).thresh95_2 = find(cumsum(s(I).explained2)>95,1,'first') ;
% % % %     s(I).thresh95 = find(cumsum(s(I).explained)>95,1,'first') ; 
% % % %     s(I).thresh90_2 = find(cumsum(s(I).explained2)>90,1,'first') ;
% % % %     s(I).thresh90 = find(cumsum(s(I).explained)>90,1,'first') ; 
% % % %     s(I).thresh99_2 = find(cumsum(s(I).explained2)>99,1,'first') ;
% % % %     s(I).thresh99 = find(cumsum(s(I).explained)>99,1,'first') ; 
% % % %     s(I).thresh100_2 = find(cumsum(s(I).explained2)>=100,1,'first') ;
% % % %     s(I).thresh100 = find(cumsum(s(I).explained)>=100,1,'first') ; 
% % % % end
% % % % figure; 
% % % % subplot(2,2,1) ; plot([s.thresh95_2],[s.thresh95],'ok');title('95')
% % % % subplot(2,2,3) ; plot([s.thresh90_2],[s.thresh90],'ok');title('90')
% % % % subplot(2,2,2) ; plot([s.thresh100_2],[s.thresh100],'ok');title('100')
% % % % subplot(2,2,4) ; plot([s.thresh99_2],[s.thresh90_2],'ok');title('99 x 90')
% % % % 
% % % % s = sortStruct(s,'thresh95');
% % % % 
% % % % %%
% % % % fh = figure('units','centimeters','position',[5 5 10 6]);
% % % % hold on; 
% % % % clrs = parula(numel(s)+2);
% % % % for I = 1:numel(s)
% % % %     txt = sprintf('S%d' , s(I).Segment ) ;
% % % %     plot( vertcat( 0 , cumsum(s(I).explained) )  ,'-' ,'DisplayName',txt,'Color',clrs(I,:));
% % % % end
% % % % xlim([0.5 25])
% % % % ylim([0 100])
% % % % legend('location','EastOutside');
% % % % xlabel('# PCs')
% % % % ylabel('% variance explained')
% % % % grid on;


