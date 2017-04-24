% As no single mutation is always lethal (Figure [LBCF1]), 
%  combinations of natural amino acids leading to low fitness can only be
%  explained by epistsis...



%% version 1

%% find the fraction of fit and unfit variants for each AAxPos 
s = struct();

fast_threshold = 0.4 ; 
slow_threshold = 0.2 ;
all_fraction_unfit = NaN(0);
all_fraction_fit = NaN(0);
all_N_unfit = NaN(0);
all_N_fit = NaN(0);
all_segments = NaN(0);
for SegN = 1:12
    fprintf('%d... ',SegN);tic; 
    T = EpistasisLBC__LoadData_FindVariableRegion_GenSparseMatrix( SegN , 'ONLY_NATLIB_FLAG',false , 'ONLY_LIB_FLAG',false , 'ONLY_MIDDLE_FLAG',true) ; 
    fprintf('%0.01f min ... ',toc/60); tic; 
    variant_mat = cell2mat( T.SparseVect );
    positions_that_vary = find(mean(variant_mat)>0) ; 
    variant_mat = variant_mat( : , positions_that_vary);

    s(SegN).fraction_fit  = NaN( 1 , size(variant_mat,2));
    s(SegN).fraction_unfit  = NaN( 1 , size(variant_mat,2));
    s(SegN).N_fit  = NaN( 1 , size(variant_mat,2));
    s(SegN).N_unfit  = NaN( 1 , size(variant_mat,2));

    for I = 1:size(variant_mat,2)
        idx = (variant_mat(:,I)==1);
        s(SegN).fraction_unfit(1,I) = mean(T.s(idx) < slow_threshold);
        s(SegN).fraction_fit(1,I) = mean(T.s(idx) > fast_threshold);
        s(SegN).N_unfit(I) = sum(T.s(idx) < slow_threshold);
        s(SegN).N_fit(I) = sum(T.s(idx) > fast_threshold);
    end
    
    all_fraction_unfit = [ all_fraction_unfit s(SegN).fraction_unfit];
    all_fraction_fit = [ all_fraction_fit s(SegN).fraction_fit];
    all_N_unfit = [ all_N_unfit s(SegN).N_unfit];
    all_N_fit = [ all_N_fit s(SegN).N_fit];
    all_segments = [ all_segments repmat(SegN,size(s(SegN).fraction_fit)) ];
    fprintf('%0.01f min . done!\n',toc/60); 

    s(SegN).T = T ; 
    s(SegN).SegN = SegN ; 
end


%% can remove some stuff
%                   aa_seq                    s       len    middle    nat_lib    nat    lib    size       aa_seq_variable        columns_that_vary     SparseMatrix       SparseVect  
for I = 1:numel(s)
    s(I).T.s = [];
    s(I).T.len = [];
    s(I).T.middle = [];
    s(I).T.SparseMatrix = [];
    s(I).T.columns_that_vary = [];
    s(I).T.size = [];
end    
save('EpistasisLBC__SupFig_ShowThatAllVariantsAreOKinSomeContext____nofilter.mat','s','-v7.3')


%% stupid pie chart :)
ph = pie( [ mean( all_N_fit==0) 1-mean( all_N_fit==0)] , {'',''})
set(ph(1),'FaceColor',[0.9 .9 .9])
set(ph(3),'FaceColor',[1 .85 .85])
title( mean( all_N_fit==0) )

%% misc old figures below here
% %% add 10 to 0s so that we can do log10()
% all_N_fit_2 = all_N_fit ; 
% all_N_unfit_2 = all_N_unfit ; 
% 
% all_N_fit_2(all_N_fit_2==0) = .1 ; 
% all_N_unfit_2(all_N_unfit_2==0) = .1 ; 
% 
% %%
% fh = figure('units','centimeters','position',[5 5 7 7 ]);
% %plot( all_fraction_fit*100 , all_fraction_unfit*100 ,'ok','MarkerFaceColor',[.7 .7 .7])
% scatter( all_fraction_fit*100 , all_fraction_unfit*100 , 20 , all_segments , 'filled')
% xlabel('% fit')
% ylabel('% unfit')
% xlim([0 100])
% ylim(xlim)
% line([0 100],[100 0],'LineStyle','--','Color',[.7 .7 .7])
% set(gca,'xtick',0:25:100)
% set(gca,'ytick',0:25:100)
% 
% 
% %%
% 
% fh = figure('units','centimeters','position',[5 5 7 7 ]);
% plot( all_N_fit_2 , all_N_unfit_2 ,'ok','MarkerFaceColor',[.7 .7 .7])
% %scatter( all_N_fit_2  , all_N_unfit_2   , 20 , all_segments , 'filled')
% xlabel('# fit')
% ylabel('# unfit')
% set(gca,'xscale','log')
% set(gca,'yscale','log')
% set(gca,'xtick',logspace(-1,7,9))
% set(gca,'ytick',logspace(-1,7,9))
% 
% xlim([.1 1e5])
% ylim([.1 1e5])
% 
% 
% set(gca,'xticklabel',[0 1 10 100 1e3 1e4 1e5 1e6 1e7 ])
% set(gca,'yticklabel',[0 1 10 100 1e3 1e4 1e5 1e6 1e7 ])
% 
% %%
% fh = figure('units','centimeters','position',[5 5 7 7 ]);
% hold on;
% histogram( log10(all_N_unfit+10)  ,50 ,'EdgeColor',[0 0 0] , 'FaceColor',[0 0 0]);
% histogram( log10(all_N_fit+10)  ,50 ,'EdgeColor','b' , 'FaceColor','b');
% legend( {'# unfit' '# fit' },'location','nw')
% set(gca,'xtick',[1.1 2 3 4])
% set(gca,'xticklabel',[0 100 1e3 1e4])
% ylabel('# of variants fit or unfit')
% axis tight;
% xlim([ min(xlim)*0.9 max(xlim)*1.01])
% xlabel('# of strains w/this AA')
% 
% %%
% fh = figure('units','centimeters','position',[5 5 7 7 ]);
% hold on;
% histogram( all_fraction_unfit * 100 ,50 ,'EdgeColor',[0 0 0] , 'FaceColor',[0 0 0]);
% histogram( all_fraction_fit  * 100 ,50 ,'EdgeColor','b' , 'FaceColor','b');
% legend( {'% unfit' '% fit' },'location','nw')
% ylabel('# of strains')
% xlim([-1 100.1])
% set(gca,'xtick',0:25:100)
% xlabel('% of strains w/this AA')
% 
% %%
% h = histogram2(  log10(all_N_fit_2) , log10(all_N_unfit_2) , 50 ,'DisplayStyle','bar3','ShowEmptyBins','off');
% h.FaceColor = 'flat';
% xlabel('# fit')
% ylabel('# unfit')
% zlabel('# of AA variants')
% set(gca,'xtick',[-1 0 1 2 3 4 5])
% set(gca,'xticklabel',[0 1 10 1e2 1e3 1e4 1e5])
% set(gca,'ytick',[-1 0 1 2 3 4 5])
% set(gca,'yticklabel',[0 1 10 1e2 1e3 1e4 1e5])
% xlim([-1 5])
% ylim(xlim)
% 
% %% for those found 0 times in fit / unfit, how common are they in the other category? 
% fh = figure('units','centimeters','position',[5 5 7 7 ]);
% hold on ;
% [f,x]=ecdf( all_N_fit( all_N_unfit==0) ) ; 
% plot(x,f,'LineWidth',3,'DisplayName','# unfit == 0')
% [f,x]=ecdf( all_N_unfit( all_N_fit==0)  ) ;
% plot(x,f,'LineWidth',3,'DisplayName','# fit == 0')
% xlim([0 1000])
% set(gca,'xscale','log')
% set(gca,'xtick',[1 2 5 10 25 100 250 1000])
% grid on ;
% ylabel('Fraction of variants')
% set(gca,'XMinorGrid','off')
% legend('location','se')
% xlabel('# with 0 in the other category')
% set(gca,'YMinorGrid','off')

