function T = EpistasisLBC__SynVariantsCAI(SegI)
%% Do synonymous variants that differ by codon bias have different expression?
%% EpistasisLBC__SynVariantsCAI.m
%% load data
%%
DataDir = '~/Develop/HIS3InterspeciesEpistasis/Data/Synonymous/';
fns = dir( [DataDir filesep 'S*_sum_*fit.csv']);
%for SegI = 1:12
metadata = regexp( fns(SegI).name , '_' ,'split');
fns(SegI).segment = str2double(metadata{1}(2:end));
T = readtable( [ fns(SegI).folder filesep fns(SegI).name ] );
T = T( cellfun( @length , T.seq) ==  mode(cellfun( @length , T.seq)) , :); % only mode seq length
T = T( ~isnan(T.s) , :);
% Calc CAI / nTE , etc
tic;
T.nTE  = cellfun(@(X)calc_tAI(X,'nTE'),T.seq) ;
T.tAI  = cellfun(@(X)calc_tAI(X,'tAI'),T.seq) ;
fprintf('Calc CAI took %0.01f sec\n' , toc);
% Correlation between CAI & fitness in different ranges of mean AA fitness?

tic;
G = grpstats(T ,'aa_seq' ,'median','DataVars','s');
G = G( G.GroupCount  > 1 , :);
D = T( ismember( T.aa_seq , G.aa_seq) ,:);
C = NaN( numel(G.aa_seq) , 8);

nTE_all = cell( height(G) , 1);
tAI_all = cell( height(G) , 1);
fit_all = cell( height(G) , 1);

for I = 1:height(G)
    idx = find(strcmp( D.aa_seq , G.aa_seq{I}));
    [s,o] = sort(D.s(idx)) ;
    nTE =  D.nTE(idx(o));
    tAI = D.tAI(idx(o));
    
    if numel(tAI)>=6
        C(I,3) = log2( mean(nTE(end-2:end)) /  mean(nTE(1:3)) );
        C(I,6) = log2( mean(tAI(end-2:end)) /  mean(tAI(1:3)) );
    end
    C(I,7) = log2( mean(nTE(end)) /  mean(nTE(1)) );
    C(I,8) = log2( mean(tAI(end)) /  mean(tAI(1)) );
    
    [ C(I,1) , C(I,2) ] = corr(s,nTE) ;
    [ C(I,4) , C(I,5) ] = corr(s,tAI) ;
    
    nTE_all{I} = nTE;
    tAI_all{I} = tAI;
    fit_all{I} = s;
end
G.nTE_C = C(:,1);
G.nTE_P = C(:,2);
G.tAI_C = C(:,4);
G.tAI_P = C(:,5);
G.nTE_l2_3 = C(:,3);
G.tAI_l2_3 = C(:,6);
G.nTE_l2_1 = C(:,7);
G.tAI_l2_1 = C(:,8);

G.tAI_all = tAI_all ;
G.fit_all = fit_all ;
G.nTE_all = nTE_all ;
fprintf('Calc corr took %0.01f sec\n' , toc);

fns(SegI).T = T ;
fns(SegI).G = G ;

SegI
s = fns(SegI) ;
save( sprintf('EpistasisLBC__SynVariantsCAI_%02d_%s.mat',SegI,char(datetime)),'s')
%end
% %% plot results
% G.s = round(G.median_s*5)/5 ;
% MS = 0;
% N = 20 ;
% min_syn_vars_for_plotting = MS ;
% Q = G( G.GroupCount>=min_syn_vars_for_plotting,:);
% %       Q = Q( Q.nTE_P<0.05 | Q.tAI_P<0.05 ,:);
% [Y,E] = discretize(Q.median_s , N);
% Q.s = Y;
%
%
% figure;
% subplot(2,2,1)
% boxplot( Q.nTE_C , Q.s ,'notch','on');line(xlim,[0 0]);
% ylabel('nTE')
% ylim([-0.5 0.5]); set(gca,'ytick',-1:.1:1);grid on;
% xlabel('Measured fitness')
%
% subplot(2,2,2)
% boxplot( Q.nTE_l2 , Q.s ,'notch','on');line(xlim,[0 0]);
% ylabel('nTE log2(FIT/unfit)')
% ylim([-0.1 0.1]); set(gca,'ytick',-1:.1:1);grid on;
% xlabel('Measured fitness')
%
% subplot(2,2,3)
% boxplot( Q.tAI_C , Q.s ,'notch','on');line(xlim,[0 0]);
% ylabel('tAI')
% ylim([-0.5 0.5]); set(gca,'ytick',-1:.1:1);grid on;
% xlabel('Measured fitness')
%
% subplot(2,2,4)
% boxplot( Q.tAI_l2 , Q.s ,'notch','on');line(xlim,[0 0]);
% ylabel('tAI log2(FIT/unfit)')
% ylim([-0.1 0.1]); set(gca,'ytick',-1:.1:1);grid on;
% xlabel('Measured fitness')
%
% %%
%
% for MS = 0
%     figname = sprintf('SynVarParSearch_%d.eps' ,MS);
%     delete(figname);
%     for N = 20
%
%         min_syn_vars_for_plotting = MS ;
%         Q = G( G.GroupCount>=min_syn_vars_for_plotting,:);
%         %    Q = Q( Q.nTE_P<0.05 | Q.tAI_P<0.05 ,:);
%         [Y,E] = discretize(Q.median_s , N);
%         Q.s = Y ;
%         GG = grpstats( Q ,'s' ,{'mean'  'median' 'std'} , 'DataVars' ,{'median_s' 'nTE_C' 'tAI_C' 'nTE_l2' 'tAI_l2'});
%         GG = GG(GG.GroupCount>1,:);
%         figure; hold on ;
%         plot(  GG.median_median_s , GG.median_nTE_C ,'.-','DisplayName','nTE C');
%         plot(  GG.median_median_s , GG.median_tAI_C ,'.-','DisplayName','tAI C');
%         ylim([0 0.2])
%         legend('location','best')
%         set(gca,'ytick',-1:0.01:1)
%         grid on ;
%         title( sprintf( 'MS=%d N=%d' , MS , N));
%         xlim([-1.8 0.75])
%         set(gca,'xtick',-5:0.2:5);
%         %  print('-dpsc2', figname,'-append');
%         %   close;
%
%     end
% end
%
%
%
