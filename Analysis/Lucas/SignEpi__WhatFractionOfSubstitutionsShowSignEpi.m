%% Figure: X % show sign epi
load('SignEpi_Pairs_04.mat')

%%
IDS1 = strcat( R.Perm , num2str(R.VarPos) , num2str(R.SegN));
IDS2 = strcat( R.SubPerm , num2str(R.SubPos) , num2str(R.SegN));
IDS = unique(vertcat(IDS1,IDS2));

%%
figure; 
hold on; 
t = [0.00001 0.0001 0.001 0.01 0.1 1 100000] ; 
for I = 1:numel(t)
    idx = R.pBon <= t(I) ;

    i1 = strcat( R.Perm(idx) , num2str(R.VarPos(idx)) , num2str(R.SegN(idx)));
    i2 = strcat( R.SubPerm(idx) , num2str(R.SubPos(idx)) , num2str(R.SegN(idx)));
    ids = unique(vertcat(i1,i2));
    Y = mean(ismember(IDS,ids));
    plot(I,Y,'ok','MarkerFaceColor',[.7 .7 .7])
end

set(gca,'xtick',1:numel(t))
set(gca,'xticklabel',t)
set(gca,'ytick',0:0.05:1)
grid on ;

%% how often w/sig bonfer

load('SignEpi_BigG_04.mat');
%%
   idx = R.pBon <= 0.01 ; 
    i1 = strcat( R.Perm(idx) , num2str(R.VarPos(idx)) , num2str(R.SegN(idx)));
    i2 = strcat( R.SubPerm(idx) , num2str(R.SubPos(idx)) , num2str(R.SegN(idx)));
    sigids = unique(vertcat(i1,i2));
    
IDS = strcat( BigG.Perm , num2str(BigG.VarPos) , num2str(BigG.SegN));

biggidx = ismember( IDS, sigids);

%%
figure;
dscatter( log10(BigG.mean_HasMinorSignFitEffect(biggidx)*100)  ,  log10( 1+BigG.sum_HasMinorSignFitEffect(biggidx)) )


%%

load('SignEpi_Ronly.mat')
Q = s(4).R; 
%%

Q10IP = Q( Q.VarPos==10 & strcmp(Q.Perm,'IP'),:);
Q2AG = Q( Q.VarPos==2 & strcmp(Q.Perm,'AG'),:);
%%
ip10I1 = cellfun( @(X) X(10)  , Q2AG.Seq1) == 'I' ;
ip10I2 = cellfun( @(X) X(10)  , Q2AG.Seq2) == 'I' ;
ip10P1 = cellfun( @(X) X(10)  , Q2AG.Seq1) == 'P' ;
ip10P2 = cellfun( @(X) X(10)  , Q2AG.Seq2) == 'P' ;

figure ; hold on; 
plot( Q2AG.Fit1 , Q2AG.Fit2,'o','Color',[.7 .7 .7]); 
plot( Q2AG.Fit1(ip10I1) , Q2AG.Fit2(ip10I1),'.'); 
%plot( Q2AG.Fit1(ip10I2) , Q2AG.Fit2(ip10I2),'.'); 
plot( Q2AG.Fit1(ip10P1) , Q2AG.Fit2(ip10P1),'.'); 
%plot( Q2AG.Fit1(ip10P2) , Q2AG.Fit2(ip10P2),'.'); 

xlabel('Fit w/A') ; ylabel('Fit w/G')
line([0 1 ] , [0 1])
