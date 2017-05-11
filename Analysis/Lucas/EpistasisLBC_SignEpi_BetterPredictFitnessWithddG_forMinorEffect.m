%% can we better predict minor sign effects w/structure (ddG) than w/fitness potnetial?
%
%

figname = 'FIp_04.eps';
delete(figname);


load('~/Desktop/HIS3scratch/SignEpi/SignEpi_Pairs_04.mat');
load('~/Desktop/HIS3scratch/SignEpi/SignEpi_BigG_04.mat');
load('/Users/lcarey/Desktop/HIS3scratch/SignEpi/SignEpi_Ronly.mat');

FitImpactThreshold = 0.4 ;

ddG = readtable('~/Develop/HIS3InterspeciesEpistasis/Analysis/Sasha/rosetta_runs/run-170503-results.csv','FileType','text','Delimiter','\t');
ddG.dG = median( table2array(ddG(:,3:5)) , 2) ;

%%    
for SegN = 1:12
    
    R = s(SegN).R ;
    
    G = BigG( BigG.SegN==SegN ,:);
    FP = readtable([ '~/Develop/HIS3InterspeciesEpistasis/Analysis/Katya/NN/residuals/S' num2str(SegN) '.csv'],'FileType','text','Delimiter',',');
    
    T  = EpistasisLBC__LoadData( SegN  ,'ONLY_NATLIB_FLAG',true, 'ONLY_MIDDLE_FLAG',true,'NO_STOP_FLAG',true);
    
    T = join( T(:,{'aa_seq' 's'}) , FP ,'Key','aa_seq');
    T = innerjoin( T  ,ddG(:,[1 7]) ,'Key','aa_seq');
   
    % Predict fitness from ddG
    % generate residuals
    x = T.fitnessPotential';
    t = T.s';
    trainFcn = 'trainscg';  % Scaled conjugate gradient backpropagation. 'trainlm'  OR 'trainbr' OR 'trainscg'
    hiddenLayerSize = 25;
    net = fitnet(hiddenLayerSize,trainFcn);
    net.divideParam.trainRatio = 70/100;
    net.divideParam.valRatio = 15/100;
    net.divideParam.testRatio = 15/100;
    [net,~] = train(net,x,t);
    T.Residuals_FitPotential = gsubtract( t , net(x) )' ;
    
    x = T.dG';
    t = T.s';
    trainFcn = 'trainscg';  % Scaled conjugate gradient backpropagation. 'trainlm'  OR 'trainbr' OR 'trainscg'
    hiddenLayerSize = 25;
    net = fitnet(hiddenLayerSize,trainFcn);
    net.divideParam.trainRatio = 70/100;
    net.divideParam.valRatio = 15/100;
    net.divideParam.testRatio = 15/100;
    [net,tr] = train(net,x,t);
    T.Residuals_dG = gsubtract( t , net(x) )' ;
    
    T.predictedMinusObserved = [];
    
    %
    R1 = innerjoin( dataset2table(R) , T ,'LeftKey','Seq1','RightKey','aa_seq');
    vn = R1.Properties.VariableNames ;
    idx = find(ismember(vn,{'fitnessPotential'  'dG'  'Residuals_FitPotential'    'Residuals_dG'}));
    for tvi = idx
        R1.Properties.VariableNames{tvi} = [ vn{tvi} '_Seq1'] ;
    end
    R0 = innerjoin( R1, T ,'LeftKey','Seq2','RightKey','aa_seq');
    vn = R0.Properties.VariableNames ;
    idx = find(ismember(vn,{'fitnessPotential' 'dG'   'Residuals_FitPotential'    'Residuals_dG'}));
    for tvi = idx
        R0.Properties.VariableNames{tvi} = [ vn{tvi} '_Seq2'] ;
    end
    
    % dG is Mut stability (more negative == more stable)
    R0.dG_log2_s1_s2 = log2( R0.dG_Seq1 ./ R0.dG_Seq2);
 
    SignChange = zeros(height(R0),1);
    SignChange( R0.FitImpact< (-1 * FitImpactThreshold)) = -1 ;
    SignChange( R0.FitImpact>FitImpactThreshold) = 1 ;
    
    for GI = 1:length(G)
        ridx = strcmp(R0.Perm,G.Perm{GI}) & R0.VarPos == G.VarPos(GI)    ;
        if numel(unique(SignChange(ridx) ))==3
            [~,p]=ttest2( R0.dG_log2_s1_s2(ridx & SignChange==-1) , R0.dG_log2_s1_s2(ridx & SignChange==1) );
            p = p * length(G) ; % Bonferroni corrrection
            
            if p < 0.1
            fh = figure('units','centimeters','position',[5 5 7 10]);
            hold on; 
            h = notBoxPlot( R0.dG_log2_s1_s2(ridx)   , SignChange(ridx));
            set( [h.data] ,'MarkerSize',2);
            set(gca,'xtick',-1:1)
             
            [a,b]=count_unique(SignChange(ridx));
%            b = b./sum(b)*100;  % convert to percent
            ylabel('log2( dG_1 / dG_2 )')
            line( [-10 10] , [0 0] ,'LineStyle','--','Color',[.7 .7 .7])
            set(gca,'xticklabel' , { sprintf('- (%d)',b(1)) ,sprintf('0 (%d)',b(2)) ,sprintf('+ (%d)',b(3)) })
%            set(gca,'xticklabel' , { sprintf('- (%0.0f%%)',b(1)) ,sprintf('0 (%0.0f%%)',b(2)) ,sprintf('+ (%0.0f%%)',b(3)) })
            if b(1)>b(3)
                MajorSign = '-';
            else
                MajorSign = '+';
            end
            
            title( sprintf( '%d %s %d %s %d p=%0.03f' ,SegN, MajorSign, GI , G.Perm{GI} , G.VarPos(GI)  , p ) );
            print('-dpsc2',figname,'-append');
            close;
            
        end
        end
    end
    
    
end

%% pick one of interest
SegN=2 ; 
GI = 4 ; %  GN 2
Q = R0(ridx,:);

%%
figure;
scatter( Q.dG_Seq1 , Q.dG_Seq2 , 30 , SignChange(ridx),'filled')
line([-560 -540] ,[-560 -540])
axis tight
xlabel('dG Seq 1')
ylabel('dG Seq 2')

figure;
scatter( Q.Fit1 , Q.Fit2 , 30 , SignChange(ridx),'filled')
line( [ min(Q.Fit1) max(Q.Fit1) ] , [ min(Q.Fit1) max(Q.Fit1) ] )
axis tight
xlabel('Fit 1')
ylabel('Fit 2')