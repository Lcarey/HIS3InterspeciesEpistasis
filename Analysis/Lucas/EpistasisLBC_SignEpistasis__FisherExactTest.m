load('~/Desktop/HIS3scratch/SignEpi/SignEpi_BigG_04.mat')
%load('~/Desktop/HIS3scratch/SignEpi/SignEpi_Ronly.mat');
%%

R = dataset();
R.SegN = NaN(1e6,1);
R.VarPos  = NaN(1e6,1);
R.Perm = cell(1e6,1);
R.SubPos =  NaN(1e6,1);
R.SubPerm = cell(1e6,1);
R.p    = NaN(1e6,1);
R.logodds = NaN(1e6,1);
R.X = cell(1e6,1);
R.GI = NaN(1e6,1);
c = 0 ;

threshold = 0.4 ;


% test for interactions w/Sign Epistasis
for GI = 1:length(BigG)
    fprintf('%d/%d (%0.0f%%)\n' , GI , length(BigG) , GI/length(BigG)*100)
    FitImpacts = BigG.FitImpact{GI} ;
    Seqs = BigG.Seq1{GI} ;
    
    if sum(abs(FitImpacts) > threshold )>0
    
    Seqs = Seqs( abs(FitImpacts) > threshold ) ;
    FitImpacts = FitImpacts( abs(FitImpacts) > threshold ) ;
    
    
    if mean(FitImpacts > threshold ) <= mean(FitImpacts < (-1 * threshold ) )
        has_minor_impact = FitImpacts > threshold ;
        has_major_impact = FitImpacts < (-1 * threshold ) ;
    else
        has_major_impact = FitImpacts > threshold ;
        has_minor_impact = FitImpacts < (-1 * threshold ) ;
    end
    
    for posI = 1:length(Seqs{1})
        aas_at_this_pos = cellfun(@(X)X(posI) , Seqs) ;
        unique_aas_at_this_position = unique(aas_at_this_pos) ;
        if numel(unique_aas_at_this_position)>1 
            permutations = nchoosek(unique_aas_at_this_position,2) ;
            for permI = 1:size(permutations,1)
                c = c + 1 ;
                idx1 = ( aas_at_this_pos==permutations(permI,1) ) ;
                idx2 = ( aas_at_this_pos==permutations(permI,2) ) ;
                X = [ sum(idx1 & has_major_impact)  sum(idx1 & has_minor_impact)  ; sum(idx2 & has_major_impact)  sum(idx2 & has_minor_impact) ] ;
                [~ , R.p(c) , stats] = fishertest( X ) ;
                R.logodds(c) = stats.OddsRatio ;
                R.SubPerm{c} = permutations(permI,:);
                R.VarPos(c) = BigG.VarPos(GI) ;
                R.SubPos(c) = posI ;
                R.SegN(c) = BigG.SegN(GI) ;
                R.Perm{c} = BigG.Perm{GI} ;
                R.X{c} = X ;
                R.GI(c) = GI ;
            end
        end
    end
    end
end
R = R( ~isnan(R.SegN) , :);
R = sortrows(R , 'p' ,'ascend');
R.pBon = R.p * length(R.p) ;
R.ReallyPositivePair = R.pBon < 0.01 ; 
R.ReallyNegativePair = R.p > 0.1 & R.logodds < 1 ; 
%% convert to abs S cer locations

[ ~ , PosTable ] = EpistasisLBC_ConvertRelativeSegPosToAbsScerPos( 1  , 1 );

R.ScerPos_WithSign = NaN(length(R),1);
R.ScerPos_PartnerSite = NaN(length(R),1);
for I = 1:length(R)
    R.ScerPos_WithSign(I) = EpistasisLBC_ConvertRelativeSegPosToAbsScerPos( R.SegN(I) , R.VarPos(I) , PosTable );
    R.ScerPos_PartnerSite(I) = EpistasisLBC_ConvertRelativeSegPosToAbsScerPos( R.SegN(I) , R.SubPos(I)   , PosTable);   
end
  


%% make sure that all node IDs are in the same order
for I = 1:length(R)
    R.Perm{I} = sort(R.Perm{I},'ascend');
    R.SubPerm{I} = sort(R.SubPerm{I},'ascend');
end

R.Xf = cellfun( @(X)X./sum(X(:)) , R.X,'UniformOutput',false);
save('~/Desktop/HIS3scratch/SignEpi/SignEpi_Pairs_04.mat','R')
writetable( dataset2table(R) ,'SignEpiPairs.xlsx')
writetable( dataset2table(R(R.pBon<0.01 & R.logodds>5 ,: ) ) ,'SignEpiPairs_ReallyInterestingPositives.xlsx')

%% generate node IDs
R.IDS_SignEpi = cell(length(R),1);
R.IDS_Partner = cell(length(R),1);
for I = 1:length(R)
    R.IDS_SignEpi{I} = sprintf('SegN=%d VarPos=%d Perm=%s' , R.SegN(I) , R.VarPos(I) , R.Perm{I} ) ;
    R.IDS_Partner{I} = sprintf('SegN=%d SubPos=%d SubPerm=%s' , R.SegN(I) , R.SubPos(I) , R.SubPerm{I}) ;
end
R.Nid_SignEpi = double(categorical(R.IDS_SignEpi)) ; 
R.Nid_Partner = double(categorical(R.IDS_Partner)) ; 
%%
writetable( dataset2table(R( R.SegN==5  , {'Nid_SignEpi' 'Nid_Partner'}))  , 'edgelist_S5_all.tab','FileType','text','Delimiter','\t','WriteVariableNames',false)
writetable( dataset2table(R( R.pBon<0.01 &   R.SegN==5  , {'Nid_SignEpi' 'Nid_Partner'}))  , 'edgelist_S5_all.tab','FileType','text','Delimiter','\t','WriteVariableNames',false)
writetable( dataset2table(R( R.pBon<0.01 &   R.SegN==5  , {'Nid_SignEpi' 'Nid_Partner'}))  , 'edgelist_S5_01.tab','FileType','text','Delimiter','\t','WriteVariableNames',false)
writetable( dataset2table(R( R.pBon<0.001 &   R.SegN==5  , {'Nid_SignEpi' 'Nid_Partner'}))  , 'edgelist_S5_001.tab','FileType','text','Delimiter','\t','WriteVariableNames',false)
writetable( dataset2table(R( R.pBon<0.01  & R.logodds>1 &   R.SegN==5  , {'Nid_SignEpi' 'Nid_Partner'}))  , 'edgelist_S5_01_1.tab','FileType','text','Delimiter','\t','WriteVariableNames',false)
writetable( dataset2table(R( R.pBon<0.01  & R.logodds>5 &   R.SegN==5  , {'Nid_SignEpi' 'Nid_Partner'}))  , 'edgelist_S5_01_5.tab','FileType','text','Delimiter','\t','WriteVariableNames',false)
writetable( dataset2table(R( R.pBon<0.01  & R.logodds>10 &   R.SegN==5  , {'Nid_SignEpi' 'Nid_Partner'}))  , 'edgelist_S5_01_10.tab','FileType','text','Delimiter','\t','WriteVariableNames',false)
writetable( dataset2table(R( R.pBon<0.001  & R.logodds>1 &   R.SegN==5  , {'Nid_SignEpi' 'Nid_Partner'}))  , 'edgelist_S5_001_1.tab','FileType','text','Delimiter','\t','WriteVariableNames',false)
writetable( dataset2table(R( R.pBon<0.001  & R.logodds>5 &   R.SegN==5  , {'Nid_SignEpi' 'Nid_Partner'}))  , 'edgelist_S5_001_5.tab','FileType','text','Delimiter','\t','WriteVariableNames',false)
writetable( dataset2table(R( R.pBon<0.001  & R.logodds>10 &   R.SegN==5  , {'Nid_SignEpi' 'Nid_Partner'}))  , 'edgelist_S5_001_10.tab','FileType','text','Delimiter','\t','WriteVariableNames',false)

%%
qqq = grpstats(R , {'SegN' 'VarPos' 'Perm'  'SubPos' 'SubPerm'} ,'min','DataVars','pBon')
    
    