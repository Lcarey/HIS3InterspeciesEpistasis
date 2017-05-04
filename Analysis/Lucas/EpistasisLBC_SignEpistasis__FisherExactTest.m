
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

threshold = 0.2 ;


% test for interactions w/Sign Epistasis
for GI = 1:length(BigG)
    FitImpacts = BigG.FitImpact{GI} ;
    Seqs = BigG.Seq1{GI} ;
    
    
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
R = R(1:c,:);
R = sortrows(R , 'p' ,'ascend');
R.pBon = R.p * length(R.p) ;
R = R( R.pBon < 0.05 ,:);
% R.X11 = cellfun(@(X)X(1,1) ,R.X);
% R.X12 = cellfun(@(X)X(1,2) ,R.X);
% R.X21 = cellfun(@(X)X(2,1) ,R.X);
% R.X22 = cellfun(@(X)X(2,2) ,R.X);

%%
hmi = find(has_minor_impact);
for I = 1:numel(hmi)
    thisidx = hmi(I);
    s1 = char(Seqs(thisidx)) ;
    d = find( cellfun( @(X) HammingDistance( s1  , X) , Seqs) == 1 )' ;
    
    Seqs( [ thisidx d])
    FitImpacts( [ thisidx d])
end

