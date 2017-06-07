%% How many nat_lib have 1,2,3,4,5,6,7 synonymous variants
% May 17, 2017
%% load data
Efitness = NaN(0);
Esegn = NaN(0);
Ensyn = NaN(0);
Lfitness = NaN(0);
Lsegn = NaN(0);
Lnsyn = NaN(0);
for SegN = 1:12
    T = readtable( [ '~/Develop/HIS3InterspeciesEpistasis/Data/S' num2str(SegN) '_scaled_info_v2.csv' ] ,'FileType','text','Delimiter','\t');
    Eidx = logical(T.nogap) & logical(T.middle) & ~logical(T.stop) & ~logical(T.nonsense) & logical(T.nat_lib);
    Efitness = vertcat( Efitness , T.s(Eidx));
    Esegn = vertcat( Esegn , repmat( SegN , sum(Eidx) , 1) );
    Ensyn = vertcat( Ensyn , T.size(Eidx));
    Lidx = logical(T.nogap) & logical(T.middle) & ~logical(T.stop) & ~logical(T.nonsense) & ~logical(T.nat_lib) & ~logical(T.lib);
    Lfitness = vertcat( Lfitness , T.s(Lidx));
    Lsegn = vertcat( Lsegn , repmat( SegN , sum(Lidx) , 1) );
    Lnsyn = vertcat( Lnsyn , T.size(Lidx));    
end

%%
h = figure('units','centimeters','position',[5 5 6 5]);hold on;
for I = 1:10
    bar( I , 100*mean( Lnsyn==I) ,'FaceColor','b','FaceAlpha',1);
    bar( I , 100*mean( Ensyn==I) ,'FaceColor','r','FaceAlpha',1);
    if I>1
            bar( I , 100*mean( Lnsyn==I) ,'FaceColor','b','FaceAlpha',1);
    end
end
set(gca,'xtick',1:100)
set(gca,'ytick',0:20:100)

xlim([0.5 I+0.5])
xlabel('# of synonymous variants')
ylabel('% of extant library')
