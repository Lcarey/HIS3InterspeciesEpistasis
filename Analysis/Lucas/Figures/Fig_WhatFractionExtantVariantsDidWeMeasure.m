% Fig_WhatFractionExtantVariantsDidWeMeasure.m
% LBC May 11, 2017
% load the multiple sequence alignment and 
%   find out what fraction of extant variation IN FUNGI did we measure? 
%
%  this is different from all extant variation, and different from nat_lib (569 species)
%
%
%
% CODE DOESN"T WORK PROPERLY!
%   can't deal with insertions in the alignment. 
DDIR = '~/Develop/HIS3InterspeciesEpistasis/Data_Small_Tables/' ;
ALIGN = readtable( [ DDIR 'HIS3_21species_Alignment.stab'],'ReadVariableNames',false,'Delimiter','\t','Format','%s%s','FileType','text');
P = readtable('~/Develop/HIS3InterspeciesEpistasis/Data_Small_Tables/positions.csv','ReadRowNames',true);
ALIGN = char( ALIGN.Var2);
%% load data
for SegN = 1:12
    T = readtable( [ '~/Develop/HIS3InterspeciesEpistasis/Data/S' num2str(SegN) '_scaled_info_v2.csv'] , 'FileType','text','Delimiter','\t');
%    T = T( logical(T.nat_lib) & logical(T.nogap) & logical(T.middle) & ~logical(T.stop) & ~logical(T.nonsense) ,:);
    T = T(  logical(T.nogap) & logical(T.middle) & ~logical(T.stop) & ~logical(T.nonsense) ,:);
    T = T( T.len == mode(T.len) ,:); % I think unnecessary
    T.SegN = repmat( SegN , height(T) , 1);
    T = T(: , {'aa_seq' 'SegN'});
    if SegN==1
        Q = T;
    else
        Q = vertcat(Q,T);
    end
end
T = Q;

%% Calculate how many of the extant AAs at each position we measured
naa = NaN( numel(ALIGN(1,:)) , 2) ; % extant ; measured ; segn
aa_measured = cell( length(naa),1);
aa_extant   = cell( length(naa),1);

for SegN = 1:12
    measured_seqs = char(T.aa_seq( T.SegN==SegN)) ; 
    positions_in_seg = str2double(regexp( regexprep( char(P{ 'positions' , ['S' num2str(SegN)]}) , '[\[\] ]','') , ',' ,'split'));
    positions_in_seg = positions_in_seg +1 ; %in file as base 0
    rel_positions_in_seg    = positions_in_seg - positions_in_seg(1) + 1 ; 
    for PosI = 1:numel(positions_in_seg)
        c = positions_in_seg(PosI) ;
        aa_extant{c} = unique(ALIGN(:,positions_in_seg(PosI)))' ;
        aa_measured{c} = unique(measured_seqs(:, rel_positions_in_seg(PosI)));
      
        naa( positions_in_seg(PosI)  , 1) = mean(ismember(aa_extant{c},aa_measured{c}));
        naa( positions_in_seg(PosI)  , 3) = SegN;
    end
end
%%
figure
bar( 1:length(naa) , naa(:,1)*100  , 1.0 ,'FaceColor','k' );
axis tight;
