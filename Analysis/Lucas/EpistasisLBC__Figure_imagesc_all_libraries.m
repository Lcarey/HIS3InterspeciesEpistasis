function files_struct = EpistasisLBC__Figure_imagesc_all_libraries()
%% EpistasisLBC__Figure_imagesc_all_libraries()
%
%  For each library, plot imagesc() of the variants
%   both sorted and not sorted by fitness
%
% LBC March 2017


%% First step: how much variability is found in our library at each position?
%% load variability at each position
WTseqsFile = '~/Develop/HIS3InterspeciesEpistasis/Data_Small_Tables/AA_variants.csv';
tmpfilename = tempname() ;
system([ 'cut -f 1,3,5  ' WTseqsFile  ' > ' tmpfilename ]);
V = readtable(tmpfilename,'FileType','text','Delimiter','\t','HeaderLines',1,'Format','%d%s%s','ReadVariableNames',false);
delete( tmpfilename)
V.Properties.VariableNames = {'LibN','Nat','Lib'};
V.Nat = cellfun( @(X)regexp(X,',','split') , V.Nat , 'UniformOutput',false);
V.Lib = cellfun( @(X)regexp(X,',','split') , V.Lib , 'UniformOutput',false);
V.PctFound = arrayfun(@(I)mean(ismember(V.Nat{I},V.Lib{I})) , 1:height(V)  )' ;
V = V( cellfun(@length,V.Lib)>1 ,:); % no single AA pos
%%
MapAA2I = containers.Map(  arrayfun(@(X){X},['A':'Z' '_' ])  , 1:27 ) ; % all AAs + stop

files_struct = dir('~/Develop/HIS3InterspeciesEpistasis/Data/S*_scaled_info.csv');

for I = 1:12
    metadata = regexp( files_struct(I).name  ,'_','split');
    files_struct(I).LibN = str2double(metadata{1}(2:end));
    T = readtable( [ files_struct(I).folder filesep files_struct(I).name] ,'FileType','text','Delimiter','\t');
    T = T( T.len == mode(T.len) , : );
    T = T( logical(T.middle) , :);
    T = T( ~logical(T.stop) , :);
    T = T( ~logical(T.nonsense) , :);
    T = T( logical(T.nat_lib) ,:);
    [ T  , ~ , ~ ] = EpistasisLBC__ShortedAAseqsOnlyPositionsThatVary( T );
    
    files_struct(I).T = T ;
    
    files_struct(I).AAmat = cell2mat( cellfun( @(X)arrayfun(@(I)MapAA2I(I), X) ,  T.aa_seq_variable,'UniformOutput',false) ) ;
    %    files_struct(I).nat_lib = logical(T.nat_lib) ;
    %    files_struct(I).lib = logical(T.lib) ;
    %    files_struct(I).middle = logical(T.middle) ;
    %    files_struct(I).stop = logical(T.stop) |  logical(T.nonsense) ;
    %    files_struct(I).lib = logical(T.lib) ;
    files_struct(I).fitness = T.s ;
    
    fprintf(' %d '  , I );
end

AAs = sort('GALMFWKQESPVICYHRNDT');
%%
N=7;
FigSize = [5 5 8 8] ;
for I = 1:numel(files_struct)
    R=7;
    C=9;
    
    fitness_bar_spot = C:C:(C*R) ;
    bar_spot = (C*R)-C+1:(C*R)-1 ;
    major_spot = 1:(C*R) ; 
    major_spot = major_spot( ~ismember( major_spot,fitness_bar_spot));
    major_spot = major_spot( ~ismember( major_spot,bar_spot)) ;
    
    n = numel(files_struct(I).fitness) ;
    
    [~,o] = sort( files_struct(I).fitness , 'descend');
    
    fh = figure('units','centimeters','position',[5 5 8 8]);
    
    sph = subplot(R,C,major_spot);
    colormap(sph,parula)
    imagesc( files_struct(I).AAmat(o,:) );
    set(gca,'ytick',[])
    set(gca,'xtick',[])
    
    title( sprintf(  'S%d N=%d' ,files_struct(I).LibN , n ) )
    
    subplot(R,C,fitness_bar_spot)
    colormap(jet);
    imagesc( files_struct(I).fitness(o) , [0 0.55]);
    set(gca,'ytick',[])
    set(gca,'xtick',[])
    
    subplot(R,C,bar_spot)
    data = V.PctFound( V.LibN==files_struct(I).LibN );
    bar( data ,'FaceColor','k' )
    axis tight;
    ylim([0 1])
    set(gca,'ytick',[])
    set(gca,'xtick',1:100)
    line( xlim , [0.25 0.25] , 'LineStyle','-','Color',[.75 .75 .75])
    line( xlim , [0.5 0.5] , 'LineStyle','-','Color',[.75 .75 .75])
    line( xlim , [0.75 0.75] , 'LineStyle','-','Color',[.75 .75 .75])
    box off ;
    
    print('-dpng2',['HeatMaps_' num2str(files_struct(I).LibN) '_FIT_NATLIB_sorted.png'],'-r600');
    close;
    
end
