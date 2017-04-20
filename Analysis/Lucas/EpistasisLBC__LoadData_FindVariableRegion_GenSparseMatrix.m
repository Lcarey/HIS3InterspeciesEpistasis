function T  = EpistasisLBC__LoadData_FindVariableRegion_GenSparseMatrix( data_file_name  ,  varargin )
%% T = EpistasisLBC__LoadData_FindVariableRegion_GenSparseMatrix( data_file_name  ,  varargin )
%  loads data_file_name (eg: 7)
%  finds all positions that vary
%  
%  generates a sparse vector w/each AA and each position that varies
%
% LBC March 2017

if ~exist('varargin','var');varargin={};end;
p = inputParser;
addRequired(p,'data_file_name');
addOptional(p,'ONLY_NATLIB_FLAG',1,@islogical)
addOptional(p,'ONLY_LIB_FLAG',1,@islogical)
addOptional(p,'ONLY_MIDDLE_FLAG',1,@islogical)
addOptional(p,'NO_STOP_FLAG',1,@islogical)


parse(p,data_file_name,varargin{:});

MapAA2I = containers.Map(  arrayfun(@(X){X},['A':'Z' '_' ])  , 1:27 ) ; % all AAs + stop

% if numel(gcp('nocreate'))==0 , 	parpool local , end
%% load data
% data_file_name can be a segment # or file name
DataDir = '~/Develop/HIS3InterspeciesEpistasis/Data/';
if isnumeric(data_file_name)
    data_file_name = [ num2str(data_file_name) '.tab'];
end
% format = 
%  aa_seq	s	len	nogap	stop	middle	nat_lib	nat	lib	size
T = readtable([ DataDir data_file_name ] , 'Format' , '%s%f%d%d%d%d%d%d%d%d%d' ,'FileType','text','Delimiter','\t');

T = T( logical(T.nogap)  , :) ; % can't handle seq's w/varying length
T = T( T.len == mode(T.len) , :); % can't handle seq's w/varying length

if p.Results.NO_STOP_FLAG
    T = T(T.nogap & ~T.stop & ~T.nonsense , :) ;
end
if p.Results.ONLY_NATLIB_FLAG
    T = T( logical(T.nat_lib) , :);
end
if p.Results.ONLY_LIB_FLAG
    T = T( logical(T.lib) , :);
end
if p.Results.ONLY_MIDDLE_FLAG
    T = T( logical(T.middle) , :);
end

T.nogap = [];
T.stop = [];
T.nonsense = [];


%% shrink sequences down to the variable part
[ T  , columns_that_vary , uniq_aas ] = EpistasisLBC__ShortedAAseqsOnlyPositionsThatVary( T );
T.columns_that_vary = repmat( uint8(find(columns_that_vary)) , height(T) , 1) ;

%% generate sparse matrix

n_positions = numel(T.aa_seq_variable{1}) ;
n_AAs = 27 ;

T.SparseMatrix = cell( height(T) , 1);
T.SparseVect = cell( height(T) , 1);

for I = 1:height(T)
    this_seq_numeric = arrayfun( @(X)MapAA2I(X) , T.aa_seq_variable{I} ) ;
    sparse_vect = zeros( 1, length(MapAA2I) * n_positions );
    idx = ( ((0:(n_positions-1)) .* n_AAs)+this_seq_numeric )  ;
    sparse_vect(idx) = 1;
    sparse_mat = reshape( sparse_vect , n_AAs , n_positions );
    T.SparseVect{I} = sparse_vect ;
    T.SparseMatrix{I} = sparse_mat ;
end

end

