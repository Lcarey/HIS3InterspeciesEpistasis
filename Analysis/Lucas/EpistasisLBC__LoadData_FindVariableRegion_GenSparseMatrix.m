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
addOptional(p,'N_variants_to_fit',999999999,@isnumeric)
addOptional(p,'RUN_LINEAR_FLAG',0,@islogical)
addOptional(p,'ONLY_NATLIB_FLAG',1,@islogical)
addOptional(p,'ONLY_LIB_FLAG',1,@islogical)
addOptional(p,'NO_STOP_FLAG',1,@islogical)
addOptional(p,'ParamID','XXXX',@ischar)


parse(p,data_file_name,varargin{:});

MapAA2I = containers.Map(  arrayfun(@(X){X},['A':'Z' '_' ])  , 1:27 ) ; % all AAs + stop


%% setup paths
addpath(genpath('~/Develop/matlab/'));
addpath(genpath('~/Develop/Matlab/'));
addpath(genpath('~/Develop/HIS3InterspeciesEpistasis/'));


% if numel(gcp('nocreate'))==0 , 	parpool local , end
%% load data
% data_file_name can be a segment # or file name
DataDir = '~/Develop/HIS3InterspeciesEpistasis/Data/';
if isnumeric(data_file_name)
    data_file_name = [ 'S' num2str(data_file_name) '_scaled_info.csv'];
end
T = readtable([ DataDir data_file_name ] );

T = T( logical(T.nogap)  , :) ; % can't handle seq's w/varying length
T = T( T.len == mode(T.len) , :); % can't handle seq's w/varying length

if p.Results.NO_STOP_FLAG
    T = T( T.middle & T.nogap & ~T.stop & ~T.nonsense , :) ;
end
if p.Results.ONLY_NATLIB_FLAG
    T = T( logical(T.nat_lib) , :);
end
if p.Results.ONLY_LIB_FLAG
    T = T( logical(T.lib) , :);
end



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

%% edit distance between all pairs matrix
%  hamming_distance_matrix = squareform(pdist( cell2mat(T.aa_seq_variable) , 'hamming')) * n_positions ; 
%%
%% Distance matrix == work in progress
% % % 
% % % 
% % % 
% % % %% compare to WT seqs
% % % WT = readtable('~/Develop/HIS3InterspeciesEpistasis/Data_Small_Tables/wt_seq.csv');
% % % 
% % % 
% % % WT_seq = WT.aa_seq{ strcmp(WT.Var1,regexprep(data_file_name,'_.*',''))} ;
% % % WT_seq_varies = WT_seq(columns_that_vary) ;
% % % 
% % % R = dataset();
% % % R.matricies = {'BLOSUM30' 'BLOSUM40' 'BLOSUM50' 'BLOSUM62' 'BLOSUM70' 'BLOSUM80' 'BLOSUM90' 'BLOSUM100'...
% % %      'PAM30' 'PAM40' 'PAM50' 'PAM60' 'PAM70' 'PAM80' 'PAM90' ...
% % %     'DAYHOFF' 'GONNET'}' ;
% % % R.corrs = NaN( numel(R.matricies),1);
% % % 
% % % for J = 1:numel(R.matricies)
% % %     tic;
% % %     mat = R.matricies{J} ; %much faster
% % %     scores = cellfun( @(I) swalign( WT_seq_varies , I , 'Alphabet', 'AA' , 'ScoringMatrix' ,mat) , T.aa_seq_variable );
% % %     R.corrs(J) =  corr(scores,T.s);
% % %     fprintf('%s\t%s\t%0.03f\t%0.2f\n' , data_file_name , R.matricies{J}, R.corrs(J) , toc );
% % % end