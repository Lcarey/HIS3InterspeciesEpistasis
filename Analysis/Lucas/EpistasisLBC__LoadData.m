function T  = EpistasisLBC__LoadData( data_file_name  ,  varargin )
%% T = EpistasisLBC__LoadData( data_file_name  ,  varargin )
%  loads data_file_name (eg: 7)
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

% if numel(gcp('nocreate'))==0 , 	parpool local , end
%% load data
% data_file_name can be a segment # or file name
DataDir = '~/Develop/HIS3InterspeciesEpistasis/Data/';
if isnumeric(data_file_name)
    data_file_name = [ num2str(data_file_name) '_all.tab'];
end
% format = 
%  aa_seq	s	len	nogap	stop	middle	nat_lib	nat	lib	size
T = readtable([ DataDir data_file_name ] , 'Format' , '%s%f%d%d%d%d%d%d%d%d%d' ,'FileType','text','Delimiter','\t');

T = T( logical(T.nogap)  , :) ; % can't handle seq's w/varying length
T = T( T.len == mode(T.len) , :); % can't handle seq's w/varying length

if p.Results.NO_STOP_FLAG
    T = T( ~T.stop & ~T.nonsense , :) ;
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
