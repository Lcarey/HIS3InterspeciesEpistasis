function EpistasisLBC__Figure_imagesc_all_libraries()
%% EpistasisLBC__Figure_imagesc_all_libraries()
%  
%  For each library, plot imagesc() of the variants 
%   both sorted and not sorted by fitness
%
% LBC March 2017

files_struct = dir('~/Develop/HIS3InterspeciesEpistasis/Data/S*_scaled_info.csv');
for I = 1:numel(files_struct)
    files_struct(I).LibN = 