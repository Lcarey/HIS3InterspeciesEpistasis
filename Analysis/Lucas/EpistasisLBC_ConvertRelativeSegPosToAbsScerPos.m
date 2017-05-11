function [ AbsScerPos  , PosTable ] = EpistasisLBC_ConvertRelativeSegPosToAbsScerPos( SegN , RelPos , PosTable )
%% [ AbsScerPos  , PosTable ] = EpistasisLBC_ConvertRelativeSegPosToAbsScerPos( SegN , RelPos , PosTable )
% 
% given a segment and realitive postion (1 - N)
% return the absolute pos according to S cer
%
% LBC May 04, 2017

if ~exist('PosTable','var')
    PosTable = LoadPosTable();
end

AbsScerPos = PosTable(SegN) + RelPos ;

end


function PosTable = LoadPosTable()
P = readtable('~/Develop/HIS3InterspeciesEpistasis/Data_Small_Tables/positions.csv','ReadRowNames',true);
PosTable = NaN(12,1);
for SegN = 1:12
    %library_variable_positions = str2double( regexp( regexprep(char(P{'positions',SegN}),'[\[\] ]','') , ',' ,'split') ) 
    PosTable(SegN) = str2double( char(P{'start_Scer',1}) ) - 1 ;
end
end