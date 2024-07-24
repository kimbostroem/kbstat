function DataLong = table2long(Data, groups, groupNames, groupLevels, levelVar)
%% Convert table to long format by grouping selected variables
%
% SYNTAX
%   DataLong = table2long(Data, groups, groupNames, groupLevels, levelVar)
%
% INPUT
%   Data        (table) Table containing data
%   groups      (cell cell char) Cell array of cell arrays of strings grouping 
%               the variables to obtain the long format 
%   groupNames  (cell char) Cell array of strings that denote each group
%   groupLevels (cell char) Cell array of strings that denote the (common) levels of each group
%   levelVar    (char) String that denotes the new variable holding the group levels
%
% OUTPUT
%   DataLong    (table) Table resulting from rearranging into long format   
%
% EXAMPLE
% DataLong = table2long(Data, {{'a_x', 'a_y', 'a_z'}, {'b_x', 'b_y', 'b_z'}}, {'A', 'B'}, {'x', 'y', 'z'}, 'Component')


DataLong = stack(Data,groups, 'NewDataVariableName',groupNames, 'IndexVariableName',levelVar);
DataLong.(levelVar) = string(categorical(DataLong.(levelVar),unique(DataLong.(levelVar)), groupLevels));

end