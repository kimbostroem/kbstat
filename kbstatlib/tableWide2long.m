function longTable = tableWide2long(wideTable, groups, groupNames, groupLevels, levelVar)
%% Convert table from wide format to long format by grouping selected variables
%
% SYNTAX
%   tableLong = tableWide2long(tableWide, groups, groupNames, groupLevels, levelVar)
%
% INPUT
%   tableWide   (table) Table containing data
%   groups      (cell cell char) Cell array of cell arrays of strings grouping 
%               the variables to obtain the long format 
%   groupNames  (cell char) Cell array of strings that denote each group
%   groupLevels (cell char) Cell array of strings that denote the (common) levels of each group
%   levelVar    (char) String that denotes the new variable holding the group levels
%
% OUTPUT
%   tableLong   (table) Table resulting from rearranging into long format   
%
% EXAMPLE
% tableLong = tableWide2long(tableWide, {{'a_x', 'a_y', 'a_z'}, {'b_x', 'b_y', 'b_z'}}, {'A', 'B'}, {'x', 'y', 'z'}, 'Component')


longTable = stack(wideTable,groups, 'NewDataVariableName', groupNames, 'IndexVariableName', levelVar);
longTable.(levelVar) = string(categorical(longTable.(levelVar), unique(longTable.(levelVar)), groupLevels));

end