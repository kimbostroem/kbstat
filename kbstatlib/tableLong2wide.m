function wideTable = tableLong2wide(longTable, factorColumns, valueColumn, idColumn)
%% Convert a long-format table to wide format with multiple factors.
%
% INPUTS:
%   longTable     - Input table in long format.
%   factorColumns - Cell array of column names (strings) to use as grouping factors for wide headers.
%   valueColumn   - The name (as a string) of the column containing the values.
%   idColumn      - The name (as a string) of the column to use as row identifiers in the wide table.
%
% OUTPUT:
%   wideTable     - Output table in wide format.
%
% EXAMPLE:
%   % Example long-format table
%   longTable = table([1; 1; 2; 2], {'A'; 'A'; 'B'; 'B'}, {'X'; 'Y'; 'X'; 'Y'}, [10; 20; 30; 40], ...
%                     'VariableNames', {'ID', 'Category', 'SubCategory', 'Value'});
%   % Convert to wide format
%   wideTable = tableLong2wide(longTable, {'Category', 'SubCategory'}, 'Value', 'ID');

% Generate unique combinations of factor levels
uniqueFactors = unique(longTable(:, factorColumns), 'rows');
nWide = size(uniqueFactors, 1);

% Remaining column headers and their factors
remHeaders = string(setdiff(longTable.Properties.VariableNames, [factorColumns, valueColumn, idColumn]));
remTypes = varfun(@class,longTable(:, remHeaders),'OutputFormat','cell');

% Create wide-format column headers from the unique factor combinations
wideHeaders = strings(1, size(uniqueFactors, 1));
for iWide = 1:nWide
    wideHeaders(iWide) = strjoin([valueColumn, string(uniqueFactors{iWide, :})], '_');
end

% Get unique identifiers for rows
uniqueIDs = unique(longTable.(idColumn));

allHeaders = [remHeaders, wideHeaders];
nHeaders = numel(allHeaders);

% Initialize the wide-format table
wideTable = table('Size', [numel(uniqueIDs), nHeaders + 1], ...
                  'VariableTypes', [{'double'}, remTypes, repmat({'double'}, 1, nWide)], ...
                  'VariableNames', [{idColumn}, cellstr(allHeaders)]);

% Assign IDs to the first column
wideTable.(idColumn) = uniqueIDs;

% Initialize all data columns with NaN
for iWide = 2:width(wideTable)
    wideTable{:, iWide} = NaN;
end

% Populate the table
for iID = 1:numel(uniqueIDs)
    id = uniqueIDs(iID);
    for iWide = 1:numel(wideHeaders)
        mask = ismember(longTable(:, factorColumns), uniqueFactors(iWide, :), 'rows');
        valueMask = mask & (longTable.(idColumn) == id);
        if any(valueMask)
            wideTable{iID, wideHeaders(iWide)} = longTable{valueMask, valueColumn};
            for iRem = 1:numel(remHeaders)
                wideTable{iID, remHeaders(iRem)} = longTable{valueMask, remHeaders(iRem)};
            end
        end
    end
end
end
