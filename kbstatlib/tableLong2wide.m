function wideTable = tableLong2wide(longTable, factorColumns, valueColumn, idColumn, newVarPrefix)
%% Convert a long-format table to wide format with multiple factors.
% If 1st argument is a file path, the table is imported from there.
%
% INPUTS:
%   longTable       - Input table in long format, or file path to such table.
%   factorColumns   - Cell array of column names (strings) to use as grouping factors for wide headers.
%   valueColumn     - The name (as a string) of the column containing the values.
%   idColumn        - The name (as a string) of the column to use as row identifiers in the wide table.
%   newVarPrefix    - Prefix for the name of the new column header
%                     OPTIONAL, default = ''
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

if nargin < 5
    newVarPrefix = '';
end

isPath = false;
if (isstring(longTable)||ischar(longTable)) && isfile(longTable)
    isPath = true;
    fpath = longTable;
    [fdir, fname, fext] = fileparts(fpath);
    longTable = readtable(fpath);
end

factorColumns = cellstr(factorColumns);

% Generate unique combinations of factor levels
uniqueFactors = unique(longTable(:, factorColumns), 'rows');
nWide = size(uniqueFactors, 1);

% Remaining column headers and their factors
remHeaders = string(setdiff(longTable.Properties.VariableNames, [factorColumns, valueColumn, idColumn]));
remTypes = varfun(@class, longTable(:, remHeaders), 'OutputFormat', 'cell');

% Create wide-format column headers from the unique factor combinations
wideHeaders = strings(1, size(uniqueFactors, 1));
for iWide = 1:nWide
    if ~isempty(newVarPrefix)
        wideHeaders(iWide) = strjoin([newVarPrefix, string(uniqueFactors{iWide, :})], '_');
    else
        wideHeaders(iWide) = strjoin(string(uniqueFactors{iWide, :}), '_');
    end
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
wideTypes = varfun(@class, wideTable, 'OutputFormat', 'cell');
for iWide = 2:width(wideTable)
    if isnumeric(wideTypes{iWide})
        wideTable{:, iWide} = NaN;
    end
end

% Populate the table
for iID = 1:numel(uniqueIDs)
    id = string(uniqueIDs(iID));
    for iWide = 1:numel(wideHeaders)
        mask = ismember(longTable(:, factorColumns), uniqueFactors(iWide, :), 'rows');
        valueMask = mask & (longTable.(idColumn) == id);
        if any(valueMask)
            entries = longTable{valueMask, valueColumn};
            nEntries = length(entries);
            if nEntries > 1
                fprintf('Subject %s has %d entries for %s -> averaging\n', string(id), nEntries, wideHeaders(iWide));
                entry = mean(entries);
            elseif nEntries == 0
                fprintf('Subject %s has missing data for %s\n', string(id), wideHeaders(iWide));
                entry = NaN;
            else
                entry = entries;
            end
            wideTable{iID, wideHeaders(iWide)} = entry;
            for iRem = 1:numel(remHeaders)
                wideTable{iID, remHeaders(iRem)} = longTable{valueMask, remHeaders(iRem)};
            end
        end
    end
end

% write table to disk if desired
if isPath
    myTable = wideTable;
    myFname = sprintf('%s_wide', fname);
    myFext = fext;
    switch myFext
        case '.xlsx'
            % save as Excel table
            outpath = fullfile(fdir, [myFname, myFext]);
            writetable(myTable, outpath, 'WriteMode', 'replacefile');
        case '.csv'
            % save as CSV table
            outpath = fullfile(fdir, [myFname, myFext]);
            writetable(myTable, outpath, 'WriteMode', 'overwrite', 'Delimiter', ';', 'QuoteStrings', 'all');
        otherwise
            error('Unknown format %s', myFext);
    end
end
end
