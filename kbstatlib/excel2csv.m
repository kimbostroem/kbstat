function outFile = excel2csv(inFile, outFile, sheet)
%% Convert Excel file to CSV file
%
% SYNTAX
%   outFile = excel2csv(inFile, outFile, sheet)
%
% INPUT
%   inFile      (char) Path to input file. Must have suffix ".xls" or 
%               ".xlsx"
%   outFile     (char) Path to output file. If no suffix is given, ".csv" 
%               is used.
%               OPTIONAL. Defaults to the same folder and file name as the
%               input file, except that the suffix is ".csv"
%   sheet       (char) Name of the Excel sheet within the file that is to 
%               be converted.
%               OPTIONAL. If no sheet or empty sheet is given, the 1st
%               sheet is converted.
%
% OUTPUT
%   outFile     (char) Path to the generated output file.
%
%
% EXAMPLES
%   excel2csv('/path/to/infile.xlsx')
%   excel2csv('/path/to/infile.xlsx', '/other/path/to/outfile.csv')
%   excel2csv('/path/to/infile.xlsx', '/other/path/to/outfile.csv', 'someSheet')
%   excel2csv('/path/to/infile.xlsx', [], 'someSheet')
%
% (c) 2022 by Predimo GmbH
% Website: http://www.predimo.com
% Authors: Kim Joris Boström
% Last revision: 2022-09-16 (KB) created file

[fdir, fname, fext] = fileparts(inFile);

if nargin < 3
    sheet = '';
end
if nargin < 2
    outFile = fullfile(fdir, [fname, '.csv']);
end

if ~(isfile(inFile) && any(strcmpi(fext, {'.xls', '.xlsx'})))
    error('Input file must be a valid Excel file with suffix ".xls" or ".xlsx"');
end

if isempty(sheet)
    T = readtable(inFile, 'Format', 'auto');
else
    T = readtable(inFile, 'Format', 'auto', 'Sheet', sheet);
end

[outdir, outname, outext] = fileparts(outFile);
% ensure that outfile suffix is ".csv"
if ~strcmpi(outext, '.csv')
    outFile = fullfile(outdir, [outname, '.csv']);
end
writetable(T, outFile, 'WriteMode', 'overwrite', 'Delimiter', ';');

end