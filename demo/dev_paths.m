fprintf('Reset paths...\n');

restoredefaultpath; % restore default path

% add paths to full model
oldPaths = strsplit(path, pathsep);
addpath(genpath('library'));
paths = strsplit(path, pathsep);
for iPath = 1:(length(paths)-length(oldPaths))
    myPath = strrep(paths{iPath},'\','\\'); % escape backslashes
    msg = sprintf('Added path %s\n', myPath);
    fprintf(msg);
end

clear oldPaths paths iPath myPath msg RESTOREDEFAULTPATH_EXECUTED