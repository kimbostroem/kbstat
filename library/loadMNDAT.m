function data = loadMNDAT(file)
[fdir, fname, fext] = fileparts(file); % parse file parts
if isempty(fext)
    fext = '.mndat';
end
fpath = fullfile(fdir, [fname, fext]);
tmpDir = fileparts(tempname);
fpaths = unzip(fpath, tmpDir); % unzip to JSON file
data = fileread(fpaths{1}); % read JSON file into memory
data = jsondecode(data); % decode JSON data
delete(fpaths{1}); % delete JSON file

end