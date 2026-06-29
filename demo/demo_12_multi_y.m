%% Demo 12 - Multiple dependent variables in one call
% Research questions often span several outcomes at once; kbstat can run the whole
% pipeline for each in a single call. Shown on the iris dataset (R base; Fisher
% 1936): the four flower measurements of 150 plants, with setosa excluded so only
% versicolor and virginica are compared. Passing a comma-separated list to
% options.y runs the full pipeline once per measurement, into per-variable
% subfolders:
%
%     {Sepal_Length, Sepal_Width, Petal_Length, Petal_Width} ~ Species
%
% The constraint demonstrates categorical row selection.

%% Initialize
thisDir = fileparts(mfilename('fullpath'));
cd(thisDir);
restoredefaultpath;
addpath(genpath(fullfile(thisDir, '..')));

%% Options
options = struct;
options.inFile     = 'data/iris.csv'; % input data file
options.outDir     = 'results/demo_12_multi_y'; % output folder (subfolders per variable)
options.y          = 'Sepal_Length, Sepal_Width, Petal_Length, Petal_Width'; % dependent variables
options.yUnits     = 'cm'; % unit label for y-axis (same for all variables)
options.x          = 'Species'; % fixed-effect factor
options.constraint = 'Species ~= "setosa"'; % exclude setosa: compare versicolor vs virginica only
options.rename     = ['Sepal_Length -> SepalLength; Sepal_Width -> SepalWidth; ' ...
                      'Petal_Length -> PetalLength; Petal_Width -> PetalWidth'];

kbstat(options);
