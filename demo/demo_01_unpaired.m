%% Demo 1 - Unpaired t-test equivalent (single two-level factor)
% The simplest equivalence: a linear model with one two-level factor reproduces
% the classical independent-samples t-test exactly. Shown on the sleep dataset
% (R base; Cushny & Peebles 1905 / Student 1908) - the extra hours of sleep
% gained by 10 patients under two soporific drugs, treated as independent groups.
%
%     extra ~ group

%% Initialize
thisDir = fileparts(mfilename('fullpath'));
cd(thisDir);
restoredefaultpath;
addpath(genpath(fullfile(thisDir, '..'))); % add the kbstat library (parent folder)

%% Options
options = struct;
options.inFile = 'data/sleep.csv'; % input data file
options.outDir = 'results/demo_01_unpaired'; % output folder
options.y      = 'extra'; % dependent variable
options.yUnits = 'h'; % unit label for y-axis
options.x      = 'group'; % fixed-effect factor (no random effects -> plain linear model)
options.rename = 'extra -> ExtraSleep; group -> DrugGroup';

kbstat(options);
