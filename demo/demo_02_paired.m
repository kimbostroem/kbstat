%% Demo 2 - Paired t-test equivalent (random intercept per subject)
% When the same subjects are measured under both conditions, accounting for the
% pairing is both more correct and more powerful. Shown on the sleep dataset
% (R base), where both drugs were measured on the same 10 patients. A random
% intercept per patient absorbs each person's baseline sleep tendency:
%
%     extra ~ group + (1 | ID)
%
% This is equivalent to a paired t-test; comparing its p-value with demo 1 shows
% how accounting for the pairing increases power.

%% Initialize
thisDir = fileparts(mfilename('fullpath'));
cd(thisDir);
restoredefaultpath;
addpath(genpath(fullfile(thisDir, '..')));

%% Options
options = struct;
options.inFile         = 'data/sleep.csv'; % input data file
options.outDir         = 'results/demo_02_paired'; % output folder
options.y              = 'extra'; % dependent variable
options.yUnits         = 'h'; % unit label for y-axis
options.x              = 'group'; % fixed-effect factor
options.id             = 'ID'; % random-effect grouping variable (subject)
options.isRandomSlopes = false; % random intercept only -> (1 | ID)
options.rename         = 'extra -> ExtraSleep; group -> DrugGroup';

kbstat(options);
