%% Demo 6 - Linear mixed model with random slopes
% Subjects often differ not just in their baseline but in how strongly a treatment
% affects them - something repeated-measures ANOVA cannot express, but random
% slopes can. Shown on the sleepstudy dataset (lme4; Belenky et al. 2003):
% reaction times of 18 subjects across a rested vs. a sleep-deprived block.
% Each subject gets their own intercept and their own Period slope:
%
%     Reaction ~ Period + (1 + Period | Subject)

%% Initialize
thisDir = fileparts(mfilename('fullpath'));
cd(thisDir);
restoredefaultpath;
addpath(genpath(fullfile(thisDir, '..')));

%% Options
options = struct;
options.inFile         = 'data/sleepstudy.csv'; % input data file
options.outDir         = 'results/demo_06_lmm_slopes'; % output folder
options.y              = 'Reaction'; % dependent variable
options.yUnits         = 'ms'; % unit label for y-axis
options.x              = 'Period'; % fixed-effect factor
options.id             = 'Subject'; % random-effect grouping variable
options.isRandomSlopes = true; % estimate random slopes ...
options.randomSlopes   = 'Period'; % ... each subject gets their own Period slope
options.rename         = 'Reaction -> ReactionTime; Period -> Day';

kbstat(options);
