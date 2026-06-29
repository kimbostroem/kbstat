%% Demo 4 - One-way repeated-measures ANOVA (random intercept)
% A within-subject factor handled by a random intercept reproduces the one-way
% repeated-measures ANOVA. Shown on the ergoStool dataset (nlme; Wretenberg et al.
% 1993): 9 subjects each rate the perceived effort (Borg scale) of rising from
% four stool types - one rating per subject per stool, balanced and unreplicated.
%
%     effort ~ Type + (1 | Subject)

%% Initialize
thisDir = fileparts(mfilename('fullpath'));
cd(thisDir);
restoredefaultpath;
addpath(genpath(fullfile(thisDir, '..')));

%% Options
options = struct;
options.inFile         = 'data/ergostool.csv'; % input data file
options.outDir         = 'results/demo_04_lmm'; % output folder
options.y              = 'effort'; % dependent variable
options.yUnits         = 'Borg'; % unit label for y-axis
options.x              = 'Type'; % within-subject fixed-effect factor
options.id             = 'Subject'; % random-effect grouping variable
options.isRandomSlopes = false; % random intercept only -> (1 | Subject)
options.rename         = 'effort -> Effort; Type -> Stool type';

kbstat(options);
