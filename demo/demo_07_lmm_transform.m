%% Demo 7 - Linear mixed model on log-transformed response
% Strictly positive, right-skewed outcomes strain the constant-variance assumption
% of an ordinary linear model; one remedy is to model them on the log scale. Shown
% on the sleepstudy dataset (lme4) - the same reaction-time data as demo 6. The
% response is fitted in log space and back-transformed for plots and tables:
%
%     log(Reaction) ~ Period + (1 | Subject)
%
% The back-transformed mean is a geometric mean (median-like), which contrasts
% instructively with the gamma GLMM of demo 8.

%% Initialize
thisDir = fileparts(mfilename('fullpath'));
cd(thisDir);
restoredefaultpath;
addpath(genpath(fullfile(thisDir, '..')));

%% Options
options = struct;
options.inFile         = 'data/sleepstudy.csv'; % input data file
options.outDir         = 'results/demo_07_lmm_transform'; % output folder
options.y              = 'Reaction'; % dependent variable
options.yUnits         = 'ms'; % unit label for y-axis
options.transform      = 'log(x)'; % fit on log scale; back-transformed for plots and tables
options.x              = 'Period'; % fixed-effect factor
options.id             = 'Subject'; % random-effect grouping variable
options.isRandomSlopes = false; % random intercept only -> (1 | Subject)
options.rename         = 'Reaction -> ReactionTime; Period -> Day';

kbstat(options);
