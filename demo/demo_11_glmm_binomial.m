%% Demo 11 - Binomial GLMM (binary outcome, logit link)
% A binary outcome has a mean that is a probability bounded between 0 and 1, so a
% gaussian model is invalid - a binomial GLMM with a logit link is the right tool.
% Shown on the bacteria dataset (MASS): presence/absence of H. influenzae in 50
% children measured at five time points under three treatments (220 observations).
% With a random intercept per child for the repeated measures:
%
%     present ~ trt + week + (1 | ID)
%
% Estimated marginal means are back-transformed from the logit scale to
% probabilities - a binary outcome that classical t-tests and ANOVA cannot handle.

%% Initialize
thisDir = fileparts(mfilename('fullpath'));
cd(thisDir);
restoredefaultpath;
addpath(genpath(fullfile(thisDir, '..')));

%% Load data and order the weeks numerically for the plot
% week is numeric, so it cannot be ordered via xOrder (numeric order strings are
% read as indices). Instead we sort the rows by week and let levelOrder = 'stable'
% follow that appearance order, giving weeks 0, 2, 4, 6, 11.
T = readtable('data/bacteria.csv', 'Format', 'auto');
T.week = str2double(string(T.week));
T = sortrows(T, 'week');

%% Options
options = struct;
options.DataRaw        = T; % pass the sorted in-memory table (no inFile)
options.outDir         = 'results/demo_11_glmm_binomial'; % output folder
options.y              = 'present'; % binary outcome: 1 = bacteria present
options.x              = 'trt, week'; % fixed-effect factors
options.id             = 'ID'; % random intercept per child
options.isRandomSlopes = false; % random intercept only -> (1 | ID)
options.distribution   = 'binomial'; % binary outcome
options.link           = 'logit'; % logit link: maps linear predictor to probability
options.levelOrder     = 'stable'; % follow the sorted data order (weeks 0, 2, 4, 6, 11)
options.rename         = 'present -> Bacteria present; trt -> Treatment; week -> Week';

kbstat(options);
