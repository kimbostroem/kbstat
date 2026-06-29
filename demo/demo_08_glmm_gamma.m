%% Demo 8 - Generalized linear mixed model with gamma distribution
% Strictly positive, right-skewed outcomes can be modelled directly with a
% distribution that fits, instead of being forced into a normal model. Shown on
% the Oats dataset (nlme; Yates 1935), a split-plot field trial: oat yield for
% three varieties under four nitrogen levels across six blocks (72 plots).
% A gamma GLMM with a log link and a random intercept per block:
%
%     yield ~ Variety * Nitrogen + (1 | Block)
%
% The gamma family captures the mean-variance relationship of skewed positive data
% directly, and the log link keeps fitted yields positive. Compare with demo 7,
% which targets similar data by log-transforming a normal model instead.

%% Initialize
thisDir = fileparts(mfilename('fullpath'));
cd(thisDir);
restoredefaultpath;
addpath(genpath(fullfile(thisDir, '..')));

%% Options
options = struct;
options.inFile         = 'data/oats.csv'; % input data file
options.outDir         = 'results/demo_08_glmm_gamma'; % output folder
options.y              = 'yield'; % dependent variable
options.yUnits         = 'qt/plot'; % unit label for y-axis (quarter-pounds per plot)
options.x              = 'Variety, Nitrogen'; % fixed-effect factors
options.id             = 'Block'; % random-effect grouping variable
options.interact       = 'Variety, Nitrogen'; % test the Variety x Nitrogen interaction
options.isRandomSlopes = false; % random intercept only -> (1 | Block)
options.distribution   = 'gamma'; % gamma GLMM (positive, right-skewed outcome)
options.link           = 'log'; % log link function
options.rename         = 'yield -> CropYield';

kbstat(options);
