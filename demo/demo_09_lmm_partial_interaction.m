%% Demo 9 - Partial interaction (test one interaction, not all)
% With three or more factors you often want to test one specific interaction
% without assuming the others - a partial interaction. Shown on the npk dataset
% (R base; Venables & Ripley 2002): pea yield under three binary treatments -
% nitrogen (N), phosphate (P), potassium (K) - across six complete blocks.
% Only the hypothesised P x K interaction is included:
%
%     yield ~ N + P*K + (1 | block)

%% Initialize
thisDir = fileparts(mfilename('fullpath'));
cd(thisDir);
restoredefaultpath;
addpath(genpath(fullfile(thisDir, '..')));

%% Options
options = struct;
options.inFile         = 'data/npk.csv'; % input data file
options.outDir         = 'results/demo_09_lmm_partial_interaction'; % output folder
options.y              = 'yield'; % dependent variable
options.yUnits         = 'lb/plot'; % unit label for y-axis
options.x              = 'N, P, K'; % three binary treatment factors
options.id             = 'block'; % random-effect grouping variable (spatial block)
options.interact       = 'P, K'; % test only the P x K interaction (not N x P or N x K)
options.isRandomSlopes = false; % random intercept only -> (1 | block)
options.rename         = ['yield -> CropYield; N -> Nitrogen; P -> Phosphate; K -> Potassium; block -> Block; ' ...
                          'N: 0 -> absent, 1 -> applied; P: 0 -> absent, 1 -> applied; K: 0 -> absent, 1 -> applied'];

kbstat(options);
