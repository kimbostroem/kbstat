%% Demo 10 - Principled outlier removal (pre-fit IQR + post-fit residual)
% Influential outliers can dominate a fit, so removing them on a principled basis
% - and seeing how much they mattered - is part of a careful analysis. Shown on
% the stackloss dataset (R base; Brownlee 1965): the percentage of ammonia lost
% during an industrial oxidation process over 21 plant runs, where observations
% 1, 3, 4, and 21 are textbook influential outliers.
%
%     stack.loss ~ Air.Flow + Water.Temp + Acid.Conc.
%
% The same model is fitted twice - once untouched, once after pre-fit IQR and
% post-fit residual outlier removal - to show how the estimates shift. Outliers
% are set to NaN rather than dropped, so excluded points stay visible and the
% data layout (hence each point's identity) is preserved under the imbalance.

%% Initialize
thisDir = fileparts(mfilename('fullpath'));
cd(thisDir);
restoredefaultpath;
addpath(genpath(fullfile(thisDir, '..')));

%% Load data and bin the continuous Air_Flow into three operating regimes
% kbstat has no automatic binning, so we derive the categorical factor here and
% pass the table in via options.DataRaw. Water_Temp and Acid_Conc_ stay numeric.
T = readtable('data/stackloss.csv', 'Format', 'auto');
T.AirFlow = discretize(T.Air_Flow, 3, 'categorical', {'low', 'medium', 'high'});

%% Common options
options = struct;
options.DataRaw   = T; % pass the in-memory table (no inFile)
options.y         = 'stack_loss'; % dependent variable
options.yUnits    = '%'; % unit label for y-axis
options.x         = 'AirFlow'; % binned operating regime (violin factor)
options.covariate = 'Water_Temp, Acid_Conc_'; % numeric covariates, in model, excluded from plots
options.rename    = 'stack_loss -> StackLoss; AirFlow -> Air flow; Water_Temp -> WaterTemp; Acid_Conc_ -> AcidConc';

%% Run 1: no outlier removal
fprintf('=== Run 1: no outlier removal ===\n');
options.outDir = 'results/demo_10_outliers/default';
options.outlierRemoval = 'none';
kbstat(options);

%% Run 2: pre-fit IQR + post-fit residual removal
fprintf('=== Run 2: pre-fit IQR + post-fit residual removal ===\n');
options.outDir = 'results/demo_10_outliers/clean';
options.outlierRemoval   = 'quartiles'; % pre-fit IQR rule
options.postOutlierMethod = 'quartiles'; % post-fit residual rule
kbstat(options);
