%% Demo 3 - Two-way ANOVA (two crossed factors + interaction)
% Two crossed factors and their interaction are the classic two-way ANOVA, which
% a linear model reproduces exactly on a balanced design. Shown here on the
% ToothGrowth dataset (R base): odontoblast length of 60 guinea pigs given vitamin
% C by two delivery methods - orange juice (OJ) and ascorbic acid (VC) - at three
% dose levels (low, medium, high), 10 animals per cell.
%
% Both factors enter as crossed between-subject effects:
%     len ~ supp * dose

%% Initialize
thisDir = fileparts(mfilename('fullpath'));
cd(thisDir);
restoredefaultpath;
addpath(genpath(fullfile(thisDir, '..'))); % add the kbstat library (parent folder)

%% Options
options = struct;
options.inFile      = 'data/toothgrowth.csv'; % input data file
options.outDir      = 'results/demo_03_twoway'; % output folder
options.y           = 'len'; % dependent variable
options.yUnits      = 'mm'; % unit label for y-axis
options.x           = 'supp, dose'; % fixed-effect factors
options.interact    = 'supp, dose'; % test the supp x dose interaction
options.rename      = 'len -> ToothLength; supp -> Supplement; dose -> Dose; supp: OJ -> orange_juice, VC -> vitamin_c';
options.xOrder2     = 'low, medium, high'; % order the dose levels (2nd factor)

kbstat(options);
