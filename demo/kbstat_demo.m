%% Demo for using the kbstat library
% For more information, see the header of kbstat.m. For most options there
% is a default, so they do not all have to be provided. Here are the
% available options for chosing the distribution and link function.
%
%       distribution    Distribution used for the GLM fit.
%                       OPTIONAL, default = 'Normal'.
%                       Possible values:
%                       'normal'	        Normal distribution
%                       'logNormal'         Normal distribution on log-values
%                       'binomial'	        Binomial distribution
%                       'poisson'	        Poisson distribution
%                       'gamma'	            Gamma distribution
%                       'inverseGaussian'	Inverse Gaussian distribution
%
%       link            Link function used for the GLM fit.
%                       OPTIONAL, default depends on chosen distribution.
%                       Possible values:
%                       'identity'	    g(mu) = mu.             Default for Normal distribution.
%                       'log'	        g(mu) = log(mu).        Default for Poisson, Gamma, and InverseGaussian.
%                       'logit'	        g(mu) = log(mu/(1-mu))  Default for Binomial distribution.
%                       'loglog'	    g(mu) = log(-log(mu))
%                       'probit'	    g(mu) = norminv(mu)
%                       'comploglog'	g(mu) = log(-log(1-mu))
%                       'reciprocal'	g(mu) = mu.^(-1)
%                       Scalar p	    g(mu) = mu.^p           Canonical for Gamma (p = -1) and InverseGaussian (p= -2)

fprintf('Initializing...\n');

% Clear workspace and close open figures
clear
close all

% Restore default path
restoredefaultpath;

% Add library and subfolders to path
addpath(genpath('../../kbstat'));

% Provide options
options = struct;
options.inFile = 'Chocolate.csv';
options.outDir = 'Results';

options.y = 'Distance';
options.yUnits = 'm';

options.x = 'Chocolate, Gender';
options.id = 'Subject';
options.within = 'Chocolate';

options.formula = 'Distance ~ Chocolate*Gender + (Chocolate*Gender|Subject)';

% The following options are optional. If not provided, defaults are taken.
% See kbstat header for more information.
%
% options.interact = 'Chocolate, Gender';
% options.distribution = 'gamma';
% options.link = 'log';
% options.fitMethod = 'MPL';
% options.outlierMethod = 'quartiles';
% options.showVarNames = 'Levels';


% Call main script
kbstat(options);