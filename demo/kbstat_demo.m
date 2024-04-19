%% Initilization

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

% clear workspace
clear
close all

% restore default path
restoredefaultpath;

% add library and subfolders to path
addpath(genpath('../../kbstat'));
resultsDir = 'Results';

options = struct;
options.inFile = 'Chocolate.csv';
options.y = 'Distance';
options.yUnits = 'm';
options.x = 'Chocolate, Gender';
options.id = 'Subject';
options.within = 'Chocolate';
options.interact = 'Chocolate, Gender';
options.distribution = 'gamma';
options.fitMethod = 'MPL';
options.posthocMethod = 'emm';
options.posthocLevel = 1;
options.preOutlierMethod = 'quartiles';
options.postOutlierMethod = 'none';
options.showVarNames = 'Levels';
kbstat(options);