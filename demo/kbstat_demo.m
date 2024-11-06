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
options = struct; % Init empty structure
options.inFile = 'Chocolate.csv'; % Input file in long format as CSV table
options.outDir = 'Results'; % Output folder relative to current working directory

% For the linear model, either provide model variables explicitly, or
% provide a Wilkinson formula instead.

% Option 1: Explicit variables (comment-out if you want to use the Wilkinson formula below)
options.y = 'Distance'; % Dependent variable of model
options.yUnits = 'm'; % Units of dependent variable
options.x = 'Chocolate, Gender'; % Fixed-effect variables
options.id = 'Subject'; % Random-effect variable

% % Option 2: Wilkinson formula (if uncommented, overwrites the model definition above)
% options.formula = 'Distance ~ Chocolate*Gender + (1|Subject)';

% The following options are entirely optional.
% See kbstat header for more information.
%
% options.interact = 'Chocolate, Gender'; % Analyse interaction between these variables
% options.distribution = 'gamma'; % Distribution family of GLM
% options.link = 'log'; % Link function of GLM
% options.fitMethod = 'MPL'; % Fitting method of GLM
% options.outlierMethod = 'quartiles'; % Remove outliers using this function 
% options.showVarNames = 'Levels'; % Write levels with capitalized first letter


% Call main script with the given options
kbstat(options);