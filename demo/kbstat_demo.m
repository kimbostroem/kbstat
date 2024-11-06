%% Demo for using the kbstat library
% For more information, see the header of kbstat.m. For most options there is 
% a default, so they do not all have to be provided. Here are the available options 
% for chosing the distribution and link function.
%%
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
%%
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
%

fprintf('Initializing...\n');

%% Clear workspace and close open figures

clear
close all

%% Restore default path

restoredefaultpath;

%% Add library and subfolders to path

addpath(genpath('../../kbstat'));

%% Provide basic options

options = struct; % Init empty structure
options.inFile = 'Chocolate.csv'; % Input file in long format as CSV table
options.outDir = 'Results'; % Output folder relative to current working directory
options.distribution = 'gamma'; % Distribution family of GLM
options.link = 'log'; % Link function of GLM
options.fitMethod = 'MPL'; % Fitting method of GLM
options.outlierMethod = 'quartiles'; % Remove outliers using this function 
options.showVarNames = 'Levels'; % Write levels with capitalized first letter

%% Alternative 1: Explicitly define variable types
% Define model variables and their roles explicitly

options.y = 'Distance'; % Dependent variable of model
options.yUnits = 'm'; % Units of dependent variable
options.x = 'Chocolate, Gender'; % Fixed-effect variables
options.interact = 'Chocolate, Gender'; % Analyse interaction between these variables
options.id = 'Subject'; % Random-effect variable
%% 
% Call main script with the given options

kbstat(options);

%% Alternative 2: Provide Wilkinson formula
% The Wilkinson formula uniquely defines all variables and their role. If it 
% is provided, then it overrides any explicitly defined variables.

options.formula = 'Distance ~ Chocolate*Gender + (1|Subject)';
options.yUnits = 'm'; % Units of dependent variable
%% 
% Call main script with the given options

kbstat(options);