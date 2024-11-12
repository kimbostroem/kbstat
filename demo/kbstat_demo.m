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

addpath(genpath('../../kbstat')); % Location of kbstat library relative to current working directory

%% Provide basic options

options = struct; % Init empty structure
options.inFile = 'Chocolate.csv'; % Relative path to input file in long format as CSV table
options.outDir = 'Results'; % Relative path to output folder
options.distribution = 'gamma'; % Distribution family of GLM
options.link = 'log'; % Link function of GLM
options.fitMethod = 'MPL'; % Fitting method of GLM
options.outlierMethod = 'quartiles'; % Remove outliers using this method 
options.showVarNames = 'Names_and_levels'; % Display names and levels, capitalize first letter of names

%% Alternative 1: Explicitly define variable types
% Define model variables and their roles explicitly

options.y = 'Distance'; % Dependent variable of model
options.yUnits = 'm'; % Units of dependent variable
options.x = 'Chocolate, Gender'; % Fixed-effect variables
options.id = 'Subject'; % Random-effect variable

%% 
% Call main script with the given options

kbstat(options);

%% Alternative 2: Provide Wilkinson formula
% The Wilkinson formula uniquely defines all variables and their role. If
% provided, it overrides any explicitly defined variables. The following
% Wilkinson formula produces the same results as Alternative 1

options.formula = 'Distance ~ Chocolate*Gender + (1|Subject)';

%% 
% Call main script with the given options

kbstat(options);