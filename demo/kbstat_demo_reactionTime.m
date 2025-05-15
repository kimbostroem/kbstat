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


%% Initialize

fprintf('Initializing...\n');

% Clear workspace and close open figures
clear
close all

% Restore default path
restoredefaultpath;

% Add library and subfolders to path
addpath(genpath('../../kbstat')); % Location of kbstat library relative to current working directory

%% Define method of analysis

% 1 = use link function, 2 = use variable transformation
analyzeMethod = 1;

%% Provide basic options

options = struct; % Init empty structure
options.inFile = 'reaction_time.csv'; % Relative path to input file in long format as CSV table
options.outDir = 'Results_rt'; % Relative path to output folder
options.fitMethod = 'MPL'; % Fitting method of GLM
options.outlierMethod = 'quartiles'; % Remove outliers using this method
options.showVarNames = 'Names_and_levels'; % Display names and levels, capitalize first letter of names

% Define model variables and their roles explicitly
options.y = 'rt'; % Dependent variable of model
options.yUnits = 'ms'; % Units of dependent variable
options.x = 'A, B'; % Fixed-effect variables
options.id = 'id'; % Random-effect variable

%% Optional: Provide explicit Wilkinson formula

% The Wilkinson formula uniquely defines all variables and their role. If
% provided, it overrides any explicitly defined variables. The following
% Wilkinson formula produces the same results as Alternative 1

options.formula = 'rt ~ A*B + (A+B|id)';

switch analyzeMethod

    case 1 % Alternative 1 (recommmended): Use link function
        fprintf('Analyze using link function\n');
        options.distribution = 'gamma'; % Distribution family of GLM
        options.link = 'log'; % Link function of GLM
    
    case 2 % Alternative 2: Use data transformation
        fprintf('Analyze using data transformation\n');
        options.transform = 'log(1+x)'; % transformation of dependent variable
        options.distribution = 'gamma'; % Distribution family of GLM
        options.link = 'log'; % Link function of GLM
end

% Perform analysis with the given options
kbstat(options);
