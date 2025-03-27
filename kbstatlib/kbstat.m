function results = kbstat(options)
%% Perform statistical analysis based on a generalized linear mixed model.
% The result is a model fit summary, diagnostic plots, a data bar plot, a
% descriptive statistics table, an ANOVA table, and a posthoc pairwise
% comparison table for selected variables. For more information see the
% README file and the demo subfolder.
%
% REMARK: This function can do multivariate analysis which can be done in
% either one of two possible ways: 1) Declare a dependent variable
% options.y as a cell array with more than 1 entry, or 2) declare
% options.multiVar as the variable whose levels encode the individual
% components of the multivariate variable indicated by options.y. By
% default, the multivariate analysis is split up into several independent
% univariate analyses, because this works in a more stable fashion. Note,
% however, that the results are not statistically corrected for multiple
% testing. To enable true multivariate analysis, set
% options.multiVariate=true. However, the posthoc test then only works
% when options.posthocMethod='ttest' instead of 'emm', due to limitations
% of the external package 'emmeans'.
%
% DISCLAIMER: This library come without any warranty and is for free use
% only. See the LICENSE file in the root folder. The library uses an
% adaptation of Bastian Bechtold's Violin plot library, available at
% https://github.com/bastibe/Violinplot-Matlab, for plotting and John
% Hartman's emmeans library, available at
% https://de.mathworks.com/matlabcentral/fileexchange/71970-emmeans, for
% post-hoc analysis.
%
% SYNTAX
%   kbstat(in, options)
%
% INPUT
%   options     (struct) Structure containing parameters needed for
%               analysis. All fields must be char or string, and if a list
%               is to be given, it must be a string of words separated by
%               comma or semicolon.
%
%       inFile          (char) Path to an Excel or CSV file. The table
%                       must be in long format, so one row per data point.
%
%       DataRaw         (char) Matlab table. must be in long format, so
%                       one row per data point.
%
%       Data            (char) Matlab table. must be in long format, so
%                       one row per data point.
%
%       y               Name of the dependent variable or comma-separated
%                       list of names of multiple dependent variables. In
%                       the latter case, by default each component is
%                       treated as a separate dependent variable and as
%                       many independent univariate analyses are performed.
%                       If "multiVariate" is set to "true", then y is
%                       treated as a vector-valued dependent variable and a
%                       multivariate analysis is performed.
%
%       yName           Display name of dependent variable appearing in
%                       tables and plots
%                       OPTIONAL, default = y
%
%       yLabel          Label for y axis in data plots.
%                       If y has more than one component, then this
%                       parameter can have either one component, in which
%                       case all units are equal to this component, or it
%                       must have the same number of components as y.
%                       OPTIONAL, default = y
%
%       yUnits          Physical units of the dependent variable.
%                       If y has more than one component, then this
%                       parameter can have either one component, in which
%                       case all units are equal to this component, or it
%                       must have the same number of components as y.
%                       OPTIONAL, default = ''.
%
%       yMult           Factor to multiply (each dimension of) y
%                       OPTIONAL, default = 1
%
%       formula         Formula in Wilkinson Notation. If given, it
%                       overrides the automatically created formula. The
%                       script tries to identify relevant variables from
%                       the given formula, including 'x', 'y', 'subject', when
%                       they are not provided.
%
%       x               Comma-separated list of categorical independent variables to be included 
%                       in the analysis as factors.
%
%       xName           Display name of independent variable(s) appearing
%                       in tables and plots
%                       OPTIONAL, default = x
%
%       subject         Name of the subject variable.
%                       OPTIONAL, default = ''.
%
%       id              Same as 'subject'.
%
%       trial           Name of the trial variable nested with subjects.
%                       OPTIONAL, default = ''.
%
%       nested          Same as 'trial'.
%
%       stimulus        Name of stimulus variable crossed with subjects.
%                       OPTIONAL, default = ''.
%
%       crossed         Same as 'stimulus'.
%
%       interact        Comma-separated list of variables whose interaction
%                       with each other is to be analyzed. Can be a subset
%                       of "x", or else its members are added to "x". When
%                       not set, all members of x are assumed to mutually
%                       interact.
%                       Example:
%                       options.subject = 'subject'
%                       options.x = 'time, dose'
%                       options.interact = 'dose, age'.
%                       OPTIONAL, default = options.x
%
%       covariate       Comma-separated list of (continuous or categorical)
%                       variables that are added as covariates to improve
%                       the model without being represented in the plots or
%                       included in the posthoc comparisons. They are,
%                       however, included in the ANOVA table.
%
%       multiVar        Name of the variable that encodes levels of a
%                       multivariate dependent variable.
%                       OPTIONAL, default = ''.
%
%       multiVarLevels  Comma-separated list of the levels of multiVar that
%                       will be taken into account. If empty or undefined,
%                       all levels are taken.
%                       OPTIONAL, default = ''.
%
%       multiVariate    Flag if, when the dependent variable has multiple
%                       components, these components should be analyzed
%                       together. If set to false, they are treated as
%                       independent variables instead and the results are
%                       statistically corrected for these multiple tests.
%                       OPTIONAL, default = false.
%
%       isRandomSlopes  Flag if random slopes should be estimated.
%                       OPTIONAL, default = true.
%
%       isRandomInteract Flag if the slopes of the interactions should also
%                       be estimated. If isRandomSlopes = false, this
%                       parameter is ignored.
%                       OPTIONAL, default = false.
%
%       randomSlopes    Comma-separated list of random slope variables
%                       OPTIONAL, default = all independent variables
%
%       fitMethod       Fit method used for the GLM fit.%
%                       Possible values:
%                       'none'                  Skip linear model fit
%                       'MPL'                   Maximum pseudo likelihood
%                       'REMPL'                 Restricted maximum pseudo likelihood
%                       'Laplace'               Maximum likelihood using Laplace approximation
%                       'ApproximateLaplace'    Maximum likelihood using approximate Laplace approximation with fixed effects profiled out
%                       OPTIONAL, default = 'MPL'.
%
%       distribution    Distribution used for the GLM fit.
%                       If y has more than one component, then this
%                       parameter can have either one component, in which
%                       case all units are equal to this component, or it
%                       must have the same number of components as y.
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
%                       If y has more than one component, then link can
%                       have either one component, in which case all units
%                       are equal to this component, or it must have the
%                       same number of components as y.
%                       Possible values:
%                       'identity'	    g(mu) = mu.             Default for Normal distribution
%                       'log'	        g(mu) = log(mu).        Default for Poisson
%                       'logit'	        g(mu) = log(mu/(1-mu))  Default for Binomial distribution
%                       'loglog'	    g(mu) = log(-log(mu))
%                       'probit'	    g(mu) = norminv(mu)
%                       'comploglog'	g(mu) = log(-log(1-mu))
%                       'reciprocal'	g(mu) = mu.^(-1)        Default for Gamma
%                       Scalar p	    g(mu) = mu.^p           Default for InverseGaussian (p= -2)
%                       'auto'          depends on chosen distribution
%                       OPTIONAL, default = 'auto'.
%
%       dummyCoding     Coding of the dummy variable.
%                       Possible values:
%                       'effects'       uses –1 to represent the last category
%                       'full'          one dummy variable for each category
%                       'reference'     first category is reference
%                       OPTIONAL, default = 'effects'
%
%       posthocMethod   Method for the posthoc pairwise comparison.
%                       Possible values:
%                       'none'      Do not perform posthoc analysis
%                       'ttest'     t-test with Holm-Bonferroni correction
%                       'utest'     Mann-Whitney u-Test (Wilcoxon ranksum test)
%                                   with Holm-Bonferroni correction
%                       'emm'       Extract contrasts from linear model fit
%                       'auto'      Perform posthoc analysis using
%                                   'emm' if univariate or multi-valued y
%                                   with multiVariate=false, 'ttest' if
%                                   distribution='normal', and 'utest'
%                                   otherwise.
%                       OPTIONAL, default = 'emm'.
%
%       posthocCorrection   Method to correct the multiple comparisons
%                       Possible values:
%                       'none'      Do not correct
%                       'holm'      Holm-Bonferroni correction
%                       'bonf',     Same as 'sidak'.
%                       'sidak'     Sidak is approx. equal to
%                                   Bonferroni for small p and cares
%                                   for p not exceeding 1.
%                       OPTIONAL, default = 'holm'.
%
%       posthocMain  Flag if also the posthoc main effects should be
%                       calculated, i.e. the comparison between one
%                       variable set to 'any'.
%                       OPTIONAL, default = false.
%
%       posthocLevel    Indiciate to which level of the independent
%                       variables the posthoc comparison should be
%                       calculated. For example, level 1 means that only
%                       level pairs of the 1st independent variable are
%                       compared, level 2 means that the level pairs of the
%                       1st and 2nd dependent variable are compared.
%                       OPTIONAL, default = 1.
%
%       correctForN     Number of test to be statistically corrected
%                       for in addition to the correction already performed
%                       within the current analysis. Most probably, this
%                       number corresponds to the number of calls of
%                       'kbstat' for one dataset.
%
%       rescale         Flag if the y-axis of each
%                       panel of the data plot is to be resized to a
%                       common scale.
%                       OPTIONAL, default = true.
%
%       outlierRemoval  Method to remove pre-fit outliers from the data.
%                       Possible values:
%                       'none'      Do not remove outliers
%                       'quartiles' Remove values outside 1.5 times the
%                                   interquartile range [.25, .75]
%                       'median'    Remove values outside more than three
%                                   scaled median absolute deviations (MAD)
%                       'mean'      Remove values outside 3 standard
%                                   deviations from the mean.
%                       OPTIONAL, default = 'none'.
%
%       postOutlierRemoval  Method to remove post-fit outliers, i.e.,
%                       outliers in the residuals.
%                       Possible values: see preOutlierMethod
%                       OPTIONAL, default = 'none'.
%
%       constraint      One or more restrictive constraints on the data before analysis.
%						Must be of the form
%						'con1 & con2 & ...'
%						where con1, con2, ... are conditions of the form
%                       'variable == value'
%                       'variable ~= value'
%                       'variable < value'
%                       'variable <= value'
%                       'variable > value'
%                       'variable >= value'
%                       where
%                       "variable" is an actual variable in the data table,
%                       and "variable" is either an existing category, in
%                       the case of categorical variables, or a numeric
%                       value. In the case of a category, the value must be put in
%						double quotes, as in 'bla = "bli"'.
%                       OPTIONAL, default = ''.
%
%       transform       Choose how to transform the dependent variable.
%                       The linear model is fit on the transformed data,
%                       but the data are plotted using the original data.
%                       OPTIONAL, default = ''.
%                       Possible values:
%                       'mean'      Divide by mean
%                       'std'       Divide by standard deviation
%                       'Z'         Z-transform to zero-centered
%                                   distribution with unit standard
%                                   derivation.
%                       'q50' or
%                       'median'
%                       'q%d' or
%                       'p%d'
%                       'q%dq%d' or
%                       'p%dp%d'
%                       'IQR', or   Rescale to interval [0,1] using upper
%                       'quartiles' and lower limits of the inter-quartile
%                                   range (IQR).
%                       'IQRmax'
%                       'MAD'
%                       'MADmax'
%                       'max'
%                       'minmax'
%                       'f(x)'      Arbitrary function of x. Examples:
%                                   'log(x)'
%                                   'atanh(x)'
%
%       outDir          Output folder for generated files, relative to the 
%                       working directory where the main script is called.
%                       If not defined or set to empty, no output files are
%                       generated. 
%                       OPTIONAL, default = '' (nothing is saved)                      
%
%       isPlot          Plot data as grouped bars with significance brackets
%
%       barType         What the bars and error bars indicate.
%                       'auto'   Bar height and error bars are defined by posthocMethod
%                       'mean'   Bar height is mean, error bars are standard deviation
%                       'meanSE' Bar height is mean, error bars are standard error
%                       'meanCI' Bar height is mean, error bars are 95% confidence interval
%                       'median' Bar height is median, error bars are
%                                interquartile range, i.e. the .25 and .75
%                                quantiles.
%                       OPTIONAL, default = 'auto'
%
%       levelOrder      The order in which the levels are displayed in the
%                       plots.
%                       'sorted'    sorted alphanumerically
%                       'stable'    sorted in the order of occurrence in
%                                   the data table
%                       'ascend'    sorted by ascending value
%                       'descend'   sorted by descending value
%                       OPTIONAL, default = 'sorted'
%
%       plotStyle       The style with which the data are plotted.
%                       Possible values:
%                       'violin'
%                       'boxplot'
%                       'bar'
%                       'prettybar'
%                       OPTIONAL, default = 'violin'
%
%       markerSize      Size of the data dots in violin plots.
%                       OPTIONAL, default = NaN (use default marker size).
%
%       title           Title for plots
%
%       showVarNames    Flag how variable names and levels are displayed in data
%                       plots.
%                       'levels'              = display levels as given in the data
%                       'names_levels'    = display names and levels as given in the data
%                       'Levels'              = display capitalized levels only
%                       'Names_levels'    = display names and levels, capitalize names
%                       'names_Levels'    = display names and levels, capitalize levels
%                       'Names_Levels'    = display names and levels, capitalize both
%                       OPTIONAL, default = 'levels'.
%
%       xOrder          Ordering of the items of the 1st x-variable on the 
%                       x-axis in data plots, based on the order of the
%                       level names defined by 'levelOrder'. 
%                       OPTIONAL, default = []. 
%                       Example: options.xOrder = '1 3 2'
%
%       xOrder<n>       Ordering of the items of the n-th x-variable (with 
%                       respect to the provided list) on the x-axis in data
%                       plots, based on the order of the level names
%                       defined by 'levelOrder'.
%                       OPTIONAL, default = []. 
%                       Example: 
%                       options.xOrder1 = '1 4 3 2'
%                       options.xOrder2 = '2 1'
%
%       plotLines      Flag if the data plots should display the median as
%                       a horizontal line in the color of the corresponding
%                       dataset.
%
%       closeFigures    Flag if the figures created during the analysis
%                       should be closed or left open.
%                       OPTIONAL, default = false
%
% OUTPUT
%       mdl             (Generalized linear mixed-effects model) Result
%                       from linear model fit
%
% EXAMPLE (for more info see demo files)
% options = struct; % Init empty structure
% options.inFile = 'Chocolate.csv'; % Input file in long format as CSV table
% options.outDir = 'Results'; % Output folder relative to current working directory
% options.y = 'Distance'; % Dependent variable of model
% options.yUnits = 'm'; % Units of dependent variable
% options.x = 'Chocolate, Gender'; % Independent variables of model
% options.id = 'Subject'; % Random variable of model
%
% Version 1.0, available at https://github.com/kimbostroem/kbstat
% (c) 2024 by Kim Joris Boström
% Website: http://www.kim-bostroem.de

%% Init

% Initialize output argument
results = struct;

%% Switch off specific warning

% The following warning occurs when distributions such as 'gamma' are
% chosen, which require the data to be strictly nonnegative, even if they
% are in fact nonnegative, making the warning pointless.
%
% Warning: The 'Reciprocal' and 'Power' links require the linear predictor
% to be non-negative. However, the model assumes that the linear predictor
% is unconstrained.

warning('off','stats:classreg:regr:lmeutils:StandardGeneralizedLinearMixedModel:BadDistLinkCombination1');

%% Get parameters

if ~isstruct(options)
    error('Input argument must be a struct');
end

% tolerance for differences between real numbers
tol = 1e-12;

% load Data
if isfield(options, 'inFile') && isfile(options.inFile)
    inFile = options.inFile;
    DataRaw = readtable(inFile, 'Format', 'auto');
elseif isfield(options, 'DataRaw')
    DataRaw = options.DataRaw;
elseif isfield(options, 'Data')
    DataRaw = options.Data;
else
    error('Please provide either a valid path to a table file (Excel or CSV) in options.inFile, or a data table in options.Data or options.DataRaw');
end

% store raw data
options.DataRaw = DataRaw;

% init 1st Data set
Data1 = DataRaw;

% all table variables
tableVars = DataRaw.Properties.VariableNames;

% convert all character cell arrays to string arrays
for iVar = 1:length(tableVars)
    tableVar = tableVars{iVar};
    if iscell(Data1.(tableVar))
        Data1.(tableVar) = string(Data1.(tableVar));
    end
end

% independent variable(s)
if isfield(options, 'x') && ~isempty(options.x)
    x = getList(options.x);
else
    fprintf('No independent variable(s) specified\n');
    x = '';
end

%% dependent variable(s)

if isfield(options, 'y') && ~isempty(options.y)
    y = cellstr(options.y);
    nY = length(y);
    for iY = 1:nY
        myY = y{iY};
        if ~ismember(myY, tableVars)
            error('Dependent variable "%s" not found in data table', myY);
        end
    end
else
    fprintf('No dependent variable(s) specified\n');
    y = '';
end

% multiVar
if isfield(options, 'multiVar') && ~isempty(options.multiVar)
    multiVar = options.multiVar;
else
    multiVar = '';
end

% multiVarLevels
if isfield(options, 'multiVarLevels') && ~isempty(options.multiVarLevels)
    multiVarLevels = getList(options.multiVarLevels);
else
    multiVarLevels = '';
end

%% more variables

% correctForN
if isfield(options, 'correctForN') && ~isempty(options.correctForN)
    correctForN = getValue(options.correctForN);
else
    correctForN = 1;
end

% transform
if isfield(options, 'transform') && ~isempty(options.transform)
    transform = options.transform;
else
    transform = '';
end

% subject variable for random effects
if isfield(options, 'subject') && ~isempty(options.subject)
    subject = options.subject;
elseif isfield(options, 'id') && ~isempty(options.id)
    subject = options.id;
else
    subject = '';
end

% trial variable nested in subjects
if isfield(options, 'trial') && ~isempty(options.trial)
    trial = options.trial;
elseif isfield(options, 'nested') && ~isempty(options.nested)
    trial = options.nested;
else
    trial = '';
end

% stimulus variable crossed with subjects
if isfield(options, 'stimulus') && ~isempty(options.stimulus)
    stimulus = options.stimulus;
elseif isfield(options, 'crossed') && ~isempty(options.crossed)
    stimulus = options.crossed;
else
    stimulus = '';
end

% interaction variables
if isfield(options, 'interact') && ~isempty(options.interact) % option provided and not empty
    interact = getList(options.interact);
    x = union(x, interact, 'stable'); % add interaction variables to independent variables
elseif isfield(options, 'interact') && isempty(options.interact) % option provided as empty
    interact = {};
else % option not provided
    interact = x;
end

% covariates
if isfield(options, 'covariate') && ~isempty(options.covariate) % option provided and not empty
    covariate = getList(options.covariate);
else % option not provided
    covariate = {};
end

% random slopes
if isfield(options, 'randomSlopes') && ~isempty(options.randomSlopes) % option provided and not empty
    randomSlopes = getList(options.randomSlopes);
elseif isfield(options, 'randomSlopes') && isempty(options.randomSlopes) % option provided as empty
    randomSlopes = {};
else % option not provided -> use all IVs except covariates
    randomSlopes = x;
end

% estimate random slopes?
if isfield(options, 'isRandomSlopes') && ~isempty(options.isRandomSlopes)
    isRandomSlopes = getValue(options.isRandomSlopes);
else
    isRandomSlopes = true;
end

% isRandomInteract
if isfield(options, 'isRandomInteract') && ~isempty(options.isRandomInteract)
    isRandomInteract = getValue(options.isRandomInteract);
else
    isRandomInteract = false;
end

% extract variables from formula, if given
randomVars = {};
if isfield(options, 'formula') && ~isempty(options.formula) % formula is given and not empty

    formula = options.formula;
    eqParts = strtrim(strsplit(formula, '~'));

    % lefthand side of formula
    if isempty(y)
        y = eqParts(1);
        nY = 1;
    end
    % righthand side of formula
    xStr = eqParts{2};

    % separate between fixed and random effects
    fixedStop = regexp(xStr, '\s*\+\s*\(.*?\)');    
    if ~isempty(fixedStop) % there is a random variable term
        xStrFixed = xStr(1:fixedStop(1));
        fixedVars = strtrim(strsplit(xStrFixed, {'+', '*', ':', '|'}));
        xStrRand = xStr(fixedStop(1):end);
        tokens = regexp(xStrRand, '\((.*?)\)', 'tokens');
        nTokens = length(tokens);
        for iToken = 1:nTokens
            token = tokens{iToken}{1};
            hits = strsplit(token, '|');
            randomVars = union(randomVars, strtrim(strsplit(hits{2}, {'+', '*', ':', '|'})), 'stable');
        end
    else % there is no random variable term
        tokens = regexp(xStr, '(.*)\s*', 'tokens');
        hits = tokens{:};
        fixedVars = strtrim(strsplit(hits{1}, {'+', '*', ':', '|'}));
        randomVars = {};
    end
    randomVars = reshape(randomVars, 1, []); % make horizontal array
   
    % x variables are all fixed-effect IVs except random variables and covariates
    if isempty(x)
        x = fixedVars;
        x = setdiff(x, randomVars, 'stable');
        x = setdiff(x, covariate, 'stable');
    end

    % set subject variable to 1st extracted random variable, if not given
    if isempty(subject)
        if ~isempty(randomVars)
            subject = randomVars{1};
        end
    end

    % set trial variable to 2nd extracted random variable, if not given
    if isempty(trial)
        if length(randomVars) > 1
            trial = randomVars{2};
        end
    end
else % no formula given or empty

    % random variable(s)
    if ~isempty(subject)
        randomVars = {subject};
    end
    if ~isempty(trial)
        randomVars = union(randomVars, {trial}, 'stable');
    end
    if ~isempty(stimulus)
        randomVars = union(randomVars, {stimulus}, 'stable');
    end

    % for safety, set formula to empty
    formula = '';
end

if nY > 1
    yVar = 'yVar';
    yVal = 'Y';
    Data2 = stack(Data1, y, 'NewDataVariableName', yVal, 'IndexVariableName', yVar);
    Data2.(yVar) = string(string(Data2.(yVar)));
    depVar = yVal; % set dependent variable to y
elseif ~isempty(multiVar)
    yVar = multiVar;
    yVal = y{1};
    if isempty(multiVarLevels)
        Data2 = Data1;
    else
        idxDesc = ismember(string(Data1.(yVar)), string(multiVarLevels));
        Data2 = Data1(idxDesc, :);
    end
    depVar = y{1}; % set dependent variable to y
    Data2.(yVar) = string(string(Data2.(yVar)));
    y = cellstr(unique(Data2.(yVar)));
    nY = length(y);
else
    Data2 = Data1;
    depVar = y{1};
    yVal = y{1};
end

% multiVariate
if isfield(options, 'multiVariate') && ~isempty(options.multiVariate)
    multiVariate = getValue(options.multiVariate);
elseif isfield(options, 'separateMulti') && ~isempty(options.separateMulti)
    multiVariate = not(getValue(options.separateMulti));
else
    multiVariate = false;
end

% fit method
if isfield(options, 'fitMethod') && ~isempty(options.fitMethod)
    fitMethod = options.fitMethod;
else
    fitMethod = 'MPL';
end

% posthoc method
if isfield(options, 'posthocMethod') && ~isempty(options.posthocMethod)
    posthocMethod = options.posthocMethod;
else
    posthocMethod = 'emm';
end

% posthoc correction
if isfield(options, 'posthocCorrection') && ~isempty(options.posthocCorrection)
    posthocCorrection = options.posthocCorrection;
else
    posthocCorrection = 'holm';
end

% posthoc main effects
if isfield(options, 'posthocMain') && ~isempty(options.posthocMain)
    posthocMain = getValue(options.posthocMain);
else
    posthocMain = false;
end

% posthoc comparison level
if isfield(options, 'posthocLevel') && ~isempty(options.posthocLevel)
    posthocLevel = getValue(options.posthocLevel);
else
    posthocLevel = 1;
end

% distribution (for GLM)
if isfield(options, 'distribution') && ~isempty(options.distribution)
    distribution = cellstr(options.distribution);
else % no parameter given
    distribution = {'normal'};
end

% link (for GLM)
if isfield(options, 'link') && ~isempty(options.link)
    % link may be numerical for power functions, e.g. -1 or -2
    % make them into cell array of strings
    if ischar(options.link) || isstring(options.link)
        link = cellstr(options.link);
    elseif isnumeric(options.link)
        link = {num2str(options.link)};
    elseif iscell(options.link)
        link = options.link;
        for iCell = 1:length(link)
            if isnumeric(link{iCell})
                link{iCell} = num2str(link{iCell});
            end
        end
    else
        error('No valid value given for ''link''');
    end
else % no parameter given
    link = {'auto'};
end

% create hash map for canonical link functions
distributions = {
    'normal'
    'gamma'
    'inverseGaussian'
    'binomial'
    'poisson'
    };
links = {
    'identity'
    'reciprocal'
    '-2'
    'logit'
    'log'
    };
canonicalLink = containers.Map(distributions,links);

% dummy coding
if isfield(options, 'dummyCoding') && ~isempty(options.dummyCoding)
    dummyCoding = options.dummyCoding;
else
    dummyCoding = 'effects';
end

% level order
if isfield(options, 'levelOrder') && ~isempty(options.levelOrder)
    levelOrder = options.levelOrder;
else
    levelOrder = 'sorted';
end
sortValues = '';
if ismember(levelOrder,{'ascend', 'descend'})
    sortValues = levelOrder;
    levelOrder = 'sorted';
end

% plot style
if isfield(options, 'plotStyle') && ~isempty(options.plotStyle)
    plotStyle = options.plotStyle;
else
    plotStyle = 'violin';
end

% plot style
if isfield(options, 'markerSize') && ~isempty(options.markerSize)
    markerSize = options.markerSize;
else
    markerSize = NaN;
end

% Flag to plot data
if isfield(options, 'isPlot') && ~isempty(options.isPlot)
    isPlot = getValue(options.isPlot);
else
    isPlot = true;
end

if isfield(options, 'barType') && ~isempty(options.barType)
    barType = lower(options.barType);
else
    barType = 'auto';
end

% flag to rescale all panel plots to the same y-scale
if isfield(options, 'rescale') && ~isempty(options.rescale)
    rescale = getValue(options.rescale);
else
    rescale = true;
end

% flag if pre-fit outliers should be removed
if isfield(options, 'outlierRemoval') && ~isempty(options.outlierRemoval)
    outlierRemoval = options.outlierRemoval;
elseif isfield(options, 'preOutlierMethod') && ~isempty(options.preOutlierMethod)
    outlierRemoval = options.preOutlierMethod;
elseif isfield(options, 'outlierMethod') && ~isempty(options.outlierMethod)
    outlierRemoval = options.outlierMethod;
else
    outlierRemoval = 'none';
end

% flag if post-fit outliers should be removed
if isfield(options, 'postOutlierMethod') && ~isempty(options.postOutlierMethod)
    postOutlierMethod = options.postOutlierMethod;
else
    postOutlierMethod = 'none';
end

% output folder
if isfield(options, 'outDir') && ~isempty(options.outDir)
    outDir = options.outDir;
    isSave = true;
else
    outDir = '';
    isSave = false;
end
if isSave && ~isfolder(outDir)
    mkdir(outDir);
end

if isfield(options, 'title') && ~isempty(options.title)
    plotTitle = options.title;
else
    plotTitle = '';
end

if isfield(options, 'showVarNames') && ~isempty(options.showVarNames)
    showVarNames = getValue(options.showVarNames);
else
    showVarNames = 'levels';
end

if isfield(options, 'xName') && ~isempty(options.xName)
    xName = getList(options.xName);
else
    xName = cellstr(x);
end

if isfield(options, 'yLabel') && ~isempty(options.yLabel)
    yLabel = cellstr(options.yLabel);
else
    yLabel = cellstr(y);
end

if isfield(options, 'yUnits')
    yUnits = cellstr(options.yUnits);
else
    yUnits = {''};
end

if isfield(options, 'yMult') && ~isempty(options.yMult)
    yMult = getValue(options.yMult);
else
    yMult = 1;
end

if isfield(options, 'xOrder') && ~isempty(options.xOrder)
    xOrder1 = getValue(options.xOrder);
elseif isfield(options, 'xOrder1') && ~isempty(options.xOrder1)
    xOrder1 = getValue(options.xOrder1);
else
    xOrder1 = [];
end

if isfield(options, 'xOrder2') && ~isempty(options.xOrder2)
    xOrder2 = getValue(options.xOrder2);
else
    xOrder2 = [];
end

if isfield(options, 'xOrder3') && ~isempty(options.xOrder3)
    xOrder3 = getValue(options.xOrder3);
else
    xOrder3 = [];
end

if isfield(options, 'xOrder4') && ~isempty(options.xOrder44)
    xOrder4 = getValue(options.xOrder4);
else
    xOrder4 = [];
end

if isfield(options, 'plotLines') && ~isempty(options.plotLines)
    plotLines = getValue(options.plotLines);
else
    plotLines = false;
end

% close figures
if isfield(options, 'closeFigures') && ~isempty(options.closeFigures)
    isCloseFigures = getValue(options.closeFigures);
elseif isSave
    isCloseFigures = true;
else
    isCloseFigures = false;
end

%% Apply constraint, if given

if isfield(options, 'constraint') && ~isempty(options.constraint)
    constraint = options.constraint;
    [conds, ops] = strsplit(constraint, {'&', '|'});
    allIdx = true(size(Data2, 1), 1);
    auxVars = {};
    for iCond = 1:length(conds)
        cond = strtrim(conds{iCond});
        [parts, splitmatches] = strsplit(cond, {'=', '==', '~=', '<', '<=', '>', '>='});
        constraintVar = strtrim(parts{1});
        constraintVal = strtrim(parts{2});
        if startsWith(constraintVal, '"')
            constraintVal = strrep(constraintVal, '"', ''); % remove double quotes
            constraintVal = string(constraintVal);
        elseif startsWith(constraintVal, '''')
            constraintVal = strrep(constraintVal, '''', ''); % remove single quotes
            constraintVal = string(constraintVal);
        else
            [num, isNum] = str2num(constraintVal);
            if isNum
                constraintVal = num;
            else
                constraintVal = string(constraintVal);
            end
        end
        % constraintVal = strrep(constraintVal, '''', ''); % remove single quotes
        % constraintVal = strrep(constraintVal, '"', ''); % remove double quotes
        % [num, isNum] = str2num(constraintVal);
        % if isNum
        %     constraintVal = num;
        % else
        %     constraintVal = string(constraintVal);
        % end
        compVar = strtrim(splitmatches{1});
        if ~ismember(constraintVar, Data2.Properties.VariableNames)
            auxVar = constraintVar;
            Data2.(constraintVar) = DataRaw.(constraintVar);
            auxVars = [auxVars; auxVar]; %#ok<AGROW>
        end
        constraintVals = Data2.(constraintVar);
        if isordinal(Data2.(constraintVar))
            constraintVals = str2double(string(constraintVals));
        elseif iscell(Data2.(constraintVar))
            constraintVals = string(constraintVals);
        end
        if ~strcmp(class(constraintVals), class(constraintVal))
            error('Comparison of ''%s'' is a comparison between ''%s'' and ''%s'', which is not possible', cond, class(constraintVals), class(constraintVal));
        end
        switch compVar
            case {'=', '=='}
                idxDesc = (constraintVals == constraintVal);
            case '~='
                idxDesc = (constraintVals ~= constraintVal);
            case '<'
                idxDesc = (constraintVals < constraintVal);
            case '<='
                idxDesc = (constraintVals <= constraintVal);
            case '>'
                idxDesc = (constraintVals > constraintVal);
            case '>='
                idxDesc = (constraintVals >= constraintVal);
            otherwise
                error('No valid comparison operator ''%s''', compVar);
        end
        if  iCond > 1 && iCond-1 <= length(ops)
            op = ops{iCond-1};
            cmd = sprintf('allIdx %s idxDesc', op);
            allIdx = eval(cmd);
        else
            allIdx = idxDesc;
        end
    end
    if any(allIdx)
        Data2 = Data2(allIdx, :);
        fprintf('Constraints removed %d of %d observations (%.2f%%)\n', sum(~allIdx), length(allIdx), sum(~allIdx)/length(allIdx)*100);
    else
        error('Constraint ''%s'' cannot be fulfilled', constraint);
    end
    switch length(unique(Data2.(constraintVar)))
        case 0
            fprintf('Variable "%s" has no levels left on constraint "%s" -> remove variable\n', constraintVar, constraintVal);
            x = setdiff(x, constraintVar, 'stable');
        case 1
            fprintf('Variable "%s" has only 1 level left on constraint "%s" -> remove variable\n', constraintVar, constraintVal);
            x = setdiff(x, constraintVar, 'stable');
    end
    for iVar = 1:length(auxVars)
        Data2 = removevars(Data2, auxVars{iVar});
    end
end

%% Create Data table

% init Data table
Data = table;

% subject variable
if ~isempty(subject)
    Data.(subject) = string(string(Data2.(subject))); % make categorical
end

% make fixed-effect IVs categorical by converting them to string
for iVar = 1:length(x)
    myVar = x{iVar};
    Data2.(myVar) = string(string(Data2.(myVar)));
end
nFactors = length(x);

% collect all independent variables
IVs = x; % start with fixed-effect IVs
IVs = union(IVs, covariate, 'stable'); % include covariates
IVs = union(IVs, randomVars, 'stable'); % include random variables
for iIV = 1:length(IVs)
    myIV = IVs{iIV};
    Data.(myIV) = Data2.(myIV);
end

% multivariate variable
if nY > 1
    Data.(yVar) = Data2.(yVar);
end

% get dependent variable
Data.(depVar) = Data2.(depVar);

% store Data table in options
options.Data = Data;

%% Create variables

% 1st factor = member variable
memberVar = x{1};
members = unique(Data.(memberVar), levelOrder);
nMembers = length(members);

% 2nd factor = group variable
if nFactors > 1
    groupVar = x{2};
    groups = unique(Data.(groupVar), levelOrder);
    nGroups = length(groups);
else
    groupVar = 'none';
    groups = string(x{1});
    nGroups = 1;
end

% 3rd factor = column variable
if nFactors > 2
    colVar = x{3};
    cols = unique(Data.(colVar), levelOrder);
    nCols = length(cols);
else
    colVar = 'none';
    cols = "";
    nCols = 1;
end

% 4th factors = row variable
if nFactors > 3
    rowVar = x{4};
    rows = unique(Data.(rowVar), levelOrder);
    nRows = length(rows);
else
    rowVar = 'none';
    rows = "";
    nRows = 1;
end

% apply potential re-ordering
if ~isempty(xOrder1)
    members = members(xOrder1);
end
if ~isempty(xOrder2)
    groups = groups(xOrder2);
end
if ~isempty(xOrder3)
    cols = cols(xOrder3);
end
if ~isempty(xOrder4)
    rows = rows(xOrder4);
end

%% Apply data transformation, if given

if ~isempty(transform)
    transVar = sprintf('%sTrans', depVar);
    for iVar = 1:nY
        if nY > 1
            idxDesc = (Data.(yVar) == y{iVar});
        else
            idxDesc = true(size(Data, 1), 1);
        end
        switch transform
            case 'mean'
                transformFcn = @(x) x/mean(x, 'omitnan');
            case 'std'
                transformFcn = @(x) x/std(x, 'omitnan');
            case 'Z'
                transformFcn = @(x) (x - mean(x, 'omitnan'))/std(x, 'omitnan');
            case {'median'}
                upperThresh = quantile(Data.(depVar)(idxDesc), 0.5);
                transformFcn = @(x) x/upperThresh;
            case 'IQR'
                [~, lowerThresh, upperThresh, ~] = isoutlier(Data.(depVar)(idxDesc), 'quartiles');
                transformFcn = @(x) x/(upperThresh - lowerThresh);
            case 'IQRmax'
                [~, ~, upperThresh, ~] = isoutlier(Data.(depVar)(idxDesc), 'quartiles');
                transformFcn = @(x) x/upperThresh;
            case 'MAD'
                [~, lowerThresh, upperThresh, ~] = isoutlier(Data.(depVar)(idxDesc), 'median');
                transformFcn = @(x) (x - lowerThresh)/(upperThresh - lowerThresh);
            case 'MADmax'
                [~, ~, upperThresh, ~] = isoutlier(Data.(depVar)(idxDesc), 'median');
                transformFcn = @(x) x/upperThresh;
            case 'max'
                upperThresh = max(Data.(depVar)(idxDesc));
                transformFcn = @(x) x/upperThresh;
            case 'minmax'
                lowerThresh = min(Data.(depVar)(idxDesc));
                upperThresh = max(Data.(depVar)(idxDesc));
                transformFcn = @(x) x/(upperThresh - lowerThresh);
            otherwise % any other expression
                tokens = regexp(transform, '[q,p](\d+)(?:[q,p]?)(\d+)?', 'tokens', 'once');
                tokens(cellfun(@isempty, tokens)) = [];
                if isscalar(tokens) % upper percentile given in the form 'q%d' or 'p%d'
                    upperThresh = prctile(Data.(depVar)(idxDesc), str2double(tokens{1}));
                    transformFcn = @(x) x/upperThresh;
                elseif length(tokens) == 2 % lower and upper percentile given in the form 'q%dq%d' or 'p%dp%d'
                    lowerThresh = prctile(Data.(depVar)(idxDesc), str2double(tokens{1}));
                    upperThresh = prctile(Data.(depVar)(idxDesc), str2double(tokens{2}));
                    transformFcn = @(x) x/(upperThresh - lowerThresh);
                else % any other function given in the form 'f(x)'
                    transformFcn = eval(sprintf('@(x) %s', transform));
                end
        end
        Data.(transVar)(idxDesc) = transformFcn(Data.(depVar)(idxDesc));
    end
else
    transVar = depVar;
end

%% Remove pre-fit outliers

nPreOutliers = 0;
nPreObs = size(Data, 1);

outlierLevel = nFactors;
if ~strcmp(outlierRemoval, 'none')

    idxOut = false(size(Data, 1), 1);

    for iVar = 1:nY

        if nY > 1
            idxDesc = (Data.(yVar) == y{iVar});
        else
            idxDesc = true(size(Data, 1), 1);
        end

        idxTest = idxDesc;

        if outlierLevel == 0 % level 0: all data
            yData = Data.(transVar)(idxTest);
            idxOut(idxTest) = isoutlier(yData, outlierRemoval);

        elseif outlierLevel == 1 % level 1: 1st dependent variable
            for iMember = 1:nMembers
                idxTest = idxDesc;
                member = members(iMember);
                idxTest = idxTest & (Data.(memberVar) == member);
                yData = Data.(transVar)(idxTest);
                idxOut(idxTest) = isoutlier(yData, outlierRemoval);
            end

        elseif outlierLevel == 2 % level 2: 2nd dependent variable, if given
            for iMember = 1:nMembers
                for iGroup = 1:nGroups
                    idxTest = idxDesc;
                    member = members(iMember);
                    idxTest = idxTest & (Data.(memberVar) == member);
                    if nGroups > 1
                        group = groups(iGroup);
                        idxTest = idxTest & (Data.(groupVar) == group);
                    end
                    yData = Data.(transVar)(idxTest);
                    idxOut(idxTest) = isoutlier(yData, outlierRemoval);
                end
            end

        elseif outlierLevel == 3 % level 3: 3rd dependent variable, if given
            for iMember = 1:nMembers
                for iGroup = 1:nGroups
                    for iCol = 1:nCols
                        idxTest = idxDesc;
                        member = members(iMember);
                        idxTest = idxTest & (Data.(memberVar) == member);
                        if nGroups > 1
                            group = groups(iGroup);
                            idxTest = idxTest & (Data.(groupVar) == group);
                        end
                        if nCols > 1
                            col = cols(iCol);
                            idxTest = idxTest & (Data.(colVar) == col);
                        end
                        yData = Data.(transVar)(idxTest);
                        idxOut(idxTest) = isoutlier(yData, outlierRemoval);
                    end
                end
            end
        elseif outlierLevel == 4 % level 4: 4th dependent variable, if given
            for iMember = 1:nMembers
                for iGroup = 1:nGroups
                    for iCol = 1:nCols
                        for iRow = 1:nRows
                            idxTest = idxDesc;
                            member = members(iMember);
                            idxTest = idxTest & (Data.(memberVar) == member);
                            if nGroups > 1
                                group = groups(iGroup);
                                idxTest = idxTest & (Data.(groupVar) == group);
                            end
                            if nCols > 1
                                col = cols(iCol);
                                idxTest = idxTest & (Data.(colVar) == col);
                            end
                            if nRows > 1
                                row = rows(iRow);
                                idxTest = idxTest & (Data.(rowVar) == row);
                            end
                            yData = Data.(transVar)(idxTest);
                            idxOut(idxTest) = isoutlier(yData, outlierRemoval);
                        end
                    end
                end
            end
        end
    end

    % remove pre-fit outliers
    nPreOutliers = sum(idxOut);
    if nPreOutliers > 0
        Data = Data(~idxOut, :);
    end
end

%% Multiply dependent variable data with provided factor(s)

for iVar = 1:nY
    % dependent variable
    myVar = y{iVar};

    % multiplication factor(s)
    if length(yMult) == nY
        myMult = yMult(iVar);
    else
        myMult = yMult(1);
    end

    % do multiplication
    if nY > 1 % multiple dependent variables or multivariate dependent variable
        idxDesc = (Data.(yVar) == myVar);
        Data.(transVar)(idxDesc) = myMult * Data.(transVar)(idxDesc);
        if myMult ~= 1
            fprintf('Multiplying %s with %g...\n', myVar, myMult);
        end
    else % one dependent variable
        Data.(transVar) = myMult * Data.(transVar);
        if myMult ~= 1
            fprintf('Multiplying %s with %g...\n', depVar, myMult);
        end
    end


end

%% Fit linear model and perform ANOVA

ticAnova = tic;

% perform nY fits only if multiVariate=false
nFits = 1;
if strcmp(fitMethod, 'none')
    nFits = 0;
elseif ~multiVariate
    nFits = nY;
end
mdls = cell(1, nFits);
anovas = cell(1, nFits);
Datasets = cell(1, nFits);
for iFit = 1:nFits % if multiVariate, this loop is left after the 1st iteration
    % dependent variable
    myVar = y{iFit};

    % yLabel
    if length(yLabel) == nY
        myLabel = yLabel{iFit};
    elseif isscalar(yLabel)
        myLabel = yLabel{1};
    end

    % create output folder
    outSubDir = sprintf('%s/%s', outDir, myVar);
    if isSave
        if ~isfolder(outSubDir)
            mkdir(outSubDir);
        end
    end

    % prepare options for output
    fields = fieldnames(options);
    fields = setdiff(fields, 'Data', 'stable'); % remove "Data" from options
    fields = setdiff(fields, 'DataRaw', 'stable'); % remove "DataRaw" from options
    paramFields = sort(fields); % sort fields alphabetically

    % add options to results
    for iField = 1:length(paramFields)
        field = paramFields{iField};
        results.options.(field) = mat2str(string(options.(field)));
    end

    % save options to file
    if isSave
        fpath = fullfile(outSubDir, 'Options.txt');
        fidOptions = fopen(fpath, 'w+');
        for iField = 1:length(paramFields)
            field = paramFields{iField};
            fprintf(fidOptions, 'options.%s = %s;\n', field, mat2str(string(options.(field))));
        end
        fclose(fidOptions);
    end

    if length(distribution) >= nFits
        myDistribution = distribution{iFit};
    else
        myDistribution = distribution{1};
    end

    if length(link) >= nFits
        myLink = link{iFit};
    else
        myLink = link{1};
    end

    if nY > 1 && ~multiVariate % separate univariate analyses of multi-valued dependent variable
        idxDesc = (Data.(yVar) == myVar);
        DataOrig = Data;
        Data = Data(idxDesc, :);
        fprintf('Performing GLMM analysis for %s...\n', myVar);
    elseif nY > 1
        fprintf('Performing multivariate GLMM analysis...\n');
    else
        fprintf('Performing GLMM analysis for %s...\n', depVar);
    end

    % open summary file for writing
    if isSave
        fidSummary = fopen(fullfile(outSubDir, 'Summary.txt'), 'w+');
    else
        fidSummary = 2; % screen
    end

    if isempty(formula)

        productTerm = strjoin(interact, '*');
        if nY > 1 && multiVariate
            productTerm = sprintf('-1 + %s:(%s)', yVar, productTerm);
        end
        xNoInteract = setdiff(x, interact, 'stable');

        sumTerm = strjoin(union(xNoInteract, covariate, 'stable'), ' + ');

        if ~isempty(subject) && length(unique(Data.(subject))) > 1
            if isRandomSlopes && ~isempty(randomSlopes)

                myRandomIntercept = '';
                if nY > 1 && multiVariate
                    myRandomEffects = strjoin(cellfun(@(x) sprintf('(%s|%s:%s)', x, yVar, subject), randomSlopes, 'UniformOutput', false), ' + ');
                else
                    % determine if random slopes or not and choose operator
                    if isRandomInteract % estimate random interactions in addition to random slopes
                        randomSlopesTerm = strjoin(randomSlopes, '*');
                    else % estimate random slopes (without interactions)
                        randomSlopesTerm = strjoin(randomSlopes, '+');
                    end

                    % determine if trial variable is given, and compose random effect terms
                    myRandomEffectTerms = {sprintf('(%s|%s)', randomSlopesTerm, subject)};
                    if ~isempty(trial)
                        myRandomEffectTerms = [myRandomEffectTerms, {sprintf('(%s|%s:%s)', randomSlopesTerm, subject, trial)}]; %#ok<AGROW>
                    end

                    % determine if stimulus variable is given, and compose random effect terms
                    if ~isempty(stimulus)
                        myRandomEffectTerms = [myRandomEffectTerms, {sprintf('(%s|%s)', randomSlopesTerm, stimulus)}]; %#ok<AGROW>
                    end

                    % compose random effects string
                    myRandomEffects = strjoin(myRandomEffectTerms, ' + ');
                end
            elseif ~isRandomSlopes && ~isempty(union(interact, covariate)) % no interaction-effect variables and no covariates
                myRandomIntercept = strjoin(cellfun(@(x) sprintf('(1|%s:%s)', x, subject), randomSlopes, 'UniformOutput', false), ' + ');
                myRandomEffects = '';
            else
                if ~isempty(trial)
                    myRandomIntercept = sprintf('(1|%s:%s)', subject, trial);
                else
                    myRandomIntercept = sprintf('(1|%s)', subject);
                end                
                myRandomEffects = '';
            end
        else
            myRandomIntercept = '';
            myRandomEffects = '';
        end

        myTerms = {};
        if ~isempty(sumTerm)
            myTerms = [myTerms, sumTerm]; %#ok<AGROW>
        end
        if ~isempty(productTerm)
            myTerms = [myTerms, productTerm]; %#ok<AGROW>
        end
        if ~isempty(myRandomIntercept)
            myTerms = [myTerms, myRandomIntercept]; %#ok<AGROW>
        end
        if ~isempty(myRandomEffects)
            myTerms = [myTerms, myRandomEffects]; %#ok<AGROW>
        end

        formulaOrig = formula;
        formula = sprintf('%s ~ %s', transVar, strjoin(myTerms, ' + '));
    else
        formulaOrig = formula;
        parts = strtrim(strsplit(formula, '~'));
        parts{1} = transVar;
        formula = strjoin(parts, ' ~ ');
    end
    fprintf('\t%s\n', formula);

    fprintf('Fitting the model...');
    ticFit = tic;

    try
        if ~isempty(myLink) && ~strcmp(myLink, 'auto') % link function given -> use it
            if ismember(myDistribution, distributions) && ~strcmp(myLink, canonicalLink(myDistribution))
                [numVal, isNum] = str2num(myLink);
                if isNum
                    myLink = numVal;
                end
                mdl = fitglme(Data, formula, ...
                    'DummyVarCoding', dummyCoding, ...
                    'FitMethod', fitMethod, ...
                    'Distribution', myDistribution, ...
                    'link', myLink, ...
                    'EBMethod', 'TrustRegion2D'); % use this EBMethod, if custom link function is not canonical link function
            else
                [numVal, isNum] = str2num(myLink);
                if isNum
                    myLink = numVal;
                end
                mdl = fitglme(Data, formula, ...
                    'DummyVarCoding', dummyCoding, ...
                    'FitMethod', fitMethod, ...
                    'Distribution', myDistribution, ...
                    'link', myLink);
            end
        else % no link function given or set to 'auto' -> use built-in default
            mdl = fitglme(Data, formula, ...
                'DummyVarCoding', dummyCoding, ...
                'FitMethod', fitMethod, ...
                'Distribution', myDistribution);
        end

    catch ME
        message = sprintf('%s', ME.message);
        fprintf('The linear model fit returned an error:\n\t%s\n', message);
        fprintf('Please try again, using fewer interactions by defining "interact" with only those independent variables whose interaction you want to investigate\n');
        return
    end
    fprintf('done in %g seconds\n', toc(ticFit));

    % remove post-fit outliers and refit model
    nPostOutliers = 0;
    mdlResiduals = residuals(mdl, 'ResidualType', 'Pearson');
    nPostObs = length(mdlResiduals);
    if ~strcmp(postOutlierMethod, 'none')
        postOutliers = isoutlier(mdlResiduals, postOutlierMethod);
        nPostOutliers = sum(postOutliers);
        msg = sprintf('Removing %d post-fit outliers from %d observations (%.1f %%%%) using ''%s''...\n', nPostOutliers, nPostObs, nPostOutliers/nPostObs*100, postOutlierMethod);
        fprintf(msg);
        fprintf(fidSummary, msg);
        if nPostOutliers > 0
            % remove outliers from Data
            Data(postOutliers, :) = [];
            % remove outliers from DataOrig
            if nY > 1 && ~multiVariate
                idxTmp = false(size(DataOrig, 1), 1);
                idxTmp(idxDesc) =  postOutliers;
                DataOrig(idxTmp, :) = [];
            end
            msg = sprintf('Re-fitting model...');
            fprintf(msg);
            fprintf(fidSummary, msg);
            try
                if ~isempty(myLink) && ~strcmp(myLink, 'auto') % link function given -> use it
                    mdl = fitglme(Data, formula, ...
                        'DummyVarCoding', dummyCoding, ...
                        'FitMethod', fitMethod, ...
                        'Distribution', myDistribution, ...
                        'link', myLink);
                else % no link function given or set to 'auto' -> use built-in default
                    mdl = fitglme(Data, formula, ...
                        'DummyVarCoding', dummyCoding, ...
                        'FitMethod', fitMethod, ...
                        'Distribution', myDistribution);
                end
            catch ME
                message = sprintf('%s', ME.message);
                fprintf('The linear model fit returned an error:\n\t%s\n', message);
                fprintf('Please try again, probably using fewer interactions by defining "interact" with only those independent variables whose interaction you want to investigate\n');
                return
            end
        end
    end

    % print results of model fit into file
    mdlOutput = formattedDisplayText(mdl, 'SuppressMarkup', true);
    fprintf(fidSummary, 'Formula:\n\t%s\n', formula);
    fprintf(fidSummary, '\t%s', mdlOutput);
    if ~strcmp(posthocMethod, 'none')
        fprintf(fidSummary, 'Performing post-hoc comparison using method ''%s''\n', posthocMethod);
    end
    if ~isempty(transform)
        fprintf(fidSummary, 'Data have been transformed using f(x) = %s\n', transform);
    end

    % close file
    if isSave
        fclose(fidSummary);
    end

    % report outlier removal to file
    if isSave
        fidOutliers = fopen(fullfile(outSubDir, 'Outliers.txt'), 'w+');
    else
        fidOutliers = 2; % screen
    end
    msg = sprintf('Removed %d pre-fit outlier(s) from %d observations (%.1f %%%%) using removal method ''%s''\n', nPreOutliers, nPreObs, nPreOutliers/nPreObs*100, outlierRemoval);
    fprintf(msg);
    fprintf(fidOutliers, msg);
    msg = sprintf('Removed %d post-fit outlier(s) from %d observations (%.1f %%%%) using removal method ''%s''\n', nPostOutliers, nPostObs, nPostOutliers/nPostObs*100, postOutlierMethod);
    fprintf(msg);
    fprintf(fidOutliers, msg);

    % close file
    if isSave
        fclose(fidOutliers);
    end

    %% Plot diagnostics

    % Since raw residuals for generalized linear mixed-effects models do not
    % have a constant variance across observations, we use the conditional
    % Pearson residuals instead.

    panelWidth = 300;
    panelHeight = 300;
    nPanels = 6;
    nPanelRows = 2;
    nPanelCols = ceil(nPanels/nPanelRows);
    figWidth = nPanelCols * panelWidth;
    figHeight = nPanelRows * panelHeight;
    figName = 'Diagnostics';
    fig = figure('Name', figName, 'Position', [0, 0, figWidth, figHeight]);

    if nY > 1 && ~multiVariate
        if strcmp(depVar, yVal)
            if isempty(plotTitle)
                sgtitle(sprintf('Diagnostics for %s', myLabel), 'interpreter', 'none', 'FontWeight', 'bold', 'FontSize', 14);
            else
                sgtitle(sprintf('Diagnostics for %s (%s)', plotTitle, myLabel), 'interpreter', 'none', 'FontWeight', 'bold', 'FontSize', 14);
            end
        else
            if isempty(plotTitle)
                sgtitle(sprintf('Diagnostics for %s (%s)', myVar, depVar), 'interpreter', 'none', 'FontWeight', 'bold', 'FontSize', 14);
            else
                sgtitle(sprintf('Diagnostics for %s (%s, %s)', plotTitle, myVar, depVar), 'interpreter', 'none', 'FontWeight', 'bold', 'FontSize', 14);
            end
        end
    elseif nY > 1
        if isempty(plotTitle)
            sgtitle(sprintf('Diagnostics for %s - multivariate Analysis', depVar), 'interpreter', 'none', 'FontWeight', 'bold', 'FontSize', 14);
        else
            sgtitle(sprintf('Diagnostics for %s (%s) - multivariate Analysis', plotTitle, depVar), 'interpreter', 'none', 'FontWeight', 'bold', 'FontSize', 14);
        end
    else
        if isempty(plotTitle)
            sgtitle(sprintf('Diagnostics for %s', myLabel), 'interpreter', 'none', 'FontWeight', 'bold', 'FontSize', 14);
        else
            sgtitle(sprintf('Diagnostics for %s (%s)', plotTitle, myLabel), 'interpreter', 'none', 'FontWeight', 'bold', 'FontSize', 14);
        end
    end

    iPanel = 0;

    % histogram
    iPanel = iPanel+1;
    subplot(nPanelRows, nPanelCols, iPanel);
    mdlResiduals = residuals(mdl, 'ResidualType', 'Pearson');
    h = histfit(mdlResiduals);
    h(1).EdgeColor = 'none';
    h(2).Color = [1 0 0];
    h(2).LineStyle = ':';
    hold on
    resFit = fitdist(mdlResiduals, 'Normal');
    xline(resFit.mu, 'Color', 'r', 'LineWidth', 2);
    ylims = ylim;
    rectangle('Position', [resFit.mu - resFit.sigma, ylims(1), 2*resFit.sigma, ylims(2)], 'FaceColor', [0 0 0], 'FaceAlpha', 0.1, 'EdgeColor', 'none');
    title('Histogram of residuals');
    xlabel('Residuals');

    % normality of residuals
    iPanel = iPanel+1;
    subplot(nPanelRows, nPanelCols, iPanel);
    plotResiduals(mdl, 'probability', 'ResidualType', 'Pearson');

    % symmetry
    iPanel = iPanel+1;
    subplot(nPanelRows, nPanelCols, iPanel);
    plotResiduals(mdl, 'symmetry', 'ResidualType', 'Pearson');

    % fitted-response
    iPanel = iPanel+1;
    subplot(nPanelRows, nPanelCols, iPanel);
    F = fitted(mdl);
    R = response(mdl);
    plot(R, F, 'rx');
    ylabel('Fitted');
    title('Fitted vs. Response');
    xlabel('Response');

    % residuals vs fitted
    iPanel = iPanel+1;
    subplot(nPanelRows, nPanelCols, iPanel);
    plotResiduals(mdl, 'fitted', 'ResidualType', 'Pearson');

    % lagged residuals
    iPanel = iPanel+1;
    subplot(nPanelRows, nPanelCols,iPanel);
    plotResiduals(mdl, 'lagged', 'ResidualType', 'Pearson');

    % save figure
    if isSave
        kbsaveFigure(fig, fullfile(outSubDir, figName), {'.pdf', '.png', '.fig'});
    end

    % close figure
    if isCloseFigures
        close(fig);
    end

    % store GLM fit
    mdls{iFit} = mdl;

    % get ANOVA table from model fit and sore
    anovaResult = anova(mdl);
    anovas{iFit} = anovaResult;

    Datasets{iFit} = Data;

    % if separate multi, restore Data and formula
    if nY > 1 && ~multiVariate
        Data = DataOrig; % set Data to original Data minus post-fit outliers
        formula = formulaOrig;
    end

    % add model fit and stuff to results
    results.(myVar).formula = formula;
    results.(myVar).mdl = mdl;
    results.(myVar).anovaResult = anovaResult;
end

%% Statistically correct ANOVA tables

if ~isempty(anovas)
    % Sidak correction
    for iFit = 1:nFits
        % anovas{iFit}.pValue = sidak_corr(anovas{iFit}.pValue, nFits); % correct for multivariate dimension
        anovas{iFit}.pValue = sidak_corr(anovas{iFit}.pValue, correctForN); % correct for custom number of tests
    end
end

%% Create ANOVA table

for iFit = 1:nFits
    % dependent variable
    myVar = y{iFit};

    % create output folder
    outSubDir = sprintf('%s/%s', outDir, myVar);
    if isSave && ~isfolder(outSubDir)
        mkdir(outSubDir);
    end

    % retrieve anova results
    anovaResult = anovas{iFit};

    % init ANOVA table
    anovaTable = table;

    anovaTable.Term = string(anovaResult.Term);
    anovaTable.DF1 = anovaResult.DF1;
    anovaTable.DF2 = anovaResult.DF2;
    anovaTable.F = anovaResult.FStat;
    anovaTable.p = anovaResult.pValue;
    anovaTable.etaSqp = f2etaSqp(anovaTable.F, anovaTable.DF1, anovaTable.DF2);
    anovaTable.effectSize = string(effprint(anovaTable.etaSqp, 'eta2'));
    anovaTable.significance = string(sigprint(anovaTable.p));

    % add data to results
    results.(myVar).Data = Datasets{iFit};
    results.(myVar).DataRaw = DataRaw;
    results.(myVar).anovaTable = anovaTable;

    % Save tables
    if isSave
        % save Data
        saveTable(Datasets{iFit}, 'Data', {'csv'}, outSubDir);
        % save raw Data table
        saveTable(DataRaw, 'DataRaw', {'csv'}, outSubDir);
        % save ANOVA table
        saveTable(anovaTable, 'Anova', {'xlsx'}, outSubDir);
    end

    % display table
    disp(anovaTable)
end

%% Plot data and make posthoc comparisons

allVars = {memberVar, groupVar, colVar, rowVar};
allVarLevels = {members, groups, cols, rows};
nPosthocLevels = min(posthocLevel, nFactors);

for iLevel = 1:nPosthocLevels

    % define what is member and what are the others
    memberVar = allVars{iLevel};
    members = allVarLevels{iLevel};
    newIdx = 1:4;
    newIdx([1, iLevel]) = newIdx([iLevel, 1]);
    groupVar = allVars{newIdx(2)};
    groups = allVarLevels{newIdx(2)};
    colVar = allVars{newIdx(3)};
    cols = allVarLevels{newIdx(3)};
    rowVar = allVars{newIdx(4)};
    rows = allVarLevels{newIdx(4)};
    nMembers = length(members);
    nGroups = length(groups);
    nCols = length(cols);
    nRows = length(rows);

    % display names
    idxDisp = strcmp(x, memberVar);
    if any(idxDisp)
        memberVarDisp = xName{idxDisp};
    else
        memberVarDisp = memberVar;
    end
    idxDisp = strcmp(x, groupVar);
    if any(idxDisp)
        groupVarDisp = xName{idxDisp};
    else
        groupVarDisp = groupVar;
    end
    idxDisp = strcmp(x, colVar);
    if any(idxDisp)
        colVarDisp = xName{idxDisp};
    else
        colVarDisp = colVar;
    end
    idxDisp = strcmp(x, rowVar);
    if any(idxDisp)
        rowVarDisp = xName{idxDisp};
    else
        rowVarDisp = rowVar;
    end

    % define posthoc comparison pairs
    pairs = nchoosek(members, 2);
    nPairs = size(pairs, 1);

    % create arrays for the estimated marginal means and confidence intervals
    emmCenter = nan(nGroups, nMembers, nRows, nCols, nY);
    emmBottom = nan(nGroups, nMembers, nRows, nCols, nY);
    emmTop = nan(nGroups, nMembers, nRows, nCols, nY);

    % create arrays for the bar data
    barCenter = nan(nGroups, nMembers, nRows, nCols, nY);
    barBottom = nan(nGroups, nMembers, nRows, nCols, nY);
    barTop = nan(nGroups, nMembers, nRows, nCols, nY);

    % create (nRows x nCols x nGroups x nMembers x nY) arrays of plot data
    bar_p = nan(nGroups, nPairs, nRows, nCols, nY);
    bar_test = nan(nGroups, nPairs, nRows, nCols, nY);
    bar_DF = nan(nGroups, nPairs, nRows, nCols, nY);
    bar_aux = nan(nGroups, nPairs, nRows, nCols, nY);
    bar_smd = nan(nGroups, nPairs, nRows, nCols, nY);
    bar_eff = nan(nGroups, nPairs, nRows, nCols, nY);
    bar_diff = nan(nGroups, nPairs, nRows, nCols, nY);
    bar_diffpct = nan(nGroups, nPairs, nRows, nCols, nY);
    bar_emm_1 = nan(nGroups, nPairs, nRows, nCols, nY);
    bar_CI95low_1 = nan(nGroups, nPairs, nRows, nCols, nY);
    bar_CI95up_1 = nan(nGroups, nPairs, nRows, nCols, nY);
    bar_emm_2 = nan(nGroups, nPairs, nRows, nCols, nY);
    bar_CI95low_2 = nan(nGroups, nPairs, nRows, nCols, nY);
    bar_CI95up_2 = nan(nGroups, nPairs, nRows, nCols, nY);

    % create (nPairs x nY) arrays of main effects data
    main_p = nan(nPairs, nRows, nCols, nY);
    main_test = nan(nPairs, nRows, nCols, nY);
    main_DF = nan(nPairs, nRows, nCols, nY);
    main_aux = nan(nPairs, nRows, nCols, nY);
    main_smd = nan(nPairs, nRows, nCols, nY);
    main_eff = nan(nPairs, nRows, nCols, nY);
    main_diff = nan(nPairs, nRows, nCols, nY);
    main_diffpct = nan(nPairs, nRows, nCols, nY);
    main_emm_1 = nan(nPairs, nRows, nCols, nY);
    main_CI95low_1 = nan(nPairs, nRows, nCols, nY);
    main_CI95up_1 = nan(nPairs, nRows, nCols, nY);
    main_emm_2 = nan(nPairs, nRows, nCols, nY);
    main_CI95low_2 = nan(nPairs, nRows, nCols, nY);
    main_CI95up_2 = nan(nPairs, nRows, nCols, nY);

    maxNValues = 0;
    for iVar = 1:nY
        myVar = y{iVar};
        if nY > 1
            if nGroups > 1
                myMaxNValues = max(cell2mat(arrayfun(@(s1,s2) sum(Data.(memberVar)==s1 & Data.(groupVar)==s2 & Data.(yVar)==myVar), repmat(members(:)',nGroups,1), repmat(groups(:),1,nMembers), 'UniformOutput', false)),[],'all');
            else
                myMaxNValues = max(cell2mat(arrayfun(@(s1,s2) sum(Data.(memberVar)==s1 & Data.(yVar)==myVar), repmat(members(:)',nGroups,1), repmat(groups(:),1,nMembers), 'UniformOutput', false)),[],'all');
            end
        else
            if nGroups > 1
                myMaxNValues = max(cell2mat(arrayfun(@(s1,s2) sum(Data.(memberVar)==s1 & Data.(groupVar)==s2), repmat(members(:)',nGroups,1), repmat(groups(:),1,nMembers), 'UniformOutput', false)),[],'all');
            else
                myMaxNValues = max(cell2mat(arrayfun(@(s1,s2) sum(Data.(memberVar)==s1), repmat(members(:)',nGroups,1), repmat(groups(:),1,nMembers), 'UniformOutput', false)),[],'all');
            end
        end
        maxNValues = max([maxNValues, myMaxNValues]);
    end
    DataPoints = nan(nGroups, nMembers, maxNValues, nRows, nCols, nY);

    for iVar = 1:nY
        % dependent variable
        myVar = y{iVar};
        mdl = mdls{iVar};

        % Calc estimated marginal means of all analyzed factors. Continuous
        % variables cannot be included, because then emmeans gives an
        % error. We only include the analyzed factors as given by 'x',
        % which may deviate from the list of all fixed effect variables
        % occurring in the formula, when explicitly specified by the user.
        emm = emmeans(mdl, reshape(x, 1, []), dummyCoding, 'unbalanced', {mdl.ResponseName, mdl.Link.Link, mdl.Link.Inverse});

        % create output folder
        outSubDir = sprintf('%s/%s', outDir, myVar);
        if isSave && ~isfolder(outSubDir)
            mkdir(outSubDir);
        end

        % posthoc method
        if strcmp(posthocMethod, 'auto')
            if nY > 1 && ~multiVariate && ~strcmp(fitMethod, 'none')
                posthocMethod = 'emm';
            elseif nY == 1 && ~strcmp(fitMethod, 'none')
                posthocMethod = 'emm';
            elseif strcmp(distribution, 'normal')
                posthocMethod = 'ttest';
            else
                posthocMethod = 'utest';
            end
        end

        % init statistics table
        StatsTable = table;

        % calc plot data and emm data
        for iRow = 1:nRows % loop over 4th x
            for iCol = 1:nCols % loop over 3rd x
                for iGroup = 1:nGroups % loop over 2nd x
                    for iMember = 1:nMembers % loop over 1st x

                        % init descriptive statistics table row
                        statsRow = table;

                        if nY > 1
                            idxDesc = (Data.(yVar) == myVar);
                            if multiVariate
                                idxEmm = (emm.table.(yVar) == myVar);
                            else
                                idxEmm = true(size(emm.table, 1), 1);
                            end
                        else
                            idxDesc = true(size(Data, 1), 1);
                            idxEmm = true(size(emm.table, 1), 1);
                        end

                        % 4th x, if given
                        if nRows > 1
                            row = rows(iRow);
                            idxDesc = idxDesc & Data.(rowVar) == row;
                            idxEmm = idxEmm & (emm.table.(rowVar) == row);
                            statsRow.(rowVarDisp) = string(row);
                        end

                        % 3rd x, if given
                        if nCols > 1
                            col = cols(iCol);
                            idxDesc = idxDesc & Data.(colVar) == col;
                            idxEmm = idxEmm & (emm.table.(colVar) == col);
                            statsRow.(colVarDisp) = string(col);
                        end

                        % 2nd x, if given
                        if nGroups > 1
                            group = groups(iGroup);
                            idxDesc = idxDesc & (Data.(groupVar) == group);
                            idxEmm = idxEmm & (emm.table.(groupVar) == group);
                            statsRow.(groupVarDisp) = string(group);
                        end

                        % 1st x
                        member = members(iMember);
                        idxDescMember = idxDesc & (Data.(memberVar) == member);
                        idxEmmMember = idxEmm & (emm.table.(memberVar) == member);
                        statsRow.(memberVarDisp) = string(member);

                        emMean = mean(emm.table.Estimated_Marginal_Mean(idxEmmMember));
                        statsRow.emMean = mdl.Link.Inverse(emMean);
                        emCI95 = mean(emm.table.CI_95_0pct(idxEmmMember, :), 1);
                        statsRow.emCI95_lower = min(mdl.Link.Inverse(emCI95));
                        statsRow.emCI95_upper = max(mdl.Link.Inverse(emCI95));

                        vals = Data.(transVar)(idxDescMember);
                        statsRow.N = length(vals(~isnan(vals)));
                        statsRow.mean = mean(vals, 'omitnan');
                        statsRow.std = std(vals, 'omitnan');
                        statsRow.SE = statsRow.std / sqrt(statsRow.N);
                        statsRow.median = median(vals, 'omitnan');
                        statsRow.q25 = quantile(vals, 0.25);
                        statsRow.q75 = quantile(vals, 0.75);
                        tci95 = tinv([0.025 0.975], statsRow.N-1); % 95% of t-Distribution
                        ci95 = statsRow.mean + tci95 * statsRow.SE;
                        statsRow.CI95_lower = ci95(1);
                        statsRow.CI95_upper = ci95(2);
                        StatsTable = [StatsTable; statsRow]; %#ok<AGROW>

                        DataPoints(iGroup, iMember, 1:length(vals), iRow, iCol, iVar) = vals;
                        emmCenter(iGroup, iMember, iRow, iCol, iVar) = statsRow.emMean;
                        emmBottom(iGroup, iMember, iRow, iCol, iVar) = statsRow.emCI95_lower;
                        emmTop(iGroup, iMember, iRow, iCol, iVar) = statsRow.emCI95_upper;

                        switch posthocMethod

                            case 'emm'
                                barCenter(iGroup, iMember, iRow, iCol, iVar) = statsRow.emMean;
                                barBottom(iGroup, iMember, iRow, iCol, iVar) = statsRow.emCI95_lower;
                                barTop(iGroup, iMember, iRow, iCol, iVar) = statsRow.emCI95_upper;

                            case 'ttest'
                                barCenter(iGroup, iMember, iRow, iCol, iVar) = statsRow.mean;
                                barBottom(iGroup, iMember, iRow, iCol, iVar) = statsRow.mean - statsRow.std;
                                barTop(iGroup, iMember, iRow, iCol, iVar) = statsRow.mean + statsRow.std;

                            case 'utest'
                                barCenter(iGroup, iMember, iRow, iCol, iVar) = statsRow.median;
                                barBottom(iGroup, iMember, iRow, iCol, iVar) = statsRow.q25;
                                barTop(iGroup, iMember, iRow, iCol, iVar) = statsRow.q75;
                        end
                    end

                    % get post-hoc p-values
                    switch posthocMethod

                        case {'ttest', 'utest'} % perform post-hoc analysis using paired t-tests or u-tests

                            % pairwise comparison
                            for iPair = 1:nPairs
                                pair = pairs(iPair, :);

                                % calc contrasts
                                L1 = idxDesc & (Data.(memberVar) == pair(1));
                                L2 = idxDesc & (Data.(memberVar) == pair(2));

                                val1 = Data.(transVar)(L1);
                                val2 = Data.(transVar)(L2);
                                switch posthocMethod
                                    case 'ttest'
                                        [~, bar_p(iGroup, iPair, iRow, iCol, iVar), ~, stats] = ttest2(val1, val2);
                                        tValue = stats.tstat;
                                        df = stats.df;
                                        sPool = stats.sd;
                                        dCohen = (mean(val2, 'omitnan') - mean(val1, 'omitnan')) / sPool;
                                        bar_test(iGroup, iPair, iRow, iCol, iVar) = tValue;
                                        bar_DF(iGroup, iPair, iRow, iCol, iVar) = df;
                                        bar_eff(iGroup, iPair, iRow, iCol, iVar) =  dCohen;
                                        bar_diff(iGroup, iPair, iRow, iCol, iVar) = mean(val2, 'omitnan') - mean(val1, 'omitnan');
                                        bar_diffpct(iGroup, iPair, iRow, iCol, iVar) = bar_diff(iGroup, iPair, iRow, iCol, iVar) / mean(val1, 'omitnan') * 100;
                                    case 'utest'
                                        [bar_p(iGroup, iPair, iRow, iCol, iVar), ~, stats] = ranksum(val1, val2);
                                        N1 = sum(~isnan(val1));
                                        N2 = sum(~isnan(val2));
                                        W = stats.ranksum;
                                        U = W - N1*(N1+1)/2;
                                        r = 1 - 2*U/(N1*N2);
                                        bar_test(iGroup, iPair, iRow, iCol, iVar) = U;
                                        bar_eff(iGroup, iPair, iRow, iCol, iVar) =  r;
                                        bar_diff(iGroup, iPair, iRow, iCol, iVar) = median(val2, 'omitnan') - median(val1, 'omitnan');
                                        bar_diffpct(iGroup, iPair, iRow, iCol, iVar) = bar_diff(iGroup, iPair, iRow, iCol, iVar) / median(val1, 'omitnan') * 100;
                                end

                                % calc main contrasts
                                if posthocMain && nGroups * nRows * nCols * nPairs > 1
                                    L1 = (Data.(memberVar) == pair(1));
                                    L2 = (Data.(memberVar) == pair(2));
                                    val1 = Data.(transVar)(L1);
                                    val2 = Data.(transVar)(L2);
                                    switch posthocMethod
                                        case 'ttest'
                                            [~, main_p(iPair, iRow, iCol, iVar), ~, stats] = ttest2(val1, val2);
                                            tValue = stats.tstat;
                                            df = stats.df;
                                            sPool = stats.sd;
                                            dCohen = (mean(val2, 'omitnan') - mean(val1, 'omitnan')) / sPool;
                                            main_test(iPair, iRow, iCol, iVar) = tValue;
                                            main_DF(iPair, iRow, iCol, iVar) = df;
                                            main_eff(iPair, iRow, iCol, iVar) =  dCohen;
                                            main_diff(iPair, iRow, iCol, iVar) =  mean(val2, 'omitnan') - mean(val1, 'omitnan');
                                            main_diffpct(iPair, iRow, iCol, iVar) = main_diff(iPair, iRow, iCol, iVar) / mean(val1, 'omitnan') * 100;
                                        case 'utest'
                                            [main_p(iPair, iRow, iCol, iVar), ~, stats] = ranksum(val1, val2);
                                            N1 = sum(~isnan(val1));
                                            N2 = sum(~isnan(val2));
                                            W = stats.ranksum;
                                            U = W - N1*(N1+1)/2;
                                            r = 1 - 2*U/(N1*N2);
                                            main_test(iPair, iRow, iCol, iVar) = U;
                                            main_eff(iPair, iRow, iCol, iVar) =  r;
                                            main_diff(iPair, iRow, iCol, iVar) =  median(val2, 'omitnan') - median(val1, 'omitnan');
                                            main_diffpct(iPair, iRow, iCol, iVar) = main_diff(iPair, iRow, iCol, iVar) / median(val1, 'omitnan') * 100;
                                    end
                                end
                            end

                        case 'emm' % perform post-hoc analysis using emmeans

                            % pairwise comparison
                            for iPair = 1:nPairs
                                pair = pairs(iPair, :);

                                % calc contrasts
                                L1 = idxEmm & (emm.table.(memberVar) == pair(1));
                                L2 = idxEmm & (emm.table.(memberVar) == pair(2));
                                L = (L1 - L2)';
                                contrasts = kbcontrasts_wald(mdl, emm, L);
                                bar_p(iGroup, iPair, iRow, iCol, iVar) = contrasts.pVal;
                                bar_test(iGroup, iPair, iRow, iCol, iVar) = contrasts.F;
                                bar_DF(iGroup, iPair, iRow, iCol, iVar) = contrasts.DF1;
                                bar_aux(iGroup, iPair, iRow, iCol, iVar) = contrasts.DF2;
                                bar_smd(iGroup, iPair, iRow, iCol, iVar) = f2smd(contrasts.F, contrasts.DF1, contrasts.DF2);
                                bar_eff(iGroup, iPair, iRow, iCol, iVar) = f2etaSqp(contrasts.F, contrasts.DF1, contrasts.DF2);
                                bar_diff(iGroup, iPair, iRow, iCol, iVar) = mean(mdl.Link.Inverse(contrasts.table.Estimated_Marginal_Mean(L2)) - mdl.Link.Inverse(contrasts.table.Estimated_Marginal_Mean(L1)));
                                bar_diffpct(iGroup, iPair, iRow, iCol, iVar) = bar_diff(iGroup, iPair, iRow, iCol, iVar) / mean(mdl.Link.Inverse(contrasts.table.Estimated_Marginal_Mean(L1))) * 100;
                                bar_emm_1(iGroup, iPair, iRow, iCol, iVar) =  emm.table.Estimated_Marginal_Mean(L1);
                                bar_CI95low_1(iGroup, iPair, iRow, iCol, iVar) =  emm.table.CI_95_0pct(L1, 1);
                                bar_CI95up_1(iGroup, iPair, iRow, iCol, iVar) =  emm.table.CI_95_0pct(L1, 2);
                                bar_emm_2(iGroup, iPair, iRow, iCol, iVar) =  emm.table.Estimated_Marginal_Mean(L2);
                                bar_CI95low_2(iGroup, iPair, iRow, iCol, iVar) =  emm.table.CI_95_0pct(L2, 1);
                                bar_CI95up_2(iGroup, iPair, iRow, iCol, iVar) =  emm.table.CI_95_0pct(L2, 2);

                                % calc main contrasts
                                if posthocMain && nGroups * nRows * nCols * nPairs > 1
                                    L1 = (emm.table.(memberVar) == pair(1));
                                    L2 = (emm.table.(memberVar) == pair(2));
                                    L = (L1 - L2)';
                                    contrasts = kbcontrasts_wald(mdl, emm, L);
                                    main_p(iPair, iRow, iCol, iVar) = contrasts.pVal;
                                    main_test(iPair, iRow, iCol, iVar) = contrasts.F;
                                    main_DF(iPair, iRow, iCol, iVar) = contrasts.DF1;
                                    main_aux(iPair, iRow, iCol, iVar) = contrasts.DF2;
                                    main_smd(iPair, iRow, iCol, iVar) = f2smd(contrasts.F, contrasts.DF1, contrasts.DF2);
                                    main_eff(iPair, iRow, iCol, iVar) = f2etaSqp(contrasts.F, contrasts.DF1, contrasts.DF2);
                                    main_diff(iPair, iRow, iCol, iVar) = mean(mdl.Link.Inverse(contrasts.table.Estimated_Marginal_Mean(L2)) - mdl.Link.Inverse(contrasts.table.Estimated_Marginal_Mean(L1)));
                                    main_diffpct(iPair, iRow, iCol, iVar) = main_diff(iPair, iRow, iCol, iVar) / mean(mdl.Link.Inverse(contrasts.table.Estimated_Marginal_Mean(L1))) * 100;
                                    main_emm_1(iPair, iRow, iCol, iVar) =  emm.table.Estimated_Marginal_Mean(L1);
                                    main_CI95low_1(iPair, iRow, iCol, iVar) =  emm.table.CI_95_0pct(L1, 1);
                                    main_CI95up_1(iPair, iRow, iCol, iVar) =  emm.table.CI_95_0pct(L1, 2);
                                    main_emm_2(iPair, iRow, iCol, iVar) =  emm.table.Estimated_Marginal_Mean(L2);
                                    main_CI95low_2(iPair, iRow, iCol, iVar) =  emm.table.CI_95_0pct(L2, 1);
                                    main_CI95up_2(iPair, iRow, iCol, iVar) =  emm.table.CI_95_0pct(L2, 2);
                                end
                            end
                    end
                end
            end
        end

        %% Save descriptive statistics

        % display and save statistics (only necessary for 1st posthoc level)
        if iLevel == 1
            % descriptive statistics
            disp(StatsTable);
            tableName = 'Statistics';
            % add table to results
            results.(myVar).(tableName) = StatsTable;
            % save table to file
            if isSave                
                saveTable(StatsTable, tableName, {'xlsx'}, outSubDir);
            end
        end
    end

    %% Statistical correction of posthoc comparisons
    % Holm-Bonferroni correction of all p-values

    if multiVariate
        nCorrs = 1;
    else
        nCorrs = nY;
    end

    bar_pCorr = bar_p;
    main_pCorr = main_p;
    for iVar = 1:nCorrs
        if multiVariate
            myBar_p = bar_p;
            myBar_pCorr = myBar_p; % create array of corrected p-values
        else
            myBar_p = bar_p(:, :, :, :, iVar);
            myBar_pCorr = myBar_p; % create array of corrected p-values
        end
        idxDesc = ~isnan(myBar_p(:)); % identify NaN-entries
        if ~isempty(idxDesc)
            sizeOrig = size(myBar_p); % store original array shape
            myBar_p = myBar_p(:); % make column vector

            % check for duplicate p-values indicating a severe problem
            for iP = 1:length(myBar_p)
                pValue = myBar_p(iP);
                % Problem: When there are very small p-values, they count as "same" in
                % certain algorithms such as setdiff, which we use to decide if there are
                % duplicate p-values. So, we forcus on p-values sufficiently greater
                % than zero.
                myBar_p_nonzero = myBar_p(myBar_p > tol);
                if length(setdiff(myBar_p_nonzero, pValue)) < length(myBar_p_nonzero) - 1
                    warning('Duplicate p-values found -> GLM fit is not good.');
                    break
                end
            end

            myBar_pCorr = myBar_pCorr(:); % make column vector
            switch posthocCorrection
                case 'none'
                    % do nothing
                case 'holm'
                    [~, myBar_pCorr(idxDesc)] = bonferroni_holm(myBar_p(idxDesc)); % correct p-values, omitting NaNs
                case {'bonf', 'sidak'}
                    nTests = sum(idxDesc);
                    myBar_pCorr(idxDesc) = sidak_corr(myBar_pCorr(idxDesc), nTests);
            end
            % myBar_pCorr(idxDesc) = sidak_corr(myBar_pCorr(idxDesc), nPosthocLevels); % correct for posthoc levels
            myBar_pCorr(idxDesc) = sidak_corr(myBar_pCorr(idxDesc), correctForN); % correct for custom number of tests
            myBar_pCorr = reshape(myBar_pCorr, sizeOrig); % bring corrected p-value array into the same shape as p-value array
            if multiVariate
                bar_pCorr = myBar_pCorr;
            else
                bar_pCorr(:, :, :, :, iVar) = myBar_pCorr;
            end
        end

        % statistical correction of main posthoc p-Values
        if posthocMain && nGroups * nRows * nCols * nPairs > 1            
            if multiVariate
                myMain_p = main_p;
                myMain_pCorr = myMain_p; % create array of corrected p-values
            else
                myMain_p = main_p(:, :, :, iVar);
                myMain_pCorr = myMain_p; % create array of corrected p-values
            end

            idxDesc = ~isnan(myMain_p(:)); % identify NaN-entries
            if ~isempty(idxDesc)
                sizeOrig = size(myMain_p); % store original array shape
                myMain_p = myMain_p(:); % make column vector
                myMain_pCorr = myMain_pCorr(:); % make column vector
                [~, myMain_pCorr(idxDesc)] = bonferroni_holm(myMain_p(idxDesc)); % correct p-values, omitting NaNs
                % myMain_pCorr(idxDesc) = sidak_corr(myMain_pCorr(idxDesc), nPosthocLevels); % correct for posthoc levels
                myMain_pCorr(idxDesc) = sidak_corr(myMain_pCorr(idxDesc), correctForN); % correct for custom number of tests
                myMain_pCorr = reshape(myMain_pCorr, sizeOrig); % restore original dimensions of p-value array
            end
            if multiVariate
                main_pCorr = myMain_pCorr;
            else
                main_pCorr(:, :, :, iVar) = myMain_pCorr;
            end
        end

    end

    %% Loop over dependent variables

    for iVar = 1:nY
        % dependent variable
        myVar = y{iVar};

        % create output folder
        outSubDir = sprintf('%s/%s', outDir, myVar);
        if isSave && ~isfolder(outSubDir)
            mkdir(outSubDir);
        end

        % yLabel
        if length(yLabel) == nY
            myLabel = yLabel{iVar};
        elseif isscalar(yLabel)
            myLabel = yLabel{1};
        else
            error('The number of dependent variables does not match the number of their labels');
        end

        % yUnits
        if length(yUnits) == nY
            myUnits = yUnits{iVar};
        elseif isscalar(yUnits)
            myUnits = yUnits{1};
        end

        % figure size
        panelWidth = nGroups * 300; % width of each panel
        panelHeight = 300; % height of each panel

        %% Plot data

        if isPlot

            figWidth = nCols * panelWidth;
            figHeight = nRows * panelHeight + 50;
            if nPosthocLevels > 1
                figName = sprintf('DataPlots_%d', iLevel);
            else
                figName = 'DataPlots';
            end
            fig = figure('Name', figName, 'Position', [0, 0, figWidth, figHeight]);
            layout = tiledlayout(nRows, nCols);
            if nY > 1 && ~multiVariate
                if isempty(plotTitle)
                    title(layout, capitalize(sprintf('%s', myLabel)), 'interpreter', 'none', 'FontWeight', 'bold', 'FontSize', 14);
                else
                    title(layout, capitalize(sprintf('%s (%s)', plotTitle, myLabel)), 'interpreter', 'none', 'FontWeight', 'bold', 'FontSize', 14);
                end
            else
                if isempty(plotTitle)
                    title(layout, capitalize(sprintf('%s', myLabel)), 'interpreter', 'none', 'FontWeight', 'bold', 'FontSize', 14);
                else
                    title(layout, capitalize(sprintf('%s (%s)', plotTitle, myLabel)), 'interpreter', 'none', 'FontWeight', 'bold', 'FontSize', 14);
                end
            end

            % prepare to display variable names and levels
            displayMemberVar = memberVarDisp;
            displayMembers = members;
            displayGroupVar = groupVarDisp;
            displayGroups = groups;
            displayColVar = colVarDisp;
            displayCols = cols;
            displayRowVar = rowVarDisp;
            displayRows = rows;
            % capitalize names and/or levels, if needed
            if contains(showVarNames, 'Names')
                displayMemberVar = string(capitalize(memberVarDisp));
                displayGroupVar = string(capitalize(groupVarDisp));
                displayColVar = string(capitalize(colVarDisp));
                displayRowVar = string(capitalize(rowVarDisp));
            end
            if contains(showVarNames, 'Levels')
                displayMembers = string(strsplit(capitalize(strjoin(cellstr(members), ', '), 'all'), ', '));
                displayGroups = string(strsplit(capitalize(strjoin(cellstr(groups), ', '), 'all'), ', '));
                displayCols = string(strsplit(capitalize(strjoin(cellstr(cols), ', '), 'all'), ', '));
                displayRows = string(strsplit(capitalize(strjoin(cellstr(rows), ', '), 'all'), ', '));
            end
            
            % plot grouped data
            for iRow = 1:nRows

                for iCol = 1:nCols

                    % % start panel subplot
                    iPanel = sub2ind([nCols, nRows], iCol, iRow);
                    % subplot(nRows, nCols, iPanel);
                    panel = tiledlayout(layout, 1, 1);
                    panel.Layout.Tile = iPanel;
                    panel.Layout.TileSpan = [1 1];
                    ax = nexttile(panel);
                    ax.XAxis.TickValues = [];
                    ax.YAxis.TickValues = [];

                    % plot panel
                    if isempty(myUnits)
                        yLabelStr = sprintf('%s', myLabel);
                    else
                        yLabelStr = sprintf('%s [%s]', myLabel, myUnits);
                    end

                    if contains(showVarNames, 'names', 'IgnoreCase', true) && contains(showVarNames, 'levels', 'IgnoreCase', true) % display variable names and levels
                        if length(cols) > 1 && length(rows) > 1
                            panelTitle = sprintf('%s = %s, %s = %s', displayRowVar, displayRows(iRow), displayColVar, displayCols(iCol));
                        elseif length(rows) > 1
                            panelTitle = sprintf('%s = %s', displayRowVar, displayRows(iRow));
                        elseif length(cols) > 1
                            panelTitle = sprintf('%s = %s', displayColVar, displayCols(iCol));
                        else
                            panelTitle = '';
                        end
                    else % display only levels
                        if length(cols) > 1 && length(rows) > 1
                            panelTitle = sprintf('%s, %s', displayRows(iRow), displayCols(iCol));
                        elseif length(rows) > 1
                            panelTitle = sprintf('%s', displayRows(iRow));
                        elseif length(cols) > 1
                            panelTitle = sprintf('%s', displayCols(iCol));
                        else
                            panelTitle = '';
                        end
                    end

                    myDataPoints = DataPoints(:, :, :, iRow, iCol, iVar);
                    myBarCenter = emmCenter(:, :, iRow, iCol, iVar);
                    myBarBottom = emmBottom(:, :, iRow, iCol, iVar);
                    myBarTop = emmTop(:, :, iRow, iCol, iVar);
                    switch barType
                        case 'auto'
                            switch posthocMethod
                                case 'emm'
                                    myBarType = 'custom';
                                case 'ttest'
                                    myBarType = 'mean';
                                case 'utest'
                                    myBarType = 'median';
                                otherwise
                                    myBarType = 'emm';
                            end
                        otherwise
                            myBarType = barType;
                    end
                    plotGroups(myDataPoints, displayMembers, displayGroups, displayMemberVar, displayGroupVar, bar_pCorr(:, :, iRow, iCol, iVar), panelTitle, yLabelStr, plotStyle, panel, showVarNames, markerSize, myBarType, myBarCenter, myBarBottom, myBarTop, plotLines, sortValues);
                end
            end

            % rescale plots to achieve the same scale for all panels
            if rescale
                axs = findobj(fig, 'type', 'axes');
                ylimits = NaN(length(axs), 2);
                for iAx = 1:length(axs)
                    if ~isempty(axs(iAx).XAxis.TickValues)
                        ylimits(iAx, :) = get(axs(iAx), 'ylim');
                    end
                end
                maxYlimits(1, 1) = min(ylimits(:, 1));
                maxYlimits(1, 2) = max(ylimits(:, 2));
                for iAx = 1:length(axs)
                    set(axs(iAx), 'ylim', maxYlimits);
                end
            end

            % save figure
            if isSave
                kbsaveFigure(fig, fullfile(outSubDir, figName), {'.pdf', '.png', '.fig'});
            end

            % close figure
            if isCloseFigures
                close(fig);
            end
        end

        %% Post-hoc table

        posthocTable = table;

        % posthoc main effects
        if posthocMain && nGroups * nRows * nCols * nPairs > 1
            for iPair = 1:nPairs
                tableRow = table;

                if nRows > 1
                    tableRow.(rowVar) = "any";
                end

                if nCols > 1
                    tableRow.(colVar) = "any";
                end

                if nGroups > 1
                    tableRow.(groupVar) = "any";
                end

                tableRow.([memberVar, '_1']) = string(pairs(iPair, 1));
                tableRow.([memberVar, '_2']) = string(pairs(iPair, 2));
                tableRow.emmCI_1 = string(sprintf('%.2f (%.2f, %.2f)', main_emm_1(iPair, iRow, iCol, iVar), main_CI95low_1(iPair, iRow, iCol, iVar), main_CI95up_1(iPair, iRow, iCol, iVar)));
                tableRow.emmCI_2 = string(sprintf('%.2f (%.2f, %.2f)', main_emm_2(iPair, iRow, iCol, iVar), main_CI95low_2(iPair, iRow, iCol, iVar), main_CI95up_2(iPair, iRow, iCol, iVar)));
                % tableRow.emm_1 = main_emm_1(iPair, iRow, iCol, iVar);
                % tableRow.CI95low_1 = main_CI95low_1(iPair, iRow, iCol, iVar);
                % tableRow.CI95up_1 = main_CI95up_1(iPair, iRow, iCol, iVar);
                % tableRow.emm_2 = main_emm_2(iPair, iRow, iCol, iVar);
                % tableRow.CI95low_2 = main_CI95low_2(iPair, iRow, iCol, iVar);
                % tableRow.CI95up_2 = main_CI95up_2(iPair, iRow, iCol, iVar);
                tableRow.p = main_p(iPair, iRow, iCol, iVar);
                tableRow.pCorr = main_pCorr(iPair, iRow, iCol, iVar);
                tableRow.diff = main_diff(iPair, iRow, iCol, iVar);
                tableRow.diffpct = main_diffpct(iPair, iRow, iCol, iVar);                
                switch posthocMethod
                    case 'ttest'
                        tableRow.t = main_test(iPair, iRow, iCol, iVar);
                        tableRow.DF = main_DF(iPair, iRow, iCol, iVar);
                        tableRow.d = main_eff(iPair, iRow, iCol, iVar);
                        tableRow.effectSize = string(effprint(main_eff(iPair, iRow, iCol, iVar), 'd'));
                    case 'utest'
                        tableRow.U = main_test(iPair, iRow, iCol, iVar);
                        tableRow.r = main_eff(iPair, iRow, iCol, iVar);
                        tableRow.effectSize = string(effprint(main_eff(iPair, iRow, iCol, iVar), 'r'));
                    case 'emm'
                        tableRow.F = main_test(iPair, iRow, iCol, iVar);
                        tableRow.DF1 = main_DF(iPair, iRow, iCol, iVar);
                        tableRow.DF2 = main_aux(iPair, iRow, iCol, iVar);
                        tableRow.SMD = main_smd(iPair, iRow, iCol, iVar);
                        tableRow.etaSqp = main_eff(iPair, iRow, iCol, iVar);
                        tableRow.effectSize = string(effprint(main_eff(iPair, iRow, iCol, iVar), 'eta2'));
                end
                tableRow.significance = string(sigprint(tableRow.pCorr));
                posthocTable = [posthocTable; tableRow]; %#ok<AGROW>
            end
        end

        % posthoc pairwise comparisons
        for iPair = 1:nPairs
            for iRow = 1:nRows
                for iCol = 1:nCols
                    for iGroup = 1:nGroups
                        tableRow = table;

                        if nRows > 1
                            tableRow.(rowVar) = string(rows(iRow));
                        end

                        if nCols > 1
                            tableRow.(colVar) = string(cols(iCol));
                        end

                        if nGroups > 1
                            tableRow.(groupVar) = string(groups(iGroup));
                        end

                        tableRow.([memberVar, '_1']) = string(pairs(iPair, 1));
                        tableRow.([memberVar, '_2']) = string(pairs(iPair, 2));
                        tableRow.emmCI_1 = string(sprintf('%.2f (%.2f, %.2f)', bar_emm_1(iGroup, iPair, iRow, iCol, iVar), bar_CI95low_1(iGroup, iPair, iRow, iCol, iVar), bar_CI95up_1(iGroup, iPair, iRow, iCol, iVar)));
                        tableRow.emmCI_2 = string(sprintf('%.2f (%.2f, %.2f)', bar_emm_2(iGroup, iPair, iRow, iCol, iVar), bar_CI95low_2(iGroup, iPair, iRow, iCol, iVar), bar_CI95up_2(iGroup, iPair, iRow, iCol, iVar)));
                        % tableRow.emm_1 = bar_emm_1(iGroup, iPair, iRow, iCol, iVar);
                        % tableRow.CI95low_1 = bar_CI95low_1(iGroup, iPair, iRow, iCol, iVar);
                        % tableRow.CI95up_1 = bar_CI95up_1(iGroup, iPair, iRow, iCol, iVar);
                        % tableRow.emm_2 = bar_emm_2(iGroup, iPair, iRow, iCol, iVar);
                        % tableRow.CI95low_2 = bar_CI95low_2(iGroup, iPair, iRow, iCol, iVar);
                        % tableRow.CI95up_2 = bar_CI95up_2(iGroup, iPair, iRow, iCol, iVar);
                        tableRow.p = bar_p(iGroup, iPair, iRow, iCol, iVar);
                        tableRow.pCorr = bar_pCorr(iGroup, iPair, iRow, iCol, iVar);
                        tableRow.diff = bar_diff(iGroup, iPair, iRow, iCol, iVar);
                        tableRow.diffpct = bar_diffpct(iGroup, iPair, iRow, iCol, iVar);
                        switch posthocMethod
                            case 'ttest'
                                tableRow.t = bar_test(iGroup, iPair, iRow, iCol, iVar);
                                tableRow.DF = bar_DF(iGroup, iPair, iRow, iCol, iVar);
                                tableRow.d = bar_eff(iGroup, iPair, iRow, iCol, iVar);
                                tableRow.effectSize = string(effprint(bar_eff(iGroup, iPair, iRow, iCol, iVar), 'd'));
                            case 'utest'
                                tableRow.U = bar_test(iGroup, iPair, iRow, iCol, iVar);
                                tableRow.r = bar_eff(iGroup, iPair, iRow, iCol, iVar);
                                tableRow.effectSize = string(effprint(bar_eff(iGroup, iPair, iRow, iCol, iVar), 'r'));
                            case 'emm'
                                tableRow.F = bar_test(iGroup, iPair, iRow, iCol, iVar);
                                tableRow.DF1 = bar_DF(iGroup, iPair, iRow, iCol, iVar);
                                tableRow.DF2 = bar_aux(iGroup, iPair, iRow, iCol, iVar);
                                tableRow.SMD = bar_smd(iGroup, iPair, iRow, iCol, iVar);
                                tableRow.etaSqp = bar_eff(iGroup, iPair, iRow, iCol, iVar);
                                tableRow.effectSize = string(effprint(bar_eff(iGroup, iPair, iRow, iCol, iVar), 'eta2'));
                        end
                        tableRow.significance = string(sigprint(tableRow.pCorr));
                        posthocTable = [posthocTable; tableRow]; %#ok<AGROW>
                    end
                end
            end
        end

        % display table
        disp(posthocTable);

        if nPosthocLevels > 1
            tableName = sprintf('Posthoc_%d', iLevel);
        else
            tableName = 'Posthoc';
        end

        % add table to results
        results.(myVar).(tableName) = posthocTable;

        % save table
        if isSave            
            saveTable(posthocTable, tableName, {'xlsx'}, outSubDir);
        end
    end
end

fprintf('DONE in %g seconds\n', toc(ticAnova));
fprintf('Results saved in %s\n', outDir);

end


%% Helper functions %%%%%%%%%%%%%%%%%%%%%%

function value = getValue(in)
% get numerical, logical, or string value from string
if ischar(in) || isstring(in)
    [numValue, isNum] = str2num(in);
    if isNum
        value = numValue;
    else
        value = in;
    end
else % bool or whatever
    value = in;
end
end

function value = getList(in)
% get numerical, logical, or string array from string
if ischar(in) || isstring(in)
    [numValue, isNum] = str2num(in);
    if isNum
        value = numValue;
    else
        value = strtrim(strsplit(in, {',', ';'}));
    end
else % bool or whatever
    value = in;
end
end