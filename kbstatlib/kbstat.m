function mdl = kbstat(options)
%% Analyze data using generalized linear mixed-model fit.
% The result is a fit summary, a diagnostic plot, a data bar plot, a
% descriptive statistics table, an ANOVA table, and a posthoc pairwise
% comparison table.
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
% options.separateMulti=false. However, the posthoc test then only works
% when options.posthocMethod='ttest' instead of 'emm', due to limitations
% of the external package 'emmeans'.
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
%                       If "separateMulti" is set to "false", then y is
%                       treated as a vector-valued dependent variable and a
%                       multivariate analysis is performed.
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
%       x               Comma-separated list of independent variables. Up
%                       to 4 independent variables are supported.
%                       String-valued variables are considered as
%                       categorical, i.e. as factors. Only factors are
%                       included in the barplot chart. Numerical variables
%                       are not considered factors, except they are
%                       included in the "catVar" list, see below.
%
%       catVar          Comma-separated list of (dependent or independent)
%                       variables that are taken as categorical, i.e. as
%                       factors, even if they have numerical values.
%                       OPTIONAL, default = '';
%
%       id              Name of the subject variable.
%                       OPTIONAL, default = ''.
%
%       within          Comma-separated list of within-subject variables.  
%                       Within-subject variables are nested within the
%                       subject variable, i.e. they vary for each subject.
%                       Variables that are not declared as within-subject,
%                       are considered between-subject, i.e. they vary only
%                       between, not within, subjects. "within" can be a
%                       subset of "x", or else its members are added to
%                       "x". Example:
%                           options.id = 'subject'
%                           options.x = 'time, age, sex'
%                           options.within = 'dose'.
%                       OPTIONAL, default = ''.
%
%       interact        Comma-separated list of variables whose interaction
%                       with each other is to be analyzed. Can be a subset
%                       of "x", or else its members are added to "x". When
%                       not set, all members of x are assumed to mutually
%                       interact. Example:
%                           options.id = 'subject'
%                           options.x = 'time, dose'
%                           options.interact = 'dose, age'.
%                       OPTIONAL, default = options.x
%
%       coVar           Comma-separated list of (continuous or categorical) 
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
%       separateMulti   Flag if, when the dependent variable has multiple
%                       components, these components should be analyzed
%                       separately. In the latter case, the results are
%                       statistically corrected for these multiple tests.
%                       OPTIONAL, default = true.
%
%       randomSlopes    Flag if random slopes should be estimated.
%                       OPTIONAL, default = true.
%
%       formula         Formula in Wilkinson Notation. If given, it
%                       overrides the automatically generated formula
%                       produced from the provided variables.
%
%       fitMethod       Fit method used for the GLM fit.%
%                       Possible values:
%                           'none'                  Skip linear model fit
%                           'MPL'                   Maximum pseudo likelihood
%                           'REMPL'                 Restricted maximum pseudo likelihood
%                           'Laplace'               Maximum likelihood using Laplace approximation
%                           'ApproximateLaplace'    Maximum likelihood using approximate Laplace approximation with fixed effects profiled out
%                       OPTIONAL, default = 'REMPL'.
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
%                           'identity'	    g(mu) = mu.             Default for Normal distribution
%                           'log'	        g(mu) = log(mu).        Default for Poisson
%                           'logit'	        g(mu) = log(mu/(1-mu))  Default for Binomial distribution
%                           'loglog'	    g(mu) = log(-log(mu))
%                           'probit'	    g(mu) = norminv(mu)
%                           'comploglog'	g(mu) = log(-log(1-mu))
%                           'reciprocal'	g(mu) = mu.^(-1)        Default for Gamma
%                           Scalar p	    g(mu) = mu.^p           Default for InverseGaussian (p= -2)
%                           'auto'          depends on chosen distribution
%                       OPTIONAL, default = 'auto'.
%
%       posthocMethod   Method for the posthoc pairwise comparison.
%                       Possible values:
%                       'none'      Do not perform posthoc analysis
%                       'ttest'     t-test with Holm-Bonferroni correction
%                       'utest'     Mann-Whitney u-Test (ranksum test) with Holm-Bonferroni correction
%                       'emm'       Extract contrasts from linear model fit
%                       'auto'      Perform posthoc analysis using
%                           'emm'   if univariate or multi-valued y with separateMulti = true
%                           'ttest' if distribution = 'normal'
%                           'utest' otherwise
%                       OPTIONAL, default = 'auto'.
%
%       posthocMainEffects  Flag if also the posthoc main effects should be
%                       calculated, i.e. the comparison between one
%                       variable set to 'any'.
%                       OPTIONAL, default = true.
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
%       outlierMethod   Method to remove post-fit outliers from the data.
%                       Possible values:
%                       'none'      Do not remove outliers
%                       'quartiles' Remove values outside 1.5 times the
%                                   interquartile range [.25, .75]
%                       'median'    Remove values outside more than three
%                                   scaled median absolute deviations (MAD)
%                       'mean'      Remove values outside 3 standard
%                                   deviations from the mean.
%                       OPTIONAL, default = 'quartiles'.
%
%       preOutlierMethod   Method to remove pre-fit outliers from the data.
%                       Possible values: see outlierMethod.
%                       OPTIONAL, default = 'quartiles'.
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
%                       OPTIONAL, default = unset.
%
%       transform       Choose how to transform the dependent variable.
%                       The linear model is fit on the transformed data,
%                       but the data are plotted using the original data.
%                       OPTIONAL, default = ''. Possible values
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
%       outDir          Output folder for generated files.
%                       OPTIONAL, defaults to the parent folder of the
%                       input file. If the input is a data table, the
%                       output folder defaults to the local folder.
%
%       isPlot          Plot data as grouped bars with significance brackets
%       errorBars       What the error bars indicate. This also defines
%                       what the bar height indicates
%                       'std'   Standard deviation. Bar height is mean
%                       'se'    Standard error. Bar height is the mean
%                       'ci95'  95% confidence interval. Bar height is the
%                               mean
%                       'q25'   25% and 75% quantiles. Bar height is the
%                               median, which is the 50% quantile
%       levelOrder      The order in which the levels are displayed in the
%                       plots.
%                       'sorted'    sorted alphanumerically (default)
%                       'stable'    sorted in the order of occurrence in
%                                   the data table
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
%       showVarNames    Flag how variables and levels are displayed in data
%                       plots.
%                       0 = dispaly only variable levels
%                       1 = display variable names and levels
%                       2 = display capitalized variable names and levels
%                       3 = do not display variable names but display
%                           capitalized levels.
%                       OPTIONAL, default = 1.
%
%       xOrder          Ordering of the items on the x axis in data plots.
%                       Overrides ordering of the level names of the 1st
%                       independent variable.
%                       OPTIONAL, default = [].
%                       Example: options.xOrder = '[1 3 2]'
%
% OUTPUT
%       mdl             (Generalized linear mixed-effects model) Result
%                       from linear model fit
%
% EXAMPLE
%   options.y = 'jointEfficiency';
%   options.yUnits = '1';
%   options.x = 'shoe, speed, sex, joint';
%   options.id = 'subject';
%   options.within = 'shoe, speed, joint';
%   options.interact = 'shoe, speed';
%   options.distribution = 'gamma';
%	options.constraint = 'speed < 2 & joint == "ankle_joint"'
%   kbstat('path/to/Data.xlsx', options);
%
% Author: Kim Joris Boström

%% Get parameters

if ~isstruct(options)
    error('Input argument must be a struct');
end

inDir = '.';
inName = '';

% load Data
if isfield(options, 'inFile') && isfile(options.inFile)
    inFile = options.inFile;
    [inDir, inName, ~] = fileparts(inFile);
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

% catVar
if isfield(options, 'catVar') && ~isempty(options.catVar)
    catVar = strtrim(strsplit(options.catVar, {',', ';'}));
else
    catVar = {};
end

% independent variable(s)
x = strtrim(strsplit(options.x, {',', ';'}));

%% dependent variable(s)

y = cellstr(options.y);
nY = length(y);
for iY = 1:nY
    myY = y{iY};
    if ~ismember(myY, tableVars)
        error('Dependent variable "%s" not found in data table', myY);
    end
end

% multiVar
if isfield(options, 'multiVar') && ~isempty(options.multiVar)
    multiVar = options.multiVar;
else
    multiVar = '';
end

% multiVarLevels
if isfield(options, 'multiVarLevels') && ~isempty(options.multiVarLevels)
    multiVarLevels = strtrim(strsplit(options.multiVarLevels, ','));
else
    multiVarLevels = '';
end

if nY > 1
    yVar = 'yVar';
    yVal = 'Y';
    Data2 = stack(Data1, y, 'NewDataVariableName', yVal, 'IndexVariableName', yVar);
    Data2.(yVar) = categorical(string(Data2.(yVar)));
    depVar = yVal; % set dependent variable to y
elseif ~isempty(multiVar)
    yVar = multiVar;
    yVal = y{1};
    if isempty(multiVarLevels)
        Data2 = Data1;
    else
        idx = ismember(string(Data1.(yVar)), string(multiVarLevels));
        Data2 = Data1(idx, :);
    end
    depVar = y{1}; % set dependent variable to y
    Data2.(yVar) = categorical(string(Data2.(yVar)));
    y = cellstr(unique(Data2.(yVar)));
    nY = length(y);
else
    Data2 = Data1;
    depVar = y{1};
    yVal = y{1};
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

% subject variable
if isfield(options, 'id') && ~isempty(options.id)
    id = options.id;
else
    id = '';
end

% within-subject variables
if isfield(options, 'within') && ~isempty(options.within) % option provided and not empty
    within = strtrim(strsplit(options.within, {',', ';'}));
    x = union(x, within, 'stable'); % add within-subject variables to independent variables
else % option not provided or empty
    within = {};
end

% interaction variables
if isfield(options, 'interact') && ~isempty(options.interact) % option provided and not empty
    interact = strtrim(strsplit(options.interact, {',', ';'}));
    x = union(x, interact, 'stable'); % add interaction variables to independent variables
elseif isfield(options, 'interact') && isempty(options.interact) % option provided as empty
    interact = {};
else % option not provided
    interact = x;
end

% covariates
if isfield(options, 'coVar') && ~isempty(options.coVar) % option provided and not empty
    coVar = strtrim(strsplit(options.coVar, {',', ';'}));
else % option not provided
    coVar = {};
end

% random slopes
if isfield(options, 'randomSlopes') && ~isempty(options.randomSlopes)
    randomSlopes = getValue(options.randomSlopes);
else
    randomSlopes = true;
end

% formula
if isfield(options, 'formula') && ~isempty(options.formula)
    formula = options.formula;
else
    formula = '';
end

% separateMulti
if isfield(options, 'separateMulti') && ~isempty(options.separateMulti)
    separateMulti = getValue(options.separateMulti);
else
    separateMulti = true;
end

% fit method
if isfield(options, 'fitMethod') && ~isempty(options.fitMethod)
    fitMethod = options.fitMethod;
else
    fitMethod = 'REMPL';
end

% posthoc method
if isfield(options, 'posthocMethod') && ~isempty(options.posthocMethod)
    posthocMethod = options.posthocMethod;
else
    posthocMethod = 'auto';
end

% posthoc main effects
if isfield(options, 'posthocMainEffects') && ~isempty(options.posthocMainEffects)
    posthocMainEffects = getValue(options.posthocMainEffects);
else
    posthocMainEffects = true;
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
    link = cellstr(options.link);
else % no parameter given
    link = {'auto'};
end

% level order
if isfield(options, 'levelOrder') && ~isempty(options.levelOrder)
    levelOrder = options.levelOrder;
else
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

if isfield(options, 'errorBars') && ~isempty(options.errorBars)
    errorBars = lower(options.errorBars);
else
    errorBars = 'std';
end

% flag to rescale all panel plots to the same y-scale
if isfield(options, 'rescale') && ~isempty(options.rescale)
    rescale = getValue(options.rescale);
else
    rescale = true;
end

% flag if post-fit outliers should be removed
if isfield(options, 'outlierMethod') && ~isempty(options.outlierMethod)
    outlierMethod = options.outlierMethod;
else
    outlierMethod = 'quartiles';
end

% flag if pre-fit outliers should be removed
if isfield(options, 'preOutlierMethod') && ~isempty(options.preOutlierMethod)
    preOutlierMethod = options.preOutlierMethod;
else
    preOutlierMethod = 'quartiles';
end

% output folder
if isfield(options, 'outDir') && ~isempty(options.outDir)
    outDir = options.outDir;
else
    outDir = fullfile(inDir, inName);
end
if ~isfolder(outDir)
    mkdir(outDir);
end

if isfield(options, 'title') && ~isempty(options.title)
    plotTitle = options.title;
else
    plotTitle = depVar;
end

if isfield(options, 'showVarNames') && ~isempty(options.showVarNames)
    showVarNames = getValue(options.showVarNames);
else
    showVarNames = 1;
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

if isfield(options, 'xOrder') && ~isempty(options.xOrder)
    xOrder = getValue(options.xOrder);
else
    xOrder = [];
end

%% Apply constraint, if given

if isfield(options, 'constraint') && ~isempty(options.constraint)
    constraint = options.constraint;
    [conds, ops] = strsplit(constraint, {'&', '|'});
    allIdx = true(size(Data2, 1), 1);
    auxVars = {};
    for iCond = 1:length(conds)
        cond = strtrim(conds{iCond});
        [parts, matches] = strsplit(cond, {'=', '==', '~=', '<', '<=', '>', '>='});
        constraintVar = strtrim(parts{1});
        constraintVal = strtrim(parts{2});
        constraintVal = strrep(constraintVal, '''', ''); % remove single quotes
        constraintVal = strrep(constraintVal, '"', ''); % remove double quotes
        [num, isNum] = str2num(constraintVal);
        if isNum
            constraintVal = num;
        else
            constraintVal = string(constraintVal);
        end
        compVar = strtrim(matches{1});
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
                idx = (constraintVals == constraintVal);
            case '~='
                idx = (constraintVals ~= constraintVal);
            case '<'
                idx = (constraintVals < constraintVal);
            case '<='
                idx = (constraintVals <= constraintVal);
            case '>'
                idx = (constraintVals > constraintVal);
            case '>='
                idx = (constraintVals >= constraintVal);
            otherwise
                error('No valid comparison operator ''%s''', compVar);
        end
        if  iCond > 1 && iCond-1 <= length(ops)
            op = ops{iCond-1};
            cmd = sprintf('allIdx %s idx', op);
            allIdx = eval(cmd);
        else
            allIdx = idx;
        end
    end
    if any(allIdx)
        Data2 = Data2(allIdx, :);
    else
        warning('Constraint ''%s'' cannot be fulfilled -> leave data unchanged', constraint);
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

%% Save options to file

% fpath = fullfile(outDir, 'options.mat');
% save(fpath, 'options');
fpath = fullfile(outDir, 'Options.txt');
fidOptions = fopen(fpath, 'w+');
fields = fieldnames(options);
fields = setdiff(fields, 'Data', 'stable'); % remove "Data" from options
fields = setdiff(fields, 'DataRaw', 'stable'); % remove "DataRaw" from options
paramFields = sort(fields); % sort fields alphabetically
for iField = 1:length(paramFields)
    field = paramFields{iField};
    fprintf(fidOptions, 'options.%s = %s;\n', field, mat2str(string(options.(field))));
end
fclose(fidOptions);

%% Create Data table

% init Data table
Data = table;

% subject variable
if ~isempty(id)
    Data.(id) = categorical(string(Data2.(id))); % make categorical
end

% make catVar variables categorical
for iVar = 1:length(catVar)
    myVar = catVar{iVar};
    Data.(myVar) = categorical(string(Data2.(myVar))); % make categorical
end

% get independent variables
IVs = union(x, coVar, 'stable');
catVars = catVar;
factors = {};
for iIV = 1:length(IVs)
    myIV = IVs{iIV};
    myLevels = unique(Data2.(myIV));
    if ismember(myIV, catVars) || ~all(isnumeric(myLevels))
        Data.(myIV) = categorical(string(Data2.(myIV)));
        catVars = union(catVars, myIV, 'stable');
        if ismember(myIV, x)
            factors = union(factors, myIV, 'stable');
        end
    else % data are not in catVars and they are numerical
        Data.(myIV) = Data2.(myIV); % keep continuous values
    end
end
nFactors = length(factors);

% multivariate variable
if nY > 1
    Data.(yVar) = Data2.(yVar);
end

% get dependent variable
Data.(depVar) = Data2.(depVar);
Data.(depVar) = Data2.(depVar);

% store Data table in options
options.Data = Data;

%% Create variables

% 1st factor = member variable
memberVar = factors{1};
members = unique(Data.(memberVar), levelOrder);
nMembers = length(members);
% apply potential re-ordering
if ~isempty(xOrder)
    members = members(xOrder);
end

% 2nd factor = group variable
if nFactors > 1
    groupVar = factors{2};
    groups = unique(Data.(groupVar), levelOrder);
    nGroups = length(groups);
else
    groupVar = '';
    groups = string(factors{1});
    nGroups = 1;
end

% 3rd factor = column variable
if nFactors > 2
    colVar = factors{3};
    cols = unique(Data.(colVar), levelOrder);
    nCols = length(cols);
else
    colVar = 'none';
    cols = "";
    nCols = 1;
end

% 4th factors = row variable
if nFactors > 3
    rowVar = factors{4};
    rows = unique(Data.(rowVar), levelOrder);
    nRows = length(rows);
else
    rowVar = 'none';
    rows = "";
    nRows = 1;
end

%% Apply data transformation, if given

if ~isempty(transform)
    transVar = sprintf('%sTrans', depVar);
    for iVar = 1:nY
        if nY > 1
            idxDep = (Data.(yVar) == y{iVar});
        else
            idxDep = true(size(Data, 1), 1);
        end
        switch transform
            case 'mean'
                transformFcn = @(x) x/mean(x, 'omitnan');
            case 'std'
                transformFcn = @(x) x/std(x, 'omitnan');
            case 'Z'
                transformFcn = @(x) (x - mean(x, 'omitnan'))/std(x, 'omitnan');
            case {'median'}
                upperThresh = quantile(Data.(depVar)(idxDep), 0.5);
                transformFcn = @(x) x/upperThresh;
            case 'IQR'
                [~, lowerThresh, upperThresh, ~] = isoutlier(Data.(depVar)(idxDep), 'quartiles');
                transformFcn = @(x) x/(upperThresh - lowerThresh);
            case 'IQRmax'
                [~, ~, upperThresh, ~] = isoutlier(Data.(depVar)(idxDep), 'quartiles');
                transformFcn = @(x) x/upperThresh;
            case 'MAD'
                [~, lowerThresh, upperThresh, ~] = isoutlier(Data.(depVar)(idxDep), 'median');
                transformFcn = @(x) (x - lowerThresh)/(upperThresh - lowerThresh);
            case 'MADmax'
                [~, ~, upperThresh, ~] = isoutlier(Data.(depVar)(idxDep), 'median');
                transformFcn = @(x) x/upperThresh;
            case 'max'
                upperThresh = max(Data.(depVar)(idxDep));
                transformFcn = @(x) x/upperThresh;
            case 'minmax'
                lowerThresh = min(Data.(depVar)(idxDep));
                upperThresh = max(Data.(depVar)(idxDep));
                transformFcn = @(x) x/(upperThresh - lowerThresh);
            otherwise % any other expression
                tokens = regexp(transform, '[q,p](\d+)(?:[q,p]?)(\d+)?', 'tokens', 'once');
                tokens(cellfun(@isempty, tokens)) = [];
                if length(tokens) == 1 % upper percentile given in the form 'q%d' or 'p%d'
                    upperThresh = prctile(Data.(depVar)(idxDep), str2double(tokens{1}));
                    transformFcn = @(x) x/upperThresh;
                elseif length(tokens) == 2 % lower and upper percentile given in the form 'q%dq%d' or 'p%dp%d'
                    lowerThresh = prctile(Data.(depVar)(idxDep), str2double(tokens{1}));
                    upperThresh = prctile(Data.(depVar)(idxDep), str2double(tokens{2}));
                    transformFcn = @(x) x/(upperThresh - lowerThresh);
                else % any other function given in the form 'f(x)'
                    transformFcn = eval(sprintf('@(x) %s', transform));
                end
        end
        Data.(transVar)(idxDep) = transformFcn(Data.(depVar)(idxDep));
    end
else
    transVar = depVar;
end

%% Remove pre-fit outliers

% open outlier file for writing
fidOutliers = fopen(fullfile(outDir, 'Outliers.txt'), 'w+');
nPreOutliers = 0;
nPreObs = size(Data, 1);

outlierLevel = nFactors;
if ~strcmp(preOutlierMethod, 'none')

    idxOut = false(size(Data, 1), 1);

    for iVar = 1:nY

        if nY > 1
            idxDep = (Data.(yVar) == y{iVar});
        else
            idxDep = true(size(Data, 1), 1);
        end

        idxTest = idxDep;

        if outlierLevel == 0 % level 0: all data
            yData = Data.(transVar)(idxTest);
            idxOut(idxTest) = isoutlier(yData, preOutlierMethod);

        elseif outlierLevel == 1 % level 1: 1st dependent variable
            for iMember = 1:nMembers
                idxTest = idxDep;
                member = members(iMember);
                idxTest = idxTest & (Data.(memberVar) == member);
                yData = Data.(transVar)(idxTest);
                idxOut(idxTest) = isoutlier(yData, preOutlierMethod);
            end

        elseif outlierLevel == 2 % level 2: 2nd dependent variable, if given
            for iMember = 1:nMembers
                for iGroup = 1:nGroups
                    idxTest = idxDep;
                    member = members(iMember);
                    idxTest = idxTest & (Data.(memberVar) == member);
                    if nGroups > 1
                        group = groups(iGroup);
                        idxTest = idxTest & (Data.(groupVar) == group);
                    end
                    yData = Data.(transVar)(idxTest);
                    idxOut(idxTest) = isoutlier(yData, preOutlierMethod);
                end
            end

        elseif outlierLevel == 3 % level 3: 3rd dependent variable, if given
            for iMember = 1:nMembers
                for iGroup = 1:nGroups
                    for iCol = 1:nCols
                        idxTest = idxDep;
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
                        idxOut(idxTest) = isoutlier(yData, preOutlierMethod);
                    end
                end
            end
        elseif outlierLevel == 4 % level 4: 4th dependent variable, if given
            for iMember = 1:nMembers
                for iGroup = 1:nGroups
                    for iCol = 1:nCols
                        for iRow = 1:nRows
                            idxTest = idxDep;
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
                            idxOut(idxTest) = isoutlier(yData, preOutlierMethod);
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

% report to file
msg = sprintf('Removed %d pre-fit outlier(s) from %d observations (%.1f %%%%) using ''%s''\n', nPreOutliers, nPreObs, nPreOutliers/nPreObs*100, preOutlierMethod);
fprintf(msg);
fprintf(fidOutliers, msg);

% close outlier file
fclose(fidOutliers);

%% Fit linear model and perform ANOVA

% perform nY fits only if separateMulti=true
nFits = 1;
if strcmp(fitMethod, 'none')
    nFits = 0;
elseif separateMulti
    nFits = nY;
end
mdls = cell(1, nFits);
anovas = cell(1, nFits);
Datasets = cell(1, nFits);
for iFit = 1:nFits % if not separateMulti, this loop is left after the 1st iteration

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

    if nY > 1 && separateMulti % separate univariate analyses of multi-valued dependent variable
        myVar = y{iFit};
        idxDep = (Data.(yVar) == myVar);
        DataOrig = Data;
        Data = Data(idxDep, :);
        outDirOrig = outDir;
        outDir = sprintf('%s/%s', outDir, myVar);
        fprintf('Performing GLMM analysis for %s...\n', myVar);
    elseif nY > 1
        fprintf('Performing multivariate GLMM analysis...\n');
    else
        fprintf('Performing GLMM analysis for %s...\n', depVar);
    end

    % create output folder, if not existing
    if ~isfolder(outDir)
        mkdir(outDir);
    end

    % open summary file for writing
    fidSummary = fopen(fullfile(outDir, 'Summary.txt'), 'w+');

    productTerm = strjoin(interact, '*');    
    if nY > 1 && ~separateMulti
        productTerm = sprintf('-1 + %s:(%s)', yVar, productTerm);
    end
    xNoInteract = setdiff(x, interact, 'stable');

    sumTerm = strjoin(union(xNoInteract, coVar, 'stable'), ' + ');
    
    if ~isempty(id) && length(unique(Data.(id))) > 1
        if randomSlopes
            randomEffect = '';
            if nY > 1 && ~separateMulti
                randomSlopes = strjoin(cellfun(@(x) sprintf('(%s|%s:%s)', x, yVar, id), union(within, coVar, 'stable'), 'UniformOutput', false), ' + ');
            else
                randomSlopes = strjoin(cellfun(@(x) sprintf('(%s|%s)', x, id), union(within, coVar, 'stable'), 'UniformOutput', false), ' + ');
            end
        else
            randomEffect = strjoin(cellfun(@(x) sprintf('(1|%s:%s)', x, id), union(within, coVar, 'stable'), 'UniformOutput', false), ' + ');
            randomSlopes = '';
        end
    else
        randomEffect = '';
        randomSlopes = '';
    end

    terms = {};
    if ~isempty(sumTerm)
        terms = [terms, sumTerm]; %#ok<AGROW>
    end
    if ~isempty(productTerm)
        terms = [terms, productTerm]; %#ok<AGROW>
    end
    if ~isempty(randomEffect)
        terms = [terms, randomEffect]; %#ok<AGROW>
    end
    if ~isempty(randomSlopes)
        terms = [terms, randomSlopes]; %#ok<AGROW>
    end

    if nY > 1 && separateMulti
        formulaOrig = formula;
    end
    if isempty(formula)
        formula = sprintf('%s ~ %s', transVar, strjoin(terms, ' + '));
    end
    fprintf('\t%s\n', formula);

    try
        if ~isempty(myLink) && ~strcmp(myLink, 'auto') % link function given -> use it
            mdl = fitglme(Data, formula, ...
                'DummyVarCoding', 'effects', ...
                'FitMethod', fitMethod, ...
                'Distribution', myDistribution, ...
                'link', myLink);
        else % no link function given or set to 'auto' -> use built-in default
            mdl = fitglme(Data, formula, ...
                'DummyVarCoding', 'effects', ...
                'FitMethod', fitMethod, ...
                'Distribution', myDistribution);
        end

    catch ME
        message = sprintf('%s', ME.message);
        fprintf('The linear model fit returned an error:\n\t%s\n', message);
        fprintf('Please try again, using fewer interactions by defining "interact" with only those independent variables whose interaction you want to investigate\n');
        return
    end

    % remove post-fit outliers and refit model
    if ~strcmp(outlierMethod, 'none')
        mdlResiduals = residuals(mdl, 'ResidualType', 'Pearson');
        mdlOutliers = isoutlier(mdlResiduals, outlierMethod);
        nOutliers = sum(mdlOutliers);
        nObs = length(mdlResiduals);
        msg = sprintf('Removing %d post-fit outliers from %d observations (%.1f %%%%) using ''%s''...\n', nOutliers, nObs, nOutliers/nObs*100, outlierMethod);
        fprintf(msg);
        fprintf(fidSummary, msg);
        if nOutliers > 0
            % remove outliers from Data
            Data(mdlOutliers, :) = [];
            % remove outliers from DataOrig
            if nY > 1 && separateMulti
                idxTmp = false(size(DataOrig, 1), 1);
                idxTmp(idxDep) =  mdlOutliers;
                DataOrig(idxTmp, :) = [];
            end
            msg = sprintf('Re-fitting model...\n');
            fprintf(msg);
            fprintf(fidSummary, msg);
            try
                if ~isempty(myLink) && ~strcmp(myLink, 'auto') % link function given -> use it
                    mdl = fitglme(Data, formula, ...
                        'DummyVarCoding', 'effects', ...
                        'FitMethod', fitMethod, ...
                        'Distribution', myDistribution, ...
                        'link', myLink);
                else % no link function given or set to 'auto' -> use built-in default
                    mdl = fitglme(Data, formula, ...
                        'DummyVarCoding', 'effects', ...
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
    fclose(fidSummary);

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

    if nY > 1 && separateMulti
        if strcmp(depVar, yVal)
            sgtitle(sprintf('Diagnostics for %s', myVar), 'interpreter', 'none', 'FontWeight', 'bold', 'FontSize', 14);
        else
            sgtitle(sprintf('Diagnostics for %s %s', myVar, depVar), 'interpreter', 'none', 'FontWeight', 'bold', 'FontSize', 14);
        end
    elseif nY > 1 && strcmp(plotTitle, yVal)
        sgtitle(sprintf('Diagnostics for multivariate Analysis'), 'interpreter', 'none', 'FontWeight', 'bold', 'FontSize', 14);
    else
        sgtitle(sprintf('Diagnostics for %s', plotTitle), 'interpreter', 'none', 'FontWeight', 'bold', 'FontSize', 14);
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
    rectangle('Position', [resFit.mu - resFit.sigma, ylims(1), 2*resFit.sigma, ylims(2)], 'FaceColor', [0 0 0 0.2], 'EdgeColor', 'none');
    title('Histogram of residuals');
    xlabel('Residuals');

    % normality of residuals
    iPanel = iPanel+1;
    subplot(nPanelRows, nPanelCols, iPanel);
    plotResiduals(mdl, 'probability', 'ResidualType', 'Pearson');
    [isNotNormal, pNormal] = lillietest(mdlResiduals);
    text(gca, 0.05,0.95,sprintf('Normality = %d (p = %.3f)', ~isNotNormal, pNormal), 'Units', 'normalized');

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
    % % plot fitted line
    % hold on
    % h = gca;
    % xData = h(1).XData;
    % yData = h(1).YData;
    % [xData, sortOrder] = sort(xData, 'ascend');
    % yData = yData(sortOrder); % Need to sort b the same way.
    % idxValid = ~isnan(yData);
    % fitCoeffs = polyfit(xData(idxValid), yData(idxValid), 1);
    % xFitted = xData(idxValid); % Evalutate the fit as the same x coordinates.
    % yFitted = polyval(fitCoeffs, xFitted);
    % plot(xFitted, yFitted, 'b-', 'LineWidth', 2);
    % % test for homoscedasticity using Arch test
    % cleanResiduals = mdlResiduals(~isnan(mdlResiduals));
    % [isHetero, pHetero] = archtest(cleanResiduals);
    % text(gca, 0.05,0.95,sprintf('Homoscedasticity = %d (p = %f)', ~isHetero, pHetero), 'Units', 'normalized');

    % lagged residuals
    iPanel = iPanel+1;
    subplot(nPanelRows, nPanelCols,iPanel);
    plotResiduals(mdl, 'lagged', 'ResidualType', 'Pearson');


    % save figure
    set(fig, 'PaperPositionMode', 'auto');
    set(fig, 'PaperUnits', 'points');
    set(fig, 'PaperSize', [figWidth figHeight]);
    print(fig, fullfile(outDir, sprintf('%s.pdf', figName)), '-fillpage', '-dpdf', sprintf('-r%.0f', 300));

    % close figure
    close(fig);

    % store GLM fit
    mdls{iFit} = mdl;

    % get ANOVA table from model fit and sore
    anovaResult = anova(mdl);
    anovas{iFit} = anovaResult;

    Datasets{iFit} = Data;

    % if separate multi, restore Data and formula
    if nY > 1 && separateMulti
        Data = DataOrig; % set Data to original Data minus post-fit outliers
        formula = formulaOrig;
        outDir = outDirOrig;
    end
end

%% Statistically correct ANOVA tables

if ~isempty(anovas)
    % Sidak correction
    for iFit = 1:nFits
        anovas{iFit}.pValue = sidak_corr(anovas{iFit}.pValue, nFits);
        anovas{iFit}.pValue = sidak_corr(anovas{iFit}.pValue, correctForN); % additionally correct for multiple tests like this one
    end
end

%% Create ANOVA table

for iFit = 1:nFits

    if nY > 1 && separateMulti % separate univariate analyses of multi-valued dependent variable
        myVar = y{iFit};
        outDirOrig = outDir;
        outDir = sprintf('%s/%s', outDir, myVar);
    end

    % retrieve GLM fit
    mdl = mdls{iFit};

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

    % save Data
    saveTable(Datasets{iFit}, 'Data', {'csv'}, outDir);

    % save raw Data table
    saveTable(DataRaw, 'DataRaw', {'csv'}, outDir);

    % save ANOVA table
    saveTable(anovaTable, 'Anova', {'xlsx'}, outDir);
    disp(anovaTable) % display table

    if nY > 1 && separateMulti % separate univariate analyses of multi-valued dependent variable
        outDir = outDirOrig;
    end
end

%% Plot data and make posthoc comparisons

outDirOrig = outDir;
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

    % define posthoc comparison pairs
    pairs = nchoosek(members, 2);
    nPairs = size(pairs, 1);

    % create (nRows x nCols x nGroups x nMembers x nY) arrays of plot data
    bar_values      = nan(nRows, nCols, nGroups, nMembers, nY);
    bar_errorTop    = nan(nRows, nCols, nGroups, nMembers, nY);
    bar_errorBottom = nan(nRows, nCols, nGroups, nMembers, nY);
    bar_p           = nan(nGroups, nPairs, nRows, nCols, nY);
    bar_test        = nan(nGroups, nPairs, nRows, nCols, nY);
    bar_DF          = nan(nGroups, nPairs, nRows, nCols, nY);
    bar_aux         = nan(nGroups, nPairs, nRows, nCols, nY);
    bar_eff         = nan(nGroups, nPairs, nRows, nCols, nY);

    % create (nPairs x nY) arrays of main effects data
    main_p      = nan(nPairs, nRows, nCols, nY);
    main_test   = nan(nPairs, nRows, nCols, nY);
    main_DF     = nan(nPairs, nRows, nCols, nY);
    main_aux    = nan(nPairs, nRows, nCols, nY);
    main_eff    = nan(nPairs, nRows, nCols, nY);

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
    violin_values = nan(nGroups, nMembers, maxNValues, nRows, nCols, nY);

    for iVar = 1:nY
        myVar = y{iVar};

        % posthoc method
        if strcmp(posthocMethod, 'auto')
            if nY > 1 && separateMulti && ~strcmp(fitMethod, 'none')
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
        Stats = table;

        % calc plot data and fill arrays
        for iRow = 1:nRows % loop over 4th x
            for iCol = 1:nCols % loop over 3rd x
                for iGroup = 1:nGroups % loop over 2nd x                    
                    for iMember = 1:nMembers % loop over 1st x

                        if nY > 1
                            idxDep = (Data.(yVar) == myVar);
                        else
                            idxDep = true(size(Data, 1), 1);
                        end

                        % init statistics table row
                        statsRow = table;

                        % 1st x
                        member = members(iMember);
                        idx = (Data.(memberVar) == member) & idxDep;

                        % 2nd x, if given
                        if nGroups > 1
                            group = groups(iGroup);
                            idx = idx & (Data.(groupVar) == group);
                        end

                        % 3rd x, if given
                        if nCols > 1
                            col = cols(iCol);
                            statsRow.(colVar) = string(col);
                            idx = idx & Data.(colVar) == col;
                        end

                        % 4th x, if given
                        if nRows > 1
                            row = rows(iRow);
                            statsRow.(rowVar) = string(row);
                            idx = idx & Data.(rowVar) == row;
                        end

                        bar_data = Data(idx, :);
                        values = bar_data.(depVar);
                        if nGroups > 1
                            statsRow.(groupVar) = string(group);
                        end
                        statsRow.(memberVar) = string(member);
                        statsRow.N = length(values(~isnan(values)));
                        statsRow.median = median(values, 'omitnan');
                        statsRow.q25 = quantile(values, 0.25);
                        statsRow.q75 = quantile(values, 0.75);
                        statsRow.mean = mean(values, 'omitnan');
                        statsRow.std = std(values, 'omitnan');
                        statsRow.SE = statsRow.std / sqrt(statsRow.N);
                        statsRow.RE = statsRow.SE / statsRow.mean * 100;
                        tci95 = tinv([0.025 0.975], statsRow.N-1); % 95% of t-Distribution
                        ci95 = statsRow.mean + tci95 * statsRow.SE;
                        statsRow.ci95_1 = ci95(1);
                        statsRow.ci95_2 = ci95(2);
                        Stats = [Stats; statsRow]; %#ok<AGROW>

                        violin_values(iGroup, iMember, 1:length(values), iRow, iCol, iVar) = values; %%%%%%%

                        % plot values
                        switch errorBars
                            case 'std'
                                bar_values(iRow, iCol, iGroup, iMember, iVar) = statsRow.mean;
                                bar_errorBottom(iRow, iCol, iGroup, iMember, iVar) = statsRow.mean - statsRow.std;
                                bar_errorTop(iRow, iCol, iGroup, iMember, iVar) = statsRow.mean + statsRow.std;
                            case 'se'
                                bar_values(iRow, iCol, iGroup, iMember, iVar) = statsRow.mean;
                                bar_errorBottom(iRow, iCol, iGroup, iMember, iVar) = statsRow.mean - statsRow.SE;
                                bar_errorTop(iRow, iCol, iGroup, iMember, iVar) = statsRow.mean + statsRow.SE;
                            case 'q25'
                                bar_values(iRow, iCol, iGroup, iMember, iVar) = statsRow.median;
                                bar_errorBottom(iRow, iCol, iGroup, iMember, iVar) = statsRow.q25;
                                bar_errorTop(iRow, iCol, iGroup, iMember, iVar) = statsRow.q75;
                            case 'ci95'
                                bar_values(iRow, iCol, iGroup, iMember, iVar) = statsRow.mean;
                                bar_errorBottom(iRow, iCol, iGroup, iMember, iVar) = statsRow.q25;
                                bar_errorTop(iRow, iCol, iGroup, iMember, iVar) = statsRow.q75;
                        end


                    end

                    % get post-hoc p-values
                    for iPair = 1:nPairs
                        pair = pairs(iPair, :);

                        switch posthocMethod

                            case {'ttest', 'utest'} % perform post-hoc analysis using paired t-tests

                                if nY > 1
                                    idxDep = (Data.(yVar) == myVar);
                                else
                                    idxDep = true(size(Data, 1), 1);
                                end

                                % init idx
                                idx = idxDep;

                                % 2nd x, if given
                                if nGroups > 1
                                    group = groups(iGroup);
                                    idx = idx & (Data.(groupVar) == group);
                                end

                                % 3rd x, if given
                                if nCols > 1
                                    col = cols(iCol);
                                    idx = idx & (Data.(colVar) == col);
                                end

                                % 4th x, if given
                                if nRows > 1
                                    row = rows(iRow);
                                    idx = idx & (Data.(rowVar) == row);
                                end

                                % calc contrasts
                                L1 = (Data.(memberVar) == pair(1) & idx);
                                L2 = (Data.(memberVar) == pair(2) & idx);
                                val1 = Data.(transVar)(L1);
                                val2 = Data.(transVar)(L2);
                                switch posthocMethod
                                    case 'ttest'
                                        [~, bar_p(iGroup, iPair, iRow, iCol, iVar), ~, stats] = ttest2(val1, val2);
                                        tValue = stats.tstat;
                                        df = stats.df;
                                        sPool = stats.sd;
                                        dCohen = (mean(val1, 'omitnan') - mean(val2, 'omitnan')) / sPool;
                                        bar_test(iGroup, iPair, iRow, iCol, iVar) = tValue;
                                        bar_DF(iGroup, iPair, iRow, iCol, iVar) = df;
                                        bar_eff(iGroup, iPair, iRow, iCol, iVar) =  dCohen;
                                    case 'utest'
                                        [bar_p(iGroup, iPair, iRow, iCol, iVar), ~, stats] = ranksum(val1, val2);
                                        N1 = sum(~isnan(val1));
                                        N2 = sum(~isnan(val2));
                                        W = stats.ranksum;
                                        U = W - N1*(N1+1)/2;
                                        r = 1 - 2*U/(N1*N2);
                                        bar_test(iGroup, iPair, iRow, iCol, iVar) = U;
                                        bar_eff(iGroup, iPair, iRow, iCol, iVar) =  r;
                                end

                                % calc main contrasts
                                if posthocMainEffects
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
                                            dCohen = (mean(val1, 'omitnan') - mean(val2, 'omitnan')) / sPool;
                                            main_test(iPair, iRow, iCol, iVar) = tValue;
                                            main_DF(iPair, iRow, iCol, iVar) = df;
                                            main_eff(iPair, iRow, iCol, iVar) =  dCohen;
                                        case 'utest'
                                            [main_p(iPair, iRow, iCol, iVar), ~, stats] = ranksum(val1, val2);
                                            N1 = sum(~isnan(val1));
                                            N2 = sum(~isnan(val2));
                                            W = stats.ranksum;
                                            U = W - N1*(N1+1)/2;
                                            r = 1 - 2*U/(N1*N2);
                                            main_test(iPair, iRow, iCol, iVar) = U;
                                            main_eff(iPair, iRow, iCol, iVar) =  r;
                                    end
                                end

                            case 'emm' % perform post-hoc analysis using emmeans

                                if nY > 1 && separateMulti
                                    mdl = mdls{iVar};
                                end

                                % Calc estimated marginal means of all
                                % factors. Continuous variables cannot be
                                % included, because then emmeans gives an
                                % error
                                emm = emmeans(mdl, reshape(catVars, 1, []), 'effects', 'unbalanced');

                                if nY > 1 && ~separateMulti
                                    idxDep = (emm.table.(yVar) == myVar);
                                else
                                    idxDep = true(size(emm.table, 1), 1);
                                end

                                % init idx
                                idx = idxDep;

                                % 2nd x, if given
                                if nGroups > 1
                                    group = groups(iGroup);
                                    idx = idx & (emm.table.(groupVar) == group);
                                end

                                % 3rd x, if given
                                if nCols > 1
                                    col = cols(iCol);
                                    idx = idx & (emm.table.(colVar) == col);
                                end

                                % 4th x, if given
                                if nRows > 1
                                    row = rows(iRow);
                                    idx = idx & (emm.table.(rowVar) == row);
                                end

                                % calc contrasts
                                L1 = idx & (emm.table.(memberVar) == pair(1));
                                L2 = idx & (emm.table.(memberVar) == pair(2));
                                L = (L1 - L2)';
                                contrasts = kbcontrasts_wald(mdl, emm, L);
                                bar_p(iGroup, iPair, iRow, iCol, iVar) = contrasts.pVal;
                                bar_test(iGroup, iPair, iRow, iCol, iVar) = contrasts.F;
                                bar_DF(iGroup, iPair, iRow, iCol, iVar) = contrasts.DF1;
                                bar_aux(iGroup, iPair, iRow, iCol, iVar) = contrasts.DF2;
                                bar_eff(iGroup, iPair, iRow, iCol, iVar) = f2etaSqp(contrasts.F, contrasts.DF1, contrasts.DF2);

                                % calc main contrasts
                                if posthocMainEffects
                                    L1 = idxDep & (emm.table.(memberVar) == pair(1));
                                    L2 = idxDep & (emm.table.(memberVar) == pair(2));
                                    L = (L1 - L2)';
                                    contrasts = kbcontrasts_wald(mdl, emm, L);
                                    main_p(iPair, iRow, iCol, iVar) = contrasts.pVal;
                                    main_test(iPair, iRow, iCol, iVar) = contrasts.F;
                                    main_DF(iPair, iRow, iCol, iVar) = contrasts.DF1;
                                    main_aux(iPair, iRow, iCol, iVar) = contrasts.DF2;
                                    main_eff(iPair, iRow, iCol, iVar) = f2etaSqp(contrasts.F, contrasts.DF1, contrasts.DF2);
                                end
                        end
                    end
                end
            end
        end

        %% Save descriptive statistics

        if nY > 1
            outDir = sprintf('%s/%s', outDirOrig, myVar);
            if ~isfolder(outDir)
                mkdir(outDir);
            end
        end

        % dispaly and save statistics (only necessary for 1st posthoc level)
        if iLevel == 1
            % display statistics
            disp(Stats);
            fileName = 'Statistics';
            saveTable(Stats, fileName, {'xlsx'}, outDir);
        end
        outDir = outDirOrig;
    end

    %% Statistical correction of posthoc comparisons
    % Holm-Bonferroni correction of all p-values

    bar_pCorr = bar_p; % create array of corrected p-values
    idx = ~isnan(bar_p(:)); % identify NaN-entries
    if ~isempty(idx)
        sizeOrig = size(bar_p); % store original array shape
        bar_p = bar_p(:); % make column vector
        bar_pCorr = bar_pCorr(:); % make column vector
        [~, bar_pCorr(idx)] = bonferroni_holm(bar_p(idx)); % correct p-values, omitting NaNs
        bar_pCorr(idx) = sidak_corr(bar_pCorr(idx), nPosthocLevels); % additionally correct for multiple sets of posthoc comparisons
        bar_pCorr(idx) = sidak_corr(bar_pCorr(idx), correctForN); % additionally correct for multiple tests like this one
        bar_p = reshape(bar_p, sizeOrig); % restore original dimensions of p-value array
        bar_pCorr = reshape(bar_pCorr, sizeOrig); % bring corrected p-value array into the same shape as p-value array
    end
    % statistical correction of main posthoc p-Values
    if posthocMainEffects
        main_pCorr = main_p;
        idx = ~isnan(main_p(:)); % identify NaN-entries
        if ~isempty(idx)
            sizeOrig = size(main_p); % store original array shape
            main_p = main_p(:); % make column vector
            main_pCorr = main_pCorr(:); % make column vector
            [~, main_pCorr(idx)] = bonferroni_holm(main_p(idx)); % correct p-values, omitting NaNs
            main_pCorr(idx) = sidak_corr(main_pCorr(idx), nPosthocLevels); % additionally correct for multiple sets of posthoc comparisons
            main_pCorr(idx) = sidak_corr(main_pCorr(idx), correctForN); % additionally correct for multiple tests like this one
            main_p = reshape(main_p, sizeOrig); % restore original dimensions of p-value array
            main_pCorr = reshape(main_pCorr, sizeOrig); % bring corrected p-value array into the same shape as p-value array
        end
    end

    %% Loop over dependent variables

    outDirOrig = outDir;
    for iVar = 1:nY

        % dependent variable
        myVar = y{iVar};

        % yLabel
        if length(yLabel) == nY
            myLabel = yLabel{iVar};
        elseif length(yLabel) == 1
            myLabel = yLabel{1};
        end

        % yUnits
        if length(yUnits) == nY
            myUnits = yUnits{iVar};
        elseif length(yUnits) == 1
            myUnits = yUnits{1};
        end


        % set output folder
        if nY > 1 && separateMulti
            outDir = sprintf('%s/%s', outDirOrig, myVar);
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
            if nY > 1 && separateMulti
                if strcmp(depVar, yVal)
                    title(layout, sprintf('Data plots for %s', myVar), 'interpreter', 'none', 'FontWeight', 'bold', 'FontSize', 14);
                else
                    title(layout, sprintf('Data plots for %s %s', myVar, depVar), 'interpreter', 'none', 'FontWeight', 'bold', 'FontSize', 14);
                end
            elseif nY > 1 && strcmp(plotTitle, yVal)
                title(layout, sprintf('Data plots for %s', myVar), 'interpreter', 'none', 'FontWeight', 'bold', 'FontSize', 14);
            else
                title(layout, sprintf('Data plots for %s', plotTitle), 'interpreter', 'none', 'FontWeight', 'bold', 'FontSize', 14);
            end

            % prepare to display variable names and levels
            switch showVarNames
                case {2, 3} % capitalize variable names and levels
                    displayMemberVar = string(capitalize(memberVar));
                    displayMembers = string(strsplit(capitalize(strjoin(cellstr(members), ', ')), ', '));
                    displayGroupVar = string(capitalize(groupVar));
                    displayGroups = string(strsplit(capitalize(strjoin(cellstr(groups), ', ')), ', '));
                    displayColVar = string(capitalize(colVar));
                    displayCols = string(strsplit(capitalize(strjoin(cellstr(cols), ', ')), ', '));
                    displayRowVar = string(capitalize(rowVar));
                    displayRows = string(strsplit(capitalize(strjoin(cellstr(rows), ', ')), ', '));
                otherwise % use original variable names and levels
                    displayMemberVar = memberVar;
                    displayMembers = members;
                    displayGroupVar = groupVar;
                    displayGroups = groups;
                    displayColVar = colVar;
                    displayCols = cols;
                    displayRowVar = rowVar;
                    displayRows = rows;
            end
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

                    switch showVarNames
                        case {1, 2} % display variable names and levels
                            if length(cols) > 1 && length(rows) > 1
                                panelTitle = sprintf('%s = %s, %s = %s', displayRowVar, displayRows(iRow), displayColVar, displayCols(iCol));
                            elseif length(rows) > 1
                                panelTitle = sprintf('%s = %s', displayRowVar, displayRows(iRow));
                            elseif length(cols) > 1
                                panelTitle = sprintf('%s = %s', displayColVar, displayCols(iCol));
                            else
                                panelTitle = '';
                            end
                        otherwise % only display variable levels
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
                    plotGroups(violin_values(:, :, :, iRow, iCol, iVar), displayMembers, displayGroups, displayMemberVar, displayGroupVar, bar_pCorr(:, :, iRow, iCol, iVar), panelTitle, yLabelStr, plotStyle, panel, showVarNames, markerSize);
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
            if nY > 1
                outDir = sprintf('%s/%s', outDirOrig, myVar);
                if ~isfolder(outDir)
                    mkdir(outDir);
                end
            end
            set(fig, 'PaperPositionMode', 'auto');
            set(fig, 'PaperUnits', 'points');
            set(fig, 'PaperSize', [figWidth figHeight]);
            print(fig, fullfile(outDir, sprintf('%s.pdf', figName)), '-fillpage', '-dpdf', sprintf('-r%.0f', 300));
            print(fig, fullfile(outDir, sprintf('%s.png', figName)), '-dpng', sprintf('-r%.0f', 300));
            saveas(fig, fullfile(outDir, sprintf('%s.fig', figName)));
            outDir = outDirOrig;

            % close figure
            close(fig);
        end

        %% Post-hoc table

        posthocTable = table;

        % posthoc main effects
        if posthocMainEffects
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
                tableRow.p = main_p(iPair, iRow, iCol, iVar);
                tableRow.pCorr = main_pCorr(iPair, iRow, iCol, iVar);
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
                        tableRow.p = bar_p(iGroup, iPair, iRow, iCol, iVar);
                        tableRow.pCorr = bar_pCorr(iGroup, iPair, iRow, iCol, iVar);
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

        % save table
        if nY > 1
            outDir = sprintf('%s/%s', outDirOrig, myVar);
            if ~isfolder(outDir)
                mkdir(outDir);
            end
        end
        if nPosthocLevels > 1
            fileName = sprintf('Posthoc_%d', iLevel);
        else
            fileName = 'Posthoc';
        end
        saveTable(posthocTable, fileName, {'xlsx'}, outDir);
        outDir = outDirOrig;
    end

end

end

%% Helper functions %%%%%%%%%%%%%%%%%%%%%%

function value = getValue(in)
% get numerical or logical value from string, if possible
if ischar(in) || isstring(in)
    [numValue, isNum] = str2num(in);
    if isNum
        value = numValue;
    else
        value = in;
    end
else
    value = in;
end
end