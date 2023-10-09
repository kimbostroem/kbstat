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
%       y               Name of the dependent variable or list of names of 
%                       dependent variables. If there is more than one
%                       component of y, by default each component is
%                       treated as a separate dependent variable and as
%                       many independent univariate analyses are performed.
%                       If "separateMulti" is set to "false", then y is
%                       treated as a vector-valued dependent variable and a
%                       multivariate analysis is performed
%
%       yUnits          Physical units of the dependent variable.
%                       OPTIONAL, default = ''.
%
%       x               List of the names of the independent variables,
%                       separated by comma or semicolon. Up to 4
%                       independent variables are supported. Variables that
%                       are integer or non-numeric, are interpreted as
%                       categorical, i.e. as factors. Only factors
%                       are included in the barplot chart.
%
%       id              Name of the subject variable.
%                       OPTIONAL, default = ''.
%
%       within          List of within-subject variables, separated by
%                       comma or semicolon.  Within-subject variables are
%                       nested within the subject variable, i.e. they vary
%                       for each subject. Variables that are not declared
%                       as within-subject, are considered between-subject,
%                       i.e. they vary only between, not within, subjects.
%                       "within" can be a subset of "x", or else its
%                       members are added to "x". Example:
%                       options.id = 'subject'
%                       options.x = 'time, age, sex'
%                       options.within = 'dose'.
%
%       interact        List of variables whose interaction with each other
%                       is to be analyzed. Can be a subset of "x", or else
%                       its members are added to "x". Example:
%                       options.id = 'subject'
%                       options.x = 'time, dose'
%                       options.interact = 'dose, age'.
%
%       randomSlopes    Flag if random slopes should be estimated.
%                       OPTIONAL, default = true.
%
%       formula         Formula in Wilkinson Notation. If given, it
%                       overrides the automatically generated formula
%                       produced from the provided variables.
%
%       fitMethod       Fit method used for the GLM fit.
%                       OPTIONAL, default = 'MPL'.
%                       Possible values:
%                       'MPL'                   Maximum pseudo likelihood
%                       'REMPL'                 Restricted maximum pseudo likelihood
%                       'Laplace'               Maximum likelihood using Laplace approximation
%                       'ApproximateLaplace'    Maximum likelihood using approximate Laplace approximation with fixed effects profiled out
%                       'none'                  Skip linear model fit
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
%                       'identity'	    g(mu) = mu.             Default for Normal distribution
%                       'log'	        g(mu) = log(mu).        Default for Poisson
%                       'logit'	        g(mu) = log(mu/(1-mu))  Default for Binomial distribution
%                       'loglog'	    g(mu) = log(-log(mu))
%                       'probit'	    g(mu) = norminv(mu)
%                       'comploglog'	g(mu) = log(-log(1-mu))
%                       'reciprocal'	g(mu) = mu.^(-1)        Default for Gamma
%                       Scalar p	    g(mu) = mu.^p           Default for InverseGaussian (p= -2)
%
%       posthocMethod   Method for the posthoc pairwise comparison.
%                       Possible values:
%                       'ttest'     t-test with Holm-Bonferroni correction
%                       'utest'     Mann-Whitney u-Test (ranksum test) with Holm-Bonferroni correction
%                       'emm'       Extract contrasts from linear model fit
%                       'none'      Do not perform posthoc analysis
%                       OPTIONAL, default is
%                           'emm'   if univariate or multi-valued y with separateMulti = true
%                           'ttest' if distribution = 'normal'
%                           'utest' otherwise
%
%       postHocMain     Flag if also the posthoc main effects should be 
%                       calculated, i.e. the comparison between one 
%                       variable set to 'any'. 
%                       OPTIONAL, default = true.
%
%       separateMulti   Flag if a multivariate dependent variable should be
%                       analyzed for each component separately.
%                       OPTIONAL, default = true.
%
%       multiVar        Name of the variable that encodes levels of a 
%                       multivariate dependent variable.
%                       OPTIONAL, default = ''.
%
%       isRescale       Flag if the y-axis of each
%                       panel of the data plot is to be resized to a
%                       common scale.
%                       OPTIONAL, default = true.
%
%       removeOutliers  Indicate which outliers should be removed from the 
%                       data. Possible values:
%                       true        Remove pre-fit and post-fit outliers
%                       false,  
%                       'none'      Do not remove outliers
%                       'pre'       Remove pre-fit outliers
%                       'post'      Remove post-fit outliers
%                       'prepost'   Remove pre-fit and post-fit outliers.
%                       OPTIONAL, default = 'prepost'.
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
%                              'Z'  Z-transform to zero-centered
%                                   distribution with unit standard 
%                                   derivation.
%                        'IQR' or   Rescale to interval [0,1] using upper 
%                      'quartiles'  and lower limits of the inter-quartile 
%                                   range (IQR).
%                           'f(x)'  Arbitrary function of x, such as
%                                       'log(x)'
%                                       'atanh(x)'
%                                       '(x-mean(x, "omitnan"))/std(x, "omitnan")'
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
%       yLabel          Label for y axis in data plots.
%                       OPTIONAL, default = y
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
%   options.interact = 'shoe, speed, joint';
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
for iFit = 1:length(tableVars)
    tableVar = tableVars{iFit};
    if iscell(Data1.(tableVar))
        Data1.(tableVar) = string(Data1.(tableVar));
    end
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
if isfield(options, 'yUnits')
    yUnits = cellstr(options.yUnits);
else
    yUnits = {''};
end

% multiVar
if isfield(options, 'multiVar') && ~isempty(options.multiVar)
    multiVar = options.multiVar;
else
    multiVar = '';
end

if nY > 1
    yVar = 'yVar';
    yVal = 'Y';
    Data2 = stack(Data1, y, 'NewDataVariableName', yVal, 'IndexVariableName', yVar);
    Data2.(yVar) = categorical(string(Data2.(yVar)));
    depVar = yVal; % set dependent variable to y
    depVarUnits = '';
elseif ~isempty(multiVar)
    yVar = multiVar;
    yVal = y{1};
    Data2 = Data1;
    depVar = y{1}; % set dependent variable to y    
    Data2.(yVar) = categorical(string(Data2.(yVar)));
    y = cellstr(unique(Data2.(yVar)));
    nY = length(y);
    depVarUnits = yUnits{1};
else
    Data2 = Data1;
    depVar = y{1};
    yVal = y{1};
    depVarUnits = yUnits{1};
end

%% more variables

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
if isfield(options, 'within') && ~isempty(options.within)
    within = strtrim(strsplit(options.within, {',', ';'}));
    x = union(x, within, 'stable'); % add within-subject variables to independent variables
else
    within = {};
end

% interaction variables
if isfield(options, 'interact') && ~isempty(options.interact)
    interact = strtrim(strsplit(options.interact, {',', ';'}));
    x = union(x, interact, 'stable'); % add interaction variables to independent variables
else
    interact = {};
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
    fitMethod = 'MPL';
end

% distribution (for GLM)
if isfield(options, 'distribution') && ~isempty(options.distribution)
    distribution = options.distribution;
else % no parameter given
    distribution = 'normal';
end

% link (for GLM)
if isfield(options, 'link') && ~isempty(options.link)
    link = options.link;
else % no parameter given
    link = '';
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
if isfield(options, 'isRescale') && ~isempty(options.isRescale)
    isRescale = getValue(options.isRescale);
else
    isRescale = true;
end

% flag if post-fit outliers should be removed
if isfield(options, 'outlierMethod') && ~isempty(options.outlierMethod)
    outlierMethod = options.outlierMethod;
else
    outlierMethod = 'quartiles';
end

% flag if pre-fit outliers should be removed
if isfield(options, 'removeOutliers') && ~isempty(options.removeOutliers)
    removeOutliers = getValue(options.removeOutliers);
else
    removeOutliers = 'prepost';
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
    yLabel = options.yLabel;
else
    yLabel = depVar;
end

if isfield(options, 'xOrder') && ~isempty(options.xOrder)
    xOrder = getValue(options.xOrder);
else
    xOrder = [];
end

fprintf('Performing Linear Model analysis for %s...\n', depVar);

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
    nLevels = length(unique(Data2.(constraintVar)));
    switch nLevels
        case 0
            fprintf('Variable "%s" has no levels left on constraint "%s" -> remove variable\n', constraintVar, constraintVal);
            x = setdiff(x, constraintVar, 'stable');
        case 1
            fprintf('Variable "%s" has only 1 level left on constraint "%s" -> remove variable\n', constraintVar, constraintVal);
            x = setdiff(x, constraintVar, 'stable');
    end
    for iFit = 1:length(auxVars)
        Data2 = removevars(Data2, auxVars{iFit});
    end
end

%% Make variables that are non-numeric or integer, categorical

% init Data table
Data = table;

% init categories
factors = {};

% get subject variable
if ~isempty(id)
    myIV = id;
    myLevels = unique(DataRaw.(id));
    if ~all(isnumeric(myLevels)) || (all(isnumeric(myLevels)) && all(mod(myLevels,1) == 0)) % levels all integers or strings -> make categorical
        Data.(myIV) = categorical(string(Data2.(myIV))); % make categorical
    else
        Data.(myIV) = Data2.(myIV); % keep continuous values
    end
end

% get independent variables
for iIV = 1:length(x)
    myIV = x{iIV};
    myLevels = unique(Data2.(myIV));
    if all(isnumeric(myLevels)) && all(mod(myLevels,1) == 0) % levels all integer -> make categorical
        Data.(myIV) = categorical(string(Data2.(myIV)));
        factors = [factors; myIV]; %#ok<AGROW>
    elseif all(isnumeric(myLevels)) % levels all numerical (but not integer) -> leave as is
        Data.(myIV) = Data2.(myIV); % keep continuous values
    else
        Data.(myIV) = categorical(Data2.(myIV)); % else -> make categorical
        factors = [factors; myIV]; %#ok<AGROW>
    end
end

% multivariate variable
if nY > 1
    Data.(yVar) = Data2.(yVar);
end

% posthoc method
if isfield(options, 'posthocMethod') && ~isempty(options.posthocMethod)
    posthocMethod = options.posthocMethod;
elseif nY > 1 && separateMulti && ~strcmp(fitMethod, 'none')
    posthocMethod = 'emm';
elseif nY == 1 && ~strcmp(fitMethod, 'none')
    posthocMethod = 'emm';
elseif strcmp(distribution, 'normal')
    posthocMethod = 'ttest';
else
    posthocMethod = 'utest';
end

% posthoc main effects
if isfield(options, 'posthocMain') && ~isempty(options.posthocMain)
    posthocMain = getValue(options.posthocMain);
else
    posthocMain = true;
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
if length(factors) > 1
    groupVar = factors{2};
    groups = unique(Data.(groupVar), levelOrder);
    nGroups = length(groups);
else
    groupVar = '';
    groups = string(factors{1});
    nGroups = 1;
end

% 3rd factor = column variable
if length(factors) > 2
    colVar = factors{3};
    cols = unique(Data.(colVar), levelOrder);
    nCols = length(cols);
else
    colVar = 'none';
    cols = "none";
    nCols = 1;
end

% 4th factors = row variable
if length(factors) > 3
    rowVar = factors{4};
    rows = unique(Data.(rowVar), levelOrder);
    nRows = length(rows);
else
    rowVar = 'none';
    rows = "none";
    nRows = 1;
end

%% Apply data transformation, if given

if ~isempty(transform)
    trnsVar = sprintf('%sTrans', depVar);
    for iFit = 1:nY
        if nY > 1
            idxDep = (Data.(yVar) == y{iFit});
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
            case {'q50', 'median'}
                upperThresh = quantile(Data.(depVar)(idxDep), 0.5);
                transformFcn = @(x) x/upperThresh;
            case 'q75'
                upperThresh = quantile(Data.(depVar)(idxDep), 0.75);
                transformFcn = @(x) x/upperThresh;
            case 'q95'
                upperThresh = quantile(Data.(depVar)(idxDep), 0.95);
                transformFcn = @(x) x/upperThresh;            
            case 'q2575'
                lowerThresh = quantile(Data.(depVar)(idxDep), 0.25);
                upperThresh = quantile(Data.(depVar)(idxDep), 0.75);
                transformFcn = @(x) x/(upperThresh - lowerThresh);            
            case 'q0595'
                lowerThresh = quantile(Data.(depVar)(idxDep), 0.05);
                upperThresh = quantile(Data.(depVar)(idxDep), 0.95);
                transformFcn = @(x) x/(upperThresh - lowerThresh);
            case 'IQRmax'
                [~, ~, upperThresh, ~] = isoutlier(Data.(depVar)(idxDep), 'quartiles');                
                transformFcn = @(x) x/upperThresh;
            case 'IQR'
                [~, lowerThresh, upperThresh, ~] = isoutlier(Data.(depVar)(idxDep), 'quartiles');                
                transformFcn = @(x) x/(upperThresh - lowerThresh);
            case 'IQRzero'
                [~, lowerThresh, upperThresh, ~] = isoutlier(Data.(depVar)(idxDep), 'quartiles');                
                transformFcn = @(x) (x - lowerThresh)/(upperThresh - lowerThresh);
            case 'MAD'
                [~, lowerThresh, upperThresh, ~] = isoutlier(Data.(depVar)(idxDep), 'median');                
                transformFcn = @(x) (x - lowerThresh)/(upperThresh - lowerThresh);
            case 'MADmax'
                [~, ~, upperThresh, ~] = isoutlier(Data.(depVar)(idxDep), 'median');                
                transformFcn = @(x) x/upperThresh;
            case '3sigma'
                [~, lowerThresh, upperThresh, ~] = isoutlier(Data.(depVar)(idxDep), 'mean');                
                transformFcn = @(x) (x - lowerThresh)/(upperThresh - lowerThresh);
            case '3sigmamax'
                [~, ~, upperThresh, ~] = isoutlier(Data.(depVar)(idxDep), 'mean');                
                transformFcn = @(x) x/upperThresh;
            case 'max'
                upperThresh = max(Data.(depVar)(idxDep));
                transformFcn = @(x) x/upperThresh;
            case 'minmax'
                lowerThresh = min(Data.(depVar)(idxDep));
                upperThresh = max(Data.(depVar)(idxDep));
                transformFcn = @(x) x/(upperThresh - lowerThresh);
            case 'minmaxzero'
                lowerThresh = min(Data.(depVar)(idxDep));
                upperThresh = max(Data.(depVar)(idxDep));
                transformFcn = @(x) (x - lowerThresh)/(upperThresh - lowerThresh);
            
            case 'qminmaxzero'
                lowerThresh = quantile(Data.(depVar)(idxDep), 0.05);
                upperThresh = quantile(Data.(depVar)(idxDep), 0.95);
                transformFcn = @(x) (x - lowerThresh)/(upperThresh - lowerThresh);
            otherwise % any other function given in the form 'f(x)'
                transformFcn = eval(sprintf('@(x) %s', transform));                
        end
        Data.(trnsVar)(idxDep) = transformFcn(Data.(depVar)(idxDep));
    end
else
    trnsVar = depVar;
end

%% Remove pre-fit outliers

outlierLevel = length(factors);
if (any(strcmp(removeOutliers, {'pre', 'prepost'})) || removeOutliers) && ~strcmp(outlierMethod, 'none')

    idxOut = false(size(Data, 1), 1);

    for iFit = 1:nY

        if nY > 1
            idxDep = (Data.(yVar) == y{iFit});
        else
            idxDep = true(size(Data, 1), 1);
        end

        idxTest = idxDep;

        if outlierLevel == 0 % level 0: all data
            yData = Data.(trnsVar)(idxTest);
            idxOut(idxTest) = isoutlier(yData, outlierMethod);

        elseif outlierLevel == 1 % level 1: 1st dependent variable
            for iMember = 1:nMembers
                idxTest = idxDep;
                member = members(iMember);
                idxTest = idxTest & (Data.(memberVar) == member);
                yData = Data.(trnsVar)(idxTest);
                idxOut(idxTest) = isoutlier(yData, outlierMethod);
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
                    yData = Data.(trnsVar)(idxTest);
                    idxOut(idxTest) = isoutlier(yData, outlierMethod);
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
                        yData = Data.(trnsVar)(idxTest);
                        idxOut(idxTest) = isoutlier(yData, outlierMethod);
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
                            yData = Data.(trnsVar)(idxTest);
                            idxOut(idxTest) = isoutlier(yData, outlierMethod);
                        end
                    end
                end
            end
        end
    end

    % remove outliers
    nPreOutliers = sum(idxOut);
    nPreObs = size(Data, 1);
    fprintf('Removed %d pre-fit outlier(s) from %d observations (%.1f %%)\n', nPreOutliers, nPreObs, nPreOutliers/nPreObs*100);
    if nPreOutliers > 0
        Data = Data(~idxOut, :);        
    end
end

%% Fit linear model and perform ANOVA

mdls = cell(1, nY);
anovas = cell(1, nY);

% perform nY fits only if separateMulti=true
nFits = 1;
if strcmp(fitMethod, 'none')
    nFits = 0;
    mdls = {};
    anovas = {};
elseif separateMulti
    nFits = nY;
end

for iFit = 1:nFits % if not separateMulti, this loop is left after the 1st iteration

    if nY > 1 && separateMulti % separate univariate analyses of multi-valued dependent variable
        myVar = y{iFit};
        idxDep = (Data.(yVar) == myVar);
        DataOrig = Data;
        Data = Data(idxDep, :);
        outDirOrig = outDir;
        outDir = sprintf('%s/%s', outDir, myVar);
    end

    % create output folder, if not existing
    if ~isfolder(outDir)
        mkdir(outDir);
    end

    productTerm = strjoin(interact, '*');
    if nY > 1 && ~separateMulti
        productTerm = sprintf('-1 + %s:(%s)', yVar, productTerm);
    end
    xNoInteract = setdiff(x, interact, 'stable');

    sumTerm = strjoin(xNoInteract, ' + ');

    if ~isempty(id) && length(unique(Data.(id))) > 1
        if randomSlopes
            randomEffect = '';
            if nY > 1 && ~separateMulti
                randomSlopes = strjoin(cellfun(@(x) sprintf('(%s|%s:%s)', x, yVar, id), within, 'UniformOutput', false), ' + ');
            else
                randomSlopes = strjoin(cellfun(@(x) sprintf('(%s|%s)', x, id), within, 'UniformOutput', false), ' + ');
            end
        else
            randomEffect = strjoin(cellfun(@(x) sprintf('(1|%s:%s)', x, id), within, 'UniformOutput', false), ' + ');
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
        formula = sprintf('%s ~ %s', trnsVar, strjoin(terms, ' + '));
    end
    fprintf('\t%s\n', formula);
   
    try
        if ~isempty(link) % link function given -> use it
            mdl = fitglme(Data, formula, ...
                'DummyVarCoding', 'effects', ...
                'FitMethod', fitMethod, ...
                'Distribution', distribution, ...
                'link', link);
        else % no link function given -> use built-in default
            mdl = fitglme(Data, formula, ...
                'DummyVarCoding', 'effects', ...
                'FitMethod', fitMethod, ...
                'Distribution', distribution);
        end

    catch ME
        message = sprintf('%s', ME.message);
        fprintf('The linear model fit returned an error:\n\t%s\n', message);
        fprintf('Please try again, using fewer interactions by defining "interact" with only those independent variables whose interaction you want to investigate\n');
        return
    end    

    % remove post-fit outliers and refit model
    if outlierMethod
        mdlResiduals = residuals(mdl, 'ResidualType', 'Pearson');
        mdlOutliers = isoutlier(mdlResiduals, 'quartiles');
        nOutliers = sum(mdlOutliers);
        nObs = length(mdlResiduals);
        if nOutliers > 0
            % remove outliers from Data
            Data(mdlOutliers, :) = [];
            % remove outliers from DataOrig
            if nY > 1 && separateMulti
                idxTmp = false(size(DataOrig, 1), 1);
                idxTmp(idxDep) =  mdlOutliers;
                DataOrig(idxTmp, :) = [];
            end
            fprintf('Removed %d post-fit outliers from %d observations (%.1f %%) and refit model...\n', nOutliers, nObs, nOutliers/nObs*100);
            try
                if ~isempty(link) % link function given -> use it
                    mdl = fitglme(Data, formula, ...
                        'DummyVarCoding', 'effects', ...
                        'FitMethod', fitMethod, ...
                        'Distribution', distribution, ...
                        'link', link);
                else % no link function given -> use built-in default
                    mdl = fitglme(Data, formula, ...
                        'DummyVarCoding', 'effects', ...
                        'FitMethod', fitMethod, ...
                        'Distribution', distribution);
                end
            catch ME
                message = sprintf('%s', ME.message);
                fprintf('The linear model fit returned an error:\n\t%s\n', message);
                fprintf('Please try again, using fewer interactions by defining "interact" with only those independent variables whose interaction you want to investigate\n');
                return
            end
        end
    end

    % print results of model fit into file
    mdlOutput = formattedDisplayText(mdl, 'SuppressMarkup', true);
    fid = fopen(fullfile(outDir, 'Summary.txt'), 'w+');
    fprintf(fid, 'Formula:\n\t%s\n', formula);
    if (any(strcmp(removeOutliers, {'pre', 'prepost'})) || removeOutliers) && ~strcmp(outlierMethod, 'none')
        fprintf(fid, 'Removed %d pre-fit outliers from %d observations (%.1f %%)\n', nPreOutliers, nPreObs, nPreOutliers/nPreObs*100);
    end
    if (any(strcmp(removeOutliers, {'post', 'prepost'})) || removeOutliers) && ~strcmp(outlierMethod, 'none')
        fprintf(fid, 'Removed %d post-fit outliers from %d observations (%.1f %%) and refitted model\n', nOutliers, nObs, nOutliers/nObs*100);
    end
    fprintf(fid, '\t%s', mdlOutput);
    fclose(fid);

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
    end
end

%% Create ANOVA table

for iFit = 1:nFits

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

    % save options
    fpath = fullfile(outDir, 'options.mat');
    save(fpath, 'options');
    fpath = fullfile(outDir, 'makeOptions.m');
    fid = fopen(fpath, 'w+');
    fields = fieldnames(options);
    fields = setdiff(fields, 'Data', 'stable'); % remove "Data" from options
    fields = setdiff(fields, 'DataRaw', 'stable'); % remove "DataRaw" from options
    paramFields = sort(fields); % sort fields alphabetically
    for iField = 1:length(paramFields)
        field = paramFields{iField};
        fprintf(fid, 'options.%s = %s;\n', field, mat2str(string(options.(field))));
    end
    fclose(fid);

    % save Data
    saveTable(Data, 'Data', {'csv'}, outDir);

    % save raw Data table
    saveTable(DataRaw, 'DataRaw', {'csv'}, outDir);

    % save ANOVA table
    saveTable(anovaTable, 'Anova', {'xlsx'}, outDir);
    disp(anovaTable) % display table
    
end

%% Calc plot data

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
main_p      = nan(nPairs, nY);
main_test   = nan(nPairs, nY);
main_DF     = nan(nPairs, nY);
main_aux    = nan(nPairs, nY);
main_eff    = nan(nPairs, nY);

maxNValues = 0;
for iFit = 1:nY
    myVar = y{iFit};
    if nY > 1
        myMaxNValues = max(cell2mat(arrayfun(@(s1,s2) sum(Data.(memberVar)==s1 & Data.(groupVar)==s2 & Data.(yVar)==myVar), repmat(members(:)',nGroups,1), repmat(groups(:),1,nMembers), 'UniformOutput', false)),[],'all');
    else
        myMaxNValues = max(cell2mat(arrayfun(@(s1,s2) sum(Data.(memberVar)==s1 & Data.(groupVar)==s2), repmat(members(:)',nGroups,1), repmat(groups(:),1,nMembers), 'UniformOutput', false)),[],'all');
    end
    maxNValues = max([maxNValues, myMaxNValues]);
end
violin_values = nan(nGroups, nMembers, maxNValues, nRows, nCols, nY);

for iFit = 1:nY
    myVar = y{iFit};

    % init statistics table
    Stats = table;    

    % calc plot data and fill arrays
    for iRow = 1:nRows
        for iCol = 1:nCols
            for iGroup = 1:nGroups
                % loop over 1st x
                for iMember = 1:nMembers

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

                    violin_values(iGroup, iMember, 1:length(values), iRow, iCol, iFit) = values; %%%%%%%

                    % plot values
                    switch errorBars
                        case 'std'
                            bar_values(iRow, iCol, iGroup, iMember, iFit) = statsRow.mean;
                            bar_errorBottom(iRow, iCol, iGroup, iMember, iFit) = statsRow.mean - statsRow.std;
                            bar_errorTop(iRow, iCol, iGroup, iMember, iFit) = statsRow.mean + statsRow.std;
                        case 'se'
                            bar_values(iRow, iCol, iGroup, iMember, iFit) = statsRow.mean;
                            bar_errorBottom(iRow, iCol, iGroup, iMember, iFit) = statsRow.mean - statsRow.SE;
                            bar_errorTop(iRow, iCol, iGroup, iMember, iFit) = statsRow.mean + statsRow.SE;
                        case 'q25'
                            bar_values(iRow, iCol, iGroup, iMember, iFit) = statsRow.median;
                            bar_errorBottom(iRow, iCol, iGroup, iMember, iFit) = statsRow.q25;
                            bar_errorTop(iRow, iCol, iGroup, iMember, iFit) = statsRow.q75;
                        case 'ci95'
                            bar_values(iRow, iCol, iGroup, iMember, iFit) = statsRow.mean;
                            bar_errorBottom(iRow, iCol, iGroup, iMember, iFit) = statsRow.q25;
                            bar_errorTop(iRow, iCol, iGroup, iMember, iFit) = statsRow.q75;
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
                            val1 = Data.(trnsVar)(L1);
                            val2 = Data.(trnsVar)(L2);
                            switch posthocMethod
                                case 'ttest'
                                    [~, bar_p(iGroup, iPair, iRow, iCol, iFit), ~, stats] = ttest2(val1, val2);
                                    tValue = stats.tstat;
                                    df = stats.df;
                                    sPool = stats.sd;
                                    dCohen = (mean(val1, 'omitnan') - mean(val2, 'omitnan')) / sPool;
                                    bar_test(iGroup, iPair, iRow, iCol, iFit) = tValue;
                                    bar_DF(iGroup, iPair, iRow, iCol, iFit) = df;
                                    bar_eff(iGroup, iPair, iRow, iCol, iFit) =  dCohen;
                                case 'utest'
                                    [bar_p(iGroup, iPair, iRow, iCol, iFit), ~, stats] = ranksum(val1, val2);
                                    N1 = sum(~isnan(val1));
                                    N2 = sum(~isnan(val2));
                                    W = stats.ranksum;
                                    U = W - N1*(N1+1)/2;
                                    r = 1 - 2*U/(N1*N2);
                                    bar_test(iGroup, iPair, iRow, iCol, iFit) = U;
                                    bar_eff(iGroup, iPair, iRow, iCol, iFit) =  r;
                            end

                            % calc main contrasts
                            if posthocMain
                                L1 = (Data.(memberVar) == pair(1));
                                L2 = (Data.(memberVar) == pair(2));
                                val1 = Data.(trnsVar)(L1);
                                val2 = Data.(trnsVar)(L2);
                                switch posthocMethod
                                    case 'ttest'
                                        [~, main_p(iPair), ~, stats] = ttest2(val1, val2);
                                        tValue = stats.tstat;
                                        df = stats.df;
                                        sPool = stats.sd;
                                        dCohen = (mean(val1, 'omitnan') - mean(val2, 'omitnan')) / sPool;
                                        main_test(iPair, iFit) = tValue;
                                        main_DF(iPair, iFit) = df;
                                        main_eff(iPair, iFit) =  dCohen;
                                    case 'utest'
                                        [main_p(iPair, iFit), ~, stats] = ranksum(val1, val2);
                                        N1 = sum(~isnan(val1));
                                        N2 = sum(~isnan(val2));
                                        W = stats.ranksum;
                                        U = W - N1*(N1+1)/2;
                                        r = 1 - 2*U/(N1*N2);
                                        main_test(iPair, iFit) = U;
                                        main_eff(iPair, iFit) =  r;
                                end
                            end

                        case 'emm' % perform post-hoc analysis using emmeans

                            if nY > 1 && separateMulti
                                mdl = mdls{iFit};
                            end

                            emm = emmeans(mdl, 'effects', 'unbalanced');

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
                            bar_p(iGroup, iPair, iRow, iCol, iFit) = contrasts.pVal;
                            bar_test(iGroup, iPair, iRow, iCol, iFit) = contrasts.F;
                            bar_DF(iGroup, iPair, iRow, iCol, iFit) = contrasts.DF1;
                            bar_aux(iGroup, iPair, iRow, iCol, iFit) = contrasts.DF2;
                            bar_eff(iGroup, iPair, iRow, iCol, iFit) = f2etaSqp(contrasts.F, contrasts.DF1, contrasts.DF2);

                            % calc main contrasts
                            if posthocMain
                                L1 = idxDep & (emm.table.(memberVar) == pair(1));
                                L2 = idxDep & (emm.table.(memberVar) == pair(2));
                                L = (L1 - L2)';
                                contrasts = kbcontrasts_wald(mdl, emm, L);
                                main_p(iPair, iFit) = contrasts.pVal;
                                main_test(iPair, iFit) = contrasts.F;
                                main_DF(iPair, iFit) = contrasts.DF1;
                                main_aux(iPair, iFit) = contrasts.DF2;
                                main_eff(iPair, iFit) = f2etaSqp(contrasts.F, contrasts.DF1, contrasts.DF2);
                            end
                    end
                end
            end
        end
    end    
end

%% Statistical correction for multiple tests

% If pairwise t-test, the posthoc comparisons must be statistically corrected
bar_pCorr = bar_p;
main_pCorr = main_p;
if any(strcmp(posthocMethod, {'ttest', 'utest', 'emm'}))
    % Holm-Bonferroni correction of all p-values so far
    sizeOrig = size(bar_p); % store original dimensions
    bar_p = bar_p(:); % make column vector
    bar_pCorr = bar_p; % create array of corrected p-values
    idx = ~isnan(bar_p); % identify NaN-entries
    [~, bar_pCorr(idx)] = bonferroni_holm(bar_p(idx)); % correct p-values, omitting NaNs
    bar_p = reshape(bar_p, sizeOrig); % restore original dimensions of p-value array
    bar_pCorr = reshape(bar_pCorr, sizeOrig); % bring corrected p-value array into the same shape as p-value array
    % statistical correction of main posthoc p-Values
    if posthocMain
        idx = ~isnan(main_pCorr); % identify NaN-entries
        [~, main_pCorr(idx)] = bonferroni_holm(main_p(idx)); % correct p-values, omitting NaNs
    end
elseif nY > 1 && separateMulti && strcmp(posthocMethod, 'emm')
    % Since 'emm' posthoc analysis evaluates each fit in one go, hence must
    % only be corrected for the number of fits.
    bar_pCorr = sidak_corr(bar_p, nFits);
    if posthocMain
        main_pCorr = sidak_corr(main_p, nFits);
    end
end

%% Plot data

outDirOrig = outDir;
for iFit = 1:nY
    myVar = y{iFit};

    % set output folder
    if nY > 1 && separateMulti        
        outDir = sprintf('%s/%s', outDirOrig, myVar);
    end

    % figure size
    panelWidth = 600; % width of each panel
    panelHeight = 300; % height of each panel

    if isPlot

        figWidth = nCols * panelWidth;
        figHeight = nRows * panelHeight + 50;
        figName = 'DataPlots';
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
                if ~isempty(yLabel)
                    yLabelStr = yLabel;
                else
                    if strcmp(depVar, yVal)
                        labelStr = myVar;
                    else
                        labelStr = depVar;
                    end
                    if isempty(depVarUnits)
                        yLabelStr = sprintf('%s', labelStr, depVarUnits);
                    else
                        yLabelStr = sprintf('%s [%s]', labelStr, depVarUnits);
                    end
                    yLabelStr = strrep(yLabelStr,'_', ' ');
                end

                switch showVarNames
                    case {1, 2} % display variable names and levels
                        if length(cols) > 1 && length(rows) > 1
                            panelTitle = sprintf('%s = %s, %s = %s', displayColVar, displayCols(iCol), displayRowVar, displayRows(iRow));
                        elseif length(rows) > 1
                            panelTitle = sprintf('%s = %s', displayRowVar, displayRows(iRow));
                        elseif length(cols) > 1
                            panelTitle = sprintf('%s = %s', displayColVar, displayCols(iCol));
                        else
                            panelTitle = '';
                        end
                    otherwise % only display variable levels
                        if length(cols) > 1 && length(rows) > 1
                            panelTitle = sprintf('%s, %s', displayCols(iCol), displayRows(iRow));
                        elseif length(rows) > 1
                            panelTitle = sprintf('%s', displayRows(iRow));
                        elseif length(cols) > 1
                            panelTitle = sprintf('%s', displayCols(iCol));
                        else
                            panelTitle = '';
                        end
                end
                plotGroups(violin_values(:, :, :, iRow, iCol, iFit), displayMembers, displayGroups, displayMemberVar, displayGroupVar, bar_pCorr(:, :, iRow, iCol, iFit), panelTitle, yLabelStr, plotStyle, panel, showVarNames);
            end
        end

        % rescale plots to achieve the same scale for all panels
        if isRescale
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

    % save descriptive statistics
    if nY > 1
        outDir = sprintf('%s/%s', outDirOrig, myVar);
        if ~isfolder(outDir)
            mkdir(outDir);
        end
    end
    fileName = 'Statistics';
    saveTable(Stats, fileName, {'xlsx'}, outDir);
    outDir = outDirOrig;

    disp(Stats) % display table

    %% Post-hoc table

    posthocTable = table;

    % posthoc main effects
    if posthocMain
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
            tableRow.p = main_p(iPair, iFit);
            tableRow.pCorr = main_pCorr(iPair, iFit);
            switch posthocMethod
                case 'ttest'
                    tableRow.t = main_test(iPair, iFit);
                    tableRow.DF = main_DF(iPair, iFit);
                    tableRow.d = main_eff(iPair, iFit);
                    tableRow.effectSize = string(effprint(main_eff(iPair, iFit), 'd'));
                case 'utest'
                    tableRow.U = main_test(iPair, iFit);
                    tableRow.r = main_eff(iPair, iFit);
                    tableRow.effectSize = string(effprint(main_eff(iPair, iFit), 'r'));
                case 'emm'
                    tableRow.F = main_test(iPair, iFit);
                    tableRow.DF1 = main_DF(iPair, iFit);
                    tableRow.DF2 = main_aux(iPair, iFit);
                    tableRow.etaSqp = main_eff(iPair, iFit);
                    tableRow.effectSize = string(effprint(main_eff(iPair, iFit), 'eta2'));
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
                    tableRow.p = bar_p(iGroup, iPair, iRow, iCol, iFit);
                    tableRow.pCorr = bar_pCorr(iGroup, iPair, iRow, iCol, iFit);
                    switch posthocMethod
                        case 'ttest'
                            tableRow.t = bar_test(iGroup, iPair, iRow, iCol, iFit);
                            tableRow.DF = bar_DF(iGroup, iPair, iRow, iCol, iFit);
                            tableRow.d = bar_eff(iGroup, iPair, iRow, iCol, iFit);
                            tableRow.effectSize = string(effprint(bar_eff(iGroup, iPair, iRow, iCol, iFit), 'd'));
                        case 'utest'
                            tableRow.U = bar_test(iGroup, iPair, iRow, iCol, iFit);
                            tableRow.r = bar_eff(iGroup, iPair, iRow, iCol, iFit);
                            tableRow.effectSize = string(effprint(bar_eff(iGroup, iPair, iRow, iCol, iFit), 'r'));
                        case 'emm'
                            tableRow.F = bar_test(iGroup, iPair, iRow, iCol, iFit);
                            tableRow.DF1 = bar_DF(iGroup, iPair, iRow, iCol, iFit);
                            tableRow.DF2 = bar_aux(iGroup, iPair, iRow, iCol, iFit);
                            tableRow.etaSqp = bar_eff(iGroup, iPair, iRow, iCol, iFit);
                            tableRow.effectSize = string(effprint(bar_eff(iGroup, iPair, iRow, iCol, iFit), 'eta2'));
                    end
                    tableRow.significance = string(sigprint(tableRow.pCorr));
                    posthocTable = [posthocTable; tableRow]; %#ok<AGROW>
                end
            end
        end
    end
    disp(posthocTable); % display table

    % save table
    if nY > 1
        outDir = sprintf('%s/%s', outDirOrig, myVar);
        if ~isfolder(outDir)
            mkdir(outDir);
        end
    end
    fileName = 'Posthoc';
    saveTable(posthocTable, fileName, {'xlsx'}, outDir);
    outDir = outDirOrig;
    
end

end

function value = getValue(in)
if ischar(in) || isstring(in)
    value = str2num(in); %#ok<ST2NM>
else
    value = in;
end
end
