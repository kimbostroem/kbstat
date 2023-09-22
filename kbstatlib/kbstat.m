function mdl = kbstat(options)
%% Analyze data using generalized linear mixed-model fit.
% The result is a fit summary, a diagnostic plot, a data bar plot, a
% descriptive statistics table, an ANOVA table, and a posthoc pairwise
% comparison table.
%
% SYNTAX
% kbstat(in, options)
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
%       y               Name of the dependent variable.
%
%       yUnits          Physical units of the dependent variable.
%                       OPTIONAL, default = '1'.
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
%                       "ttest"     t-test plus Holm-Bonferroni correction
%                       "emm"       Extract contrasts from linear model fit
%                       "none"      Do not perform posthoc analysis
%                       OPTIONAL, default = 'emm'.
%
%       isRescale       Flag if the y-axis of each
%                       panel of the data plot is to be resized to a
%                       common scale.
%                       OPTIONAL, default = false.
%
%       removeOutliers  Flag if outliers should be removed
%                       from the data before analysis.
%                       OPTIONAL, default = true.
%
%       outlierThreshold If "removeOutliers" is set, outliers are detected
%                       as points lying outside the interquartile range
%                       (IQR) multiplied by outlierThreshold.
%                       OPTIONAL, default = 1.5.
%
%       outlierLevel    If "removeOutliers" is set, outliers are detected
%                       by constraining the dataset to values of the the
%                       n-th independent variable, where n = outlierLevel.
%                       OPTIONAL, default = length(x), i.e. number of
%                       independent variables.
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
%       applyAbs        Flag if "abs" is to be applied to the values of the
%                       dependent variable before doing anything else.
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
nObsRaw = size(DataRaw, 1);

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

yVar = 'yVar';
if nY > 1
    Data2 = stack(Data1, y, 'NewDataVariableName', 'Y', 'IndexVariableName', yVar);
    Data2.multiVar = categorical(string(Data2.multiVar));
    depVar = 'Y'; % set dependent variable to y
else
    Data2 = Data1;
    depVar = y{1};
    depVarUnits = yUnits{1};
end

%% more variables

% subject variable
if isfield(options, 'applyAbs') && ~isempty(options.applyAbs)
    applyAbs = getValue(options.applyAbs);
else
    applyAbs = false;
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

% posthoc method
if isfield(options, 'posthocMethod') && ~isempty(options.posthocMethod)
    posthocMethod = options.posthocMethod;
else
    posthocMethod = 'emm';
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
    isRescale = false;
end

% flag if outliers should be removed
if isfield(options, 'removeOutliers') && ~isempty(options.removeOutliers)
    removeOutliers = getValue(options.removeOutliers);
else
    removeOutliers = true;
end
nOutliers = 0;

% threshold factor for removing outliers
if isfield(options, 'outlierThreshold') && ~isempty(options.outlierThreshold)
    outlierThreshold = getValue(options.outlierThreshold);
else
    outlierThreshold = 1.5;
end

% level for removing outliers
if isfield(options, 'outlierLevel') && ~isempty(options.outlierLevel)
    outlierLevel = getValue(options.outlierLevel);
else
    outlierLevel = length(x);
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
    for iVar = 1:length(auxVars)
        Data2 = removevars(Data2, auxVars{iVar});
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

% get dependent variable
Data.(depVar) = Data2.(depVar);
Data.(depVar) = Data2.(depVar);
if applyAbs
    Data.(depVar) = abs(Data.(depVar));
end


% store Data table in options
options.Data = Data;

% % remove duplicate rows (can be caused by constraints on a redundant higher-level table
% Data = unique(Data,'rows');

%% Create output folder, if not existing
if ~isfolder(outDir)
    mkdir(outDir);
end

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

if removeOutliers

    idxOut = false(size(Data.(depVar), 1), 1);
    idxTest = true(size(Data.(depVar), 1), 1);

    if outlierLevel == 0 % level 0: all data
        yData = Data.(depVar)(idxTest);
        idxOut = getOutliers(yData, outlierThreshold);

    elseif outlierLevel == 1 % level 1: 1st dependent variable
        for iMember = 1:nMembers
            idxTest = true(size(Data.(depVar), 1), 1);
            member = members(iMember);
            idxTest = idxTest & (Data.(memberVar) == member);
            yData = Data.(depVar)(idxTest);
            idxOut(idxTest) = getOutliers(yData, outlierThreshold);
        end

    elseif outlierLevel == 2 % level 2: 2nd dependent variable, if given
        for iMember = 1:nMembers
            for iGroup = 1:nGroups
                idxTest = true(size(Data.(depVar), 1), 1);
                member = members(iMember);
                idxTest = idxTest & (Data.(memberVar) == member);
                if nGroups > 1
                    group = groups(iGroup);
                    idxTest = idxTest & (Data.(groupVar) == group);
                end
                yData = Data.(depVar)(idxTest);
                idxOut(idxTest) = getOutliers(yData, outlierThreshold);
            end
        end

    elseif outlierLevel == 3 % level 3: 3rd dependent variable, if given
        for iMember = 1:nMembers
            for iGroup = 1:nGroups
                for iCol = 1:nCols
                    idxTest = true(size(Data.(depVar), 1), 1);
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
                    yData = Data.(depVar)(idxTest);
                    idxOut(idxTest) = getOutliers(yData, outlierThreshold);
                end
            end
        end
    elseif outlierLevel == 4 % level 4: 4th dependent variable, if given
        for iMember = 1:nMembers
            for iGroup = 1:nGroups
                for iCol = 1:nCols
                    for iRow = 1:nRows
                        idxTest = true(size(Data.(depVar), 1), 1);
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
                        yData = Data.(depVar)(idxTest);
                        idxOut(idxTest) = getOutliers(yData, outlierThreshold);
                    end
                end
            end
        end
    end

    % remove outliers
    nOutliers = sum(idxOut);
    nObsRaw = size(Data, 1);
    if nOutliers > 0
        Data = Data(~idxOut, :);
        fprintf('Removed %d outlier(s) from %d observations (%.1f %%)\n', nOutliers, nObsRaw, nOutliers/nObsRaw*100);
    end
end




% save Data table
saveTable(Data, 'Data', {'csv'}, outDir);

%% Fit linear model and perform ANOVA

anovaTable = table; % init ANOVA table

productTerm = strjoin(interact, ' * ');
xNoInteract = [{'1'}, setdiff(x, interact, 'stable')];
sumTerm = strjoin(xNoInteract, ' + ');
if ~isempty(id) && length(unique(Data.(id))) > 1
    if randomSlopes
        randomEffect = '';
        randomSlopes = strjoin(cellfun(@(x) sprintf('(%s|%s)', x, id), within, 'UniformOutput', false), ' + ');
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
    terms = [terms, sumTerm];
end
if ~isempty(productTerm)
    terms = [terms, productTerm];
end
if ~isempty(randomEffect)
    terms = [terms, randomEffect];
end
if ~isempty(randomSlopes)
    terms = [terms, randomSlopes];
end

if isempty(formula)
    formula = sprintf('%s ~ %s', depVar, strjoin(terms, ' + '));
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


% get ANOVA table from model fit
results = anova(mdl);

% print results of model fit into file
mdlOutput = formattedDisplayText(mdl);
fid = fopen(fullfile(outDir, 'Summary.txt'), 'w+');
fprintf(fid, 'Formula:\n\t%s\n', formula);
fprintf(fid, 'Removed %d outliers from %d observations (%.1f %%)\n', nOutliers, nObsRaw, nOutliers/nObsRaw*100);
fprintf(fid, '\t%s', mdlOutput);
fclose(fid);

anovaTable.Term = string(results.Term(2:end));
anovaTable.DF1 = results.DF1(2:end);
anovaTable.DF2 = results.DF2(2:end);
anovaTable.F = results.FStat(2:end);
anovaTable.p = results.pValue(2:end);
anovaTable.etaPSquare = f2etaSqp(anovaTable.F, anovaTable.DF1, anovaTable.DF2);
anovaTable.effectSize = string(etaprint(anovaTable.etaPSquare));
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

% save raw Data table
saveTable(DataRaw, 'DataRaw', {'csv'}, outDir);

% save ANOVA table
saveTable(anovaTable, 'Anova', {'xlsx'}, outDir);
disp(anovaTable) % display table

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

sgtitle(sprintf('Diagnostics for %s', plotTitle), 'interpreter', 'none', 'FontWeight', 'bold', 'FontSize', 14);
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

% probability
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
set(fig, 'PaperPositionMode', 'auto');
set(fig, 'PaperUnits', 'points');
set(fig, 'PaperSize', [figWidth figHeight]);
print(fig, fullfile(outDir, sprintf('%s.pdf', figName)), '-fillpage', '-dpdf', sprintf('-r%.0f', 300));


%% Calc plot data

% figure size
panelWidth = 600; % width of each panel
panelHeight = 300; % height of each panel

% define posthoc comparison pairs
pairs = nchoosek(members, 2);
nPairs = size(pairs, 1);

% init statistics table
Stats = table;

% create (nRows x nCols x nGroups x nMembers) arrays of plot data
bar_values = nan(nRows, nCols, nGroups, nMembers);
bar_errorTop = nan(nRows, nCols, nGroups, nMembers);
bar_errorBottom = nan(nRows, nCols, nGroups, nMembers);
bar_p = nan(nGroups, nPairs, nRows, nCols);
main_p = nan(nPairs);

maxNValues = max(cell2mat(arrayfun(@(s1,s2) sum(Data.(memberVar)==s1 & Data.(groupVar)==s2), repmat(members(:)',nGroups,1), repmat(groups(:),1,nMembers), 'UniformOutput', false)),[],'all');
violin_values = nan(nGroups, nMembers, maxNValues, nRows, nCols);


% calc plot data and fill arrays
for iRow = 1:nRows

    for iCol = 1:nCols

        for iGroup = 1:nGroups

            % loop over 1st x
            for iMember = 1:nMembers

                % init statistics table row
                statsRow = table;

                % 1st x
                member = members(iMember);
                idx = (Data.(memberVar) == member);

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

                violin_values(iGroup, iMember, 1:length(values), iRow, iCol) = values; %%%%%%%

                % plot values
                switch errorBars
                    case 'std'
                        bar_values(iRow, iCol, iGroup, iMember) = statsRow.mean;
                        bar_errorBottom(iRow, iCol, iGroup, iMember) = statsRow.mean - statsRow.std;
                        bar_errorTop(iRow, iCol, iGroup, iMember) = statsRow.mean + statsRow.std;
                    case 'se'
                        bar_values(iRow, iCol, iGroup, iMember) = statsRow.mean;
                        bar_errorBottom(iRow, iCol, iGroup, iMember) = statsRow.mean - statsRow.SE;
                        bar_errorTop(iRow, iCol, iGroup, iMember) = statsRow.mean + statsRow.SE;
                    case 'q25'
                        bar_values(iRow, iCol, iGroup, iMember) = statsRow.median;
                        bar_errorBottom(iRow, iCol, iGroup, iMember) = statsRow.q25;
                        bar_errorTop(iRow, iCol, iGroup, iMember) = statsRow.q75;
                    case 'ci95'
                        bar_values(iRow, iCol, iGroup, iMember) = statsRow.mean;
                        bar_errorBottom(iRow, iCol, iGroup, iMember) = statsRow.q25;
                        bar_errorTop(iRow, iCol, iGroup, iMember) = statsRow.q75;
                end


            end

            % get post-hoc p-values
            for iPair = 1:nPairs
                pair = pairs(iPair, :);

                switch posthocMethod

                    case 'ttest'

                        % init idx
                        idx = true(size(Data, 1), 1);

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
                        val1 = Data.(depVar)(L1);
                        val2 = Data.(depVar)(L2);
                        [~, bar_p(iGroup, iPair, iRow, iCol)] = ttest2(val1, val2);

                        % calc main contrasts
                        L1 = (Data.(memberVar) == pair(1));
                        L2 = (Data.(memberVar) == pair(2));
                        val1 = Data.(depVar)(L1);
                        val2 = Data.(depVar)(L2);
                        [~, main_p(iPair)] = ttest2(val1, val2);

                    case 'emm'

                        % perform post-hoc analysis using emm
                        emm = emmeans(mdl, 'effects', 'unbalanced');

                        % init idx
                        idx = true(size(emm.table, 1), 1);

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
                        contrasts = contrasts_wald(mdl, emm, L);
                        bar_p(iGroup, iPair, iRow, iCol) = contrasts.pVal;

                        % calc main contrasts
                        L1 = (emm.table.(memberVar) == pair(1));
                        L2 = (emm.table.(memberVar) == pair(2));
                        L = (L1 - L2)';
                        contrasts = contrasts_wald(mdl, emm, L);
                        main_p(iPair) = contrasts.pVal;
                end
            end
        end
    end
end

% If pairwise t-test, the posthoc comparisons must be statistically corrected
bar_pCorr = bar_p;
main_pCorr = main_p;
if strcmp(posthocMethod, 'ttest')
    sizeOrig = size(bar_p); % vstore original dimensions
    bar_p = bar_p(:); % make column vector
    bar_pCorr = bar_p; % create array of corrected p-values
    idx = ~isnan(bar_p); % identify NaN-entries
    [~, bar_pCorr(idx)] = bonferroni_holm(bar_p(idx)); % correct p-values, omitting NaNs
    bar_p = reshape(bar_p, sizeOrig); % restore original dimensions of p-value array
    bar_pCorr = reshape(bar_pCorr, sizeOrig); % bring corrected p-value array into the same shape as p-value array
    % statistical correction of main posthoc p-Values
    idx = ~isnan(main_pCorr); % identify NaN-entries
    [~, main_pCorr(idx)] = bonferroni_holm(main_p(idx)); % correct p-values, omitting NaNs
end

%% Plot data

if isPlot

    figWidth = nCols * panelWidth;
    figHeight = nRows * panelHeight + 50;
    figName = 'ViolinPlot';
    fig = figure('Name', figName, 'Position', [0, 0, figWidth, figHeight]);
    layout = tiledlayout(nRows, nCols);
    title(layout, sprintf('Data plots for %s', plotTitle), 'interpreter', 'none', 'FontWeight', 'bold', 'FontSize', 14);

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
            if isempty(depVarUnits)
                ylabelStr = sprintf('%s', yLabel, depVarUnits);
            else
                ylabelStr = sprintf('%s [%s]', yLabel, depVarUnits);
            end
            ylabelStr = strrep(ylabelStr,'_', ' ');
            
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
            plotViolinGroups(violin_values(:, :, :, iRow, iCol), displayMembers, displayGroups, displayMemberVar, displayGroupVar, bar_pCorr(:, :, iRow, iCol), panelTitle, ylabelStr, panel, showVarNames);
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
    set(fig, 'PaperPositionMode', 'auto');
    set(fig, 'PaperUnits', 'points');
    set(fig, 'PaperSize', [figWidth figHeight]);
    print(fig, fullfile(outDir, sprintf('%s.pdf', figName)), '-fillpage', '-dpdf', sprintf('-r%.0f', 300));
    print(fig, fullfile(outDir, sprintf('%s.png', figName)), '-dpng', sprintf('-r%.0f', 300));
    saveas(fig, fullfile(outDir, sprintf('%s.fig', figName)));

end

% save descriptive statistics
saveTable(Stats, 'Statistics', {'xlsx'}, outDir);
disp(Stats) % display table

%% Post-hoc table

posthocTable = table;
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
    tableRow.p = main_p(iPair);
    tableRow.pCorr = main_pCorr(iPair);
    tableRow.significance = string(sigprint(tableRow.pCorr));
    posthocTable = [posthocTable; tableRow]; %#ok<AGROW>
end
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
                tableRow.p = bar_p(iGroup, iPair, iRow, iCol);
                tableRow.pCorr = bar_pCorr(iGroup, iPair, iRow, iCol);
                tableRow.significance = string(sigprint(tableRow.pCorr));
                posthocTable = [posthocTable; tableRow]; %#ok<AGROW>
            end
        end
    end
end

% save table
saveTable(posthocTable, 'Posthoc', {'xlsx'}, outDir);
disp(posthocTable); % display table

% close all figures
close all

end

function value = getValue(in)
if ischar(in) || isstring(in)
    value = str2num(in); %#ok<ST2NM>
else
    value = in;
end
end
