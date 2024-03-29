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
%   in          (char or table) Either the path to an Excel or CSV file,
%               or a Matlab table containing the data. The table must be in
%               long format, so one row per data point.
%
%   options     (struct) Structure containing parameters needed for 
%               analysis. All fields must be char or string, and if a list
%               is to be given, it must be a string of words separated by
%               comma or semicolon.
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
%                       'identity'	    g(mu) = mu.             Default for Normal distribution.
%                       'log'	        g(mu) = log(mu).        Default for Poisson, Gamma, and InverseGaussian.
%                       'logit'	        g(mu) = log(mu/(1-mu))  Default for Binomial distribution.
%                       'loglog'	    g(mu) = log(-log(mu))
%                       'probit'	    g(mu) = norminv(mu)
%                       'comploglog'	g(mu) = log(-log(1-mu))
%                       'reciprocal'	g(mu) = mu.^(-1)
%                       Scalar p	    g(mu) = mu.^p           Canonical for Gamma (p = -1) and InverseGaussian (p= -2)
%
%       posthocMethod   Method for the posthoc pairwise comparison.
%                       Possible values:
%                       "ttest"     t-test plus Holm-Bonferroni correction
%                       "emm"       Extract contrasts from linear model fit
%                       "none"      Do not perform posthoc analysis
%                       OPTIONAL, default = 'emm'.
%
%       legendLocation  Location of the legend in each panel of
%                       the data plot.
%                       OPTIONAL, default = 'best'.
%
%       isRescale       Flag if the y-axis of each
%                       panel of the data plot is to be resized to a
%                       common scale.
%                       OPTIONAL, default = false.
%
%       removeOutliers  Flag if outliers should be removed
%                       from the data before analysis
%
%       thresholdFactor If "removeOutliers" is set, define the factor to
%                       multiply the outlier criterion.
%                       OPTIONAL, default = 3
%
%       constraint      Constrain the data before analysis. Must be
%                       of the form
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
%                       value.
%                       OPTIONAL, default = unset
%
%       outDir          Output folder for generated files.
%                       OPTIONAL, defaults to the parent folder of the
%                       input file. If the input is a data table, the
%                       output folder defaults to the local folder.
%
%       isPlot          Plot data as grouped bars with significance brackets
%
% OUTPUT
%   mdl                 (Generalized linear mixed-effects model) Result from linear model fit
%
% EXAMPLE
%   options.y = 'jointEfficiency';
%   options.yUnits = '1';
%   options.x = 'shoe, speed, sex, joint';
%   options.id = 'subject';
%   options.within = 'shoe, speed, joint';
%   options.interact = 'shoe, speed, joint';
%   options.distribution = 'gamma';
%   options.removeOutliers = 'true';
%   kbstat('path/to/Data.xlsx', options);
%
% (c) 2022 by Predimo GmbH
% Website: http://www.predimo.com
% Author: Kim Boström
% Last revision: 2022-09-17 (KB) Added documentation

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
nObservationsRaw = size(DataRaw, 1);

% init Data table
Data = table;

% all table variables
tableVars = DataRaw.Properties.VariableNames;

% independent variable(s)
x = strtrim(strsplit(options.x, {',', ';'}));

% dependent variable
y = options.y;
if ~ismember(y, tableVars)
    error('Dependent variable "%s" not found in data table', y);
end
if isfield(options, 'yUnits')
    yUnits = options.yUnits;
else
    yUnits = '1';
end
responseVariable = y; % set response variable to y

% subject variable
if isfield(options, 'id')
    id = options.id;
else
    id = '';
end

% within-subject variables
if isfield(options, 'within')
    within = strtrim(strsplit(options.within, {',', ';'}));
    x = union(x, within, 'stable'); % add within-subject variables to independent variables
else
    within = {};
end

% interaction variables
if isfield(options, 'interact')
    interact = strtrim(strsplit(options.interact, {',', ';'}));
    x = union(x, interact, 'stable'); % add interaction variables to independent variables
else
    interact = {};
end

% posthoc method
if isfield(options, 'posthocMethod')
    posthocMethod = options.posthocMethod;
else
    posthocMethod = 'emm';
end

% fit method
if isfield(options, 'fitMethod')
    fitMethod = options.fitMethod;
else
    fitMethod = 'MPL';
end

% distribution (for GLM)
if isfield(options, 'distribution')
    distribution = options.distribution;
else % no parameter given
    distribution = 'normal';
end

% link (for GLM)
if isfield(options, 'link')
    link = options.link;
else % no parameter given
    switch lower(distribution)
        case 'normal'
            link = 'identity';
        case 'lognormal'
            link = 'identity';
            responseVariable = sprintf('log%s', y); % re-define response variable to log of y
            Data.(responseVariable) = log(DataRaw.(y)); % create new column with log'd data
            distribution = 'normal'; % re-set distribution to 'normal'
            fprintf('Using normal-distributed logarithm of y "%s"\n', y); % tell user what's happening
        case 'binomial'
            link = 'logit';
        case 'poisson'
            link = 'log';
        case 'gamma'
            link = 'log';
        case 'inversegaussian'
            link = 'log';
        otherwise
            error('Unknown distribution "%s"', distribution);
    end
end

% Flag to plot data
if isfield(options, 'isPlot')
    isPlot = getValue(options.isPlot);
else
    isPlot = true;
end

% legend location
if isfield(options, 'legendLocation')
    legendLocation = options.legendLocation;
else
    legendLocation = 'Best';
end

% flag to rescale all panel plots to the same y-scale
if isfield(options, 'isRescale')
    isRescale = getValue(options.isRescale);
else
    isRescale = false;
end

% flag if outliers should be removed
if isfield(options, 'removeOutliers')
    removeOutliers = getValue(options.removeOutliers);
else
    removeOutliers = false;
end

if isfield(options, 'thresholdFactor')
    thresholdFactor = getValue(options.thresholdFactor);
else
    thresholdFactor = 4;
end

% output folder
if isfield(options, 'outDir')
    outDir = options.outDir;
else
    outDir = fullfile(inDir, inName);
end
if ~isfolder(outDir)
    mkdir(outDir);
end

fprintf('Performing Linear Model analysis for %s...\n', y);


%% Make variables that are non-numeric or integer, categorical

% init categories
factors = {};

% get subject variable
if ~isempty(id)
    myIV = id;
    myLevels = unique(DataRaw.(id));
    if ~all(isnumeric(myLevels)) || (all(isnumeric(myLevels)) && all(mod(myLevels,1) == 0))
        Data.(myIV) = categorical(DataRaw.(myIV)); % make values categorical
    else
        Data.(myIV) = DataRaw.(myIV); % adopt continuous values    
    end
end

% get independent variables
for iIV = 1:length(x)
    myIV = x{iIV};
    myLevels = unique(DataRaw.(myIV));
    if ~all(isnumeric(myLevels)) || (all(isnumeric(myLevels)) && all(mod(myLevels,1) == 0))
        Data.(myIV) = categorical(DataRaw.(myIV)); % make values categorical
        factors = [factors; myIV]; %#ok<AGROW>
    else
        Data.(myIV) = DataRaw.(myIV); % adopt continuous values
    end
end

% get response variable
Data.(y) = DataRaw.(y);

% store Data table in options
options.Data = Data;

%% Apply constraint, if given

if isfield(options, 'constraint')
    myConstraint = options.constraint;
    parts = strsplit(myConstraint);
    myVar = parts{1};
    cmd = sprintf('Data.%s', myConstraint);
    idx = eval(cmd);
    if any(idx)
        Data = Data(idx, :);
    else
        warning('Constraint %s cannot be fulfilled -> leave data unchanged', cmd);
    end
    nLevels = length(unique(Data.(myVar)));
    switch nLevels
        case 0
            warning('x "%s" has no levels left on constraint "%s" -> remove x');
            x = setdiff(x, myVar);
        case 1
            warning('x "%s" has only 1 level left on constraint "%s" -> remove x');
            x = setdiff(x, myVar);
    end
end


%% Create output folder, if not existing
if ~isfolder(outDir)
    mkdir(outDir);
end

%% Create variables

% 1st factor = member variable
memberVar = factors{1};
members = unique(Data.(memberVar), 'stable');
nMembers = length(members);

% 2nd factor = group variable
if length(factors) > 1
    groupVar = factors{2};
    groups = unique(Data.(groupVar), 'stable');
    nGroups = length(groups);
else
    groups = string(factors{1});
    nGroups = 1;
end

% 3rd factor = column variable
if length(factors) > 2
    colVar = factors{3};
    cols = unique(Data.(colVar), 'stable');
    nCols = length(cols);
else
    colVar = 'none';
    cols = "none";
    nCols = 1;
end

% 4th factors = row variable
if length(factors) > 3
    rowVar = factors{4};
    rows = unique(Data.(rowVar), 'stable');
    nRows = length(rows);
else
    rowVar = 'none';
    rows = "none";
    nRows = 1;
end


%% Fit linear model and perform ANOVA

anovaTable = table; % init ANOVA table

productTerm = strjoin(interact, ' * ');
xNoInteract = [{'1'}, setdiff(x, interact, 'stable')];
sumTerm = strjoin(xNoInteract, ' + ');
if ~isempty(id) && length(unique(Data.(id))) > 1
    randomEffect = sprintf('(1|%s)', id);
    randomSlopes = strjoin(cellfun(@(x) sprintf('(1|%s:%s)', x, id), within, 'UniformOutput', false), ' + ');
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

formula = sprintf('%s ~ %s', responseVariable, strjoin(terms, ' + '));
fprintf('\t%s\n', formula);

try
    mdl = fitglme(Data, formula, ...
        'DummyVarCoding', 'effects', ...
        'FitMethod', fitMethod, ...
        'Distribution', distribution, ...
        'Link', link);
catch ME
    message = sprintf('%s', ME.message);
    fprintf('The linear model fit returned an error:\n\t%s\n', message);
    fprintf('Please try again, using fewer interactions by defining "interact" with only those independent variables whose interaction you want to investigate\n');
    return
end

if removeOutliers
    mdlResiduals = residuals(mdl, 'ResidualType', 'Pearson');
    mdlOutliers = isoutlier(mdlResiduals, 'ThresholdFactor', thresholdFactor);
    nOutliers = sum(mdlOutliers);
    fprintf('Removed %d outlier(s) from %d observations (%.1f %%)\n', nOutliers, nObservationsRaw, nOutliers/nObservationsRaw*100);
    if nOutliers > 0
        fprintf('Re-fitting model...\n');
        Data = Data(~mdlOutliers,:);
        % re-fit linear model
        try
            mdl = fitglme(Data, formula, ...
                'DummyVarCoding', 'effects', ...
                'FitMethod', fitMethod, ...
                'Distribution', distribution, ...
                'Link', link);
        catch ME
            message = sprintf('%s', ME.message);
            fprintf('The linear model fit returned an error after removing outliers:\n\t%s\n', message);
            fprintf('Please try again, using a higher "thresholdFactor" and/or fewer interactions by defining "interact" with only those independent variables whose interaction you want to investigate\n');
            return
        end
    end
else
    nOutliers = 0;
end

% get ANOVA table from model fit
results = anova(mdl);

% print results of model fit into file
mdlOutput = formattedDisplayText(mdl);
fid = fopen(fullfile(outDir, 'Summary.txt'), 'w+');
fprintf(fid, 'Formula:\n\t%s\n', formula);
fprintf(fid, 'Removed %d outliers from %d observations (%.1f %%)\n', nOutliers, nObservationsRaw, nOutliers/nObservationsRaw*100);
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
fields = setdiff(fields, 'Data'); % remove "Data" from options
fields = setdiff(fields, 'DataRaw'); % remove "DataRaw" from options
paramFields = sort(fields); % sort fields alphabetically
for iField = 1:length(paramFields)
    field = paramFields{iField};
    fprintf(fid, 'options.%s = %s;\n', field, mat2str(string(options.(field))));
end
fclose(fid);

% save raw Data table
fpath = fullfile(outDir, 'DataRaw.xlsx');
writetable(DataRaw, fpath, 'WriteMode', 'replacefile');

% save Data table
fpath = fullfile(outDir, 'Data.xlsx');
writetable(Data, fpath, 'WriteMode', 'replacefile');

% save ANOVA table
fpath = fullfile(outDir, 'Anova.xlsx');
writetable(anovaTable, fpath, 'WriteMode', 'replacefile');
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
set(fig, 'renderer', 'painters');
print(fig, fullfile(outDir, sprintf('%s.pdf', figName)), '-fillpage', '-dpdf', sprintf('-r%.0f', 300));


%% Calc plot data

% figure size
panelWidth = 600; % width of each panel
panelHeight = 300; % height of each panel

% define posthoc comparison pairs
pairs = nchoosek(members, 2);
pairIdxs = nchoosek(1:nMembers, 2);
nPairs = size(pairs, 1);

% init statistics table
Stats = table;

% create (nRows x nCols x nGroups x nMembers) arrays of plot data
bar_values = nan(nRows, nCols, nGroups, nMembers);
bar_errorTop = nan(nRows, nCols, nGroups, nMembers);
bar_errorBottom = nan(nRows, nCols, nGroups, nMembers);
bar_p = nan(nRows, nCols, nGroups, nPairs);
main_p = nan(nPairs);

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
                values = bar_data.(y);
                bar_values(iRow, iCol, iGroup, iMember) = median(values, 'omitnan');
                bar_errorBottom(iRow, iCol, iGroup, iMember) = quantile(values, 0.25);
                bar_errorTop(iRow, iCol, iGroup, iMember) = quantile(values, 0.75);
                if nGroups > 1
                    statsRow.(groupVar) = string(group);
                end
                statsRow.(memberVar) = string(member);
                statsRow.N = length(values(~isnan(values)));
                statsRow.median = bar_values(iRow, iCol, iGroup, iMember);
                statsRow.dQ25 = bar_errorBottom(iRow, iCol, iGroup, iMember) - statsRow.median;
                statsRow.dQ75 = bar_errorTop(iRow, iCol, iGroup, iMember) - statsRow.median;
                statsRow.mean = mean(values, 'omitnan');
                statsRow.std = std(values, 'omitnan');
                statsRow.SE = statsRow.std / sqrt(statsRow.N);
                statsRow.RE = statsRow.SE / statsRow.mean * 100;
                Stats = [Stats; statsRow]; %#ok<AGROW>
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
                        val1 = Data.(responseVariable)(L1);
                        val2 = Data.(responseVariable)(L2);
                        [~, bar_p(iRow, iCol, iGroup, iPair)] = ttest2(val1, val2);
                        
                        % calc main contrasts
                        L1 = (Data.(memberVar) == pair(1));
                        L2 = (Data.(memberVar) == pair(2));
                        val1 = Data.(responseVariable)(L1);
                        val2 = Data.(responseVariable)(L2);
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
                        bar_p(iRow, iCol, iGroup, iPair) = contrasts.pVal;

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

% statistical correction of posthoc p-values
sizeOrig = size(bar_p); % vstore original dimensions
bar_p = bar_p(:); % make column vector
bar_pCorr = bar_p; % create array of corrected p-values
idx = ~isnan(bar_p); % identify NaN-entries
[~, bar_pCorr(idx)] = bonferroni_holm(bar_p(idx)); % correct p-values, omitting NaNs
bar_p = reshape(bar_p, sizeOrig); % restore original dimensions of p-value array
bar_pCorr = reshape(bar_pCorr, sizeOrig); % bring corrected p-value array into the same shape as p-value array
% statistical correction of main posthoc p-Values
main_pCorr = main_p; % create array of corrected p-values
idx = ~isnan(main_pCorr); % identify NaN-entries
[~, main_pCorr(idx)] = bonferroni_holm(main_p(idx)); % correct p-values, omitting NaNs

%% Plot data

if isPlot

    figWidth = nCols * panelWidth;
    figHeight = nRows * panelHeight;
    figName = 'DataPlot';
    fig = figure('Name', figName, 'Position', [0, 0, figWidth, figHeight]);
    for iRow = 1:nRows

        for iCol = 1:nCols

            % start panel subplot
            iPanel = sub2ind([nCols, nRows], iCol, iRow);
            subplot(nRows, nCols, iPanel);

            % plot panel
            myBarValues = reshape(bar_values(iRow, iCol, :, :), nGroups, nMembers);
            myBarErrorsBottom = reshape(bar_errorBottom(iRow, iCol, :, :), nGroups, nMembers); 
            myBarErrorsTop = reshape(bar_errorTop(iRow, iCol, :, :), nGroups, nMembers);
            barPositions = plotBarGroups(myBarValues, members, groups, myBarErrorsBottom, myBarErrorsTop, 'none');
            if nGroups > 1
                leg = legend(gca,'Location', legendLocation);
                set(leg,'AutoUpdate','off');
            end
            for iGroup = 1:nGroups
                for iPair = 1:nPairs
                    pairIdx = pairIdxs(iPair, :);
                    if bar_pCorr(iRow, iCol, iGroup, iPair) < 0.05
                        sigstar({barPositions(iGroup, pairIdx)}, bar_pCorr(iRow, iCol, iGroup, iPair));
                    end
                end
            end

            titleAdds = {};

            % 3rd x, if given
            if nCols > 1
                col = cols(iCol);
                titleAdds = [titleAdds, cellstr(sprintf('%s = %s', colVar, col))]; %#ok<AGROW>
            end

            % 4th x, if given
            if nRows > 1
                row = rows(iRow);
                titleAdds = [titleAdds, cellstr(sprintf('%s = %s', rowVar, row))]; %#ok<AGROW>
            end

            titleAddStr = strjoin(titleAdds, ', ');
            if ~isempty(titleAddStr)
                titleStr = sprintf('%s (%s)', y, titleAddStr);
            else
                titleStr = sprintf('%s', y);
            end

            title(titleStr, 'Interpreter', 'none');

            % y-axis label
            ylabel(sprintf('%s [%s]', y, yUnits), 'Interpreter', 'none');
        end
    end

    % rescale plots to achieve the same scale for all panels
    if isRescale
        axs = findobj(fig, 'type', 'axes');
        ylimits = NaN(length(axs), 2);
        for iAx = 1:length(axs)
            ylimits(iAx, :) = get(axs(iAx), 'ylim');
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
    set(fig, 'renderer', 'painters');
    print(fig, fullfile(outDir, sprintf('%s.pdf', figName)), '-fillpage', '-dpdf', sprintf('-r%.0f', 300));
    print(fig, fullfile(outDir, sprintf('%s.png', figName)), '-dpng', sprintf('-r%.0f', 300));
    saveas(fig, fullfile(outDir, sprintf('%s.fig', figName)));

end

% save descriptive statistics
fpath = fullfile(outDir, 'Statistics.xlsx');
writetable(Stats, fpath, 'WriteMode', 'replacefile');
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
                tableRow.p = bar_p(iRow, iCol, iGroup, iPair);
                tableRow.pCorr = bar_pCorr(iRow, iCol, iGroup, iPair);
                tableRow.significance = string(sigprint(tableRow.pCorr));
                posthocTable = [posthocTable; tableRow]; %#ok<AGROW>
            end
        end
    end    
end

% save table
fpath = fullfile(outDir, 'Posthoc.xlsx');
writetable(posthocTable, fpath, 'WriteMode', 'replacefile');
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
