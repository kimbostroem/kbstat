function plotGroups(values, members, groups, memberName, groupName, bar_pCorr, plotTitle, ylabelStr, plotStyle, parent, showVarNames, markerSize, barType, barCenter, barBottom, barTop, plotLines, sortValues, groupSize, Outliers, connectPoints)

nArgs = 21;

if nargin < nArgs
    connectPoints = false;
end

if nargin < nArgs-1
    Outliers = [];
end

if nargin < nArgs-2
    groupSize = [];
end

if nargin < nArgs-3
    sortValues = '';
end

if nargin < nArgs-4
    plotLines = false;
end

if nargin < nArgs-5
    barBottom = [];
end

if nargin < nArgs-6
    barTop = [];
end

if nargin < nArgs-7
    barCenter = [];
end

if nargin < nArgs-8
    barType = 'median';
end

if nargin < nArgs-9
    markerSize = [];
end

if nargin < nArgs-10
    showVarNames = 1;
end
if nargin < nArgs-11
    parent = [];
end

if nargin < nArgs-12
    plotStyle = 'violin';
end

if contains(showVarNames, 'names', 'IgnoreCase', true) && contains(showVarNames, 'levels', 'IgnoreCase', true) % display variable names and levels
    members = strcat(memberName, {' = '}, cellstr(members));
    groups = strcat(groupName, {' = '}, cellstr(groups));
end

members = string(members);
groups = string(groups);
nMembers = max([1, length(members)]);
nGroups = max([1, length(groups)]);

if isempty(markerSize) || isnan(markerSize)
    msize = @(x) max([2,24 - 1.5*floor(x/10)]);
    myValues = values;
    myValues(isnan(myValues)) = [];
    nPoints = numel(myValues)/nMembers/nGroups;
    markerSize = msize(nPoints);
end

colors = lines(nMembers);

if ~isempty(parent)
    htl = tiledlayout(parent, 1,nGroups, 'TileSpacing','none');
else
    htl = tiledlayout(1,nGroups, 'TileSpacing','none');
end
pairIdxs = nchoosek(1:nMembers, 2);
nPairs = size(pairIdxs, 1);
barPositions = nan(nGroups, nMembers);
for iGroup = 1:nGroups
    for iMember = 1:nMembers
        barPositions(iGroup, iMember) = iMember;
    end
end
hnt = gobjects(nGroups, 1);
for iGroup = 1:nGroups
    hnt(iGroup) = nexttile(htl, iGroup);
    switch barType
        case 'custom'
            % take over arrays provided by user
        case 'median'
            barCenter(iGroup, :) = squeeze(quantile(values(iGroup,:,:), 0.5, 3));
            barBottom(iGroup, :) = squeeze(quantile(values(iGroup,:,:), 0.25, 3));
            barTop(iGroup, :) = squeeze(quantile(values(iGroup,:,:), 0.75, 3));
        case 'mean'
            barCenter(iGroup, :) = mean(squeeze(values(iGroup,:,:)), 2, 'omitnan');
            barBottom(iGroup, :) = barCenter(iGroup, :) - std(squeeze(values(iGroup,:,:)), 0, 2, 'omitnan');
            barTop(iGroup, :) = barCenter(iGroup, :) + std(squeeze(values(iGroup,:,:)), 0, 2, 'omitnan');
        case 'meanSE'
            barCenter(iGroup, :) = mean(squeeze(values(iGroup,:,:)), 2, 'omitnan');
            barBottom(iGroup, :) = barCenter(iGroup, :) - std(squeeze(values(iGroup,:,:)), 0, 2, 'omitnan') / size(squeeze(values(iGroup,:,:)), 2);
            barTop(iGroup, :) = barCenter(iGroup, :) + std(squeeze(values(iGroup,:,:)), 0, 2, 'omitnan') / size(squeeze(values(iGroup,:,:)), 2);
        case 'meanCI'
            barCenter(iGroup, :) = mean(squeeze(values(iGroup,:,:)), 2, 'omitnan');
            barBottom(iGroup, :) = NaN(1, nMembers);
            barTop(iGroup, :) = NaN(1, nMembers);
            for iMember = 1:nMembers
                pd = fitdist(squeeze(values(iGroup,iMember,:)),'Normal');
                ci = paramci(pd);
                barBottom(iGroup, iMember) = ci(1,1);
                barTop(iGroup, iMember) = ci(2,1);
            end
    end

    switch sortValues
        case {'ascend', 'descend'}
            [barCenter(iGroup, :), sortIdx] = sort(barCenter(iGroup, :), sortValues);
            barBottom(iGroup, :) = barBottom(iGroup, sortIdx);
            barTop(iGroup, :) = barTop(iGroup, sortIdx);
            values(iGroup,:,:) = values(iGroup, sortIdx, :);
            barPositions(iGroup, sortIdx) = barPositions(iGroup, :);
            sortedMembers = members(sortIdx);
            sortedColors = colors(sortIdx,:);
        otherwise
            sortedMembers = members;
            sortedColors = colors;
            % do nothing
    end

    switch plotStyle

        case 'violin'

            % data to be plotted
            xyData = squeeze(values(iGroup,:,:))';

            dataStyle = 'scatter';
            if markerSize == 0
                dataStyle = 'none';
            end

            switch barType

                case {'auto', 'custom'}                    
                    violinAxes = violinplot(xyData, [], 'MedianMarkerSize', 48, 'BoxColor', 0.2*[1 1 1], 'MarkerSize', markerSize, 'ViolinColor', sortedColors, 'ShowBox', false, 'ShowMedian', false, 'DataStyle', dataStyle);

                    % plot bar inside violins
                    hold on
                    for iMember = 1:nMembers
                        pos = iMember;
                        % Bar ends
                        plot([pos pos], [barBottom(iGroup, iMember) barTop(iGroup, iMember)], 'LineWidth', 2, 'Color', 0.2*[1 1 1]);
                        % Center
                        scatter(pos, barCenter(iGroup, iMember), 48, [1 1 1], 'filled', 'MarkerEdgeColor', 0.2*[1 1 1]);
                    end

                case 'mean'
                    violinAxes = violinplot(xyData, [], 'MedianMarkerSize', 48, 'BoxColor', 0.2*[1 1 1], 'MarkerSize', markerSize, 'ViolinColor', sortedColors, 'ShowMean', true);

                otherwise
                    violinAxes = violinplot(xyData, [], 'MedianMarkerSize', 48, 'BoxColor', 0.2*[1 1 1], 'MarkerSize', markerSize, 'ViolinColor', sortedColors);
            end

            if connectPoints && markerSize > 0
                violins = {violinAxes.ScatterPlot};
                nViolins = length(violins);
                for iViolin = 1:nViolins-1
                    nLines = length(violins{iViolin}.XData);
                    for iLine = 1:nLines
                        xData = [violins{iViolin}.XData(iLine), violins{iViolin+1}.XData(iLine)];
                        yData = [violins{iViolin}.YData(iLine), violins{iViolin+1}.YData(iLine)];
                        plot(xData, yData, 'Color', [0 0 0 0.2])
                    end
                end
            end

        case 'boxplot'

            % data to be plotted
            xyData = squeeze(values(iGroup,:,:))';
            
            % do not plot outliers here, as they are plotted separately later. Whiskers extend to the entire range
            boxplot(xyData, 'Colors', sortedColors, 'Symbol', '', 'Whisker', Inf);

        case 'boxchart'
            hold on
            if ~ismember(barType, {'auto', 'median'}); warning('This barType is not supported for plotStyle boxchart. The default barType is used instead'); end
            switch barType
                case 'median'
                    for iMember = 1:nMembers
                        yData = squeeze(values(iGroup,iMember,:));
                        boxchart(iMember*ones(numel(yData),1), yData, 'notch', 'on', 'MarkerStyle', '.', 'JitterOutliers','on', 'MarkerSize', markerSize, ...
                            'BoxFaceColor',sortedColors(iMember,:), 'MarkerColor',sortedColors(iMember,:));
                    end
                    xticks(1:nMembers)

                otherwise
                    paramArray = struct;
                    for iMember = 1:nMembers                        
                        if ~isempty(Outliers)
                            paramArray(iMember).outliers = squeeze(Outliers(iGroup,iMember,:)); % pre-identified outliers
                        end
                        cValues = squeeze(values(iGroup,iMember,:));
                        paramArray(iMember).topWhisker = max(cValues); % non-outlier max
                        paramArray(iMember).topBox = quantile(cValues, 0.75); % 75%
                        paramArray(iMember).topNotch = barTop(iGroup, iMember); % CI top
                        paramArray(iMember).center = barCenter(iGroup, iMember); % EMM
                        paramArray(iMember).bottomNotch = barBottom(iGroup, iMember); % CI bottom
                        paramArray(iMember).bottomBox = quantile(cValues, 0.25); % 25%
                        paramArray(iMember).bottomWhisker = min(cValues); % non-outlier min
                    end
                    kbboxchart(paramArray);
                    % xticks(1:nMembers)
            end

        case 'bar'
            hbar = bar(1:nMembers, barCenter, 'LineStyle', 'none', 'FaceColor', 'flat');
            for iMember = 1:nMembers
                hbar.CData(iMember,:) = sortedColors(iMember,:);
                hold on
                errorbar(iMember, barCenter(iGroup, iMember), barCenter(iGroup, iMember) - barBottom(iGroup, iMember), barTop(iGroup, iMember) - barCenter(iGroup, iMember), 'LineStyle', 'none', 'Color', 0.6*[1 1 1], 'LineWidth', 2);
            end
    end

    % plot outliers
    if ~isempty(Outliers)
        hold on
        for iMember = 1:nMembers
            cOutliers = squeeze(Outliers(iGroup,iMember,:));
            if sum(~isnan(cOutliers)) == 1
                jitterstrength = 0;
            else
                jitterstrength = 1;
            end
            scatter(iMember+rand(numel(cOutliers),1)*jitterstrength, cOutliers, markerSize, sortedColors(iMember,:), '*');
        end
    end

    % plot horizontal lines at values if desired
    if plotLines && size(barCenter, 1) == nGroups && size(barCenter, 2) == nMembers
        for iMember = 1:nMembers
            yline(barCenter(iGroup, iMember), '-', 'Color', sortedColors(iMember,:));
        end
    end

    % xaxis labeling
    if strlength(groupName) > 0 && ~strcmp(groupName, 'none')
        title(groups(iGroup), 'interpreter', 'none');
    end
    xticklabels(sortedMembers);

    % yaxis labeling
    if iGroup == 1
        ylabel(ylabelStr, 'interpreter', 'none');
    end
    if iGroup > 1
        hnt(iGroup).YAxis.TickValues = [];
    end
    set(gca,'TickLabelInterpreter','none');

    % significance
    for iPair = 1:nPairs
        pairIdx = pairIdxs(iPair, :);
        if bar_pCorr(iGroup, iPair) < 0.05
            sigstar({barPositions(iGroup, pairIdx)}, bar_pCorr(iGroup, iPair));
        end
    end
end
linkaxes(hnt, 'y');

% Plot number of subjects if desired
if ~isempty(groupSize)
    NperGroupMember = sum(~isnan(values),3);
    NperGroup = sum(NperGroupMember,2);
    for iGroup = 1:nGroups
        text(hnt(iGroup), .02, .96, sprintf('N=%g', groupSize(iGroup)), 'Fontsize', 12, 'Units', 'Normalized'),
    end
end

% add super Title
if ~isempty(plotTitle)
    title(htl, plotTitle, 'interpreter', 'none')
end

end