function plotGroups(values, members, groups, memberName, groupName, bar_pCorr, plotTitle, ylabelStr, plotStyle, parent, showVarNames, markerSize, barType, barCenter, barBottom, barTop, plotLines, sortValues)

if nargin < 18
    sortValues = '';
end

if nargin < 17
    plotLines = false;
end

if nargin < 16
    barBottom = [];
end

if nargin < 15
    barTop = [];
end

if nargin < 14
    barCenter = [];
end

if nargin < 13
    barType = 'median';
end

if nargin < 12 || isnan(markerSize)
    markerSize = 2;
end

if nargin < 11
    showVarNames = 1;
end
if nargin < 10
    parent = [];
end

if nargin < 9
    plotStyle = 'violin';
end

if contains(showVarNames, 'names', 'IgnoreCase', true) && contains(showVarNames, 'levels', 'IgnoreCase', true) % display variable names and levels
    members = strcat(memberName, {' = '}, cellstr(members));
    groups = strcat(groupName, {' = '}, cellstr(groups));
end

members = string(members);
groups = string(groups);
nMembers = length(members);
nGroups = length(groups);

ogColors = lines(nMembers);

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
            colors = ogColors(sortIdx,:);
        otherwise
            sortedMembers = members;
            colors = ogColors;
            % do nothing
    end


    switch lower(plotStyle)

         case 'violin'
            switch barType

                case {'auto', 'custom'}
                    violinplot(squeeze(values(iGroup,:,:))', [], 'MedianMarkerSize', 48, 'BoxColor', 0.2*[1 1 1], 'MarkerSize', markerSize, 'ViolinColor', colors, 'ShowBox', false, 'ShowMedian', false);

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
                    violinplot(squeeze(values(iGroup,:,:))', [], 'MedianMarkerSize', 48, 'BoxColor', 0.2*[1 1 1], 'MarkerSize', markerSize, 'ViolinColor', colors, 'ShowMean', true);

                otherwise
                    violinplot(squeeze(values(iGroup,:,:))', [], 'MedianMarkerSize', 48, 'BoxColor', 0.2*[1 1 1], 'MarkerSize', markerSize, 'ViolinColor', colors);
            end

        case 'boxplot'
            boxplot(squeeze(values(iGroup,:,:))', 'Colors', colors);

        case 'boxchart'
            hold on
            for iMember = 1:nMembers
                yData = squeeze(values(iGroup,iMember,:));
                boxchart(iMember*ones(numel(yData),1), yData, 'notch', 'on', 'MarkerStyle', '.', 'JitterOutliers','on', 'MarkerSize', markerSize, ...
                    'BoxFaceColor',colors(iMember,:), 'MarkerColor',colors(iMember,:));
            end
            xticks(1:nMembers)

        case 'bar'
            hbar = bar(1:nMembers, barCenter, 'LineStyle', 'none', 'FaceColor', 'flat');
            for iMember = 1:nMembers
                hbar.CData(iMember,:) = colors(iMember,:);
                hold on
                errorbar(iMember, barCenter(iGroup, iMember), barCenter(iGroup, iMember) - barBottom(iGroup, iMember), barTop(iGroup, iMember) - barCenter(iGroup, iMember), 'LineStyle', 'none', 'Color', 0.6*[1 1 1], 'LineWidth', 2);
            end
    end



    % plot horizontal lines at values if desired
    if plotLines && size(barCenter, 1) == nGroups && size(barCenter, 2) == nMembers
        % colors = lines;
        for iMember = 1:nMembers
            yline(barCenter(iGroup, iMember), '-', 'Color', colors(iMember,:));
        end
    end

    if ~isempty(groupName)
        title(groups(iGroup), 'interpreter', 'none');
    end
    xticklabels(sortedMembers);

    for iPair = 1:nPairs
        pairIdx = pairIdxs(iPair, :);
        if bar_pCorr(iGroup, iPair) < 0.05
            sigstar({barPositions(iGroup, pairIdx)}, bar_pCorr(iGroup, iPair));
        end
    end

    if iGroup == 1
        ylabel(ylabelStr, 'interpreter', 'none');
    end    
    if iGroup > 1
        hnt(iGroup).YAxis.TickValues = [];
    end
    set(gca,'TickLabelInterpreter','none');
end
linkaxes(hnt, 'y');

if ~isempty(plotTitle)
    title(htl, plotTitle, 'interpreter', 'none')
end

end