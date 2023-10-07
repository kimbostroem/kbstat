function [barPositions, ylimits] = plotGroups(values, members, groups, memberName, groupName, bar_pCorr, plotTitle, ylabelStr, plotStyle, parent, showVarNames)

if nargin < 11
    showVarNames = 1;
end
if nargin < 10
    parent = [];
end

if nargin < 9
    plotStyle = 'violin';
end

switch showVarNames
    case {1, 2} % display variable names and levels
    members = strcat(memberName, {' = '}, cellstr(members));
    groups = strcat(groupName, {' = '}, cellstr(groups));
end

members = string(members);
groups = string(groups);
nMembers = length(members);
nGroups = length(groups);

if ~isempty(parent)
    htl = tiledlayout(parent, 1,2, 'TileSpacing','none');
else
    htl = tiledlayout(1,2, 'TileSpacing','none');
end
barPositions = [1:nGroups; 1:nMembers];
pairIdxs = nchoosek(1:nMembers, 2);
nPairs = size(pairIdxs, 1);
hnt = gobjects(nGroups, 1);
for iGroup = 1:nGroups
    hnt(iGroup) = nexttile(htl, iGroup);
    switch lower(plotStyle)
        case 'violin'
            violinplot(squeeze(values(iGroup,:,:))', [], 'MedianMarkerSize', 48, 'BoxColor', 0.2*[1 1 1], 'MarkerSize', 24);
        case 'boxplot'
            boxplot(squeeze(values(iGroup,:,:))', 'Colors', lines(nMembers));
        case 'bar'
            barValues = squeeze(quantile(values(iGroup,:,:), 0.5, 3));
            errorBottom = squeeze(quantile(values(iGroup,:,:), 0.25, 3));
            errorTop = squeeze(quantile(values(iGroup,:,:), 0.75, 3));
            colors = lines(nMembers);
            hbar = bar(1:nMembers, barValues, 'LineStyle', 'none', 'FaceColor', 'flat');
            for iMember = 1:nMembers
                hbar.CData(iMember,:) = colors(iMember,:);
                hold on
                errorbar(iMember, barValues(iMember), barValues(iMember)-errorBottom(iMember), errorTop(iMember)-barValues(iMember), 'LineStyle', 'none', 'Color', 0.6*[1 1 1], 'LineWidth', 2);
            end
    end   

    title(groups(iGroup))
    xticklabels(members);

    for iPair = 1:nPairs
        pairIdx = pairIdxs(iPair, :);
        if bar_pCorr(iGroup, iPair) < 0.05
            sigstar({barPositions(iGroup, pairIdx)}, bar_pCorr(iGroup, iPair));
        end
    end
    if iGroup==1; ylabel(ylabelStr); end
end
linkaxes(hnt)
hnt(2).YAxis.TickValues = [];
if ~isempty(plotTitle)
    title(htl, plotTitle)
end
hold on;

set(gca, 'YGrid', 'on', 'XGrid', 'off')
set(gca,'TickLabelInterpreter','none');

ylimits = ylim;

end