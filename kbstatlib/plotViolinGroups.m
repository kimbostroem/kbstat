function [barPositions, ylimits] = plotViolinGroups(values, members, groups, memberName, groupName, bar_pCorr, plotTitle, ylabelStr, layout)

if nargin < 9
    layout = [];
end

members = strcat(memberName, {' = '}, char(members));
groups = strcat(groupName, {' = '}, char(groups));
members = string(members);
groups = string(groups);
nMembers = length(members);
nGroups = length(groups);

if ~isempty(layout)
    htl = tiledlayout(layout, 1,2, 'TileSpacing','none');
else
    htl = tiledlayout(1,2, 'TileSpacing','none');
end
barPositions = [1:3; 1:3];
pairIdxs = nchoosek(1:nMembers, 2);
nPairs = size(pairIdxs, 1);
hnt = gobjects(nGroups, 1);
for iGroup = 1:nGroups
    hnt(iGroup) = nexttile(htl, iGroup);
    violinplot(squeeze(values(iGroup,:,:))', [], 'MedianMarkerSize', 48, 'BoxColor', 0.2*[1 1 1]);
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
title(htl, plotTitle)

hold on;

set(gca, 'YGrid', 'on', 'XGrid', 'off')
set(gca,'TickLabelInterpreter','none');

ylimits = ylim;

end