function [barPositions, ylimits] = plotGroups(values, members, groups, memberName, groupName, bar_pCorr, plotTitle, ylabelStr, plotStyle, parent, showVarNames, markerSize, emms, errorBars)

if nargin < 14
    errorBars = 'iqr';
end

if nargin < 13
    emms = [];
end

if nargin < 12 || isnan(markerSize)
    markerSize = 24;
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

switch showVarNames
    case {1, 3, 'names_and_levels', 'Names_and_Levels'} % display variable names and levels
    members = strcat(memberName, {' = '}, cellstr(members));
    groups = strcat(groupName, {' = '}, cellstr(groups));
end

members = string(members);
groups = string(groups);
nMembers = length(members);
nGroups = length(groups);

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
    switch errorBars
        case 'iqr'
            barValues = squeeze(quantile(values(iGroup,:,:), 0.5, 3));
            errorBottom = squeeze(quantile(values(iGroup,:,:), 0.25, 3));
            errorTop = squeeze(quantile(values(iGroup,:,:), 0.75, 3));
        case 'std'
            barValues = mean(squeeze(values(iGroup,:,:)), 2, 'omitnan');
            errorBottom = barValues - std(squeeze(values(iGroup,:,:)), 0, 2, 'omitnan');
            errorTop = barValues + std(squeeze(values(iGroup,:,:)), 0, 2, 'omitnan');
        case 'se'
            barValues = mean(squeeze(values(iGroup,:,:)), 2, 'omitnan');
            errorBottom = barValues - std(squeeze(values(iGroup,:,:)), 0, 2, 'omitnan') / size(squeeze(values(iGroup,:,:)), 2);
            errorTop = barValues + std(squeeze(values(iGroup,:,:)), 0, 2, 'omitnan') / size(squeeze(values(iGroup,:,:)), 2);
        case 'ci95'
            barValues = mean(squeeze(values(iGroup,:,:)), 2, 'omitnan');
            errorBottom = NaN(nMembers, 1);
            errorTop = NaN(nMembers, 1);
            for iCol = 1:length(barValues)
                pd = fitdist(squeeze(values(iGroup,iCol,:)),'Normal');
                ci = paramci(pd);
                errorBottom(iCol) = ci(1,1);
                errorTop(iCol) = ci(2,1);
            end
    end
    switch lower(plotStyle)
        case 'violin'
            violinplot(squeeze(values(iGroup,:,:))', [], 'MedianMarkerSize', 48, 'BoxColor', 0.2*[1 1 1], 'MarkerSize', markerSize);
        case 'boxplot'
            boxplot(squeeze(values(iGroup,:,:))', 'Colors', lines(nMembers));
        case 'bar'
            colors = lines(nMembers);
            hbar = bar(1:nMembers, barValues, 'LineStyle', 'none', 'FaceColor', 'flat');
            for iMember = 1:nMembers
                hbar.CData(iMember,:) = colors(iMember,:);
                hold on
                errorbar(iMember, barValues(iMember), barValues(iMember)-errorBottom(iMember), errorTop(iMember)-barValues(iMember), 'LineStyle', 'none', 'Color', 0.6*[1 1 1], 'LineWidth', 2);
            end
    end

    % plot horizontal lines at values
    if size(emms, 1) == nGroups && size(emms, 2) == nMembers
        colors = lines;
        for iMember = 1:nMembers
            yline(emms(iGroup, iMember), '-', 'Color', colors(iMember,:));
        end
    end

    if ~isempty(groupName)
        title(groups(iGroup), 'interpreter', 'none');
    end
    xticklabels(members);

    for iPair = 1:nPairs
        pairIdx = pairIdxs(iPair, :);
        if bar_pCorr(iGroup, iPair) < 0.05
            sigstar({barPositions(iGroup, pairIdx)}, bar_pCorr(iGroup, iPair));
        end
    end

    if iGroup == 1
        ylabel(ylabelStr, 'interpreter', 'none'); 
    end    
    linkaxes(hnt);
    if iGroup > 1
        hnt(iGroup).YAxis.TickValues = [];
    end
    set(gca,'TickLabelInterpreter','none');
end

if ~isempty(plotTitle)
    title(htl, plotTitle, 'interpreter', 'none')
end

ylimits = ylim;

end