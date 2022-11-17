function [barPositions, ylimits] = plotBarGroups(values, members, groups, errorBottom, errorTop, formatSpec)

if nargin < 6
    formatSpec = '%.3f';
end
if nargin < 5
    errorTop = [];
end

if nargin < 4
    errorBottom = [];
end

if ~iscategorical(members)
    members = reordercats(categorical(members), members);
end
if ~iscategorical(groups)
    groups = reordercats(categorical(groups), groups);
end
nMembers = length(members);
nGroups = length(groups);
barPosAmp = 0.3682626 - 0.3298725 * exp(-0.407004 * (nMembers-1)); % position amplitude
barPosInc = 2 * barPosAmp / (nMembers-1);
bar(1:nGroups, values, 'LineStyle','none');
hold on;
barPositions = NaN(nGroups, nMembers);
for iGroup = 1:nGroups
    for iMember = 1:nMembers
        barPos = iGroup;
        if nMembers~=1
            barPos = barPos - barPosAmp + (iMember-1) * barPosInc;
        end
        barPositions(iGroup, iMember) = barPos;
        if ~isempty(errorBottom) && ~isempty(errorTop)
            errorbar(barPos, values(iGroup, iMember), values(iGroup, iMember) - errorBottom(iGroup, iMember), errorTop(iGroup, iMember) - values(iGroup, iMember), 'LineStyle', 'none', 'Color', 0.6*[1 1 1], 'LineWidth', 2);
        elseif ~isempty(errorBottom)
            errorbar(barPos, values(iGroup, iMember), values(iGroup, iMember) - errorBottom(iGroup, iMember), 'LineStyle', 'none');
        end
        if ~isempty(formatSpec) && ~strcmp(formatSpec, 'none')
            text(barPos, values(iGroup, iMember), num2str(values(iGroup, iMember), formatSpec), 'vert', 'bottom', 'horiz', 'center');
        end
    end    
end
if nGroups > 1
    xticklabels(groups);
    legend(members, 'Location', 'best', 'Interpreter', 'none');
else
    xticks(barPositions(1,:));
    xticklabels(members);
end
set(gca, 'YGrid', 'on', 'XGrid', 'off')
set(gca,'TickLabelInterpreter','none');

ylimits = ylim;


end