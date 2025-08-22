function kbsaveFigure(fig, filePath, fileTypes, figWidth, figHeight)

[fdir, fname, ~] = fileparts(filePath);
filePath = fullfile(fdir, fname);
fileTypes = cellstr(fileTypes); % ensure cell array
figPos = get(fig, 'Position');
if nargin < 4
    figWidth = figPos(3);
    figHeight = figPos(4);
end

fig.Position(3:4) = [figWidth, figHeight];
set(fig, 'PaperPositionMode', 'auto');
set(fig, 'PaperUnits', 'points');
set(fig, 'PaperSize', [figWidth figHeight]);
set(fig, 'renderer', 'painters');
printResolution = 300; % 72 = monitor resolution; 300 = printer resolution

for i = 1:length(fileTypes)
    fileType = fileTypes{i};
    switch lower(fileType)
        case {'fig', '.fig'}
            savefig(fig, sprintf('%s.fig', filePath));
        case {'png', '.png'}
            print(fig, sprintf('%s', filePath), '-dpng', sprintf('-r%.0f', printResolution));
        case {'pdf', '.pdf'}
            print(fig, sprintf('%s', filePath), '-dpdf', sprintf('-r%.0f', printResolution), '-bestfit');
        case {'eps', '.eps'}
            print(fig, sprintf('%s', filePath), '-deps', sprintf('-r%.0f', printResolution));
    end
end

end