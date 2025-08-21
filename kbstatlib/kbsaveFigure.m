function kbsaveFigure(fig, filePath, fileTypes, width, height)

[fdir, fname, ~] = fileparts(filePath);
filePath = fullfile(fdir, fname);
fileTypes = cellstr(fileTypes); % ensure cell array
figPos = get(fig, 'Position');
if nargin < 4
    width = figPos(3);
    height = figPos(4);
end

% scrsz = get(0, 'ScreenSize');
set(fig, 'Position', [1, 1, width, height]);
% set(fig, 'Position', [figPos(1), figPos(2), width, height]);
set(fig, 'PaperPositionMode', 'auto');
set(fig, 'PaperUnits', 'points');
set(fig, 'PaperSize', [width height]);
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