function kbboxplot(paramArray, varargin)
% KBBOXPLOT - Draws multiple stylized boxplots with explicitly defined parameters,
%             mimicking the style of MATLAB's boxchart.
%
% Inputs:
%   paramArray - struct array, each with fields:
%       bottomWhisker, bottomBox, bottomNotch, median,
%       topNotch, topBox, topWhisker, outliers (optional)
%
% Optional Name-Value Parameters:
%   'XPositions'    - x-axis positions of the boxes
%   'BoxWidth'      - width of the box (default: 0.5)
%   'NotchAlpha'    - notch fill transparency (default: 0.7)
%   'Colors'        - Nx3 matrix of RGB values
%   'Annotations'   - a struct array with fields:
%       .x1         - index or x-position of the first box
%       .x2         - index or x-position of the second box
%       .y          - y-position for the top of the bracket (optional)
%       .label      - significance label (e.g., '*', '**', 'ns')


% --- Parse input ---
N = numel(paramArray);
p = inputParser;
addParameter(p, 'XPositions', 1:N);
addParameter(p, 'BoxWidth', 0.5);
addParameter(p, 'NotchAlpha', 0.5);
addParameter(p, 'Colors', lines);
addParameter(p, 'Annotations', struct);
parse(p, varargin{:});

xPositions = p.Results.XPositions;
boxWidth = p.Results.BoxWidth;
notchAlpha = p.Results.NotchAlpha;
colors = p.Results.Colors;

hold on;

for i = 1:N
    pStruct = paramArray(i);
    xPos = xPositions(i);
    color = colors(i, :);

    % Validate fields
    requiredFields = {'bottomWhisker','bottomBox','bottomNotch','median','topNotch','topBox','topWhisker'};
    for f = requiredFields
        if ~isfield(pStruct, f{1})
            error(['Missing field: ', f{1}, ' in paramArray(', num2str(i), ')']);
        end
    end

    % Coordinates for box and notch
    halfBox = boxWidth / 2;
    quarterBox = boxWidth / 4;

    % --- Outline ---
    
    xOutline = [
        xPos - halfBox
        xPos - halfBox
        xPos - quarterBox
        xPos - halfBox
        xPos - halfBox
        xPos + halfBox
        xPos + halfBox
        xPos + quarterBox
        xPos + halfBox
        xPos + halfBox
        ];

    yOutline = [
        pStruct.bottomBox
        pStruct.bottomNotch
        pStruct.median
        pStruct.topNotch
        pStruct.topBox
        pStruct.topBox
        pStruct.topNotch
        pStruct.median
        pStruct.bottomNotch
        pStruct.bottomBox
        ];

    patch(xOutline, yOutline, color, 'EdgeColor', color, 'FaceColor', 'white');

    % --- Notch polygon ---

    xNotch = [
        xPos - halfBox
        xPos - quarterBox
        xPos - halfBox
        xPos + halfBox
        xPos + quarterBox
        xPos + halfBox
        ];

    yNotch = [
        pStruct.bottomNotch
        pStruct.median
        pStruct.topNotch
        pStruct.topNotch
        pStruct.median
        pStruct.bottomNotch
        ];

    patch(xNotch, yNotch, color, 'EdgeColor', 'none', 'FaceAlpha', notchAlpha);

    % --- Median line ---
    line([xPos - quarterBox, xPos + quarterBox], [pStruct.median, pStruct.median], ...
        'Color', color, 'LineWidth', 2);

    % --- Whiskers ---
    line([xPos, xPos], [pStruct.bottomWhisker, pStruct.bottomBox], 'Color', 'k', 'LineWidth', 1);
    line([xPos, xPos], [pStruct.topBox, pStruct.topWhisker], 'Color', 'k', 'LineWidth', 1);

    capWidth = boxWidth * 0.3;
    line([xPos - capWidth/2, xPos + capWidth/2], [pStruct.bottomWhisker, pStruct.bottomWhisker], 'Color', 'k', 'LineWidth', 1);
    line([xPos - capWidth/2, xPos + capWidth/2], [pStruct.topWhisker, pStruct.topWhisker], 'Color', 'k', 'LineWidth', 1);

    % --- Outliers ---
    if isfield(pStruct, 'outliers')
        plot(xPos + zeros(size(pStruct.outliers)), pStruct.outliers, 'ko', ...
            'MarkeredgeColor', 'none', 'MarkerFaceColor', color, 'MarkerSize', 6);
    end
end

% --- Annotations ---
if any(strcmpi(varargin, 'Annotations'))
    A = varargin{find(strcmpi(varargin, 'Annotations')) + 1};

    for j = 1:length(A)
        a = A(j);
        
        % Determine x coordinates (accept index or value)
        if all(a.x1 == floor(a.x1)) && a.x1 <= N
            x1 = xPositions(a.x1);
        else
            x1 = a.x1;
        end
        if all(a.x2 == floor(a.x2)) && a.x2 <= N
            x2 = xPositions(a.x2);
        else
            x2 = a.x2;
        end

        % Determine y-coordinate of the bracket
        if isfield(a, 'y') && ~isempty(a.y)
            y = a.y;
        else
            % Auto: max whisker + margin
            idx1 = find(xPositions == x1, 1);
            idx2 = find(xPositions == x2, 1);
            y1 = paramArray(idx1).topWhisker;
            y2 = paramArray(idx2).topWhisker;
            y = max([y1, y2]) + 0.05 * range(ylim);
        end

        % Bracket height
        h = 0.02 * range(ylim);

        % Draw bracket
        plot([x1, x1, x2, x2], [y, y+h, y+h, y], 'k', 'LineWidth', 1);

        % Label
        if isfield(a, 'label')
            text(mean([x1, x2]), y+h + 0.01 * range(ylim), a.label, ...
                'HorizontalAlignment', 'center', ...
                'VerticalAlignment', 'bottom', ...
                'FontWeight', 'bold');
        end
    end
end


end
