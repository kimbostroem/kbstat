paramArray(1).bottomWhisker = 1.5;
paramArray(1).bottomBox = 2.5;
paramArray(1).bottomNotch = 3;
paramArray(1).median = 3.5;
paramArray(1).topNotch = 4;
paramArray(1).topBox = 4.5;
paramArray(1).topWhisker = 5.2;
paramArray(1).outliers = [];

paramArray(2).bottomWhisker = 1.8;
paramArray(2).bottomBox = 2.7;
paramArray(2).bottomNotch = 3.2;
paramArray(2).median = 3.6;
paramArray(2).topNotch = 4;
paramArray(2).topBox = 4.4;
paramArray(2).topWhisker = 5;
paramArray(2).outliers = [1.9, 5.5];

paramArray(3).bottomWhisker = 2.2;
paramArray(3).bottomBox = 2.9;
paramArray(3).bottomNotch = 3.5;
paramArray(3).median = 3.7;
paramArray(3).topNotch = 4.2;
paramArray(3).topBox = 4.8;
paramArray(3).topWhisker = 5.5;
paramArray(3).outliers = [2.0, 6.5, 6.1];

annotations(1).x1 = 2;
annotations(1).x2 = 3;
annotations(1).label = '**';

figure;
% kbboxplot(paramArray, ...
%     'XPositions', [1 2], ...
%     'BoxWidth', 0.5, ...
%     'NotchAlpha', 0.5, ...
%     'Colors', lines, ...
%     'Annotations', annotations);

kbboxplot(paramArray, 'Annotations', annotations);

xticks([1 2]);
xticklabels({'Control', 'Treatment'});
ylabel('Value');
title('Custom Boxplots');
