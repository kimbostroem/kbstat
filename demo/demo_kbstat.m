options = struct;
options.inFile = 'Chocolate_LM4.xlsx';
options.y = 'Distance';
options.yUnits = 'm';
options.x = 'Chocolate, Gender';
options.id = 'Subject';
options.within = 'Chocolate';
options.interact = 'Chocolate, Gender';
options.distribution = 'gamma';
options.posthocMethod = 'emm';
options.removeOutliers = 'true';
kbstat(options);