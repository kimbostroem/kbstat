%% Run all kbstat demo scripts and report pass/fail
% Executes every demo_*.m in this folder in turn, catching errors so one failure
% does not stop the rest. Each demo writes its output under demo/results/.

thisDir = fileparts(mfilename('fullpath'));
cd(thisDir);

demos = dir(fullfile(thisDir, 'demo_*.m'));
fprintf('=== kbstat demo runner ===\n\n');

failed = {};
for iDemo = 1:numel(demos)
    name = demos(iDemo).name;
    fprintf('  %-40s', name);
    try
        run(fullfile(thisDir, name));
        close all force;
        fprintf('OK\n');
    catch ME
        fprintf('FAILED (%s)\n', ME.message);
        failed{end+1} = name; %#ok<AGROW>
    end
    cd(thisDir); % demos cd into their own folder; restore for the next one
end

fprintf('\n');
if isempty(failed)
    fprintf('=== All %d demos passed ===\n', numel(demos));
else
    fprintf('=== %d of %d demo(s) FAILED ===\n', numel(failed), numel(demos));
    for iFail = 1:numel(failed)
        fprintf('  %s\n', failed{iFail});
    end
end
