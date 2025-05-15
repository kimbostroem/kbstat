% Set random seed
rng(123);

% Meta parameters
nSubjects = 30;
nA = 2;
nB = 3;
nTrials = 3;
shape = 2;

% Define mean reaction times for each A-B condition
% Rows = A level, Columns = B level
mu = [300, 400, 500;  % A = 1
      350, 450, 600]; % A = 2

% Preallocate arrays
totalRows = nSubjects * nA * nB * nTrials;
id = strings(totalRows, 1);
A = strings(totalRows, 1);
B = strings(totalRows, 1);
trial = zeros(totalRows, 1);
rt = zeros(totalRows, 1);

row = 1;

for iSubject = 1:nSubjects
    mySubject = sprintf("P%02d", iSubject);
    for iA = 1:nA
        for iB = 1:nB
            scale = mu(iA, iB) / shape;
            for iTrial = 1:nTrials
                myRt = gamrnd(shape, scale);
                id(row) = mySubject;
                A(row) = sprintf('%d', iA);
                B(row) = sprintf('%d', iB);
                trial(row) = iTrial;
                rt(row) = myRt;
                row = row + 1;
            end
        end
    end
end

% Log-transform
rt_log = log(rt);

% Create table
DataTable = table(id, A, B, trial, rt, rt_log);

% Save master CSV
writetable(DataTable, 'reaction_time.csv');