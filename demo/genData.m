% Set random seed
rng(123);

% Parameters
n_participants = 30;
levels_A = 2;
levels_B = 3;
n_trials = 3;

shape_param = 2;

% Define mean reaction times for each A-B condition
% Rows = A level, Columns = B level
mu = [300, 400, 500;   % A = 1
       350, 450, 600]; % A = 2

% Preallocate arrays
total_rows = n_participants * levels_A * levels_B * n_trials;
participant_id = strings(total_rows, 1);
A = zeros(total_rows, 1);
B = zeros(total_rows, 1);
trial = zeros(total_rows, 1);
reaction_time = zeros(total_rows, 1);

row = 1;

for p = 1:n_participants
    pid = sprintf("P%02d", p);
    for a = 1:levels_A
        for b = 1:levels_B
            scale_param = mu(a, b) / shape_param;
            for t = 1:n_trials
                rt = gamrnd(shape_param, scale_param);

                participant_id(row) = pid;
                A(row) = a;
                B(row) = b;
                trial(row) = t;
                reaction_time(row) = rt;

                row = row + 1;
            end
        end
    end
end

% Log-transform
log_rt = log(reaction_time);

% Create table
T = table(participant_id, A, B, trial, reaction_time, log_rt);

% Save master CSV
writetable(T, 'reaction_time.csv');

% % Optional: save separate files per condition (A & B)
% for a = 1:levels_A
%     for b = 1:levels_B
%         T_ab = T(T.A == a & T.B == b, :);
%         filename = sprintf('reaction_times_A%d_B%d.csv', a, b);
%         writetable(T_ab, filename);
%     end
% end
% 
% % Optional: visualize histogram
% figure;
% histogram(reaction_time, 30);
% xlabel('Reaction Time (ms)');
% ylabel('Frequency');
% title('Simulated Reaction Times: 2×3×3 Design');
