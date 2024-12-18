function [optimal_shift, min_mean_distance, aligned_traj1] = alignTrajectories(traj1, traj2, max_initial_shifts)
% ALIGNTRAJECTORIES Finds the optimal alignment of two 3D trajectories
% by minimizing the mean point-wise Euclidean distance through iterative shifts.
% 
% INPUT:
%   traj1 - Nx3 array representing the first trajectory (3D points).
%   traj2 - Mx3 array representing the second trajectory (3D points).
%   max_initial_shifts - Maximum number of initial shifts to test in each direction (default: 4).
% 
% OUTPUT:
%   optimal_shift - The optimal shift applied to align traj1 with traj2.
%   min_mean_distance - The smallest mean Euclidean distance between aligned trajectories.
%   aligned_traj1 - The aligned version of traj1.

if nargin < 3
    max_initial_shifts = 4; % Default number of initial shifts
end

% Initialize parameters
shift_step = 1; % Initial step size for shifts
current_shift = initialShift(traj1, traj2);
direction = 1; % 1 for forward, -1 for backward

min_mean_distance = inf;
optimal_shift = current_shift;
iteration = 0;

while true
    iteration = iteration + 1;
    shifted_distances = zeros(max_initial_shifts, 1);

    % Test shifts in the current direction
    for shift_idx = 1:max_initial_shifts
        shift_amount = current_shift + direction * (shift_idx * shift_step);

        % Apply the shift and calculate distances
        shifted_traj1 = applyShift(traj1, shift_amount);
        distances = computeMeanDistances(shifted_traj1, traj2);
        shifted_distances(shift_idx) = distances;
    end

    % Find the minimum distance for this direction
    [local_min_distance, local_min_idx] = min(shifted_distances);
    local_min_shift = current_shift + direction * (local_min_idx * shift_step);

    if local_min_distance < min_mean_distance
        % Update optimal values if improvement is found
        min_mean_distance = local_min_distance;
        optimal_shift = local_min_shift;
        current_shift = local_min_shift;
    elseif all(diff(shifted_distances) > 0)
        % If distances increase in this direction, reverse direction
        direction = -direction;
    elseif all(diff(shifted_distances) < 0)
        % If distances decrease with maximum tested shifts, increase range
        max_initial_shifts = max_initial_shifts + 1;
    else
        % If no improvement in either direction, stop
        break;
    end
end

% Apply the optimal shift
aligned_traj1 = applyShift(traj1, optimal_shift);
end

function shifted_traj = applyShift(traj, shift_amount)
% Helper function to apply a temporal shift to a trajectory.
% Positive shift truncates the beginning, negative shift truncates the end.
if shift_amount > 0
    shifted_traj = traj(1+shift_amount:end, :);
elseif shift_amount < 0
    shifted_traj = traj(1:end+shift_amount, :);
else
    shifted_traj = traj;
end
end

function mean_distance = computeMeanDistances(traj1, traj2)
% Helper function to compute the mean Euclidean distance between two trajectories.
len = min(size(traj1, 1), size(traj2, 1));
traj1 = traj1(1:len, :);
traj2 = traj2(1:len, :);

mean_distance = mean(sqrt(sum((traj1 - traj2).^2, 2)));
end

function initial_shift = initialShift(traj1, traj2)
% Helper function to compute the initial shift that aligns the first point of traj1
% with the closest point in traj2.
first_point = traj1(1, :);
distances = sqrt(sum((traj2 - first_point).^2, 2));
[~, closest_idx] = min(distances);
initial_shift = closest_idx - 1;
end
