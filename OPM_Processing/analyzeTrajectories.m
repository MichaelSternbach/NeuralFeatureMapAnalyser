function [pairing, mean_distance, directions_set1, directions_set2] = analyzeTrajectories(set1, set2, min_length)
% ANALYZETRAJECTORIES Computes pairwise distances between trajectories in two sets,
% creates a checkerboard plot, finds an optimal pairing, visualizes results, 
% and calculates trajectory direction distributions.
% 
% INPUT:
%   set1 - Cell array of trajectories (Nx3 matrices).
%   set2 - Cell array of trajectories (Mx3 matrices).
%   min_length - Minimum trajectory length to include in the analysis.
% 
% OUTPUT:
%   pairing - Optimal pairing between trajectories in set1 and set2.
%   mean_distance - Mean distance of the optimal pairing.
%   directions_set1 - Distribution of trajectory directions (angles) for set1.
%   directions_set2 - Distribution of trajectory directions (angles) for set2.

% Step 1: Filter trajectories based on length
set1 = set1(cellfun(@(x) size(x, 1) >= min_length, set1));
set2 = set2(cellfun(@(x) size(x, 1) >= min_length, set2));

% Ensure non-empty sets after filtering
if isempty(set1) || isempty(set2)
    error('One or both trajectory sets are empty after filtering by length.');
end

% Step 2: Plot all trajectories in 3D color-coded by set
figure;
hold on;
for i = 1:numel(set1)
    traj = set1{i};
    plot3(traj(:, 1), traj(:, 2), traj(:, 3), 'b-', 'LineWidth', 1.5);
end
for j = 1:numel(set2)
    traj = set2{j};
    plot3(traj(:, 1), traj(:, 2), traj(:, 3), 'r--', 'LineWidth', 1.5);
end
hold off;
title('3D Visualization of All Trajectories');
xlabel('X'); ylabel('Y'); zlabel('Z');
grid on;
legend('Set 1', 'Set 2');

% Step 3: Compute pairwise distances
num1 = numel(set1);
num2 = numel(set2);
distance_matrix = zeros(num1, num2);

for i = 1:num1
    for j = 1:num2
        distance_matrix(i, j) = computeMeanDistance(set1{i}, set2{j});
    end
end

% Step 4: Create a checkerboard plot
figure;
imagesc(distance_matrix);
colorbar;
title('Checkerboard Plot of Pairwise Distances');
xlabel('Set 2 Trajectories');
ylabel('Set 1 Trajectories');

% Step 5: Find optimal pairing
if num1 <= num2
    [pairing, mean_distance] = findOptimalPairing(distance_matrix, num1);
else
    [pairing, mean_distance] = findOptimalPairing(distance_matrix', num2);
    pairing = pairing(:, [2, 1]); % Reverse order for consistency
end

% Step 6: Visualize paired trajectories in 3D
figure;
hold on;
for k = 1:size(pairing, 1)
    traj1 = set1{pairing(k, 1)};
    traj2 = set2{pairing(k, 2)};
    
    % Plot trajectory 1
    plot3(traj1(:, 1), traj1(:, 2), traj1(:, 3), 'b-', 'LineWidth', 1.5);
    
    % Plot trajectory 2
    plot3(traj2(:, 1), traj2(:, 2), traj2(:, 3), 'r--', 'LineWidth', 1.5);
    
    % Draw a line between the first points of the pair
    plot3([traj1(1, 1), traj2(1, 1)], [traj1(1, 2), traj2(1, 2)], [traj1(1, 3), traj2(1, 3)], 'k:', 'LineWidth', 1.2);
end
hold off;
title('3D Visualization of Paired Trajectories');
xlabel('X'); ylabel('Y'); zlabel('Z');
grid on;
legend('Set 1 Trajectory', 'Set 2 Trajectory', 'Pair Connection');

% Step 7: Calculate trajectory direction distributions
[directions_set1, angles_set1] = calculateTrajectoryDirections(set1);
[directions_set2, angles_set2] = calculateTrajectoryDirections(set2);

% Plot distributions
figure;
subplot(1, 2, 1);
histogram(angles_set1, 'Normalization', 'probability');
title('Trajectory Directions (Set 1)');
xlabel('Angle (rad)'); ylabel('Probability');

subplot(1, 2, 2);
histogram(angles_set2, 'Normalization', 'probability');
title('Trajectory Directions (Set 2)');
xlabel('Angle (rad)'); ylabel('Probability');

end

function mean_distance = computeMeanDistance(traj1, traj2)
% COMPUTEMEANDISTANCE Computes the mean Euclidean distance between two trajectories.
len = min(size(traj1, 1), size(traj2, 1));
traj1 = traj1(1:len, :);
traj2 = traj2(1:len, :);

mean_distance = mean(sqrt(sum((traj1 - traj2).^2, 2)));
end

function [pairing, mean_distance] = findOptimalPairing(distance_matrix, num_pairs)
% FINDOPTIMALPAIRING Finds the pairing that minimizes mean distance.
% Uses a greedy algorithm to pair the smaller set with the larger set.

[~, idx] = sort(distance_matrix(:));
[row, col] = ind2sub(size(distance_matrix), idx);

pairing = zeros(num_pairs, 2);
selected_rows = false(size(distance_matrix, 1), 1);
selected_cols = false(size(distance_matrix, 2), 1);

count = 0;
total_distance = 0;
for k = 1:numel(row)
    if count == num_pairs
        break;
    end
    if ~selected_rows(row(k)) && ~selected_cols(col(k))
        count = count + 1;
        pairing(count, :) = [row(k), col(k)];
        total_distance = total_distance + distance_matrix(row(k), col(k));
        selected_rows(row(k)) = true;
        selected_cols(col(k)) = true;
    end
end

mean_distance = total_distance / num_pairs;
end

function [directions, angles] = calculateTrajectoryDirections(set)
% CALCULATETRAJECTORYDIRECTIONS Calculates the direction vectors and angles of trajectories.
% INPUT:
%   set - Cell array of trajectories (Nx3 matrices).
% OUTPUT:
%   directions - Nx3 matrix of mean direction vectors for each trajectory.
%   angles - Vector of angles (in radians) of the direction vectors.

directions = zeros(numel(set), 3);
angles = zeros(numel(set), 1);

for i = 1:numel(set)
    traj = set{i};
    mean_vector = traj(end, :) - traj(1, :);
    directions(i, :) = mean_vector / norm(mean_vector); % Normalize the direction vector
    angles(i) = atan2(norm(cross([1, 0, 0], directions(i, :))), dot([1, 0, 0], directions(i, :))); % Angle w.r.t. x-axis
end

end
