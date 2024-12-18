function modified_set = switchZandX(set)
% SWITCHZANDX Switches the Z and X coordinates for all trajectories in a set.
% INPUT:
%   set - Cell array of trajectories (Nx3 matrices).
% OUTPUT:
%   modified_set - Cell array with Z and X coordinates swapped for each trajectory.

modified_set = cellfun(@(traj) traj(:, [3, 2, 1]), set, 'UniformOutput', false);
end