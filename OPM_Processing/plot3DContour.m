function plot3DContour(data, targetValue, roiMask)
    % Function to plot a 3D scatter plot of voxels within a specified
    % tolerance (epsilon) of a target value, considering 2*pi circular wrapping,
    % restricted by a 3D boolean ROI.
    %
    % Inputs:
    %   data        - A 3D array (e.g., matrix of size [M, N, P])
    %   targetValue - The value to match (scalar) in the range [0, 2*pi]
    %   epsilon     - The tolerance for matching (scalar)
    %   roiMask     - A 3D boolean array of the same size as 'data' specifying
    %                 the region of interest (ROI), where true values define the ROI
    epsilon = 0.01;
    % Check if the ROI mask is the same size as the data
    if ~isequal(size(data), size(roiMask))
        error('The size of roiMask must match the size of data.');
    end

    % Normalize data values to the range [0, 2*pi]
    data = mod(data, 2*pi);
    targetValue = mod(targetValue, 2*pi);

    % Find the indices of voxels that match the target value within epsilon,
    % accounting for circular wrapping at 0 and 2*pi
    matchingVoxels = (abs(data - targetValue) <= epsilon | ...
                     abs(data - (targetValue + 2*pi)) <= epsilon | ...
                     abs(data - (targetValue - 2*pi)) <= epsilon) & roiMask;

    % Get the (x, y, z) coordinates of the matching voxels
    [x, y, z] = ind2sub(size(data), find(matchingVoxels));

    % Get the corresponding data values for coloring
    values = data(matchingVoxels);

    % Create a 3D scatter plot of the matching voxels
    figure;
    scatter3(x, y, z, 36, values, 'filled'); % 36 is the marker size
    colorbar;
    xlabel('X');
    ylabel('Y');
    zlabel('Z');
    title(['3D Plot of Circular Values Matching: ', num2str(targetValue), ...
           ' \pm ', num2str(epsilon)]);
    grid on;
    axis equal;
end

