function [new_array, new_ROI] = resample_to_uniform_density(input_array, pixel_size, ROI, target_pixel_size, interp_method)
    % Function to resample a 2D, 3D, or 4D array with varying pixel densities to a uniform pixel density
    % Additionally, resamples a 3D or 2D ROI (Region of Interest) if provided, without interpolation (binary).
    %
    % input_array: 2D, 3D or 4D array with dimensions [X, Y] or [X, Y, Z] or [X, Y, Z, T]
    % pixel_size: 1x2 or 1x3 vector [dx, dy] or [dx, dy, dz] specifying the size of each pixel in each spatial dimension
    % ROI: (optional) 2D or 3D Region of Interest with dimensions [X, Y] or [X, Y, Z]. Default is a 2D or 3D array of ones.
    % target_pixel_size: (optional) scalar specifying the desired uniform pixel size for all spatial dimensions.
    %                    
    % interp_method: (optional) interpolation method (e.g., 'linear', 'nearest', 'cubic'). Default is 'linear'.
    %
    % Returns:
    % new_array: Resampled input_array
    % new_ROI: Resampled ROI (if provided), otherwise returns empty.
    
    % Set default for ROI if not specified (full ROI of ones)
    if nargin < 3 || isempty(ROI)
        if length(size(input_array)) == 3
            ROI = ones(size(input_array, 1), size(input_array, 2), size(input_array, 3));
        else
            ROI = ones(size(input_array, 1), size(input_array, 2));
        end
    end
    
    % Check if input array is 2D, 3D, or 4D
    input_dims = ndims(input_array);
    
    if input_dims == 2
        % If 2D, add a singleton dimension for consistency
        input_array = reshape(input_array, size(input_array, 1), size(input_array, 2), 1, 1);
        is_2D_input = true;
    elseif input_dims == 3
        is_2D_input = false;
    elseif input_dims == 4
        is_2D_input = false;
    else
        error('Input array must be either 2D, 3D, or 4D.');
    end
    
    % Check if ROI matches dimensions [X, Y, Z] or [X, Y]
    if input_dims == 3 || input_dims == 4
        if ~isequal(size(ROI), size(input_array(:,:,:,1)))
            error('ROI must have dimensions [X, Y, Z] matching input array.');
        end
    elseif input_dims == 2
        if ~isequal(size(ROI), size(input_array(:,:,1)))
            error('ROI must have dimensions [X, Y] matching input array.');
        end
    end
    
    % Set default for target_pixel_size if not specified
    if nargin < 4 || isempty(target_pixel_size)
        % Calculate target pixel size as the greatest common divisor of pixel sizes
        target_pixel_size = gcd_multiple(pixel_size);
    end
    
    % Set default for interpolation method if not specified
    if nargin < 5
        interp_method = 'linear';  % Default interpolation method
    end
    
    % Determine if the third dimension corresponds to Z (depth) or T (trials)
    if length(pixel_size) == 3
        % 3D spatial data: the third dimension is Z (depth)
        is_3D_spatial = true;
    else
        % 2D spatial data: the third dimension is T (trials)
        is_3D_spatial = false;
    end
    
    [X, Y, Z, T] = size(input_array);
    
    % Set target sizes based on the desired pixel size
    if is_3D_spatial
        target_size_x = round(X * (pixel_size(1) / target_pixel_size));
        target_size_y = round(Y * (pixel_size(2) / target_pixel_size));
        target_size_z = round(Z * (pixel_size(3) / target_pixel_size));
    else
        target_size_x = round(X * (pixel_size(1) / target_pixel_size));
        target_size_y = round(Y * (pixel_size(2) / target_pixel_size));
        target_size_z = 1; % No Z dimension in 2D input
    end
    
    % Create original coordinate grids for resampling
    orig_grid_x = linspace(1, X, X) * pixel_size(1);
    orig_grid_y = linspace(1, Y, Y) * pixel_size(2);
    
    if is_3D_spatial
        orig_grid_z = linspace(1, Z, Z) * pixel_size(3);
    end
    
    % Create new coordinate grids for uniform sampling
    new_grid_x = linspace(orig_grid_x(1), orig_grid_x(end), target_size_x);
    new_grid_y = linspace(orig_grid_y(1), orig_grid_y(end), target_size_y);
    
    if is_3D_spatial
        new_grid_z = linspace(orig_grid_z(1), orig_grid_z(end), target_size_z);
    end
    
    % Initialize new array
    if is_3D_spatial
        new_array = zeros(target_size_x, target_size_y, target_size_z, T);
    else
        new_array = zeros(target_size_x, target_size_y, T);
    end
    
    % Perform interpolation for each trial (or single volume in 3D/4D case)
    for t = 1:T
        % Extract the data for the t-th trial or volume
        data = input_array(:,:,:,t);
        
        % Use 'interpn' to perform multi-dimensional interpolation
        if is_3D_spatial
            [Xg, Yg, Zg] = ndgrid(orig_grid_x, orig_grid_y, orig_grid_z);
            [new_Xg, new_Yg, new_Zg] = ndgrid(new_grid_x, new_grid_y, new_grid_z);
            
            % Interpolate data for the current trial/volume with specified interpolation method
            new_array(:,:,:,t) = interpn(Xg, Yg, Zg, data, new_Xg, new_Yg, new_Zg, interp_method);
        else
            [Xg, Yg] = ndgrid(orig_grid_x, orig_grid_y);
            [new_Xg, new_Yg] = ndgrid(new_grid_x, new_grid_y);
            
            % Interpolate data for the current trial/volume with specified interpolation method
            new_array(:,:,t) = interpn(Xg, Yg, data, new_Xg, new_Yg, interp_method);
        end
    end
    
    % If the input was 2D, remove the singleton dimension in the output
    if is_2D_input
        new_array = squeeze(new_array);
    end
    
    % Resample the ROI similarly to the data (using nearest-neighbor for binary data)
    if is_3D_spatial
        new_ROI = round(interpn(Xg, Yg, Zg, double(ROI), new_Xg, new_Yg, new_Zg, 'nearest'));
    else
        new_ROI = round(interpn(Xg, Yg, double(ROI), new_Xg, new_Yg, 'nearest'));
    end
end

function gcd_value = gcd_multiple(values)
    % Calculate the greatest common divisor (GCD) for an array of values
    gcd_value = values(1);
    for i = 2:numel(values)
        gcd_value = gcd(gcd_value, values(i));
    end
end
