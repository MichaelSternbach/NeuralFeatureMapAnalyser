function map = filter_map3D(map, ROI, pixels_per_mm, highpass_cutoff, lowpass_cutoff, rise,cut_ROI, padd_size)
    % filter_map: Applies Fermi highpass and lowpass filters to a 3D map.
    %
    % Inputs:
    % - map: 3D data array to be filtered.
    % - ROI: binary mask defining the region of interest
    % - cut_ROI (optional): whether to restrict filtering to the ROI (default: true)
    % - padd_size (optional): padding size for each dimension; defaults to the largest dimension
    % - pixels_per_mm: pixels per mm, used to calculate field_size
    % - highpass_cutoff: highpass cutoff frequency (optional)
    % - lowpass_cutoff: lowpass cutoff frequency (optional)
    % - rise: steepness of the Fermi filter transition (scalar)

    if nargin < 7
        cut_ROI = true;
    end

    % Set padd_size to the largest dimension of map if not provided
    if nargin < 8 || isempty(padd_size)
        max_dim = max(size(map));
        padd_size = [max_dim, max_dim, max_dim];
    end

    % Calculate pixels and field_size based on map dimensions and pixels_per_mm
    pixels = size(map);  % Number of pixels in each dimension
    field_size = pixels ./ pixels_per_mm;  % Field size in mm for each dimension

    % Calculate filter parameters using the calculate_filter function
    filter_params = calculate_filter(padd_size, field_size, pixels, highpass_cutoff, lowpass_cutoff, rise);

    % Check if input map is a vector; convert if necessary
    vector_flag = false;
    if isvector(map)
        vector_flag = true;
        if length(map) ~= sum(ROI(:))
            error('Map input is a vector and does not match the size of the defined ROI');
        else
            tmp = zeros(size(ROI));
            tmp(ROI) = map;
            map = tmp;
        end
    end

    % Zero-pad the map and ROI
    map_padd = zeros(padd_size);
    map_padd(1:size(map,1), 1:size(map,2), 1:size(map,3)) = map;
    ROI_padd = zeros(padd_size);
    ROI_padd(1:size(ROI,1), 1:size(ROI,2), 1:size(ROI,3)) = double(ROI);

    % Apply Fermi Highpass Filter if specified
    if ~isempty(filter_params.filter_highpass)
        if cut_ROI
            map_padd(ROI_padd == 0) = 0;
            map_padd = map_padd - ...
                fftshift(ifftn(fftn(fftshift(map_padd)) .* filter_params.filter_highpass)) ...
                ./ fftshift(ifftn(fftn(fftshift(ROI_padd)) .* filter_params.filter_highpass));
            map_padd(ROI_padd == 0) = 0;
        else
            map_padd = map_padd - ...
                fftshift(ifftn(fftn(fftshift(map_padd)) .* filter_params.filter_highpass));
        end
    end

    % Apply Fermi Lowpass Filter if specified
    if ~isempty(filter_params.filter_lowpass)
        if cut_ROI
            map_padd(ROI_padd == 0) = 0;
            map_padd = fftshift(ifftn(fftn(fftshift(map_padd)) .* filter_params.filter_lowpass)) ...
                ./ fftshift(ifftn(fftn(fftshift(ROI_padd)) .* filter_params.filter_lowpass));
            map_padd(ROI_padd == 0) = 0;
        else
            map_padd = fftshift(ifftn(fftn(fftshift(map_padd)) .* filter_params.filter_lowpass));
        end
    end

    % Remove padding and normalize
    map = map_padd(1:size(map,1), 1:size(map,2), 1:size(map,3));
    map = (map - mean(map(ROI))) / std(map(ROI));

    % Return a vector if the input was a vector
    if vector_flag
        map = map(ROI);
    end
end


function filter_params = calculate_filter(padd_size, field_size, pixels, highpass_cutoff, lowpass_cutoff, rise)
    % calculate_filter: Calculates highpass and lowpass 3D Fermi filters
    %
    % Inputs:
    % - padd_size: padding size for each dimension [padd_size_y, padd_size_x, padd_size_z]
    % - field_size: field size in mm for each dimension [field_size_y, field_size_x, field_size_z]
    % - pixels: number of pixels for each dimension [pixels_y, pixels_x, pixels_z]
    % - highpass_cutoff: highpass cutoff frequency (scalar, optional)
    % - lowpass_cutoff: lowpass cutoff frequency (scalar, optional)
    % - rise: steepness of the Fermi filter transition (scalar)
    %
    % Output:
    % - filter_params: struct containing highpass and lowpass 3D Fermi filter kernels

    % Set padd_size to the largest dimension if it is empty
    if isempty(padd_size)
        max_dim = max([field_size(1), field_size(2), field_size(3)]);
        padd_size = [max_dim, max_dim, max_dim];
    end

    % Create a 3D grid for calculating distances from the center
    [xx, yy, zz] = ndgrid(-padd_size(1)/2:padd_size(1)/2-1, ...
                          -padd_size(2)/2:padd_size(2)/2-1, ...
                          -padd_size(3)/2:padd_size(3)/2-1);
    dist = sqrt(xx.^2 + yy.^2 + zz.^2);

    % Compute the padding size in mm for each dimension
    padd_in_mm = field_size .* padd_size ./ pixels;

    % Initialize output struct
    filter_params = struct('filter_highpass', [], 'filter_lowpass', []);

    % Highpass Filter
    if ~isempty(highpass_cutoff)
        hp = max(padd_in_mm) / highpass_cutoff;
        rise_hp = rise * hp;
        filter_params.filter_highpass = fftshift(1 ./ (1 + exp((dist - hp) ./ rise_hp)));
    end

    % Lowpass Filter
    if ~isempty(lowpass_cutoff)
        lp = max(padd_in_mm) / lowpass_cutoff;
        rise_lp = rise * lp;
        filter_params.filter_lowpass = fftshift(1 ./ (1 + exp((dist - lp) ./ rise_lp)));
    end
end

