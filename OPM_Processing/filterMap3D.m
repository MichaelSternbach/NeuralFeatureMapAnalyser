function z_filtered = filterMap3D(z, cutoff_length,steepness)
    % This function filters a 3D complex dataset `z` using a Fermi filter.
    % The cutoff frequency is translated to a spatial cutoff length.
    %
    % Inputs:
    % - z: the input 3D complex data array
    % - cutoff_length: the desired cutoff length in the spatial domain
    %
    % Output:
    % - z_filtered: the filtered 3D complex data

    if nargin < 3
        steepness = 2; % Define steepness of the filter
    end
    
    %% Step 1: Zero-padding to the next power of two
    padded_size = max(size(z));%2.^nextpow2(max(size(z)));  % Compute the next power of two for each dimension
    
    padded_data = padarray(z, [padded_size padded_size padded_size] - size(z), 0, 'post');  % Zero-pad the data
    
    %% Step 2: FFT of the padded complex data
    data_fft = fftn(padded_data);
    
    %% Step 3: Convert cutoff length to cutoff frequency
    % The cutoff frequency corresponds to 1/cutoff_length, but we need to
    % normalize it according to the size of the padded array.
    voxel_size = 1;  % Assuming isotropic voxel size of 1 unit, adjust if needed
    cutoff_frequency = padded_size(1) / (voxel_size * cutoff_length);  % Convert length to frequency
    
    %% Step 4: Fermi filter design
    lp=padded_size/cutoff_length;
    rise_lp = steepness*lp;
    filter = fermi_filter(padded_size, cutoff_frequency, rise_lp);
    
    %% Step 5: Apply the filter in the frequency domain (element-wise multiplication)
    data_fft_filtered = data_fft .* filter;
    
    %% Step 6: Inverse FFT to return to the spatial domain
    data_filtered_padded = ifftn(data_fft_filtered);
    
    %% Step 7: Remove padding
    %z_filtered = data_filtered_padded(1:size(z,1), 1:size(z,2), 1:size(z,3));
    z_filtered = data_filtered_padded;
end

function filter_3d = fermi_filter(shape, cutoff_frequency, steepness)
    % Create a Fermi filter in 3D frequency space
    %
    % Inputs:
    % - shape: the size of the 3D array
    % - cutoff_frequency: the cutoff frequency for the filter (normalized)
    % - steepness: the steepness of the Fermi filter transition
    %
    % Output:
    % - filter_3d: the Fermi filter in 3D

    % Generate the 3D frequency grid
    [kx, ky, kz] = ndgrid(1:shape, 1:shape, 1:shape);
    
    % Shift to center frequencies
    center = ceil([shape shape shape] / 2);
    kx = kx - center(1);
    ky = ky - center(2);
    kz = kz - center(3);
    
    % Compute the radial frequency in the grid
    radius = sqrt(kx.^2 + ky.^2 + kz.^2);
   
    % Apply the Fermi filter formula: 1 / (1 + exp((r - cutoff) / steepness))
    filter_3d = 1./(1 + exp((radius - cutoff_frequency)./steepness));
end
