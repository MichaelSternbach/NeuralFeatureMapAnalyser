function snr_values = calculate_snr(data, average_over_dimensions, output_in_db)
%CALCULATE_SNR Calculate the Signal-to-Noise Ratio (SNR) for n-dimensional arrays.
%   snr_values = CALCULATE_SNR(data, average_over_dimensions, output_in_db) computes the SNR
%   for a dataset where the last dimension represents the set iterations.
%   Optionally averages over all other dimensions if specified and provides
%   the option to output the SNR in dB or as a linear value.
%
%   Inputs:
%       data - n-dimensional array, where the last dimension represents
%              the set iterations. Supports complex values.
%       average_over_dimensions - logical flag (true/false). If true, the
%              mean and variance are calculated by averaging over all
%              dimensions, including the last.
%              Default is false.
%       output_in_db - logical flag (true/false). If true, the SNR is
%              returned in dB. If false, the linear SNR is returned.
%              Default is false.
%
%   Output:
%       snr_values - If average_over_dimensions is false, this is an
%                    n-1 dimensional array containing the SNR values
%                    for each entry along the other dimensions.
%                    If true, a scalar value representing the averaged SNR.

% Ensure the input is valid
if nargin < 1
    error('Input data is required.');
end

if nargin < 2
    average_over_dimensions = true;
end

if nargin < 3
    output_in_db = false;
end

% Get the size of the input data
input_size = size(data);

% Number of iterations (last dimension)
n_iterations = input_size(end);

if n_iterations < 2
    error('The last dimension (iterations) must have at least 2 elements.');
end

if average_over_dimensions
    % Calculate mean and variance over all dimensions, including the last
    signal_power = mean(abs(data(:)))^2;
    noise_power = var(abs(data(:)), 0);
else
    % Reshape the data to separate the last dimension
    reshaped_data = reshape(data, [], n_iterations);

    % Calculate the mean (signal) and standard deviation (noise) along the last dimension
    signal_power = mean(abs(reshaped_data), 2).^2;
    noise_power = var(abs(reshaped_data), 0, 2);
end

% Compute the linear SNR
snr_linear = signal_power ./ noise_power;

if output_in_db
    % Convert to dB
    snr_values = 10 * log10(snr_linear);
else
    snr_values = snr_linear;
end

if ~average_over_dimensions
    % Reshape the result to match the original dimensions (excluding the last dimension)
    out_size = input_size(1:end-1);
    snr_values = reshape(snr_values, out_size);
end
end
