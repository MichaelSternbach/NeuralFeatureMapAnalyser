function [bootstrapSamples, jackknifeSamples] = resample_bootstrap_jackknife(data, nBoot)
% resample_bootstrap_jackknife - Generate bootstrap and jackknife samples from a vector.
%
% Syntax:
%   [bootstrapSamples, jackknifeSamples] = resample_bootstrap_jackknife(data, nBoot)
%
% Inputs:
%   data   - A vector of data values.
%   nBoot  - Number of bootstrap samples to generate.
%
% Outputs:
%   bootstrapSamples - A matrix of bootstrap samples (length(data) x nBoot)
%   jackknifeSamples - A matrix of jackknife samples ((length(data)-1) x length(data))

    if ~isvector(data)
        error('Input "data" must be a vector.');
    end

    data = data(:); % Ensure column vector
    n = length(data);

    % Generate bootstrap samples
    bootstrapSamples = zeros(n, nBoot);
    for i = 1:nBoot
        indices = randi(n, n, 1);
        bootstrapSamples(:, i) = data(indices);
    end

    % Generate jackknife samples
    jackknifeSamples = zeros(n-1, n);
    for i = 1:n
        jackknifeSamples(:, i) = data([1:i-1, i+1:end]);
    end
end