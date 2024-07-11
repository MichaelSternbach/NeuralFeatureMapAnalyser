function [chi2Stat, pValue, numBins] = compareDistributionsChiSquared(data1, data2, numBins)
    % Compare whether two datasets follow the same underlying distribution
    % using a chi-squared test
    %
    % Inputs:
    %   data1 - First dataset (vector)
    %   data2 - Second dataset (vector)
    %   numBins - Number of bins (integer) or method ('sturges', 'rice', 'scott', 'fd')
    %
    % Outputs:
    %   chi2Stat - Chi-squared statistic
    %   pValue - P-value of the chi-squared test
    %   numBins - Number of bins used
    
    % Check the number of input arguments
    if nargin == 2
        % Default to Freedman-Diaconis rule if numBins is not provided
        disp('Number of bins not provided, using Freedman-Diaconis rule.');
        numBins = 'fd';
    end
    
    % Combine data
    combinedData = [data1; data2];
    
    % Determine number of bins if numBins is a method name
    if ischar(numBins)
        numBins = calculateNumBins(combinedData, numBins);
    end
    
    % Create bins
    [counts, edges] = histcounts(combinedData, numBins);
    
    % Get observed frequencies for each dataset
    obs1 = histcounts(data1, edges);
    obs2 = histcounts(data2, edges);
    
    % Total observations
    n1 = length(data1);
    n2 = length(data2);
    nTotal = n1 + n2;
    
    % Calculate expected frequencies under the null hypothesis
    expected1 = (counts * (n1 / nTotal));
    expected2 = (counts * (n2 / nTotal));
    
    % Combine observed and expected frequencies
    observed = [obs1; obs2];
    expected = [expected1; expected2];
    
    % Compute the chi-squared statistic
    chi2Stat = sum(((observed - expected).^2) ./ expected, 'all');
    
    % Degrees of freedom
    df = numBins - 1;
    
    % Calculate p-value
    pValue = 1 - chi2cdf(chi2Stat, df);
end
