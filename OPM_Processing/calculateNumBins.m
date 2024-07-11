function numBins = calculateNumBins(data, method)
    % Calculate the optimal number of bins for a dataset based on the specified method
    %
    % Inputs:
    %   data - Dataset (vector)
    %   method - Method for calculating number of bins ('sturges', 'rice', 'scott', 'fd')
    %
    % Output:
    %   numBins - Number of bins
    
    if nargin < 2
        method = 'fd';
        disp('Method is not provided, using Freedman-Diaconis rule.');
    end

    N = length(data);
    switch lower(method)
        case 'sturges'
            % Sturges' Formula
            numBins = ceil(log2(N) + 1);
        case 'rice'
            % Rice Rule
            numBins = ceil(2 * N^(1/3));
        case 'scott'
            % Scott's Normal Reference Rule
            sigma = std(data);
            binWidth = 3.5 * sigma * N^(-1/3);
            numBins = ceil((max(data) - min(data)) / binWidth);
        case 'fd'
            % Freedman-Diaconis Rule
            IQR = iqr(data);
            binWidth = 2 * IQR * N^(-1/3);
            numBins = ceil((max(data) - min(data)) / binWidth);
        otherwise
            error('Unknown method for determining number of bins.');
    end
end
