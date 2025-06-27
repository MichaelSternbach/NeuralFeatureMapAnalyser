function weightedhistogram(data, weights, varargin)
% WEIGHTEDHISTOGRAM Plots a weighted histogram outline (counts or pdf)
%
%   weightedhistogram(DATA, WEIGHTS) plots a histogram outline of DATA weighted by WEIGHTS.
%
%   Name-value options:
%       'BinEdges'  - custom bin edges
%       'BinWidth'  - bin width
%       'NumBins'   - number of bins
%       'Normalize' - 'count' (default) or 'pdf'
%       'Label'     - string for legend
%       ...         - any plot() styling arguments

    arguments
        data (:,1) double
        weights (:,1) double
    end
    arguments (Repeating)
        varargin
    end

    % Parse binning and other named parameters
    p = inputParser;
    addParameter(p, 'BinEdges', []);
    addParameter(p, 'BinWidth', []);
    addParameter(p, 'NumBins', []);
    addParameter(p, 'Normalize', 'count'); % or 'pdf'
    addParameter(p, 'Label', '');
    parse(p, varargin{:});

    binEdges = p.Results.BinEdges;
    binWidth = p.Results.BinWidth;
    numBins = p.Results.NumBins;
    normalize = lower(p.Results.Normalize);
    label = p.Results.Label;

    % Remove parsed parameters from varargin
    plotArgs = varargin;
    toRemove = {'BinEdges', 'BinWidth', 'NumBins', 'Normalize', 'Label'};
    for i = 1:length(toRemove)
        idx = find(strcmpi(plotArgs, toRemove{i}));
        if ~isempty(idx)
            plotArgs(idx:idx+1) = [];
        end
    end

    % Determine bin edges
    if ~isempty(binEdges)
        edges = binEdges;
    elseif ~isempty(binWidth)
        minD = min(data);
        maxD = max(data);
        edges = minD:binWidth:(maxD + binWidth);
    elseif ~isempty(numBins)
        edges = linspace(min(data), max(data), numBins + 1);
    else
        [~, edges] = histcounts(data, 'BinMethod', 'auto');
    end

    % Bin data
    binIdx = discretize(data, edges);
    valid = ~isnan(binIdx);
    binIdx = binIdx(valid);
    weights = weights(valid);

    % Weighted counts per bin
    counts = accumarray(binIdx, weights, [length(edges)-1, 1], @sum, 0);

    % Normalize if requested
    binWidths = diff(edges);
    if strcmp(normalize, 'pdf')
        counts = counts ./ (sum(weights) * binWidths');
    end

    % Build x and y for stairs
    x = repelem(edges, 2);
    x = x(2:end-1);
    y = repelem(counts', 2);
    y = [0, y, 0];
    x = [x(1), x, x(end)];

    % Plot
    if ~isempty(label)
        plot(x, y, 'DisplayName', label, plotArgs{:});
    else
        plot(x, y, plotArgs{:});
    end
end
