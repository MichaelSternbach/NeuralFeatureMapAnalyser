
function n = orientation_hist(op, mask, bins, weighted)
%ORIENTATION_HIST Calculate an orientation preference histogram.
%
% n = orientation_hist(op, mask, bins, weighted)
%
% Given an orientation preference map in complex form, calculates a histogram
% of the proportion of pixels with preferred orientation falling into each
% bin. A binary mask can be specified to only include a certain portion of a
% map. The bins used by default are -5:10:185, where the first and last bins
% are used to correctly handle angles wrapping around. If weighted is true then
% a weighted histogram is calculated using the tuning strength of pixel.
%
% To generate weighted histograms, the histwc function must be available:
% http://www.mathworks.com/matlabcentral/fileexchange/42493-generate-weighted-histogram/content/histwc/histwc/histwc.m

% Default inputs
if nargin < 2
    mask = true(size(op));
end
if nargin < 3
    bins = -5:10:185;
end
if nargin < 4
    weighted = false;
end

% Mask
op(~mask) = nan;

% Calculate orientation distribution
thetas = wrapTo2Pi(angle(op(:))+pi)/2*180/pi;
if weighted
    a = histwc(thetas, abs(op(:)), bins);
else
    a = histc(thetas, bins);
end
a(1) = a(1) + a(end-1);
a = a(1:end-2);
n = a/sum(a);
© 2020 GitHub, Inc.