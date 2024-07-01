
function [pinw, windno] = locate_pinwheels(asd)
%LOCATE_PINWHEELS Locate pinwheels in an OP map.
%
% [pinw, windno] = locate_pinwheels(map)
%
% Finds the location and winding number of all pinwheels in the given OP map
% in complex form.
% 
% Notes:
%
% - The sign of a pinwheel is defined positive if the orientation increases
%   clockwise and negative if anticlockwise. This seems to be the adopted
%   convenion in neuroscience, but unfortunately it is the opposite of the
%   one used in physics!
%
% - The positive (pinwp) and negative (pinwn) pinwheel locations can be
%   obtained as follows:
%     pinwp = pinw(windno > 0,:);
%     pinwn = pinw(windno < 0,:);
%   In principle, the winding number should be either 1 or -1 for all
%   pinwheels (and 0 for non-pinwheel points), but it might happen that it
%   is ±2 or higher for a pinwheel. I don't think that has ever happened to
%   anyone with OR maps, but keep an eye on it!
%
% Copyright (c) 2002 by Miguel A. Carreira-Perpinan

% Contour integral of radius one around each pixel (except the border pixels)
% angles = mod(degtorad(asd), 2*pi);
maps_c= exp(1i*asd/256*2*pi);
angles = mod(angle(maps_c), 2*pi);

pinwq = nan(size(asd));
for j = 2:size(asd,1)-10
    for k = 2:size(asd,2)-1
        indicies = [j-1,k-1; j-1,k; j-1,k+1; j,k+1; ...
                    j+1,k+1; j+1,k; j+1,k-1; j,k-1; ...
                    j-1,k-1];
        L = zeros(length(indicies),1);
        for p = 1:length(indicies)
            L(p) = angles(indicies(p,1),indicies(p,2));
        end
        if any(isnan(L))
            pinwq(j,k) = 0;
        else
            pinwq(j,k) = sum(diff(my_unwrap(L)));
        end
    end
end
pinwq = abs(pinwq) > 0.1;

% Join clusters of pinwheels and reduce to their centre of mass

[tmp, np] = bwlabel(pinwq(:,:,1), 8);
tmp = regionprops(tmp, 'Centroid');
pinw = reshape([tmp.Centroid], 2, np)';

% Calculate winding numbers
if nargout == 2
    if isempty(pinw)
        windno = [];
        return;
    end
    windno = zeros(length(pinw),1);
    for l = 1:size(pinw,1)
        j = round(pinw(l,2));
        k = round(pinw(l,1));
        indicies = [j-1,k-1; j-1,k; j-1,k+1; j,k+1; ...
                    j+1,k+1; j+1,k; j+1,k-1; j,k-1; ...
                    j-1,k-1];
        L = zeros(length(indicies),1);
        for p = 1:length(indicies)
            L(p) = angles(indicies(p,1),indicies(p,2));
        end
        windno(l) = sum(diff(unwrap(L)))/(2*pi);
    end

    % Check for some possible errors
    if any(abs(round(windno)-windno) > 10^(-5))
        warning('Some winding numbers are not exactly integers.');
    end
    if any(abs(windno)-1>eps*10)
        warning('Some winding numbers are different from ±1.');
    end
end


function q = my_unwrap(p)
% A simplified version of unwrap because p is always a column vector here.

% Unwrap each column of p
q = p;
% Find NaN's and Inf's
indf = find(isfinite(p(:)));
% Unwrap finite data (skip non finite entries)
q(indf) = LocalUnwrap(p(indf));

function p = LocalUnwrap(p)
%LocalUnwrap   Unwraps column vector of phase values.

m = length(p);

% Unwrap phase angles.  Algorithm minimizes the incremental phase variation 
% by constraining it to the range [-pi,pi]
dp = diff(p,1,1);                % Incremental phase variations
dps = mod(dp+pi,2*pi) - pi;      % Equivalent phase variations in [-pi,pi)
dps(dps==-pi & dp>0,:) = pi;     % Preserve variation sign for pi vs. -pi
dp_corr = dps - dp;              % Incremental phase corrections
dp_corr(abs(dp)<pi,:) = 0;       % Ignore correction when incr. variation is < CUTOFF

% Integrate corrections and add to P to produce smoothed phase values
p(2:m,:) = p(2:m,:) + cumsum(dp_corr,1);
