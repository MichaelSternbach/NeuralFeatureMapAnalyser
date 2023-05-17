function [phi, r] = oiCalcDRMap(maps, dirs, mode, varargin),
%OICALCDRMAP Calculates the direction preference map.
%
%  [PHI,R] = OICALCDRMAP(MAPS, DIRS, MODE[, opt1, val1, ...]) calculates
%  the direction preference map from the supplied single condition maps.
%
%  MAPS is a PxQxR matrix containing the single condition maps, where:
%
%    P is the height of the maps (in pixels)
%    Q is the width of the maps (in pixels)
%    R is the number of maps/directions.
%
%  DIRS is a 1xR or Rx1 vector containing the stimulus direction for each
%  map (0 < dirs < 2*pi).
%
%  MODE is one of:
%    'vector' - the preferred direction is computed by vector averaging. This
%               is mathematically equivalent to fitting a cosine function
%               to the tuning data. This is the default mode of operation.
%    'model'  - the preferred direction is computed by fitting a sum of 
%               two circular normal (von Mises) distributions to the tuning
%               data.
%
%  For details see Swindale et al. (2003) Cereb. Cortex 13:225-238. 
%
%  See also OISCALEMAPS.

% check input and output arguments
assert(nargin >= 2, ...
  'Too few input arguments. Type "help oiCalcDRMap" for usage information.');

assert(~iscell(maps), ...
  'The input argument MAPS does not appear to be a set of single condition maps. Type "help oiCalcDRMap" for usage information.');

if (isempty(maps)),
  return;
end

assert(size(maps,3) == length(dirs), ...
  'Dimension mismatch! Check that length(DIRS) == size(MAPS,3). Type "help oiCalcDRMap" for usage information.');

if nargin < 3,
  mode = 'vector';
end

numConds = size(maps,3);

switch lower(mode),
  case 'vector',
    % calculate direction preference by vector averaging
    for idx = 1:numConds,
      x(:,:,idx) = maps(:,:,idx) .* cos(dirs(idx));
      y(:,:,idx) = maps(:,:,idx) .* sin(dirs(idx));
    end

    x = sum(x,3);
    y = sum(y,3);

    phi = atan2(y,x); % atan2 returns angles between -pi and +pi
    phi = phi+(sign(phi) < 0)*(2*pi); % 0 < phi < 2*pi

    r = sqrt(x.^2 + y.^2);
  case 'model',
    % calculate direction preference by fitting a model tuning curve

    [m,n] = size(maps(:,:,1));
    
    x = 180*dirs/pi; % radians to degrees
    for i = 1:m,
      for j = 1:n,
        y = squeeze(maps(i,j,:));
        [mx,k] = max(y(:));

        p0 = [mx, x(k), 2.0, 0.8*mx, x(k)+180, 2.0];
        
        p(i,j,:) = fitVonMises(x', y, p0, 'funcType', 'double');
        
        if p(i,j,1) >= p(i,j,4),
          phi(i,j) = p(i,j,2)*pi/180;
          r(i,j) = p(i,j,1);
        else,
          phi(i,j) = p(i,j,5)*pi/180;
          r(i,j) = p(i,j,4);
        end
      end
    end
end
