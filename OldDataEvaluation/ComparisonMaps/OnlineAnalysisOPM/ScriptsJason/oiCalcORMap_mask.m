function [phi, r] = oiCalcORMap(maps, theta, mode, grnImg, varargin),
%OICALCORMAP Calculates the orientation preference map.
%
%  [PHI,R] = OICALCORMAP(MAPS, DIRS, MODE[, opt1, val1, ...]) calculates
%  the orientation preference map from the supplied single condition maps.
%
%  MAPS is a PxQxR matrix containing the single condition maps, where:
%
%    P is the height of the maps (in pixels)
%    Q is the width of the maps (in pixels)
%    R is the number of maps/directions.
%
%  THETA is a 1xR or Rx1 vector containing the stimulus orientation for each
%  map (0 < theta < pi).
%
%  MODE is one of:
%    'vector' - the preferred orientation is computed by vector averaging. This
%               is mathematically equivalent to fitting a cosine function
%               to the tuning data. This is the default mode of operation.
%    'model'  - the preferred orientation is computed by fitting a circular
%               normal (von Mises) distributions to the tuning data.
%
%  For details see Swindale et al. (2003) Cereb. Cortex 13:225-238. 
%
%  See also OISCALEMAPS.

% check input and output arguments
assert(nargin >= 2, ...
  'Too few input arguments. Type "help oiCalcORMap" for usage information.');

assert(~iscell(maps), ...
  'The input argument MAPS does not appear to be a set of single condition maps. Type "help oiCalcORMap" for usage information.');

if (isempty(maps)),
  return;
end

assert(size(maps,3) == length(theta), ...
  'Dimension mismatch! Check that length(THETA) == size(MAPS,3). Type "help oiCalcORMap" for usage information.');

if nargin < 3,
  mode = 'vector';
end

numConds = size(maps,3);

switch lower(mode),
  case 'vector',
    % calculate orientation preference by vector averaging
    for idx = 1:numConds,
      x(:,:,idx) = maps(:,:,idx) .* cos(2*theta(idx));
      y(:,:,idx) = maps(:,:,idx) .* sin(2*theta(idx));
    end

    x = sum(x,3);
    y = sum(y,3);

    phi = atan2(y,x); % atan2 returns angles between -pi and +pi
    phi = phi+(sign(phi) < 0)*(2*pi); % 0 < phi < 2*pi

    phi = 0.5*phi;
    
    phi = phi.*grnImg(1:end,1:end);
    
    r = sqrt(x.^2 + y.^2);
    

    
  case 'model',
    % calculate orientation preference by fitting a model tuning curve

    [m,n] = size(maps(:,:,1));
    
    x = 180*(2*theta)/pi; % radians to degrees
    for i = 1:m,
      for j = 1:n,
        y = squeeze(maps(i,j,:));
        [mx,k] = max(y(:));

        p0 = [mx, x(k), 0.59];
        
        p(i,j,:) = fitVonMises(x', y, p0, 'funcType', 'single');
        
        phi(i,j) = 0.5*p(i,j,2)*pi/180; % degrees to radians
        r(i,j) = p(i,j,1);
      end
    end
end
