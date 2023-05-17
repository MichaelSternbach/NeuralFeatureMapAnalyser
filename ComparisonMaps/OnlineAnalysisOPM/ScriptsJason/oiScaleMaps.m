function [maps, b, hmax] = oiScaleMaps(maps),
%OISCALEMAPS Rescales single condition maps.
%
%  MAPS = OISCALEMAPS(MAPS) rescales the responses at each map location
%  according to:
%
%      A'(k) = 100 * (A(k) - b)/Hmax
%
%  where A(k) is the response amplitude at location (i,j) for stimulus
%             condition k,
%        b    is the minimum of all response values at location (i,j)
%             (i.e., the baseline or stimulus non-specific response level), and
%        Hmax is the maximum value of H across the whole map, with H being
%             defined at the difference between the maximum and minimum
%             response levels at location (i,j).
%
%  [MAPS,B] = OISCALEMAPS(MAPS) returns the stimulus non-specific baseline
%  for each point in the map.
%
%  [MAPS,B,HMAX] = OISCALEMAPS(MAPS) returns the maximum value of H across
%  the whole map.
%
%  MAPS is a PxQxR matrix containing the single condition maps, where:
%
%    P is the height of the maps (in pixels)
%    Q is the width of the maps (in pixels)
%    R is the number of maps/directions.
%
%  For further details see Swindale et al. (2003) Cereb. Cortex 13:225-238. 

% 12/5/2011 - Shaun L. Cloherty <s.cloherty@ieee.org>

% check input and output arguments
assert(nargin == 1, ...
  'Too few input arguments. Type "help oiScaleMaps" for usage information.');

assert(~iscell(maps), ...
  'The input argument MAPS does not appear to be a set of single condition maps. Type "help oiScaleMaps" for usage information.');

if (isempty(maps)),
  return;
end

numConds = size(maps,3);

b = min(maps,[],3);
h = range(maps,3);
hmax = max(h(:));

for i = 1:numConds,
  maps(:,:,i) = 100*(maps(:,:,i) - b)/hmax;
end
