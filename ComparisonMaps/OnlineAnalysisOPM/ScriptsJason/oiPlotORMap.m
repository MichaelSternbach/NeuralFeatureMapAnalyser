function h = oiPlotORMap(images, varargin),
%OIPLOTORMAP Plot the orientation preference map.
%
%  OIPLOTORMAP(IMAGES, [, opt1, val1, ...]) plots the orientation
%  preference map based on the supplied single condition maps.
%
%  Available options are:
%   none

% FIXME: add options as appropriate
%   'xxx'       - a vector specifying something...

% $Id: oiPlotORMap.m,v 1.2 2009-06-09 00:51:59 shaunc Exp $

% % check input and output arguments
% error(nargchk(1, 4, nargin));

assert(nargin >= 1, ...
  'Too few input arguments. Type "help oiPlotORMap" for usage information.');

assert(~iscell(images), ...
  'The input argument "images" does not appear to be a set of image sequences. Type "help oiPlotORMap" for usage information.');

if (isempty(images)),
  return;
end

% default options...
normalise = false;
polar = true;

if nargin > 1,
  if isstruct(varargin{1}),
    opts = varargin{1};
  else,
    if ischar(varargin{1}),
      % parse the options strings
      opts = parseOpts({'normalise',false,'polar',true}, varargin{:});
    else,
      error('Too many input arguments, and/or no invalid options found. Type "help oiPlotORMap" for usage information.');
    end
  end
     
  % process the options structure
  if isfield(opts,'normalise'), ...
      normalise = opts.normalise; end
  if isfield(opts,'polar'), ...
      polar = opts.polar; end
end

numConds = size(images,3);

%figure;

% calculate orientation preference by vector addition, see Bonhoeffer, T.,
% and Grinvald, A. (1991). Iso-orientation domains in cat visual cortex
% are arranged in pinwheel-like patterns. Nature 353, 429–431.

if normalise,
  % FIXME: need to think about this!
  for idx = 1:numConds,
    tmp = images(:,:,idx);
    images(:,:,idx) = (tmp - mean(tmp(:))) / std(tmp(:));
  end
end

%images = abs(images);
for idx = 1:numConds,
  x(:,:,idx) = images(:,:,idx) .* cos((idx-1)*2*(pi/numConds));
  y(:,:,idx) = images(:,:,idx) .* sin((idx-1)*2*(pi/numConds));
end

x = sum(x,3);
y = sum(y,3);

% pref = atan2(y,x) + pi; % atan2 returns angles between -pi and +pi
pref = atan2(y,x); % atan2 returns angles between -pi and +pi
pref = pref+(sign(pref) < 0)*(2*pi); % 0 < pref < 2*pi

pref = 256*pref/(2*pi);
pref = double(pref);

subimage(zeros(size(pref)),gray(2));
hold on;

h = subimage(pref, hsv(256));

if polar,
  mag = sqrt(x.^2 + y.^2);
%   mag(mag < 0.1) = 0.0;
  set(h,'AlphaData', min(sqrt(mag/max(mag(:)/6)), 1));
end
hold off;

h = pref;

