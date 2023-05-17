function  h = oiPlotRegionalMean(images,varargin)
%OIPLOTREGIONALMEAN Plot the mean pixel intensity within subregions of an image sequence.
%
%  M = OIPLOTREGIONALMEAN(IMAGES, [, opt1, val1, ...]) plots the mean pixel
%  intensity within subregions of the supplied image sequence(s).
%
%  IMAGES is a PxQxR matrix containing a single image sequence where:
%
%    P is the height of a single frame (in pixels)
%    Q is the width os a single frame (in pixels)
%    R is the number of frames in the sequence
%
%  OIPLOTREGIONALMEAN returns a matrix of handles to lineseries objects, one
%  handle per plotted line.
%
%  Alternatively, IMAGES may be an MxN cell array of image sequences where:
%
%    M is the number of conditions
%    N is the number of trials
%
%  in which case, OIPLOTREGIONALMEAN returns a cell array of matricies
%  containing handles to lineseries objects. 
%
%  Available options are:
%   'subRegions' - a 2x1 or 1x2 vector specifying the number of subregions
%                  [rows, cols] within which to calculate the mean intensities
%                  (the default is to plot the mean intensities for 24
%                  subregions; 6 rows, 4 columns),
%   'showGrid'   - toggle display of the grid indicating region boundaries
%                  (default: true)
%   'refImage'   - the reference image over which to plot the regional
%                  means (default: images(:,:,1))

assert(nargin >= 1, ...
  'Too few input arguments. Type "help oiPlotRegionalMean" for usage information.');

if (isempty(images)),
  return;
end

% default options...
subRegions = [6, 4];
showGrid = true;

if iscell(images),
  refImage = images{1,1}(:,:,1);
else
  refImage = images(:,:,1);
end
  
if nargin > 1,
  if isstruct(varargin{1}),
    opts = varargin{1};
  else,
    if ischar(varargin{1}),
      % parse the options strings
      opts = parseOpts({'subRegions',subRegions,'showGrid',showGrid,'refImage',refImage}, varargin{:});
    else,
      error('Too many input arguments, and/or no invalid options found. Type "help oiPlotRegionalMean" for usage information.');
    end
  end
    
  % process the options structure
  if isfield(opts,'subRegions'), ...
      subRegions = opts.subRegions; end
  if isfield(opts,'showGrid'), ...
      showGrid = opts.showGrid; end
  if isfield(opts,'refImage'), ...
      refImage = opts.refImage; end
end

% calculate regional mean(s)
mn = oiCalcRegionalMean(images,'subRegions',subRegions);

% calculate subwindows
[hght,wdth] = size(refImage);
dh = (hght-1)/subRegions(1);
dw = (wdth-1)/subRegions(2);

% calculate indices of the top-left corner of each subregion
iH = floor([0:subRegions(1)]*dh)+1;
iW = floor([0:subRegions(2)]*dw)+1;

if iscell(mn),
  numFrames = size(mn{1,1},3);
else,
  numFrames = size(mn,3);
end

% scale regional means... vscale controls the minimum plot scale
vscale = 0.05;

mx = max(max(cellfun(@(x) max(abs(x(:))), mn)),vscale);
mn = cellfun(@(x) x/mx * dh/2, mn, 'UniformOutput', false);

imagesc(refImage,[min(refImage(:)), max(refImage(:))]);

colormap(gray);
hold on;

if showGrid,
  for i = 2:length(iW)-1,
    h(i-1) = line(iW(i)*[1,1], get(gca,'YLim'));
  end
  for j = 2:length(iH)-1,
    h(i-1+j-1) = line(get(gca,'XLim'),iH(j)*[1,1]);
  end
  set(h,'Color',0.5*[1,1,1],'LineStyle',':');
end

% this controls the horizontal 'separation' between plots...
hbuf = 0.2;

xc = ([0:(numFrames-1)]/(numFrames-1) * (1-hbuf) + hbuf/2)*dw; % x coordinate
iH = iH + dh/2;

% plot y=0 in gray
for i = 1:length(iH)-1,
  for j = 1:length(iW)-1,
    yc = zeros(size(xc));
    plot(xc+iW(j),yc+iH(i),'Color',0.5*[1,1,1]);
  end
end

% plot away...
if iscell(mn),
  [numConds, numTrials] = size(mn);
  h = cell(size(mn));
  for i = 1:numConds,
    colourSpec = (i/numConds)*[0, 0, 0.75] + [0, 0, 0.25];
    for j = 1:numTrials,
      h{i,j} = plotmn(mn{i,j},colourSpec);
    end
  end
else,
  h = plotmn(mn,'b');
end

% use a nested function to do the plotting... I hate nested functions, but
% nested functions retain access to variables in the containing function,
% which is convenient here    
function hh = plotmn(mn,colourSpec),
  for ii = 1:length(iH)-1,
    for jj = 1:length(iW)-1,
      yc = squeeze(-mn(ii,jj,:)); % -ve - 'image' origin is at the top-left, not the bottom-left
      hh(ii,jj) = plot(xc+iW(jj),yc+iH(ii),'Color',colourSpec);
    end
  end
end % end nested function

hold off;

axis image;

end % end containing function - keep this!
