function mn = oiCalcRegionalMean(images,varargin)
%OICALCREGIONALMEAN Calculate mean pixel intensity within subregions of an image sequence.
%
%  M = OICACLREGIONALMEAN(IMAGES, [, opt1, val1, ...]) calculates the mean pixel intensity
%  within subregions of the supplied image sequence.
%
%  IMAGES is a PxQxR matrix containing a single image sequence where:
%
%    P is the height of a single frame (in pixels)
%    Q is the width os a single frame (in pixels)
%    R is the number of frames in the sequence.
%
%  OICALCREGIONALMEAN returns a matrix containing the mean time series in
%  each region.
%
%  Alternatively, IMAGES may be an MxN cell array of image sequences where:
%
%    M is the number of conditions
%    N is the number of trials
%
%  in which case, OICALCREGIONALMEAN returns a cell array of matricies
%  containing the mean time series in each region within each image sequence.
%
%  Available options are:
%   'subRegions' - a 2x1 or 1x2 vector specifying the number of subregions
%                  [rows, cols] within which to calculate the mean intensities
%                  (the default is to plot the mean intensities for the
%                  whole image; 1 rows, 1 columns),

assert(nargin >= 1, ...
  'Too few input arguments. Type "help oiCalcRegionalMean" for usage information.');

% default options...
subRegions = [1, 1];

if nargin > 1,
  if isstruct(varargin{1}),
    opts = varargin{1};
  else,
    if ischar(varargin{1}),
      % parse the options strings
      opts = parseOpts({'subRegions',subRegions}, varargin{:});
    else,
      error('Too many input arguments, and/or no invalid options found. Type "help oiCalcRegionalMean" for usage information.');
    end
  end
    
  % process the options structure
  if isfield(opts,'subRegions'), ...
      subRegions = opts.subRegions; end
end


if iscell(images),
  % call oiCalcRegionalMean( ) recursively for each image sequence

  [numConds, numTrials] = size(images);

  mn = cell(size(images));
  for i = 1:numConds,
    for j = 1:numTrials,
      mn{i,j} = oiCalcRegionalMean(images{i,j},'subRegions',subRegions);
    end
  end
  
  return
end

% calculate subwindows
[hght,wdth,numFrames] = size(images);
dh = (hght-1)/subRegions(1);
dw = (wdth-1)/subRegions(2);

% calculate indices of the top-left corner of each subregion
iH = floor([0:subRegions(1)]*dh)+1;
iW = floor([0:subRegions(2)]*dw)+1;

mn = zeros([subRegions numFrames]);
for i = 1:subRegions(1),
  for j = 1:subRegions(2),
    tmp = mean(images(iH(i):iH(i+1),iW(j):iW(j+1),:),1);
    tmp = squeeze(mean(tmp,2));
    mn(i,j,:) = tmp;
  end
end
