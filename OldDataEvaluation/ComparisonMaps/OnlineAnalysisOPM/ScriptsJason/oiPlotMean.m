function h = oiPlotMontage(images, varargin),
%OIPLOTMEAN Plot the mean pixel intensity of an image sequence.
%
%  OIPLOTMEAN(IMAGES, [, opt1, val1, ...]) plots the mean pixel intensity of
%  the supplied image sequence as a time series (one value for each image).
%
%  Available options are:
%   'roi' - a vector specifying the region of interest [x,y,wdth,hght] over
%           which to calculate the mean intensity (the default is to plot
%           the mean intensity for the whole frame),

% $Id: oiPlotMean.m,v 1.1 2008-04-30 04:56:07 shaunc Exp $

assert(nargin >= 1, ...
  'Too few input arguments. Type "help oiPlotMean" for usage information.');

assert(~iscell(images), ...
  'The input argument "images" does not appear to be an image sequence. Type "help oiPlotMontage" for usage information.');

if (isempty(images)),
  return;
end

% default options...
roi = [];

if nargin > 1,
  if isstruct(varargin(1)),
    opts = varargin(1);
  else,
    if ischar(varargin(1)),
      % parse the options strings
      opts = parseOpts({'roi',[]}, varargin{:});
    else,
      error('Too many input arguments, and/or no invalid options found. Type "help oiPlotMean" for usage information.');
    end
  end
    
  % process the options structure
  if isfield(opts,'roi'), ...
      roi = opts.roi; end
end

[height,width,numFrames] = size(images);

if isempty(roi),
  roi = [1, 1, width, height];
end

refWin = [1:5];
sigWin = [31:35];

ref = images(:,:,refWin);
sig = images(:,:,sigWin);
tmp_ref = mean(mean(ref(roi(2):(roi(2)+roi(4)-1),roi(1):(roi(1)+roi(3)-1),:)));
m_r = squeeze(tmp_ref);
ref_mean = mean(m_r);

tmp = mean(mean(images(roi(2):(roi(2)+roi(4)-1),roi(1):(roi(1)+roi(3)-1),:)))
m= squeeze(tmp);

ratio = (m-ref_mean)/ref_mean*0.001;

figure;
plot(0:numFrames-1,-ratio,'LineWidth',1,'color',[0 0.5 0]);

xlabel('Frame #');
ylabel('Mean Pixel Value');
