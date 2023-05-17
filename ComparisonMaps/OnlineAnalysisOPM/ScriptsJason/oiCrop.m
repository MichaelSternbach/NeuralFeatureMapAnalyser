function img = oiCrop(images, xx, yy),
%OICROP Crop the supplied image sequence(s).
%
%  IMG = OICROP(IMAGES, X, Y) Crops (spatially) the supplied image 
%  sequence(s).
%
%  IMAGES is a PxQxR matrix containing a single image sequence where
%
%    P is the height of a single frame (in pixels)
%    Q is the width os a single frame (in pixels)
%    R is the number of frames acquired in the sequence
%
%  X and Y are 1x2 (or 2x1) vectors specifying the upper and lower (X), and the
%  left and right (Y) bounds of the cropped region.
%
%  Alternatively, IMAGES may be an MxN cell array of image sequences where:
%
%    M is the number of conditions
%    N is the number of trials
%
%  in which case, binning is calculated for each individual image
%  sequence. 

% $Id: oiCrop.m,v 1.1 2008-04-30 04:56:07 shaunc Exp $

assert(nargin == 3, ...
  'Too many or too few input arguments. Type "help oiCrop" for usage information.');

if iscell(images),
  % call oiCrop( ) recursively for each image sequence

  [numConds, numTrials] = size(images);

  img = cell(size(images));
  for i = 1:numConds,
    for j = 1:numTrials,
      img{i,j} = oiCrop(images{i,j}, xx, yy);
    end
  end
  
  return
end

% crop spatially
img = images(xx(1):xx(2),yy(1):yy(2),:);

% FIXME: crop in time?
