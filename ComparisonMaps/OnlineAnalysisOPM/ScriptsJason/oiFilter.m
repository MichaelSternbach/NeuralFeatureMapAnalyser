function img = oiFilter(images, h),
%OIFILTER Spatial filtering of images.
%
%  IMG = OIFILTER(IMAGES, H) Performs spatial filtering of the supplied
%  image sequence(s) using the two-dimensional filter H.
%
%  IMAGES is a PxQxR matrix containing a single image sequence where
%
%    P is the height of a single frame (in pixels)
%    Q is the width of a single frame (in pixels)
%    R is the number of frames in the sequence
%
%  Filtering is performed by two-dimensional convolution. The output images
%  are the same size as the input images. To minimise artifacts at the edges
%  of the images, values outside the bounds of the input images are computed
%  by mirroring the image across the image border.
%
%  Alternatively, IMAGES may be an MxN cell array of image sequences where:
%
%    M is the number of conditions
%    N is the number of trials
%
%  in which case, filtering is applied to each individual image sequence. 
%
%  See also FSPECIAL, IMFILTER.

% 12/5/2011 - Shaun L. Cloherty <s.cloherty@ieee.org>

assert(nargin == 2, ...
  'Too many or too few input arguments. Type "help oiFilter" for usage information.');

if iscell(images),
  % call oiFilter( ) recursively for each image sequence
  [numConds, numTrials] = size(images);

  img = cell(size(images));
  for i = 1:numConds,
    for j = 1:numTrials,
      img{i,j} = oiFilter(images{i,j}, h);
    end
  end
  
  return
end

for i = 1:size(images,3),
  img(:,:,i) = imfilter(images(:,:,i),h,'symmetric','conv');
end
