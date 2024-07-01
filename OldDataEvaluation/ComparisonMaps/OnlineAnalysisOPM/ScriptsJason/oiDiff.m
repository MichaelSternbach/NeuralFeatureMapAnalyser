function dimg = oiDiff(images, ref);
%OIDIFF Calcuate difference image sequence.
%
%  DIMG = OIDIFF(IMAGES, REF) Calculates the difference image sequence
%  by subtracting the reference image from all images in the sequence.
%
%  IMAGES is a PxQxR matrix containing a single image sequence where
%
%    P is the height of a single frame (in pixels)
%    Q is the width os a single frame (in pixels)
%    R is the number of frames acquired in the sequence
%
%  REF is either a PxQ matrix containing the reference image, or an Nx1
%  (or 1xN) vector of frame numbers to be averaged to obtain the reference
%  image.
%
%  Alternatively, IMAGES may be an MxN cell array of image sequences where:
%
%    M is the number of conditions
%    N is the number of trials
%
%  in which case, the difference is calculated for each individual image
%  sequence. 

assert(nargin >= 2, ...
  'Too few input arguments. Type "help oiDiff" for usage information.');

if iscell(images),
  % call oiDiff( ) recursively for each image sequence

  [numConds, numTrials] = size(images);

  dimg = cell(size(images));
  for i = 1:numConds,
    for j = 1:numTrials,
      dimg{i,j} = oiDiff(images{i,j}, ref);
    end
  end
  
  return
end

if isvector(ref),
%   ref = int16(mean(images(:,:,ref),3));
  ref = single(mean(images(:,:,ref),3));
end

for i = 1:size(images,3),
%   dimg(:,:,i) = images(:,:,i) - ref; % Delta R
  dimg(:,:,i) = (single(images(:,:,i)) - ref)./ref; % (Delta R)/R
end
