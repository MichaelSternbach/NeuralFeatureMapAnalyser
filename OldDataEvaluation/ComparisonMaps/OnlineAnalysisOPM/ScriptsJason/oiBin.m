function img = oiBin(images, sbin, varargin),
%OIBIN Spatial binning of images.
%
%  IMG = OIBIN(IMAGES, SBIN) Performs spatial binning of the supplied
%  image sequence(s).
%
%  IMG = OIBIN(IMAGES, SBIN, TBIN) Performs binning in both space and time. To
%  perform binning in time only, pass an empty matrix as SBIN.
%
%  IMAGES is a PxQxR matrix containing a single image sequence where
%
%    P is the height of a single frame (in pixels)
%    Q is the width os a single frame (in pixels)
%    R is the number of frames acquired in the sequence
%
%  SBIN is a 1x2 (or 2x1) vector specifying the spatial binning factors
%  [width,height] (default: [], no spatial binning)
%
%  If present, TBIN is an integer specifying the temporal binning factor
%  (default: 1, no temporal binning).
%
%  Alternatively, IMAGES may be an MxN cell array of image sequences where:
%
%    M is the number of conditions
%    N is the number of trials
%
%  in which case, binning is calculated for each individual image
%  sequence. 

% $Id: oiBin.m,v 1.1 2008-04-30 04:56:07 shaunc Exp $

assert(nargin >= 2, ...
  'Too few input arguments. Type "help oiBin" for usage information.');

if iscell(images),
  % call oiBin( ) recursively for each image sequence

  [numConds, numTrials] = size(images);

  img = cell(size(images));
  for i = 1:numConds,
    for j = 1:numTrials,
      img{i,j} = oiBin(images{i,j}, sbin);
    end
  end
  
  return
end

if (isempty(sbin)),
  sbin = [1, 1];
end

% bin vertically
if sbin(1) > 1,
  newsize = floor(size(images,1)/sbin(1));
  img = images(1:sbin(1):sbin(1)*newsize,:,:);
  for i = 2:sbin(1),
    img = img + images(i:sbin(1):sbin(1)*newsize,:,:);
  end
  img = img / sbin(1);
end

% bin horizontally
if sbin(2) > 1,
  if sbin(1) > 1,
    images = img;
  end
  newsize = floor(size(images,2)/sbin(2));
  img = images(:,1:sbin(2):sbin(2)*newsize,:);
  for i = 2:sbin(2),
    img = img + images(:,i:sbin(2):sbin(2)*newsize,:);
  end
  img = img / sbin(2);
end

% bin in time?
if nargin > 2,
  tbin = varargin{1};
  if tbin > 1,
    if (sum(sbin) > 2)
      images = img;
    end
    newsize = floor(size(images,3)/tbin);
    img = images(:,:,1:tbin:tbin*newsize);
    for i = 2:tbin,
      img = img + images(:,:,i:tbin:tbin*newsize);
    end
    img = img / tbin;
  end
end
