function aimg = oiAve(images)
%OIAVE Calcuate average image sequence.
%
%  AIMG = OIAVE(IMAGES) Calculates the average image sequence by
%  summing across trials and dividing by the number of trials.
%
%  IMAGES is an MxN cell array of PxQxR matricies where:
%
%    M is the number of conditions
%    N is the number of trials
%
%    P is the height of a single frame (in pixels)
%    Q is the width os a single frame (in pixels)
%    R is the number of frames acquired for each image sequence

% $Id: oiAve.m,v 1.1 2008-04-30 04:56:07 shaunc Exp $

assert(nargin >= 1, ...
  'Too few input arguments. Type "help oiAve" for usage information.');

assert(iscell(images), ...
  'The input argument "images" does not appear to be a cell array. Type "help oiAve" for usage information.');

% [numConds, numTrials] = size(images);
numConds = size(images,1);

aimg = images(:,1); % the first trial...
for i = 1:numConds,
  numTrials = sum(cellfun(@(x) ~isempty(x), images(i,:)));
  
  for j = 2:numTrials,
    aimg{i,1} = aimg{i,1} + images{i,j};
  end
  aimg{i,1} = aimg{i,1}./numTrials;
end
