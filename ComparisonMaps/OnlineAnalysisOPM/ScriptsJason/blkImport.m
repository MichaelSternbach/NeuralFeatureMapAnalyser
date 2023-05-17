function data = blkImport(fid),
%BLKImport Load the imaging data from the supplied file handle(s).
%
%  DATA = BLKIMPORT(FID) reads and returns the image sequence(s) for *all*
%  stimulus Ids from the supplied file handle(s).
%
%  FID a file handle or a vector of file handles for the Vdaq .BLK file(s).
%
%  DATA is an MxN cell array of PxQxR matricies where
%
%    M is the number of conditions
%    N is the number of trials
%
%    P is the height of a single frame (in pixels)
%    Q is the width of a single frame (in pixels)
%    R is the number of frames acquired for each stimulus
%
%    So, DATA{i,j}(:,:,k) contains the k'th frame for the j'th trial of
%    the i'th stimulus Id.

% 27/6/2012 - Shaun L Cloherty <s.cloherty@ieee.org>

data = cell([1,1]);
for i = 1:length(fid),
  stimIds = blkStimIds(fid(i));
  
  data_ = blkGetStim(fid(i));
  
  numConds = length(stimIds);
  
  cnt = ones([numConds,1]);
  for j = 1:numConds,
    if size(data,1) >= stimIds(j),
      cnt(j) = sum(cellfun(@(x) ~isempty(x), data(stimIds(j),:))) + 1;
    end
    data(stimIds(j),cnt(j)) = data_(j);
  end
end
