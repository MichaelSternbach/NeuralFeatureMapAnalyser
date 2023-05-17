function ids = blkStimIds(fid),
%BLKSTIMIDS Returns a list of stimulus Ids.
%
%  H = BLKSTIMIDS(FID) reads and returns the list of stimulus Ids
%  from the specified file handle.

% $Id: blkStimIds.m,v 1.1 2008-04-30 04:56:07 shaunc Exp $

h = blkFileHeader(fid);

if isempty(h),
  ids = [];
  return;
end;

ids = h.stimIds;
