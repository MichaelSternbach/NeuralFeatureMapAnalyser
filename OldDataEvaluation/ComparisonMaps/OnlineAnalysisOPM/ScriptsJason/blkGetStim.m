function data = blkGetStim(fid, stimId),
%BLKGETSTIM Load the imaging data for the given stimulus condition.
%
%  DATA = BLKGETSTIM(FID, STIMID) reads and returns the image sequence(s)
%  from the specified file handle.
%
%  FID the file handle for the Vdaq .blk file.
%
%  STIMID the stimulus Ids for which data should be loaded. If STIMID is
%  absent, image sequences for *all* available stimulus Ids are returned.

% $Id: blkGetStim.m,v 1.1 2008/04/30 04:56:07 shaunc Exp $

h = blkFileHeader(fid);

if isempty(h),
  data = {}; % empty cell array
  return;
end;

if nargin < 2,
  stimId = h.stimIds;
end
%idx = find(h.stimIds == stimId);
idx = arrayfun(@(x)(find(h.stimIds == x)), stimId, 'UniformOutput', true);

numStim = length(idx);

fseek(fid, h.headerLen, 'bof');
ptr = ftell(fid);

for i = 1:numStim,
  fseek(fid, ptr+(idx(i)-1)*h.stimSize, 'bof');
  for j = 1:h.framesPerStim,
%    tmp = fread(fid, h.frameWidth*h.frameHeight, 'uint16');
    tmp = int16(fread(fid, h.frameWidth*h.frameHeight, h.dataPrecision));
%    tmp = single(fread(fid, h.frameWidth*h.frameHeight, h.dataPrecision));
      
    data{i,1}(:,:,j) = reshape(tmp, h.frameWidth, h.frameHeight)';
  end
end
