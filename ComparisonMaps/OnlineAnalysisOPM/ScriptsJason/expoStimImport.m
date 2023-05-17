function [conds, idx] = expoStimImport(fid),

% $Id: blkImport.m,v 1.1 2008/04/30 04:56:07 shaunc Exp $

fseek(fid, 0, 'bof');
jnk = fgetl(fid); % column headings are the first line of the file

pat = '(?<contrast>\d+)\s+(?<direction>\d+)\s+\d+\s+(?<index>\d+)\s+.*';

contrast = []; direction = []; index = [];
while ~feof(fid),
  ln = fgetl(fid);
  
  tmp = regexp(ln, pat, 'names');
    
  contrast = [contrast; str2double(tmp.contrast)];
  direction = [direction; str2double(tmp.direction)];
  index = [index; str2double(tmp.index)];
end

[index, i] = sort(index);
contrast = contrast(i);
direction = direction(i);

contrast = 2*contrast - 1;

tmp = contrast.*direction; % directions, -ve's are the control/blank condition

conds = unique(tmp); 

idx = zeros(size(index));
for i = 1:length(idx),
  idx(i) = find(conds == tmp(i));
end

conds = [conds(2:end); NaN]'; % NaN is the control/blank condition
idx = idx - 1;
idx(find(idx == 0)) = length(conds);
