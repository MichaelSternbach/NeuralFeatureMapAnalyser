function [conds, idx] = sbStimImport(fid),

% $Id: $

tmp = [];
while ~feof(fid), ln = fgetl(fid); tmp = [tmp; str2num(ln)]; end
tmp = tmp + 1;

conds = unique(tmp); 

idx = zeros(size(tmp));
for i = 1:length(idx),
  idx(i) = find(conds == tmp(i));
end

conds = [conds(2:end); NaN]'; % NaN is the control/blank condition
idx = idx - 1;
idx(find(idx == 0)) = length(conds);
