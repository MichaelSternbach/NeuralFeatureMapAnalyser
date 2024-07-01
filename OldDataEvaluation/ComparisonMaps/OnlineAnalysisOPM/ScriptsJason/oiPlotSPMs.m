function h = oiPlotSPMs(Z, thresh, grnImg),

error(nargchk(3,3,nargin));

% find max/min across all conditions
mn = min(cellfun(@(x) min(x(:)), Z));
mx = max(cellfun(@(x) max(x(:)), Z));

%%
tmp = round(max(abs([mx, mn]))+0.5);
mn = sign(mn)*tmp;
mx = sign(mx)*tmp;
%%

[m,n] = size(Z{1,1});

for condId = 1:length(Z),
  h(condId) = figure;
  
  % plot the green image
%   subplot(1,3,1);
  
  base = ind2rgb(uint16(grnImg), gray(256));
%   subimage(base);
  
  % 1px = 24um (sensor pixels are binned 2x2)
%   axis image
%   set(gca,'XTick',[0:20.83:n]);
%   set(gca,'XTickLabel',cellfun(@(x) sprintf('%.1f',x), num2cell(0.024*[0:20.83:n]), 'UniformOutput', false));
%   set(gca,'YTick',[0:20.83:m]);
%   set(gca,'YTickLabel',cellfun(@(x) sprintf('%.1f',x), num2cell(0.024*[0:20.83:m]), 'UniformOutput', false));
%   
%   xlabel('mm');
%   ylabel('mm');
  
  % plot the composite (i.e., green image with overlayed SPM{Z})
%   subplot(1,3,2);
  subplot(1,2,2);
 
  tmp = Z{condId,1}; 
  tmp = 256*(tmp - mn)/(mx-mn); % FIXME: think about scaling, the Z scores are approximately normally distributed...?
  
%  base = ind2rgb(uint16(grnImg), gray(256));
  overlay = ind2rgb(uint16(tmp), jet(256)); %circshift(hsv(256),[0, 128]));

  alphaMap = repmat(abs(Z{condId,1}) > thresh,[1,1,3]);
  tmp = (1-alphaMap).*base + alphaMap.*overlay;

  subimage(tmp);

  % 1px = 24um (sensor pixels are binned 2x2)
  axis image
  set(gca,'XTick',[0:20.83:n]);
  set(gca,'XTickLabel',cellfun(@(x) sprintf('%.1f',x), num2cell(0.024*[0:20.83:n]), 'UniformOutput', false));
  set(gca,'YTick',[0:20.83:m]);
  set(gca,'YTickLabel',cellfun(@(x) sprintf('%.1f',x), num2cell(0.024*[0:20.83:m]), 'UniformOutput', false));
   
  xlabel('mm');
  ylabel('mm');
  
  % plot SPM{Z}
%   subplot(1,3,3);
  subplot(1,2,1);
  
  subimage(overlay);
  
  % 1px = 24um (sensor pixels are binned 2x2)
  axis image
  set(gca,'XTick',[0:20.83:n]);
  set(gca,'XTickLabel',cellfun(@(x) sprintf('%.1f',x), num2cell(0.024*[0:20.83:n]), 'UniformOutput', false));
  set(gca,'YTick',[0:20.83:m]);
  set(gca,'YTickLabel',cellfun(@(x) sprintf('%.1f',x), num2cell(0.024*[0:20.83:m]), 'UniformOutput', false));
  
  xlabel('mm');
  ylabel('mm');
end