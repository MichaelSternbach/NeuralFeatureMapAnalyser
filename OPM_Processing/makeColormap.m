function cm = makeColormap(type,num_bins)

if nargin < 2
    num_bins = 12;
end
if nargin < 1
    type = 'default';
end

try
    cm = feval(type,num_bins);
catch
    error('Colormap not defined')
end

if size(cm,1) ~= num_bins
    cm = interpolateColorMap(cm,num_bins,false);
end

end

function cm = interpolateColorMap(cm,num_bins,isCircular)

if ~isCircular
    cm = interp1(linspace(0,1,size(cm,1)),cm,linspace(0,1,num_bins));
else
    cm = [cm;cm(1,:)];
    cm = interp1(linspace(0,1,size(cm,1)),cm,linspace(0,1,num_bins+1));
    cm(end,:) = [];
end

end

function cm = orientation(num_bins)

cm=[0         1.0000         0;
    0.1569    1.0000         0;
    0.5490    1.0000         0;
    0.8235    1.0000         0;
    1.0000    1.0000         0;
    1.0000    0.8235         0;
    1.0000    0.5882         0;
    1.0000    0.3529         0;
    1.0000         0         0;
    1.0000         0         0.4784;
    0.7059         0         1.0000;
    0.3922         0         1.0000;
    0              0         1.0000;
    0              0.3137    1.0000;
    0              0.5490    0.5882;
    0              0.7843    0.1961];

cm = circshift(cm(end:-1:1,:),ceil(size(cm,1)*0.78),1);

cm = interpolateColorMap(cm,num_bins,true);

end


function cm = selectivity(num_bins)

% intensity for selectivity
min_val=0.3;
slope=1.2; % 1.2
center=1.8; % 1.5

mm = (tanh(slope*linspace(-2.5,2.5,num_bins)+center)+1)/2;
mm = (mm + min_val/(1-min_val) )/( 1 + min_val/(1-min_val));

cm = repmat(mm(:),[1 3]);

end


function cm = blueyellowred(~)
% color map from http://geog.uoregon.edu/datagraphics/color_scales.htm
cm = ...
    [41   10  216;
    38   77  255;
    63  160  255;
    114  217  255;
    170  247  255;
    224  255  255;
    255  255  191;
    255  224  153;
    255  173  114;
    247  109  94;
    216   38  50;
    165    0  33]/255;
end

function cm = bluered(~)
% color map from http://colorbrewer.org.
cm =...
    [103	0       31;
    178     24      43;
    214     96      77;
    244     165     130;
    253     219     199;
    209     229     240;
    146     197     222;
    67	    147     195;
    33      102     172;
    5       48      97]/255;

cm = flipud(cm);

end

function cm = bluewhitered(~)
cm = ...
    [0        0         0.5000;
    0         0.1667    0.6667;
    0         0.3333    0.8333;
    0         0.5000    1.0000;
    0.3333    0.6667    1.0000;
    0.6667    0.8333    1.0000;
    1.0000    1.0000    1.0000;
    1.0000    1.0000    1.0000;
    1.0000    0.6000    0.6000;
    1.0000    0.2000    0.2000;
    0.9000    0         0;
    0.7000    0         0;
    0.5000    0         0];
end

