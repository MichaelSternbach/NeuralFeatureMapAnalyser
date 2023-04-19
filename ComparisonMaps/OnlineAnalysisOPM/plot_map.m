function varargout=plot_map(data,ROI,ref_sel,black_roi,smoothing,color_type)
% [h,fig2plot]=plot_map(data,ROI,ref_sel,black_roi,smoothing,color_type)
% Function to make a color and intensity colored orientation preference
% display
% input = polar data, Region of interest
% mapPolar(data2D,ROI)
% created by chepe@nld.ds.mpg.de

% Other functions used:
% - makeColormap
% - smooth_pf, negative_smoothing, imresize


%% ---------------------  read parameters

data(isnan(data)) = 0;

% if empty ROI, use all image
if ~exist('ROI','var') || isempty(ROI)
    ROI=true(size(data));
end

% determine polar scaling for selectivity
if ~exist('ref_sel','var')
    ref_sel=3*sqrt(mean(abs(data(ROI)).^2));
elseif isempty(ref_sel)
    ref_sel=3*sqrt(mean(abs(data(ROI)).^2));
end

if ~exist('black_roi','var')
    black_roi = 0;
end

% increase size of display for exporting figures
if exist('smoothing','var') && ~isempty(smoothing)
    if smoothing>0
        for times=1:smoothing
            data=smooth_pf(data);
        end
    elseif smoothing<0
        for times=1:abs(smoothing)
            data = negative_smoothing(data);
        end
    end
    ROI = imresize(ROI,2^abs(smoothing));
end

%-  make color maps
if ~exist('color_type','var')
    color_type = 'std';
end
switch color_type
    case 'std'
        cm = makeColormap('orientation',16);        
    case 'interp'
        cm = makeColormap('orientation',24);
    case 'circ'
        cm = makeColormap('circular',50);        
end
mm = makeColormap('selectivity',64);
%% --------------------------- Asign color values

% -- ORIENTATION

% normalize from 0-1
ori=(angle(data))/(2*pi);ori(ori<0)=1+ori(ori<0);
ori=ori(:);

% find to which color interval it belongs
intervals=linspace(0,1,2*size(cm,1)+1);
[~,ori] = histc(ori,intervals(2:2:end));
ori = ori+1;
ori(ori==size(cm,1)+1)=1;

% assign
anglePlot = (cm(ori(:),:));

% -- SELECTIVITY

if ref_sel==0
    
    % if polar plot is not required
    magnitudePlot = 0.9*ones(size(anglePlot));
    
else
    
    % threshold
    sel=abs(data)/ref_sel;
    sel(sel>1)=1;
    
    % find intensity match
    sel=ceil(size(mm,1)*sel(:));
    sel(sel==0)=1;
    
    % assign
    magnitudePlot = (mm(sel(:),:));
end

%% --------------------------- Combine

if ismatrix(data)
    magnitudePlot = reshape(magnitudePlot, [size(data,1),size(data,2),3]);
    anglePlot = reshape(anglePlot, [size(data,1),size(data,2),3]);
    fig2plot=anglePlot.*magnitudePlot;
    %fig2plot=repmat(double(ROI),[1 1 3]).*anglePlot.*magnitudePlot;
    if black_roi == 1
        % black ROI
        fig2plot(repmat(~ROI,[1 1 3])) = 0;
    elseif black_roi == -1
        % white ROI
        fig2plot(repmat(~ROI,[1 1 3])) = 1;
    end
else
    % in case of volume image, return vector
    fig2plot = magnitudePlot.*anglePlot;
end

% if only the colormap is required, return without plotting
if nargout==2 || ~ismatrix(data)
    varargout{1} = [];
    varargout{2} = fig2plot;   
    return
end

%% --------------------------- Make figure

h = image(fig2plot);
set(gca,'ydir','reverse')
axis image

% set(gca,'xtick',5:5:size(data,2))
% set(gca,'ytick',5:5:size(data,1))

% set(gca,'xticklabel',[])
% set(gca,'yticklabel',[])

if nargout==1
    varargout{1} = h;    
end

end
