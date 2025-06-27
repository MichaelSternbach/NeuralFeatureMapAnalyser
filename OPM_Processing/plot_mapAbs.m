function cbh = plot_mapAbs(map,Title,maxMap,minMap,ROI,ax)
    %% argin
    if nargin==1
        Title = '';
    end
%     if nargin <= 2
%         maxMap = max(map,[],'all');
%         minMap = min(map,[],'all');
%     end
    if (nargin <= 4)
        ROI = ones(size(map));
    end
    if (ROI==1)
        ROI = ones(size(map));
    end
    
%     %% norm map
%     a = (map-minMap)/(maxMap-minMap);
    a = map;   
    
    
    %% plot map
    if nargin < 6
        ax = axes;
    end
    if (nargin <= 2) || (isnan(maxMap)||isnan(minMap))
        imagesc(ax,a);
    else
        a(~ROI)=minMap;
        imagesc(ax,a,[minMap maxMap]);
    end
    %% details Plot
    %colormap jet;
    %colormap turbo;
    m=100;
    cm_viridis=viridis(m);
    colormap(ax, cm_viridis);%gray,'jet'
    %colormap hot;
    
    set(gca,'xtick',[])
    set(gca,'ytick',[])
    title(ax,Title)%, 'interpreter', 'latex'
    %pbaspect([1 1 1])
    
    hold on
    %colorbar('Ticks',[0,0.5,3])
    cbh = colorbar(ax) ; %Create Colorbar
%     axis equal
    axis image
%     cbh.Ticks = linspace(0, 1, NTicks) ; %Create 8 ticks from zero to 1
%     cbh.TickLabels = num2cell(linspace(minMap,maxMap,NTicks)) ;
    %cbh.TickLabels = num2cell(linspace(-1,1,NTicks)) ;
end
