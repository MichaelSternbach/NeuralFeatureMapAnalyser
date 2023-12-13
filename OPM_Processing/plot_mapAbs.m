function plot_mapAbs(map,Title,maxMap,minMap,ROI,ax)
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
    if nargin <= 2
        imagesc(ax,a);
    else
        a(find(~ROI))=minMap;
        imagesc(ax,a,[minMap maxMap]);
    end
    %% details Plot
    %colormap jet;
    %colormap turbo;
    colormap(ax, 'jet');%gray
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
