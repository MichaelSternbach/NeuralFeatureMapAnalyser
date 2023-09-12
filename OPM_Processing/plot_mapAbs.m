function plot_mapAbs(map,Title,maxMap,minMap,ROI,NTicks)
    %% argin
    if nargin == 2
        maxMap = max(map,[],'all');
        minMap = min(map,[],'all');
    end
    if nargin <= 4
        ROI = ones(size(map));
    end
    if nargin <= 5
        NTicks = 8;
    end
    
    %% norm map
    a = (map-minMap)/(maxMap-minMap);
    a=a.*ROI;
    
    %% plot map
    imagesc(a); 
    
    %% details Plot
    %colormap jet;
    colormap turbo;
    %colormap gray;
    %colormap hot;
    
    set(gca,'xtick',[])
    set(gca,'ytick',[])
    title(Title)%, 'interpreter', 'latex'
    %pbaspect([1 1 1])
    
    hold on
    %colorbar('Ticks',[0,0.5,3])
    cbh = colorbar ; %Create Colorbar
%     axis equal
    axis image
    cbh.Ticks = linspace(0, 1, NTicks) ; %Create 8 ticks from zero to 1
    cbh.TickLabels = num2cell(linspace(minMap,maxMap,NTicks)) ;
    %cbh.TickLabels = num2cell(linspace(-1,1,NTicks)) ;
end
