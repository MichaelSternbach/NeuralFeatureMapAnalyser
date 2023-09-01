function plot_mapAbs(map,Title,maxMap,minMap,ROI)
    %figure;
    if nargin == 4
        ROI = ones(size(map));
    elseif nargin == 2
        maxMap = max(map,[],'all');
        minMap = min(map,[],'all');
    end
    a = (map-minMap)/(maxMap-minMap);
    a=a.*ROI;
    imagesc(a); 
    
    %colormap jet;
    colormap turbo;
    %colormap gray;
    %colormap hot;
    
%     set(gca,'xtick',[])
%     set(gca,'ytick',[])
    title(Title, 'interpreter', 'latex')
    %pbaspect([1 1 1])
    
    hold on
    %colorbar('Ticks',[0,0.5,3])
    cbh = colorbar ; %Create Colorbar
    NTicks = 8;
%     axis equal
    axis image
    cbh.Ticks = linspace(0, 1, NTicks) ; %Create 8 ticks from zero to 1
    cbh.TickLabels = num2cell(linspace(minMap,maxMap,NTicks)) ;
    %cbh.TickLabels = num2cell(linspace(-1,1,NTicks)) ;
end
