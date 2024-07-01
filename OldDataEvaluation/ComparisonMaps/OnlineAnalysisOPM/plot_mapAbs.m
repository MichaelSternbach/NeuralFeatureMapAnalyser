function plot_mapAbs(map,Title,maxMap,minMap)
    %figure;
    a = map/(maxMap-minMap);
    imagesc(a); 
    %colormap jet;
    colormap turbo;
    %colormap gray;
    %colormap hot;
    
    %set(gca,'xtick',[])
    %set(gca,'ytick',[])
    title(Title)
    %pbaspect([1 1 1])
    hold on
%     colormap(gray)
    %colorbar('Ticks',[0,0.5,3])
    cbh = colorbar ; %Create Colorbar
    NTicks = 8;
%     axis equal
    set(gca,'YDir','normal')
    set(gca,'xdir','reverse')
    axis image
    cbh.Ticks = linspace(0, 1, NTicks) ; %Create 8 ticks from zero to 1
    cbh.TickLabels = num2cell(linspace(minMap,maxMap,NTicks)) ;
    %cbh.TickLabels = num2cell(linspace(-1,1,NTicks)) ;
end
