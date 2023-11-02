function plotContourCI(contourValue,Map,ROI,CI_Map,CIType,unit)
    plot_map(Map,ROI,0,1)
    hold on
    contour(CI_Map,[contourValue,contourValue],'black')%X,Y,
    %pbaspect([1 1 1])
    title(['Countour line uncertainty ' CIType ' ' num2str(contourValue) unit])
end