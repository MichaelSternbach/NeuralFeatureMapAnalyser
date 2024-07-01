function plotContourAngleDelta(AngleDelta,AllMaps,ConfidenceBorder)
    addpath ~/Cloud/git/vone/MatlabCode/PatchAnalysis
    plot_map(AllMaps(:,:,1))
    hold on
    contour(ConfidenceBorder,[AngleDelta,AngleDelta],'white')%X,Y,
    %pbaspect([1 1 1])
    title(['Countour line uncertainty orientation ',num2str(AngleDelta),'°'])
end
