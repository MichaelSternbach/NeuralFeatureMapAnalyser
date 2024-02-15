function PlotBloodVessels(BloodVesselImg,ROI,transp)
    colormap gray;
    map = BloodVesselImg/max(BloodVesselImg,[],'all');
    map(~ROI)=min(map,[],'all');
    im = imagesc(map);
    im.AlphaData = transp;
    %axis equal
    axis image
end