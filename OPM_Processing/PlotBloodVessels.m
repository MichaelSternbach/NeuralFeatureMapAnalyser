function PlotBloodVessels(BloodVesselImg,ROI,transp)
    colormap gray;
    im = imagesc(ROI.*BloodVesselImg/max(BloodVesselImg,[],'all'));
    im.AlphaData = transp;
    %axis equal
    axis image
end