function PlotBloodVessels(BloodVesselImg,ROI,transp)
    colormap gray;
    im = imagesc(ROI.*BloodVesselImg/max(BloodVesselImg,[],'all'));
    im.AlphaData = transp;
    %axis equal
    set(gca,'YDir','normal')
    axis image
end