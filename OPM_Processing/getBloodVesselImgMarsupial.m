function BloodVesselImg = getBloodVesselImgMarsupial(data_path,ID)
    fig = openfig([data_path ID '.fig']);
    BloodVesselImg = getimage(fig);
    close(fig)
end

