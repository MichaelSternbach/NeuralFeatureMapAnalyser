function BloodVesselImg = getBloodVesselImg(data_path,ID)
    fig = openfig([data_path ID '.fig']);
    BloodVesselImg = getimage(fig);
end