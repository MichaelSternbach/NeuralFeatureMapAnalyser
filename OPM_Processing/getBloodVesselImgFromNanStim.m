function BloodVesselImg = getBloodVesselImgFromNanStim(data,stim_order)
    NanStim = find(isnan(stim_order));
    BloodVesselImg = mean(data(:,:,NanStim,:),4);
end