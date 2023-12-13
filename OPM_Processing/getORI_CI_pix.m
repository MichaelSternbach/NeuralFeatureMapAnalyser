function ORI_CI_pix = getORI_CI_pix(ORI_Stats,ROI)
    ORI_CI = getORI_CI(ORI_Stats);
    ORI_CI_pix = ORI_CI(ROI);
end