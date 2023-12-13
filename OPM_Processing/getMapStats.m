function [ORI_CI_pix,PwProb,PwStd] = getMapStats(orientation_stats,pinwheel_stats)

    %% get ORI CI
    ORI_CI_pix = getORI_CI_pix(orientation_stats,ROI);

    %% get Pw Prob
    PwProb = pinwheel_stats.probability;

    %% get PwStd
    PwStd = getPwStd(pinwheel_stats);

end