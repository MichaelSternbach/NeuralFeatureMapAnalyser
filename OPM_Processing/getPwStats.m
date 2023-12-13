function [PwProb,PwStd] = getPwStats(pinwheel_stats)

    %% get Pw Prob
    PwProb = pinwheel_stats.probability;

    %% get PwStd
    PwStd = getPwStd(pinwheel_stats);

end