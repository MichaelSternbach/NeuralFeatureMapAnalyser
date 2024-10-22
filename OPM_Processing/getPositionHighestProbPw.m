function [x_high_prob_pw,y_high_prob_pw]=getPositionHighestProbPw(pinwheel_stats)
    [~,high_prob_pw]=max(pinwheel_stats.probability);
    x_high_prob_pw = round(pinwheel_stats.x(high_prob_pw,1));
    y_high_prob_pw = round(pinwheel_stats.y(high_prob_pw,1));
end

