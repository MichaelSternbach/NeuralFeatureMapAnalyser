function selectivities_pw = getSelectivitiesHighestProbPw(pinwheel_stats,data_obj)
    %% find position of the pinwheel with the highest Bootstrap Probability
    [x_high_prob_pw,y_high_prob_pw]=getPositionHighestProbPw(pinwheel_stats);
    %% calculate the selectivity at that position in all bootstrap samples 
    selectivities_pw = getSelectivitiesPw(x_high_prob_pw,y_high_prob_pw,data_obj);

end

function [x_high_prob_pw,y_high_prob_pw]=getPositionHighestProbPw(pinwheel_stats)
    [~,high_prob_pw]=max(pinwheel_stats.probability);
    x_high_prob_pw = round(pinwheel_stats.x(high_prob_pw,1));
    y_high_prob_pw = round(pinwheel_stats.y(high_prob_pw,1));
end

function selectivity_pw = getSelectivitiesPw(x,y,data_obj)
    selectivity_pw = zeros([1 size(data_obj.samples_array,3)]);
    for ii = 2:size(data_obj.samples_array,3)
        z = data_obj.filter_map(data_obj.read_map(ii));
        selectivity_pw(ii) = abs(z(x,y))^2;
    end
end