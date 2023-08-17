function CA_sizes = getPinwheelConfidenceAreaSizes(pinwheel_stats,field_size_pix)
    CA_sizes = zeros([1 getN_PW(pinwheel_stats)]);
    for i_pw = 1:getN_PW(pinwheel_stats)
        CA_sizes(i_pw) = getCA_size(pinwheel_stats.x(i_pw,:),pinwheel_stats.y(i_pw,:),field_size_pix);
    end
end

function CA_size = getCA_size(PWx,PWy,field_size_pix)
    CI = points_confidence_region(PWx,PWy,field_size_pix,'hull');
    CA_size = sum(CI,'all');
end


function N_PW = getN_PW(pinwheel_stats)
    N_PW = size(pinwheel_stats.x,1);
end