function selectivities_pw = SumSelectivitiesPw(pinwheel_stats,data_obj)
    
    selectivities_pw = zeros([1 size(data_obj.samples_array,3)]);
    for ii = 2:size(data_obj.samples_array,3)
        S = 0;
        z = data_obj.filter_map(data_obj.read_map(ii));
        z = (z-mean(z,'all'))./std(z);
        for jj = 1:size(pinwheel_stats.x,1)
            x = round(pinwheel_stats.x(jj,1));
            y = round(pinwheel_stats.y(jj,1));
            S = S + abs(z(x,y))^2;
        end
        selectivities_pw(ii) = S/size(pinwheel_stats.x,1);
    end
    
end
