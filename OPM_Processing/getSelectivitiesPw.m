function selectivity_pw = getSelectivitiesPw(x,y,data_obj)
    selectivity_pw = zeros([1 size(data_obj.samples_array,3)]);
    for ii = 2:size(data_obj.samples_array,3)
        z = data_obj.filter_map(data_obj.read_map(ii));
        selectivity_pw(ii) = abs(z(y,x))^2;
    end
end