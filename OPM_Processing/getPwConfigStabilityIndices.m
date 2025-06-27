function StabilityIndicies=getPwConfigStabilityIndices(data_obj,num_bootstrapsamples)
    %% get randomize data obj
    data_rand = randomizeData(data_obj.data);
    data_obj_rand =  data_handle_corrected(data_obj.info,data_rand,data_obj.ROI);

    %% find pinwheel positions x,y
    z = data_obj.filter_map(data_obj.read_map());
    ROI = data_obj.ROI;
    [~,~,~,PWxList,PWyList,~,~] = find_pinwheels(z,0,ROI);

    %% get stability index data obj
    StabilityIndicies.data = calStabilityIndicies(data_obj,num_bootstrapsamples,PWxList,PWyList);

    %% get stability index data obj
    StabilityIndicies.data_rand = calStabilityIndicies(data_obj_rand,num_bootstrapsamples,PWxList,PWyList);
end

function StabilityIndicies=calStabilityIndicies(data_obj,num_bootstrapsamples,PWxList,PWyList)
    
    %% set bootstrap samples
    data_obj.prepare_samples_array(num_bootstrapsamples);
    
    %% get stability index data obj
    z = data_obj.filter_map(data_obj.read_map());
    ROI = data_obj.ROI;

    % %% find pinwheel positions x,y
    % [~,~,~,PWxList,PWyList,~,~] = find_pinwheels(z,0,ROI);

    %% get selectivities for all pinwheels
    selectivity_data_pw = zeros([num_bootstrapsamples-1,length(PWxList)]);
    selectivity_data_all = zeros([num_bootstrapsamples-1,sum(ROI(:))]);
    for ii = 2:num_bootstrapsamples
        z = data_obj.filter_map(data_obj.read_map(ii));
        selectivity_data_pw(ii,:) = abs(z(sub2ind(size(z),round(PWyList),round(PWxList))));
        selectivity_data_all(ii,:) = abs(z(ROI));
    end

    %% get full variance
    V = nanvar(selectivity_data_all,0,'all');

    %% calc stability indicies
    StabilityIndicies = (V - nanvar(selectivity_data_pw,0,1))/V;


    
end