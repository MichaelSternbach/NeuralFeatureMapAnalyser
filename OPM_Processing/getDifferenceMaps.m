function [DiffMaps,ROI]= getDifferenceMaps(data_obj,scale)
    num_boot_samples = size(data_obj.samples_array,3);
    MeanMap = data_obj.filter_map(data_obj.read_map(1));

    if scale == 1
        DiffMaps = zeros([size(MeanMap) num_boot_samples-1]);
    else
        DiffMaps = zeros([size(imresize(MeanMap,scale)) num_boot_samples-1]);
    end

    for ii = 1:num_boot_samples-1
        if scale == 1
            DiffMaps(:,:,ii) = data_obj.filter_map(data_obj.read_map(ii+1))-MeanMap;
        else
            DiffMaps(:,:,ii) = imresize(data_obj.filter_map(data_obj.read_map(ii+1))-MeanMap,scale);
        end
    end

    if scale == 1
        ROI = data_obj.ROI;
    else
        ROI = imresize(data_obj.ROI,scale);
    end
end