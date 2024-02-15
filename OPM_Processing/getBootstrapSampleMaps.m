
function [BottstapSampleMaps,MeanMap,ROI]= getBootstrapSampleMaps(data_obj,scale,DoFilter)
    num_boot_samples = size(data_obj.samples_array,3);
    
    %% convert scale if given in pixel
    scale = convert_scale(scale,data_obj);
        
    
    if scale == 1
        if DoFilter
            MeanMap = data_obj.filter_map(data_obj.read_map(1));
        else
            MeanMap = data_obj.normalize_map(1);
        end
    else
        if DoFilter
            MeanMap = imresize(data_obj.filter_map(data_obj.read_map(1)),scale);
        else
            MeanMap = imresize(data_obj.normalize_map(1),scale);
        end
    end
    
    BottstapSampleMaps = zeros([size(MeanMap) num_boot_samples-1]);


    for ii = 1:num_boot_samples-1
        if scale == 1
            if DoFilter
                BottstapSampleMaps(:,:,ii) = data_obj.filter_map(data_obj.read_map(ii+1));
            else
                BottstapSampleMaps(:,:,ii) = data_obj.normalize_map(ii+1);
            end
        else
            if DoFilter
                BottstapSampleMaps(:,:,ii) = imresize(data_obj.filter_map(data_obj.read_map(ii+1)),scale);
            else
                BottstapSampleMaps(:,:,ii) = imresize(data_obj.normalize_map(ii+1),scale);
            end
        end
    end

    if scale == 1
        ROI = data_obj.ROI;
    else
        ROI = imresize(data_obj.ROI,scale);
    end
end