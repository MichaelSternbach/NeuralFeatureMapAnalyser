function SNR = calcSNR_OPM_Data(data_obj,doFilter,direction_map)

    if nargin < 2
        doFilter = false;
    end
    if nargin <3
    	direction_map = false;
    end

    %% get mean map and noise
    [DiffMaps,MeanMap,ROI]= getDifferenceMaps(data_obj,1,'complex',doFilter,direction_map);
    %[BottstapSampleMaps,MeanMap,ROI]= getBootstrapSampleMaps(data_obj,scale,doFilter,direction_map);

    %% get power MeanMap
    signal_power = mean(MeanMap(ROI).*conj(MeanMap(ROI)),'all');
    noise_power = mean(DiffMaps(ROI).*conj(DiffMaps(ROI)),'all');

    %% get SNR
    SNR = signal_power/noise_power;
end
