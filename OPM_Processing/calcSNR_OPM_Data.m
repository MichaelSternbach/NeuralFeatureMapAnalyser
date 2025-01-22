function SNR = calcSNR_OPM_Data(data_obj,dofilter)

    if nargin < 2
        dofilter = true;
    end

    %% get mean map and noise
    [DiffMaps,MeanMap,ROI]= getDifferenceMaps(data_obj,1,'complex',dofilter);

    %% get power MeanMap
    signal_power = mean(MeanMap(ROI).*conj(MeanMap(ROI)),'all');
    noise_power = mean(DiffMaps(ROI).*conj(DiffMaps(ROI)),'all');

    %% get SNR
    SNR = signal_power/noise_power;
end