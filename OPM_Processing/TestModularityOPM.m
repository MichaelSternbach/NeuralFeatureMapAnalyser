function [power_profiles,mean_abs_squared]=TestModularityOPM(data_obj,profile_range_mm,Jackknife)
        
    %% Input parameter
    if nargin < 2
        profile_range_mm = 0.1:0.1:1.5;
    end
    if nargin <3
        Jackknife = false; 
    end

    %% determine power at peak for bootstrap samples
    BS = cell([1 size(data_obj.samples_array,3)]);
    mean_abs_squared = zeros([1 size(data_obj.samples_array,3)]);

    for ii = 1:size(data_obj.samples_array,3)
        z = data_obj.read_map(ii);
        mean_abs_squared(ii) = mean(abs(z).^2,'all');
        BS{ii} = define_filter_settings(data_obj.info,data_obj.ROI,z,profile_range_mm);
    end
    power_profiles.BS=BS;
    
    if Jackknife
        %% determine power at peak for Jackknife
        data_obj.prepare_jackknife_samples;
        JK = cell([1 data_obj.data_parameters.num_blocks]);

        parfor ii=1:data_obj.data_parameters.num_blocks
            z = data_obj.read_map(ii);
            JK{ii} = define_filter_settings(data_obj.info,data_obj.ROI,z,profile_range_mm,profile_step_mm);
        end
        power_profiles.JK=JK;
    end
end


function [value_peak,ii_peak]=findMaxPeak(x)
    %% determine peaks
    [values_peaks,ii_peaks]=findPeaks(x);
    
    %% find max peak
    [value_peak,ii]=max(values_peaks);
    ii_peak = ii_peaks(ii);

end

function [values_peaks,ii_peaks]=findPeaks(x)
    values_peaks = [];
    ii_peaks = [];
    for ii = 2:length(x)-1
        if isPeak(x(ii-1),x(ii),x(ii+1))
            ii_peaks = [ii_peaks ii];
            values_peaks = [values_peaks x(ii)];
        end
    end
end

function b = isPeak(x1,x2,x3)
    b = x1<x2 && x3<x2;
end
