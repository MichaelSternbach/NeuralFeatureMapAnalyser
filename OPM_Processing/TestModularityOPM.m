function TestModularityOPM(data_obj,data_obj_rand,profile_range_mm,profile_step_mm)
        
    %% Input parameter
    if nargin < 2
        profile_range_mm = [0.1 1.5]
    end
    if nargin < 3
        profile_step_mm = 0.1;
    end
    
    %% get powerspectrum mean
    power_profile = define_filter_settings(data_obj.info,data_obj.ROI,data_obj.data,profile_range_mm,profile_step_mm);
    
    %% determine peak
    [peak_value,ii_peak]=findMaxPeak(power_profile.values);
    peak_position_mm = power_profile.scale_mm(ii_peak);
        
    %% determine power at peak for bootstrap samples
    peak_value_BS = zeros([1 length(data_obj.samples_array)]);
    for ii = 2:length(data_obj.samples_array)
        power_profile = define_filter_settings(data_obj.info,data_obj.ROI,data_obj.data,[peak_position_mm peak_position_mm],peak_position_mm);
        peak_value_BS(ii) = power_profile.values;
    end
    
    %% determine power at peak for shuffeld
    peak_value_Rand = zeros([1 length(data_obj_rand.samples_array)]);
    for ii = 1:length(data_obj_rand.samples_array)
        power_profile = define_filter_settings(data_obj_rand.info,data_obj_rand.ROI,data_obj_rand.data,[peak_position_mm peak_position_mm],peak_position_mm);
        peak_value_Rand(ii) = power_profile.values;
    end
    
    %% determine CI spectrum
    
    %% determine distribution peak
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
