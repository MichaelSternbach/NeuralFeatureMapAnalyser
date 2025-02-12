function [peak_values,peak_values_jk,power_profiles,peak_position_mm,mean_abs_squared]=TestModularityOPM(data_obj,profile_range_mm,DoPlot,Jackknife,FullCurve)
        
    %% Input parameter
    if nargin < 2
        profile_range_mm = 0.1:0.1:1.5;
    end
    if nargin < 3
        DoPlot = true;
    end
    if nargin <4
        Jackknife = false; 
    end
    if nargin < 5
        FullCurve = false;
    end
   
    %% set peak
    if length(profile_range_mm) == 2
        
        %% get powerspectrum mean
        z = data_obj.read_map();
        % z = z./mean(abs(z));
        power_profile = define_filter_settings(data_obj.info,data_obj.ROI,z,profile_range_mm);

        %% determine peak
        [peak_value,ii_peak]=findMaxPeak(power_profile.values);
        if ~isnan(ii_peak)
            peak_position_mm = power_profile.scale_mm(ii_peak);
        else
            peak_position_mm = nan;
        end
%         %% plot peak
%         figure; plot(power_profile.scale_mm,power_profile.values)
%         hold on; plot([peak_position_mm peak_position_mm],[min(power_profile.values) max(power_profile.values)])

    elseif isscalar(profile_range_mm)
        peak_position_mm = profile_range_mm;
        ii_peak = 1;
    elseif length(profile_range_mm) > 2
        peak_position_mm = nan;
        peak_value = nan;
        ii_peak = nan;
    end
    
    %% set profile range
    if ~FullCurve||isscalar(profile_range_mm)
        profile_step_mm = peak_position_mm;
        profile_range_mm = [peak_position_mm peak_position_mm];
    end

    %% determine power at peak for bootstrap samples
    peak_values = zeros([1 size(data_obj.samples_array,3)]);
    BS = cell([1 size(data_obj.samples_array,3)]);
    mean_abs_squared = zeros([1 size(data_obj.samples_array,3)]);

    parfor ii = 1:size(data_obj.samples_array,3)
        z = data_obj.read_map(ii);
        mean_abs_squared(ii) = mean(abs(z).^2,'all');
        % z = z./mean(abs(z));
        BS{ii} = define_filter_settings(data_obj.info,data_obj.ROI,z,profile_range_mm);
        if ~isnan(ii_peak)
            peak_values(ii) = BS{ii}.values(ii_peak);
        end
    end
    power_profiles.BS=BS;
    
    if Jackknife
        %% determine power at peak for Jackknife
        data_obj.prepare_jackknife_samples;
        peak_values_jk = zeros([1 data_obj.data_parameters.num_blocks]);
        JK = cell([1 data_obj.data_parameters.num_blocks]);

        parfor ii=1:data_obj.data_parameters.num_blocks
            z = data_obj.read_map(ii);
            % z = z./mean(abs(z));
            JK{ii} = define_filter_settings(data_obj.info,data_obj.ROI,z,profile_range_mm,profile_step_mm);
            if ~isnan(ii_peak)
                peak_values_jk(ii) = JK{ii}.values(ii_peak);
        
            end
        end
        power_profiles.JK=JK;
    else
        peak_values_jk = nan;
    end
   
        
    % elseif length(profile_range_mm) == 1
    % 
    %     %% determine power at peak for bootstrap samples
    %     peak_values = zeros([1 size(data_obj.samples_array,3)]);
    %     for ii = 1:size(data_obj.samples_array,3)
    %         power_profile = define_filter_settings(data_obj.info,data_obj.ROI,data_obj.read_map(ii),[profile_range_mm profile_range_mm],profile_range_mm);
    %         peak_values(ii) = power_profile.values;
    %     end
    % 
    %     if Jackknife
    % 
    %         %% determine power at peak for Jackknife
    %         data_obj.prepare_jackknife_samples;
    %         peak_values_jk = zeros([1 data_obj.data_parameters.num_blocks]);
    %         for ii=1:data_obj.data_parameters.num_blocks
    %             power_profile = define_filter_settings(data_obj.info,data_obj.ROI,data_obj.read_map(ii),[profile_range_mm profile_range_mm],profile_range_mm);
    %             peak_values_jk(ii) = power_profile.values;
    %         end
    %     end
    % end
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
