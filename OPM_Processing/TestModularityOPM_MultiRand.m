function [power_profiles,mean_abs_squared,mean_abs_squared2] =TestModularityOPM_MultiRand(data_obj,profile_range_mm,profile_step_mm,N_seeds)


    %% determine power at peak for bootstrap samples
    power_profiles= cell([1 size(data_obj.samples_array,3)]);
    mean_abs_squared = 0;
    mean_abs_squared2 = 0;
    for ii = 1:N_seeds
        rng(ii);
        data = randomizeData(data_obj.data);
        data_obj_rand =  data_handle_corrected(data_obj.info,data,data_obj.ROI);
        %% test1
        z = data_obj_rand.read_map();
        mean_abs_squared = mean_abs_squared + mean(abs(z).^2,'all');
        
        %% test2
        z = makeMap(data,data_obj.info.stim_order);
        mean_abs_squared2 = mean_abs_squared2 + mean(abs(z).^2,'all');

        %% calc power profile
        power_profiles{ii} = define_filter_settings(data_obj.info,data_obj.ROI,z,profile_range_mm,profile_step_mm);
    end
    mean_abs_squared = mean_abs_squared/N_seeds;
    mean_abs_squared2 = mean_abs_squared2/N_seeds;
    

   
        
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



function z = makeMap(data,stim_order)
    data = mean(data,4);
    
    first = true;
    for ii = 1:size(data,3)
        if ~isnan(stim_order(ii))
            if first
                z = data(:,:,ii)*exp(2i*pi*stim_order(ii)/180);
                first = false;
            else
                z = z + data(:,:,ii)*exp(2i*pi*stim_order(ii)/180);
            end
        end
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
