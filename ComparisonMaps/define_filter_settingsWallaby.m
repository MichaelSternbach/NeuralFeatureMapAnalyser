function define_filter_settingsWallaby(data_set,set_ID,data_obj)
% Finding appropiate filter settings is crucial in the data analysis.
% This function calculates the radial profile of the power spectrum of the 
% maps in each trial to use as a reference to set initial filter cut-offs. 
% Since the cutt-offs are defined in mm, the function calculates the
% profile scale in mm instead of k (wavenumber).
% The results are all plotted in a figure and saved as power_profile 
% in exp_info.mat.
% When the map is of good quality, the power is low for small mm scales,
% has a few large peaks corresponding to the typical scale of the map around
% 0.6-1.0 mm, decreases afterwards and starts fluctuating depending on the
% larger structures of the layout.
% A good highpass is where the first peaks have dropped and before the power
% rises again.
% A good lowpass is where the first peak starts raising for all maps.
% The lowpass is crucial for the pinwheel position, so this initial
% estimate will be refined in a second step later on.
% The defined cut-offs should be saved in make_info_files.m under settings:
% e.g.
% data_info.settings.lowpass_mm = 0.47;
% data_info.settings.highpass_mm = 1.5;
%%

N_Samples = 10;
% parameters for profile
% profile_range_mm = [0.0001 0.1];
% profile_step_mm = 0.0001;
% profile_scale_mm = profile_range_mm(1):profile_step_mm:profile_range_mm(2);

profile_range_mm = [0.01 2];
profile_step_mm = 0.01;
profile_scale_mm = profile_range_mm(1):profile_step_mm:profile_range_mm(2);

% read data files   %%% read meta-data file
[data_info,data_path] = info_handle(data_set,set_ID);
ROI=ones(data_info.field_size_pix);%load([data_path,'/exp_info.mat'],'ROI')


%load([data_path,data_info.ID,'.mat'],'data');
data_obj.prepare_samples_array(N_Samples)
% map = data_obj.read_map();%make_map(data,data_info.stim_order,ROI,true);
% map(~ROI) = 0;

% trials_to_use = find(data_info.protocol.blocks>0);
% 
% % make space for data
for ii_sample = 1:N_Samples
    power_profile(ii_sample).scale_mm = [];
    power_profile(ii_sample).values = [];
end

% get profile from trials
% for trial = trials_to_use
%     
%     load([data_path,data_info.ID,'.mat'],'data');
%     map = make_map(data,data_info.stim_order,ROI,true);
%     map(~ROI) = 0;
%     
%     power_profile(trial).values = calculate_power_profile(map,profile_scale_mm*data_info.pix_per_mm);
%     power_profile(trial).scale_mm = profile_scale_mm;
% end

for ii_sample = 1:N_Samples
    map = data_obj.read_map(ii_sample);%make_map(data,data_info.stim_order,ROI,true);
    map(~ROI) = 0;
    power_profile(ii_sample).values = calculate_power_profile(map,profile_scale_mm*data_info.pixels_per_mm);
    power_profile(ii_sample).scale_mm = profile_scale_mm;
end
% make figure
figure
hold on

% for trial = trials_to_use
%     plot(power_profile(trial).scale_mm,power_profile(trial).values);
% end
for ii_sample = 1:N_Samples
    plot(power_profile(ii_sample).scale_mm,power_profile(ii_sample).values);
end

xlabel('Scale in mm')
ylabel('Power')
xlim([0 max(profile_scale_mm)])
set(gca,'fontsize',15)

% save
save([data_path,'/exp_info.mat'],'power_profile','-append')
end

function profile = calculate_power_profile(map,profile_scale_pixels)
% This function makes a radial average of the power spectrum of a map to
% get a spectrum profile.
%%

% initiate variable
profile = zeros(length(profile_scale_pixels),1);

% calculate padding
pad_size_pixels = 2^nextpow2(max(size(map)));

% calculate power spectrum
PS = abs(fftshift(fft2(map,pad_size_pixels,pad_size_pixels))).^2;
PS = PS/mean(PS(:));

% create function
[kx,ky] = meshgrid(-pad_size_pixels/2:pad_size_pixels/2-1,-pad_size_pixels/2:pad_size_pixels/2-1);
PowerFun = scatteredInterpolant(kx(:),ky(:),PS(:));

% calculate radial profile
theta=linspace(0,2*pi,101);theta(end) = [];
for ind = 1:length(profile_scale_pixels)
    r = pad_size_pixels/profile_scale_pixels(ind);
    profile(ind) = mean(PowerFun(bsxfun(@times,r,cos(theta)),bsxfun(@times,r,sin(theta))));
end

end
