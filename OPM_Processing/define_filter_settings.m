function power_profile = define_filter_settings(data_info,ROI,data,profile_range_mm,profile_step_mm)
%(experiment_num)
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

% parameters for profile
% profile_range_mm = [0.1 5];
% profile_step_mm = 0.01;

if length(profile_range_mm)>2 || nargin < 5
    profile_scale_mm = profile_range_mm;
else
    profile_scale_mm = profile_range_mm(1):profile_step_mm:profile_range_mm(2);
end

stim_order = data_info.stim_order;
% read meta-data file
%[data_info,data_path] = info_handle(experiment_num);
%load([data_path,'/exp_info.mat'],'ROI')
%trials_to_use = find(data_info.protocol.blocks>0);

% make space for data
% for trial = 1:length(data_info.protocol.blocks)
%     power_profile(trial).scale_mm = [];
%     power_profile(trial).values = [];
% end

power_profile.scale_mm = [];
power_profile.values = [];

% get profile from trials
% for trial = trials_to_use
%     
%     load([data_path,'Processed_2/trial_',num2str(trial),'.mat'],'data');
%     map = make_map(data,data_info.stim_order,ROI,true);
%     map(~ROI) = 0;
%     
%     power_profile(trial).values = calculate_power_profile(map,profile_scale_mm*data_info.pix_per_mm);
%     power_profile(trial).scale_mm = profile_scale_mm;
% end
%load([data_path,data_info.ID,'.mat'],'data');
if ndims(data)< 3
    map = data;
else
    map = make_map(data,stim_order,ROI,false);
end
    %map(~ROI) = 0;


power_profile.values = calculate_power_profile(map,profile_scale_mm*data_info.pix_per_mm);
power_profile.scale_mm = profile_scale_mm;



%% normalize powerspectrum
power = mean(abs(map).^2);
power_profile.k_mm_inv = 1./power_profile.scale_mm;
scale_mm_fine = linspace(min(power_profile.k_mm_inv), max(power_profile.k_mm_inv), 1000);  % Create a fine grid over the range of x
values_fine = interp1(power_profile.k_mm_inv, power_profile.values, scale_mm_fine, 'linear');  % Linear interpolation on the fine grid
integral_power = trapz(scale_mm_fine, values_fine); 

power_profile.values_kspace = power_profile.values/integral_power*power;



%% calculate power spectrum



% make figure
% figure
% hold on

% for trial = trials_to_use
%     plot(power_profile(trial).scale_mm,power_profile(trial).values);
% end
% plot(power_profile.scale_mm,power_profile.values);
% xlabel('Scale in mm')
% ylabel('Power')
% xlim([0 max(profile_scale_mm)])
% set(gca,'fontsize',15)

% save
%save([data_path,'/exp_info.mat'],'power_profile','-append')
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
%PS = PS/mean(PS(:));

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
