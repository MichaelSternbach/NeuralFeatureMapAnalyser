function filters = find_lowpass(data_obj,data_info,highpass_mm,lowpass_mm,lowpass_cutoffs,average_w,local_w)
%    case {'mouse lemur'}
%        lowpass_mm = 0.2;
%         highpass_mm = 1.2;
%         lowpass_cutoffs = linspace(0.1, 0.79,70);


filters.design_analysis.lowpass = lowpass_mm;
filters.design_analysis.highpass = highpass_mm;
%filters.design_analysis.pw_density = data_info.design.pw_dens;

%% get lowpass with estimated plateau in global density
disp('Obtaining lowpass with global plateau fitting...')

%ROI=ones(data_info.field_size_pix);



% get pinwheel density for filter settings
pw_density_global = zeros(size(lowpass_cutoffs));
%z_base = data_obj.filter_map(data_obj.read_map(base))
for ii = 1:length(lowpass_cutoffs)
    %z_filtered = filter_map(data_obj.read_map,data_info.pixels_per_mm,lowpass_cutoffs(ii),highpass_mm);
    data_obj.set_filter_parameters('lowpass',lowpass_cutoffs(ii))
    z_filtered =data_obj.filter_map(data_obj.read_map);
    pw_density_global(ii) = get_pinwheel_density(z_filtered,local_w.*data_obj.ROI);%ROI,
end

% set axis and fit plateau to global data
x = lowpass_cutoffs.*data_info.pixels_per_mm./average_w;
y = pw_density_global;

% figure
% plot(x,y)

min_length = 0.15*data_info.pixels_per_mm./average_w;
fit_range = [0.2 1.0];

[pw_dens, center_x, plateau] = fit_piecewise_linear_Manuel(x,y,fit_range,min_length,true);

global_lp = center_x/(data_info.pixels_per_mm./average_w);
global_plt = plateau/(data_info.pixels_per_mm./average_w);

% add to data
filters.global_plateau.lowpass = global_lp;
filters.global_plateau.highpass = highpass_mm;

filters.global_plateau.lowpass_vs_density = [lowpass_cutoffs(:) pw_density_global(:)];

filters.global_plateau.pw_density = pw_dens;
filters.global_plateau.plateau = global_plt;
end

