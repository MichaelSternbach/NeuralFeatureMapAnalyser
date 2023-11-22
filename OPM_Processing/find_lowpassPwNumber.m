function filters = find_lowpassPwNumber(data_obj,data_info,highpass_mm,lowpass_mm,lowpass_cutoffs)
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
pw_number = zeros(size(lowpass_cutoffs));
%z_base = data_obj.filter_map(data_obj.read_map(base))
z = data_obj.read_map;
for ii = 1:length(lowpass_cutoffs)
    %z_filtered = filter_map(data_obj.read_map,data_info.pixels_per_mm,lowpass_cutoffs(ii),highpass_mm);
    data_obj.set_filter_parameters('lowpass',lowpass_cutoffs(ii))
    z_filtered =data_obj.filter_map(z);
    %pw_density_global(ii) = get_pinwheel_density(z_filtered,local_w.*data_obj.ROI);%ROI,
    [count,~,~,~,~,~,~] = find_pinwheels(z_filtered,0,data_obj.ROI,0);
    pw_number(ii)=count;
end
data_obj.set_filter_parameters('lowpass',lowpass_mm)
% set axis and fit plateau to global data
x = lowpass_cutoffs;%.*data_info.pixels_per_mm./average_w;
y = pw_number;

% figure
% plot(x,y)

% min_length = 0.15;%*data_info.pixels_per_mm./average_w;
% fit_range = [0.2 1.0];

% [pw_dens, center_x, plateau] = fit_piecewise_linear_Manuel(x,y,fit_range,min_length,false);

% global_lp = center_x;%/(data_info.pixels_per_mm./average_w);
% global_plt = plateau;%/(data_info.pixels_per_mm./average_w);

% % add to data
% filters.global_plateau.lowpass = global_lp;
% filters.global_plateau.highpass = highpass_mm;

filters.global_plateau.lowpass_vs_density = [lowpass_cutoffs(:) pw_number(:)];

% filters.global_plateau.pw_density = pw_dens;
% filters.global_plateau.plateau = global_plt;
end

